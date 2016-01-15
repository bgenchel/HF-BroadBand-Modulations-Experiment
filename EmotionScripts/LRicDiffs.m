% Compares power spectrum between two conditions at two different CHANNELS
%
%  LRicDiffs(datset,datpath,events,chans,frqlim);
%
% INPUTS:
% datset -- [string or cell array of strings] if single dataset, then will 
%           look in same dataset for 'events', which must be indices in this case
% datpath -- [string] full directory path where 'datset' can be found
% events -- [cell array of strings or integers] can be strings, if taking
%           data from 2 different datasets, or event indices. Data will be
%           taken from between these events for spectral decomposition.
% chans -- [vector] two channel indices from which to take EEG data
% gdcomps -- [vector] of IC indices corresponding to dipolar brain components
% artifact -- [vector] of ICs to remove from data before spectral decomp
% frqlim -- [minfrq maxfrq] in Hz for plotting purposes
% freq -- [integer in Hz] frequency to find IC contributions to 'chans' spectra
%


function [contrib,speccomp,cics,p,fbins] = LRicDiffs(datset,datpath,events,Lics,Rics,frqlim,freq);
    
    
    gdcomps = sort([Lics,Rics]);
    if iscell(datset)
        for ds = 1:length(datset)
            EEG = pop_loadset(datset{ds}, datpath);

            if iscell(events{ds})
                for ev = 1:length(EEG.event)
                    if strcmp(EEG.event(ev).type,events{ds}{1})
                        evlat(1,1) = EEG.event(ev).latency;
                    elseif strcmp(EEG.event(ev).type,events{ds}{2})
                        evlat(1,2) = EEG.event(ev).latency;
                    end;
                end;
                if isempty(events{ds}{2})
                    evlat(1,2) = size(EEG.data,2); % if event #2 is empty take end of data
                end;
            end;
            EEG = pop_select(EEG,'point',evlat);
            
            %%%%%%%%%%  Divide continuous data into 2 sec epochs:
            overlap = 4;
            x = EEG.event(1).type;
            for evtm = EEG.srate*2:EEG.srate/overlap:size(EEG.data,2)-EEG.srate*2 
                % create events to make overlapping 2 sec epochs
                EEG.event(end+1) =  EEG.event(1);% appends events to the end
                EEG.event(end).latency = evtm;
                if ischar(x)
                    EEG.event(end).type = 'fake';% for string event codes  
                else
                    EEG.event(end).type = 1000;% for string event codes 
                end;
            end;
            if ischar(x)
                EEG = pop_epoch( EEG,{'fake'} , [-1 1]);
            else
                EEG = pop_epoch( EEG,{1000} , [-1 1]);
            end;        
            EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
            EEG = eeg_checkset(EEG);
            %%%%%%%  Run automatic rejection on newly epoched data:
            %-------------------------------------------------------
            % probability and kurtosis are based on standard deviations, 
            % so iterative rejection is useful
            
            fprintf('\nRunning auto-rejection protocol...\n');
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)] ,-1000,1000,EEG.xmin,EEG.xmax,0,0);
            numrej = length(find(EEG.reject.rejthresh));  % count number of epochs marked
            if numrej > 0
                EEG = pop_rejepoch( EEG,EEG.reject.rejthresh,0); % actually reject high prob epochs
                fprintf('\nRe-baselining after large amplitude artifact removed...\n');
                EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
            end;
            %--------------------------------------------
            %% run probability test
            startlim = 5; % start probability threshold that will be up'ed if too many trials found
            maxrej = 40;
            EEG = pop_jointprob(EEG,0,gdcomps ,startlim,startlim,0,0);% calculate component probabilities
            numrej = length(find(EEG.reject.icarejjp));  % count number of epochs marked
            if numrej < maxrej
                EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0); % actually reject high prob epochs
            else
                fprintf('Re-adjusting probability limits and running again...*************************\n');
                startlim = startlim + .5;
            end;                    
            repeat = 1; maxiter = 0; 
            while repeat == 1 % keep running probability until there are no epochs above threshold
                if numrej > 0
                    EEG = pop_jointprob(EEG,0,gdcomps ,startlim,startlim,0,0);                    
                    numrej = length(find(EEG.reject.icarejjp));
                    if numrej < maxrej
                        EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0);
                    else
                        startlim = startlim + .5; EEG.reject.icarejjp = [];EEG.reject.rejjpE = [];
                        fprintf('Re-adjusting probability limits and running again...*************************\n');
                    end;                    
                else
                    if startlim > 5 & maxiter < 8 % don't decrease and startover more than 8 times
                        fprintf('Decreasing probability limits for final pruning...###########################\n');
                        startlim = startlim - .5; numrej = 1; maxiter = maxiter+1; % repeat process back to 5 stds
                    else
                        if maxiter > 8
                            maxrej = 100; % go through last round with a high threshold
                        else                
                            repeat = 0;
                        end;              
                    end;            
                end; 
            end;       
            % run kurtosis check
            EEG = pop_rejkurt(EEG,0,gdcomps ,6,6,0,0);        
            numrej = length(find(EEG.reject.icarejkurt));  % count number of epochs marked
            if numrej > 0
                EEG = pop_rejepoch( EEG,EEG.reject.icarejkurt,0);
            end;
            %--------------------------------------------        
            EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
            %%% collect RMS IC acts first
            EEG.icaact = [];EEG = eeg_checkset(EEG);
            multfact = sqrt(mean(EEG.icawinv(:,gdcomps).^2));
            for cp = 1:length(multfact)
                allics{ds}(cp,:,:) = EEG.icaact(cp,:,:)*multfact(cp); 
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        end;
    else % for two time stretches in one dataset
% $$$         EEG = pop_loadset(datset, fullpaths{nx});
% $$$         % insert option to get two time stretches from single datset (requires event #s)
% $$$         ALLEEG(1) = pop_select(EEG,'point',events{1});
% $$$         ALLEEG(2) = pop_select(EEG,'point',events{2});
% $$$         
% $$$         alldat{ds}(1,:,:) = EEG.data(chans(1),:,:);
% $$$         alldat{ds}(2,:,:) = EEG.data(chans(2),:,:);

        %%>>>>> Not finished!!
    end;
    %%%%%%%%%%%%%%%%  Begin decompostion of ICs epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nCalculating pwr for each 2 second epoch for all gdcomps...\n');
    for ds = 1:length(allics)
        for cp = 1:size(allics,1)
            eplength = size(allics{ds},2);
            for ic = 1:length(gdcomps)
                for epoch = 1:size(allics{ds},3)
                    [icpwr(epoch,:,ic) frqs] = pwelch(allics{ds}(ic,:,epoch),EEG.srate,EEG.srate/2,EEG.srate*4,EEG.srate); 
                    icpwr(epoch,:,ic) = 10*log10(icpwr(epoch,:,ic));
                end;   
            end;
            mnicpwr{ds} = icpwr; % makes an ic x freqs matrix
        end;
    end;

    
    ntrials(1,1) = size(allics{1},3);
    ntrials(1,2) = size(allics{2},3);
    dectrials = min(ntrials);

    
    fbins = {[2 5],[5 7],[7 9],[9 11],[11 15],[15 18],[18 22],[22 32],[32 42]};
    fr = find(frqs < frqlim(2));
    for cp1 = 1:length(gdcomps)-1
        for cp2 = cp1+1:length(gdcomps)
            for f = 1:length(fbins)
                fr = find(frqs > fbins{f}(1) & frqs < fbins{f}(2));
                statmat = zeros(0,2); % make a trials*datasets x left/right matrix
                for ds = 1:length(mnicpwr)
                    spwr = sort(mean(mnicpwr{ds}(:,fr,cp1),2));
                    tpvec(:,1) = spwr(round([.5:dectrials]*size(mnicpwr{ds},1)/dectrials));
                    spwr = sort(mean(mnicpwr{ds}(:,fr,cp2),2));
                    tpvec(:,2) = spwr(round([.5:dectrials]*size(mnicpwr{ds},1)/dectrials));
                    statmat(end+1:end+size(tpvec,1),:) = tpvec;
                end;
                [p(f,:),table,stats] = anova2(statmat,dectrials);  close
            end; 
            allps{cp1,cp2} = p;
        end;
        fprintf('.');
    end;
    
    
    %%  find ICs with sig interaction with emotion
    frrange =  find(frqs > 9 & frqs < 11); % alpha
    alpha = .00000000001; pl = 1; clear sigpairs
    for cp1 = 1:length(gdcomps)-1
        for cp2 = cp1+1:length(gdcomps)
            if allps{cp1,cp2}(4,3) < alpha % 4th bin, interaction term
                sigpairs(pl,:) = [gdcomps(cp1) gdcomps(cp2)]; pl = pl+1;
            end;            
% $$$             manysig = 0;
% $$$             for f = frrange(1):frrange(end)
% $$$                 if allps{cp1,cp2}(f,3) < alpha
% $$$                     manysig = manysig+1;
% $$$                 end;
% $$$             end;
% $$$             if manysig > round(3*(length(frrange)/4))
% $$$                 sigpairs(pl,:) = [gdcomps(cp1) gdcomps(cp2)]; pl = pl+1;
% $$$             end;
        end;
    end;
    %%  plot power diffs between datasets
    figure;
    fr =  find(frqs > 9 & frqs < 11); % alpha
    for ic = 1:length(gdcomps)
        onepwr1 = mean(mnicpwr{1}(:,fr,ic),2);
        onepwr1 = mean(onepwr1,1);
        onepwr2 = mean(mnicpwr{2}(:,fr,ic),2);
        onepwr2 = mean(onepwr2,1);
        diff = onepwr2 - onepwr1;
        ph = bar(ic,diff); hold on;
        if diff < 0
            set(ph,'facecolor','b');
        else
            set(ph,'facecolor','r');
        end;
    end;
    set(gca,'fontsize',16);
    set(gca,'xtick',[1:length(gdcomps)]);
    set(gca,'xticklabel',gdcomps);
    xlabel('Independent Components');
    ylabel(['Alpha power difference (',datset{2}(1:end-4),' - ',datset{1}(1:end-4),') (dB)']);
    set(gca,'xlim',[0 length(gdcomps)+1]);
        
    

    figure; row = 3; col = 3; pl = 1;
    lw = 2; cols = lines(14);
    for c1 = 1:length(Lics)   
        leg = cell(1,0);              
        cp1 = find(gdcomps == Lics(c1));
        sbplot(row,col,pl); pl = pl+1;
        ph = plot(frqs,mean(mnicpwr{1}(:,:,cp1),1),'b','linewidth',lw); hold on;
        leg{end+1} = [datset{1},' IC ',int2str(gdcomps(cp1))];
        ph = plot(frqs,mean(mnicpwr{2}(:,:,cp1),1),'r','linewidth',lw); hold on;
        leg{end+1} = [datset{2},' IC ',int2str(gdcomps(cp1))];
        for c2 = 1:length(Rics)            
            cp2 = find(gdcomps == Rics(c2));
            ph = plot(frqs,mean(mnicpwr{1}(:,:,cp2),1),'c','linewidth',lw); hold on;
            %leg{end+1} = [datset{1},' IC ',int2str(gdcomps(cp2))];
            ph = plot(frqs,mean(mnicpwr{2}(:,:,cp2),1),'m','linewidth',lw); hold on;
            %leg{end+1} = [datset{2},' IC ',int2str(gdcomps(cp2))];
            set(gca,'xlim',[frqlim(1) frqlim(2)]);
        end;
        title(['IC: ',int2str(gdcomps(cp1)),'vs all R']);
        legend(leg);                set(gca,'fontsize',10);
    end;
    figure; ds = 1;
    lw = 2; cols = lines(length(gdcomps));
        leg = cell(1,0);              
    for ic = 1:length(gdcomps)          
        ph = plot(frqs,mean(mnicpwr{ds}(:,:,ic),1),'b','linewidth',lw); hold on;
        set(ph,'color',cols(ic,:));
        leg{end+1} = [' IC ',int2str(gdcomps(ic))];        
    end;
    set(gca,'xlim',[frqlim(1) frqlim(2)]);
    legend(leg);                set(gca,'fontsize',16);
    title(datset{ds});
 
    figure;  row = round(sqrt(length(gdcomps))); col = ceil(sqrt(length(gdcomps))); 
    for cp1= 1:length(gdcomps)
        sbplot(row,col,cp1)
        leg = cell(1,0);          
        ph = plot(frqs,mean(mnicpwr{1}(:,:,cp1),1),'b','linewidth',lw); hold on;
        leg{end+1} = [datset{1},' IC ',int2str(gdcomps(cp1))];
        ph = plot(frqs,mean(mnicpwr{2}(:,:,cp1),1),'r','linewidth',lw); hold on;
        leg{end+1} = [datset{2},' IC ',int2str(gdcomps(cp1))];
        set(gca,'xlim',[frqlim(1) frqlim(2)]);
        %set(gca,'ylim',[-35 8]);
        title(['IC ',int2str(gdcomps(cp1))]);
        %legend(leg);            
    end;
    
    
