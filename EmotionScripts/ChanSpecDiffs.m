% Compares power spectrum between two conditions at two different CHANNELS
%
%  ChanSpecDiffs(datset,datpath,events,chans,frqlim);
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


function [contrib,speccomp,cics,p,fbins] = ChanSpecDiffs(datset,datpath,events,chans,gdcomps,artifact,frqlim,freq);
    
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
            %%% collect 'comps' first
            EEG.icaact = [];EEG = eeg_checkset(EEG);
            multfact = sqrt(mean(EEG.icawinv.^2));
            for ch = 1:length(chans)
                multfact = EEG.icawinv(chans(ch),:);                
                for cp = 1:length(multfact)
                    allics{ds,ch}(cp,:,:) = EEG.icaact(cp,:,:)*multfact(cp); 
                end;
            end;

            %EEG = pop_select(EEG,'nochannel',[151 152]); % take out heart for hf45
            
% $$$             for cp = 1:size(EEG.icawinv,2)
% $$$                 [pvaf(1,cp),pvafs(cp,:),vars] = eeg_pvaf(EEG,cp,'artcomps',artifact,'fraction',1,'chans',chans,'plot','off');
% $$$             end;
% $$$             for ch = 1:size(pvafs,2)
% $$$                 gdvaf = pvafs(gdcomps,ch);
% $$$                 figure; row = round(sqrt(length(gdcomps))); col = ceil(sqrt(length(gdcomps)));
% $$$                 for cp = 1:length(gdcomps)
% $$$                     sbplot(row,col,cp)
% $$$                     topoplot(EEG.icawinv(:,gdcomps(cp)),EEG.chanlocs);
% $$$                     title(num2str(round(gdvaf(cp)*100)/100));
% $$$                 end;
% $$$                 textsc(['Percent Variance Accounted for at Channel ',int2str(chans(ch))],'title');
% $$$             end;
            % for first channel:
            %[spectra(1,:,ds),freqs,speccomp{ds}(:,:,1),contrib(1,:,ds),specstd] = spectopo(EEG.data(chans(1),:,:),size(EEG.data,2), EEG.srate, 'freqrange',[3 25],'nfft',EEG.srate,'winsize',EEG.srate,'overlap',EEG.srate/4, 'plot','off', 'icacomps',[1:size(EEG.icawinv,2)],'weights',EEG.icaweights(:,chans(1)), 'icawinv',EEG.icawinv(chans(1),:) ,'chanlocs',EEG.chanlocs, 'freq', freq);
            % for second channel:
            %[spectra(2,:,ds),freqs,speccomp{ds}(:,:,2),contrib(2,:,ds),specstd] = spectopo(EEG.data(chans(2),:,:),size(EEG.data,2), EEG.srate, 'freqrange',[3 25],'nfft',EEG.srate,'winsize',EEG.srate,'overlap',EEG.srate/4, 'plot','off', 'icacomps', [1:size(EEG.icawinv,2)],'weights',EEG.icaweights(:,chans(2)), 'icawinv',EEG.icawinv(chans(2),:) ,'chanlocs',EEG.chanlocs, 'freq', freq);
            
            % now remove artifacts and take channels
           %if ~isempty(artifact)
           %     EEG = pop_subcomp( EEG,artifact, 0);
           % end;        
            alldat{ds}(1,:,:) = EEG.data(chans(1),:,:);
            alldat{ds}(2,:,:) = EEG.data(chans(2),:,:);
            
           if ~isempty(artifact)
               EEG = pop_subcomp( EEG,artifact, 0);
            end;        
            alldatclean{ds}(1,:,:) = EEG.data(chans(1),:,:);
            alldatclean{ds}(2,:,:) = EEG.data(chans(2),:,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        end;
    else
        EEG = pop_loadset(datset, fullpaths{nx});
        % insert option to get two time stretches from single datset (requires event #s)
        ALLEEG(1) = pop_select(EEG,'point',events{1});
        ALLEEG(2) = pop_select(EEG,'point',events{2});
        
        alldat{ds}(1,:,:) = EEG.data(chans(1),:,:);
        alldat{ds}(2,:,:) = EEG.data(chans(2),:,:);

        %%>>>>> Not finished!!
    end;
            
    
    %%%%%%%%%%%%%%%%  Begin FFT of single epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nCalculating pwr for each 2 second epoch for each selected channel...\n');
    for ds = 1:length(alldat)
        eplength = size(alldat{ds},2);
        for ch = 1:size(alldat{ds},1)
            for epoch = 1:size(alldat{ds},3)
                [pwr(epoch,:,ch) frqs] = pwelch(alldatclean{ds}(ch,:,epoch),EEG.srate,EEG.srate/2,EEG.srate*6,EEG.srate); 
                pwr(epoch,:,ch) = 10*log10(pwr(epoch,:,ch));% convert to log 
            end;   
        end;
        mnpwr{ds} = pwr;    
        mnpwr{ds}= pwr;
    end;
    
    figure;  lw = 2;
    set(gca,'fontsize',16);
    plot(frqs,mean(mnpwr{1}(:,:,1),1),'b','linewidth',lw); hold on;
    plot(frqs,mean(mnpwr{1}(:,:,2),1),'c','linewidth',lw);
    plot(frqs,mean(mnpwr{2}(:,:,1),1),'r','linewidth',lw);
    plot(frqs,mean(mnpwr{2}(:,:,2),1),'m','linewidth',lw);
    set(gca,'xlim',[frqlim(1) frqlim(2)]);
    legend({[datset{1}(1:end-4),',Chan ',int2str(chans(1))], [datset{1},',Chan ',int2str(chans(2))],[datset{2},',Chan ',int2str(chans(1))], [datset{2},',Chan ',int2str(chans(2))]})
    title(['Mean Clean Power Spectra from Subj ',datpath(end-4:end-1)]);

    ntrials(1,1) = size(mnpwr{1},1);
    ntrials(1,2) = size(mnpwr{2},1);
    dectrials = min(ntrials);
    
    fbins = {[2 5],[5 7],[7 9],[9 11],[11 15],[15 18],[18 22],[22 32],[32 42]};
    for f = 1:length(fbins)
        fr = find(frqs > fbins{f}(1) & frqs < fbins{f}(2));
        statmat = zeros(0,2); % make a trials*datasets x left/right matrix
        for ds = 1:length(mnpwr)
            for ch = 1:size(mnpwr{ds},3)
                spwr = sort(mean(mnpwr{ds}(:,fr,ch),2));
                tpvec(:,ch) = spwr(round([.5:dectrials]*size(mnpwr{ds},1)/dectrials));
            end;
            statmat(end+1:end+size(tpvec,1),:) = tpvec;
        end;        
        [p(f,:),table,stats] = anova2(statmat,dectrials); 
        close
    end;
    
    keyboard
% $$$     for ds = 1:size(contrib,3)
% $$$         for ch = 1:size(contrib,1)
% $$$             [val contribic] = sort(contrib(ch,gdcomps,ds));
% $$$             contribic = contribic(end-3:end);
% $$$             contribvals = val(end-3:end);            
% $$$             cics(ch,:,ds) = gdcomps(contribic);
% $$$             cvals(ch,:,ds) = contribvals;
% $$$         end;
% $$$     end;
% $$$     combics1 = union(cics(1,:,1),cics(1,:,2));
% $$$     combics2 = union(cics(2,:,1),cics(2,:,2));
% $$$     combics = union(combics1,combics2);
    %%%%%%%%%%%%%%%%  Begin decompostion of ICs epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nCalculating pwr for each 2 second epoch for all gdcomps...\n');
    for ds = 1:2%size(cics,3)
        for ch = 1:size(allics,2)
            eplength = size(allics{ds,ch},2);
            for icc = 1:length(combics)
                ic = combics(icc);
                for epoch = 1:size(allics{ds},3)
                    [icpwr(epoch,:,icc) frqs] = pwelch(allics{ds,ch}(ic,:,epoch),EEG.srate,EEG.srate/2,EEG.srate*6,EEG.srate); 
                    icpwr(epoch,:,icc) = 10*log10(icpwr(epoch,:,icc));% convert to log 
                end;   
            end;
            mnicpwr{ds,ch} = squeeze(mean(icpwr,1))'; % makes an ic x freqs matrix   
        end;
    end;
    
    
    lw = 2; cols = lines(size(mnicpwr{2,1},1));
    figure; row = 2; col = 2; pl = 1;
    for ds = 1:size(mnicpwr,1)
        for ch = 1:2%size(cics,1)            
            sbplot(row,col,pl); leg = cell(1,0);pl = pl+1;
            set(gca,'fontsize',14);
           for icc = 1:20%length(combics)
                ic = icc;
                ph = plot(frqs,mnicpwr{ds,ch}(ic,:),'b','linewidth',lw); hold on;
                set(ph,'color',cols(ic,:));
                leg{end+1} = ['IC ',int2str(combics(ic))];
            end;
            ph = plot(frqs,mean(mnpwr{ds}(:,:,ch),1),'k','linewidth',lw+1); hold on;
            set(gca,'xlim',[frqlim(1) frqlim(2)]);
            leg{end+1} = ['Chan'];
            if pl == 5
                legend(leg)
            end;
            title([datset{ds}(1:end-4),', chan ',int2str(chans(ch))]);
        end;
    end;        
    textsc(['Subj ',datpath(end-4:end-1),':  ICs contributing to 10 Hz power'],'title');
    set(gcf,'PaperOrientation','landscape'); set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    
    lw = 2; cols = lines(size(mnicpwr{2,1},1));
    figure; row = 2; col = 2; pl = 1;
    for ch = 1:2%size(cics,1)            
        sbplot(row,col,pl); leg = cell(1,0);pl = pl+1;
        set(gca,'fontsize',14);
        ph = plot(frqs,mean(mnpwr{2}(:,:,ch),1)-mean(mnpwr{1}(:,:,ch),1),'k','linewidth',lw+1); hold on;
        set(gca,'xlim',[frqlim(1) frqlim(2)]);
        leg{end+1} = ['Chan'];
        
        for icc = 1:10%length(combics)
            ic = icc;
            ph = plot(frqs,mnicpwr{2,ch}(ic,:)-mnicpwr{1,ch}(ic,:),'b','linewidth',lw); hold on;
            set(ph,'color',cols(ic,:));
            leg{end+1} = ['IC ',int2str(combics(ic))];
        end;
        title([datset{ds}(1:end-4),', chan ',int2str(chans(ch))]);
    end;
    sbplot(2,2,3);   
    for icc = 1:length(combics)
        ic = icc;
        ph = plot(frqs,mnicpwr{2,ch}(ic,:)-mnicpwr{1,ch}(ic,:),'b','linewidth',lw); hold on;
        set(ph,'color',cols(ic,:));
        leg{end+1} = ['IC ',int2str(combics(ic))];
    end;
    legend(leg)
    textsc(['Subj ',datpath(end-4:end-1),':  Happy-Sad Power Spectra'],'title');
    set(gcf,'PaperOrientation','landscape'); set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    % Plot ic pwr diffs
    figure;row = 6; col = 6;
    for ic = 1:36%length(combics)
    sbplot(row,col,ic);   
        %ph = plot(frqs,mnicpwr{2,ch}(ic,:)-mnicpwr{1,ch}(ic,:),'b','linewidth',lw); hold on;
        ph = plot(frqs,mnicpwr{2,1}(ic,:)-mnicpwr{1,1}(ic,:),'g','linewidth',lw); hold on;
        ph = plot(frqs,mnicpwr{2,2}(ic,:)-mnicpwr{1,2}(ic,:),'g','linewidth',lw); hold on;
        set(ph,'color',[1 .6 0]); 
        set(gca,'xlim',[3 45]);set(gca,'xticklabel',[]);
        plot([get(gca,'xlim')],[0 0],'k-'); title(['IC ',int2str(ic)]);
        %set(gca,'ylim',[-2 2]);
        %set(ph,'color',cols(ic,:));
        %leg{end+1} = ['IC ',int2str(combics(ic))];
    end;

    % plot on IC of interest at a time
    fr = find(frqs>3&frqs<45);
    ic = 16; 
    figure; 
    ph = plot(frqs(fr),mnicpwr{1,1}(ic,fr),'b-','linewidth',2.5); hold on;
    ph = plot(frqs(fr),mnicpwr{1,2}(ic,fr),'c-','linewidth',2.5);
    ph = plot(frqs(fr),mnicpwr{2,1}(ic,fr),'r-','linewidth',2.5);
    ph = plot(frqs(fr),mnicpwr{2,2}(ic,fr),'m-','linewidth',2.5);
    set(gca,'fontsize',20);
    set(gca,'xlim',[3 45]); set(gca,'box','off');
        
    alfr = [61,62,63];
    emo = 2;    
    HmnalphaL = mean(mnicpwr{emo,1}(ic,alfr)); % left
    HmnalphaR = mean(mnicpwr{emo,2}(ic,alfr)); % right
    emo = 1;    
    GmnalphaL = mean(mnicpwr{emo,1}(ic,alfr));% left
    GmnalphaR = mean(mnicpwr{emo,2}(ic,alfr));% right
    
