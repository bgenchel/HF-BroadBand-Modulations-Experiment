% plot all spectral factor envelopes on a single component mean spectrum
% Takes the 95th percentile of trials to exclude outliers...
%
% [freqout] = PlotSpecFacEnv(datset,savedat,fullpath,templates,comps,trialset,side,frqlim,mps,perc,mnfac);
% 
% INPUTS:
% datset -- [string] name of a .set file used to plot scalp maps for each component
% savedat -- [string] name of data file with parameters and name of .wts and .sph files
% fullpath -- [string] full data path where datset and savedat can be found
% templates -- [vector] of factor indices to plot for each component
% comps -- [vector] of component indices to plot (must be a subset of 'origlist'
% trialset -- [number OR vector] If not [], will plot MEAN back-proj of only selected trials(vector of trial indices) OR number of 'dstrial' block to plot. (set 'mnfac' to 1)
% side -- [ 1 | -1 | [] ] 1 will plot only positive-weighted trials, -1: negative trials. 
%         [] plots pos/neg back-projection.  
%         OR enter a vector of -1s and 1, one for each modulator. This will turn off envelope. 
%         Use 'mnfac' == 1 to plot mean of selected trials (will override using mean of all trials).
% frqlim -- [minfr maxfr] in Hz to plot
% mps -- [0 | 1] if 1, then will plot component scalp maps, otherwise only spectra
% perc -- [decimal] percent (from 0-1) of data to plot (total data, pca reduced and factor back-projections)
% mnfac -- [0 | 1] if 1, will plot the mean factor back-projections instead of the outer envelopes

function [freqout] = PlotSpecFacEnv(datset,savedat,fullpath,templates,comps,trialset,side,frqlim,mps,perc,mnfac);
    
    rawbkgd = [.89 .89 .89]; % .89 .89 .89 color of raw data background, [] to turn off
    pcabkgd = [.75 .75 .75];%  .75 .75 .75color of pca reduced data background
    lnwdth = 2.2;%2.2
    
    cols(1,:) = [0 .8 .8];
    cols(2,:) = [1 0 1];
    cols(3,:) = [1 0 0];
    cols(4,:) = [0 0 1];
    cols(5,:) = [0 1 0];
    cols(6,:) = [1 .4 0]; % orange
    cols(7,:) = [.6 0 .6]; % purple
    cols(8,:) = [1 .4 .4]; %light red
    cols(9,:) = [.4 .6 0]; % army green
    cols(10,:) = [1 .8 0]; % mustard
    
    % or, if you run out of the above:
    %cols = colorcube(64);
    
%cols(2,:) = [1.0000    0.5625         0]; % changes yellow/green to orange for 5 IM plots
    %%% this part is for colorcube colors:----------------------
    %cols = colorcube(64);
    %useset = [2,3,4,5,8,9,10,11,13,17,19,21,23,25,27,30,32,37,41,44,54];
    %useset = shuffle(useset);
    %if length(useset) < length(templates)
    %    useset = [1:64];
    %    useset = shuffle(useset);
    %end;
    %newset = useset([1:length(templates)]);
    %cols = cols(newset,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if exist('mps') == 0
        mps = 0;
    end;
    
    if mps == 1
        EEG = pop_loadset(datset,fullpath); 
    end;
    
    s = load([fullpath,savedat,'.mat']); 
    if size(s.freqs,2) < 2
        s.freqs = s.freqs';
    end;    
    if isempty(comps)
        comps = s.complist;
    end;

    figure; 
    if mps == 1
        row = ceil(sqrt(length(comps)*3));
        col = ceil(sqrt(length(comps)*3));
        rep = 1;
        while rep == 1
            if col ==3 | col==6 | col ==9 | col==12
                rep = 0;
            else                
                col = col+1; % 
                if length(comps)*3 <= (row-1)* col
                    row = row -1;
                end;
            end;
        end;   
    else    
        row = round(sqrt(length(comps)));
        col = ceil(sqrt(length(comps)));
    end;
    
    fprintf('\n Loading spectral decomposition weights and data...');
    sph=floatread([fullpath,savedat,'.sph'],[s.numrows s.numrows],[],0); 
    wts=floatread([fullpath,savedat,'.wts'],[s.pcs s.numrows],[],0); 
    data = floatread([fullpath,savedat,'.fdt'],[s.numrows s.numframes],[],0);    
    origdat = floatread([fullpath,savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);    
    ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws 
    speceig = floatread([fullpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  
    winv = specwts;    
    
    if ~isempty(trialset)
        if length(trialset) == 1 % if specified by 'dstrials'
            newset = [sum(s.dstrials(1:trialset-1))+1:sum(s.dstrials(1:trialset))];
            trialset = newset;
        end;
    end;
    if isempty(templates)
        templates = [1:size(activations,1)];
        lnwdth = lnwdth*.5;% if all, decrease linewidth
        cols = hsv(length(templates)); % determine colors for IM traces
    end;
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(2));
    freqs = s.freqs(fr);
    pl = 1; 

    for cp = 1:length(comps)
        rcp = find(comps(cp) == s.complist);
        if mps == 1
            sbplot(row,col,pl)
            topoplot(EEG.icawinv(:,comps(cp)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off');
            title(int2str(comps(cp)));
            pl = pl+1;
        end;        
        plotdata = origdat(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp) + repmat(s.meanpwr(rcp,:),[size(origdat,1) 1]);
        plotdata = plotdata(:,fr);
        %sbplot(row,col,[pl pl+1]);         
         sbplot(row,col,pl);         
       %%%%%%%%%%%%%%  Plot full data %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        plotdata2 = zeros(2,length(freqs));
        [tpsort sortidx] = sort(plotdata,1);% sort in trial dimension
        if perc < 1 % shave the outermost trials
            cutends = round(size(tpsort,1)-(size(tpsort,1)*perc));
            tpsort(1:cutends,:) = [];
            tpsort(end-cutends:end,:) = [];
        end;
        mdn = median(tpsort);
        if side == 1 % plot only positive deviations
            plotdata2 = [tpsort(end,:);mdn];
        elseif side == -1 %plot only negative deviations
            plotdata2 = [mdn;tpsort(1,:)];
        else % if side=[] or length(side) > 1
            plotdata2 = [tpsort(end,:);tpsort(1,:)];
        end;
  
        if ~isempty(rawbkgd)
            if strcmp(s.freqscale,'linear')
                f1 = [freqs,[freqs(end):-1:1]];
                f2 = [plotdata2(1,:),[plotdata2(2,end:-1:1)]];
                ph = fill(f1,f2,rawbkgd);hold on;
            elseif strcmp(s.freqscale,'quad')                
                lfreqs = sqrt(freqs); % find linear freqs
                f1 = [lfreqs,[lfreqs(end:-1:1)]];
                f2 = [plotdata2(1,:),[plotdata2(2,end:-1:1)]];
                ph = fill(f1,f2,rawbkgd);hold on;
                %%  --- quad-spaced freqs
                labelx = [4,10,20,40,70,100,150,200,250]; % xtick labels
                labelx = labelx(find(labelx > freqs(1) & labelx < freqs(end)));
                for f = 1:length(labelx)
                    [val mch(f)] = min(abs(freqs-labelx(f)));
                end;
                realx = lfreqs(mch);        
                set(gca,'xtick',realx); % appropriately label
                set(gca,'xticklabel',labelx);
            elseif strcmp(s.freqscale,'log')
                lfreqs = exp(freqs); % find linear freqs
                f1 = [lfreqs,[lfreqs(length(lfreqs):-1:1)]];
                f2 = [plotdata2(1,:),[plotdata2(2,end:-1:1)]];
                ph = fill(f1,f2,rawbkgd);hold on;
                %%  --- log-spaced freqs
                labelx = [4,10,20,40,70,100,150,200,250]; % xtick labels
                labelx = labelx(find(labelx > freqs(1) & labelx < freqs(end)));
                for f = 1:length(labelx)
                    [val mch(f)] = min(abs(freqs-labelx(f)));
                end;
                realx = lfreqs(mch);        
                set(gca,'xtick',realx); % appropriately label
                set(gca,'xticklabel',labelx);
            end;
            set(ph,'edgecolor',rawbkgd); %set(ph,'edgealpha',0);       
            minl = min(plotdata2(:));
            maxl = max(plotdata2(:));
            clear plotdata plotdata2 envdata
        else
            minl = 0; maxl = 0;
        end;
        %%%  Plot PCA-reduced data back-projection  %%%%%%%%%%%%%%%%%%%%     
        pcaproj = winv*activations(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);        
        pcaproj = pcaproj + repmat(s.meanpwr(rcp,:),[size(pcaproj,1) 1]);
        pcaproj=pcaproj(:,fr);
        pcaproj2 = zeros(2,length(freqs));
        [tpsort sortidx] = sort(pcaproj,1);% sort in pc dimension
        if perc < 1 % shave the outermost trials
            cutends = round(size(tpsort,1)-(size(tpsort,1)*perc));
            tpsort(1:cutends,:) = [];
            tpsort(end-cutends:end,:) = [];
        end;
        mdn = median(tpsort);
        if side == 1 % plot only positive deviations
            pcaproj2 = [tpsort(end,:);mdn];
        elseif side == -1 %plot only negative deviations
            pcaproj2 = [mdn;tpsort(1,:)];
        else % if side=[] or length(side) > 1
            pcaproj2 = [tpsort(end,:);tpsort(1,:)];
        end;
        if ~isempty(pcabkgd) % only plot if given a color
        if strcmp(s.freqscale,'linear')
            f1 = [freqs,[freqs(end):-1:1]];
            f2 = [pcaproj2(1,:),[pcaproj2(2,end:-1:1)]];
            ph = fill(f1,f2,pcabkgd);hold on;
        elseif strcmp(s.freqscale,'quad')
            lfreqs = sqrt(freqs); % find linear freqs
            f1 = [lfreqs,[lfreqs(end:-1:1)]];
            f2 = [pcaproj2(1,:),[pcaproj2(2,end:-1:1)]];
            ph = fill(f1,f2,pcabkgd);hold on;
            %%  --- quad-spaced freqs
            labelx = [4,10,20,40,70,100,150,200,250]; % xtick labels
            labelx = labelx(find(labelx > freqs(1) & labelx < freqs(end)));
            for f = 1:length(labelx)
                [val mch(f)] = min(abs(freqs-labelx(f)));
            end;
            realx = lfreqs(mch);        
            set(gca,'xtick',realx); % appropriately label
            set(gca,'xticklabel',labelx);
        elseif strcmp(s.freqscale,'log')
            lfreqs = exp(freqs); % find linear freqs
            f1 = [lfreqs,[lfreqs(length(lfreqs):-1:1)]];
            f2 = [pcaproj2(1,:),[pcaproj2(2,end:-1:1)]];
            ph = fill(f1,f2,pcabkgd);hold on;
            %%  --- log-spaced freqs
            labelx = [4,10,20,40,70,100,150,200,250]; % xtick labels
            labelx = labelx(find(labelx > freqs(1) & labelx < freqs(end)));
            for f = 1:length(labelx)
                [val mch(f)] = min(abs(freqs-labelx(f)));
            end;
            realx = lfreqs(mch);        
            set(gca,'xtick',realx); % appropriately label
            set(gca,'xticklabel',labelx);
        end;
        set(ph,'edgecolor',pcabkgd);%set(ph,'edgealpha',0);    
        xlabel('Frequency (Hz)'); ylabel('Power (dB)');

        minl2 = min(pcaproj2(:));
        maxl2 = max(pcaproj2(:)); 
        if minl2 < minl
            minl = minl2;
        end;
        if maxl2 > maxl
            maxl = maxl2;
        end;  
        clear pcaproj pcaproj2
        end; % to conditional plotting
        %%%%%%%%%%%%%%  Plot IM back-projections   %%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~isempty(trialset)
            imwinv = winv(trialset,:); % cut down to specified trials
        else
            imwinv = winv;
        end;    
        for tpp = 1:length(templates)            
            tp = templates(tpp);
            backproj = imwinv(:,tp)*activations(tp,:); % backproj curr IM

            onebkprj = backproj(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
            [tpsort sortidx] = sort(imwinv(:,tp)); % sort trials by curr IM
            onebkprj = onebkprj(sortidx,:);
            onebkprj =  onebkprj + repmat(s.meanpwr(rcp,:),[size(onebkprj,1) 1]);
            if perc < 1 % to shave off outer trials
                cutends = round(length(tpsort)-(length(tpsort)*perc));
                tpsort(1:cutends) = [];onebkprj(1:cutends,:) = [];
                tpsort(end-cutends:end) = [];onebkprj(end-cutends:end,:) = [];
            end;
            if side == 1 % only pos proj
                envdata2(1,:) = onebkprj(end,:); 
            elseif side < 0% only neg proj
                envdata2(1,:) = onebkprj(1,:); 
            else % side=[] or length(side) > 1
                envdata2 = [onebkprj(end,:);onebkprj(1,:)]; 
            end;
            if mnfac == 1 % to plot mean instead of env
                envdata2 = mean(onebkprj,1);
            end;
            if length(templates) < 2 & isempty(side)% do red/blue if only one
                snglcols = {'r','b'};
                for e = 1:size(envdata2,1) 
                    if strcmp(s.freqscale,'linear')
                        ph = plot(freqs,envdata2(e,fr),'k-','linewidth',lnwdth); hold on;
                    elseif strcmp(s.freqscale,'quad')
                        [ph realx labelx] = quadplot(freqs,envdata2(e,fr),lnwdth,'k'); hold on;
                    elseif strcmp(s.freqscale,'log')
                        [ph realx labelx] = logplot(freqs,envdata2(e,fr),lnwdth,'k'); hold on;
                    end;
                    set(ph,'color',snglcols{e});
                end;
            else                
                if strcmp(s.freqscale,'linear')
                    ph = plot(freqs,envdata2(:,fr),'k-','linewidth',lnwdth); hold on;
                elseif strcmp(s.freqscale,'quad')
                    [ph realx labelx] = quadplot(freqs,envdata2(:,fr),lnwdth,'k'); hold on;
                elseif strcmp(s.freqscale,'log')
                    [ph realx labelx] = logplot(freqs,envdata2(:,fr),lnwdth,'k'); hold on;
                end;
                set(ph,'color',cols(tpp,:));
            end;
            collmaxes(1,tpp) = max(envdata2(:));
            collmins(1,tpp) = min(envdata2(:));       
            
            clear plotprj envdata2
        end;
        
        
        
        minl2 = min(collmins); % for use without full back-proj
        maxl2 = max(collmaxes); 
        if minl2 < minl
            minl = minl2;
        end;
        if maxl2 > maxl
            maxl = maxl2;
        end;     
        maxl = maxl + maxl*.02;
        minl = minl - minl*.02;
        set(gca,'ylim',[minl maxl]); set(gca,'box','off');

        minl = 10*(floor(minl/10)); % get to the closest 10
        maxl = 10*(ceil(maxl/10)); % get to the closest 10
        
        set(gca,'ytick',[minl:5:maxl]);        set(gca,'yticklabel',[minl:5:maxl]);
        set(gca,'ticklength',[.02 .02]);
        %title(['Cp ',int2str(comps(cp))]);
        
        % plot mean power:----
        if strcmp(s.freqscale,'linear')
            ph = plot(freqs,s.meanpwr(rcp,fr),'k-','linewidth',lnwdth+.1); hold on;
        elseif strcmp(s.freqscale,'quad')
            [ph realx labelx] = quadplot(freqs,s.meanpwr(rcp,fr),lnwdth+.1,'k'); 
        elseif strcmp(s.freqscale,'log')
            [ph realx labelx] = logplot(freqs,s.meanpwr(rcp,fr),lnwdth+.1,'k'); hold on;
        end;      
        fprintf('\nComp %s of %s done.',int2str(cp),int2str(length(comps)));
        %set(gca,'xscale','log'); %%%% ******  added this for log freq axis
        pl = pl+2;% advance for scalp map and spectra
    end;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    set(gcf,'color','w');
    axcopy
    freqout = s.freqs(fr);
