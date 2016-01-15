% plots results from SpecCoModAnal.m
%
% SpecCoModPlot(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone,largeone,auxpath);
%
% INPUTS:
% datset: dataset with EEG.icawinv to plot scalp maps
% datpath: full directory path where dataset can be found
% plotcomps: list of components to plot (must be in original 'complist')(default: all)
% factors: list of factor numbers to plot (default: all)
% savedat: filename that float files were saved as (same as input to SpecCoModAnal.m)
% frqlim: [minfreq maxfreq] to plot
% superpos -- ['ims', 'ics' or 'off'] 'ics' will plot a single axis with all comps 
%             superimposed.Not available when 'justone' is used.
%             'ims' will show each IC scalp map with super imposed IM templates.
%             'off' will plot all ICs and IMs plotted separately. 
% meanvar -- [0 | 1] if 1, then will plot mean spectra and trial backprojections, else will plot templates (only valid for 'justone' plots), otherwise plots templates
% justone: (integer) factor number to plot alone. [] to plot all. 
%          Plots mean with spectral variations on top
% largeone -- {'on' or 'off'} plots just the largest IC template or back-projection
%              (if meanvar == 1) for each IM with associated scalp map.
%             If 'largeone' == 'on', this trumps all other options
% auxpath -- [string] if .mat is saved in a different folder from .set
%             this string should be the one pointing to the .mat
%
% Plotting options:
% 1) regular view with scalp maps on top row, corresponding IM templates below 
%    can plot any number of IMs up to 'mxims' per page.
% 2) plot one axis of ONE IM with templates superimposed (no scalp maps)('superpos')
% 3) plot only one IM with scalp maps, but templates (or back-projections, 
%    if 'meanvar'==1) arranged to fill a page. 
% 4) plot only the largest template for each requested IM with assoc. scalp map


function SpecCoModPlot(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone,largeone,auxpath);
    
    mxims = 15; % number of factors to plot each page (15 usually)
    scalelog = 'off'; % if 'log' will scale x-axis in log
    lnwdth = 2.5;
    nlim = []; mlim = [];
    EEG = pop_loadset(datset, datpath);
    fprintf('\nLoading subject %s data...\n',savedat);
    if exist('auxpath')
        if ~isempty(auxpath)
            datpath = auxpath;
        end;
    end;
    if ~exist('largeone')
        largeone = 'off';
    end;
    s = load([datpath,savedat,'.mat']);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph=floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);        
    icamatall = floatread([datpath,savedat,'.fdt'],[s.pcs s.numframes],[],0);    
    ws = wts*sph;    
    winv = pinv(ws);
    activations = ws*icamatall;    clear wts sph ws icamatall
    speceig = floatread([datpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  
    winv = specwts;
    
    if isempty(plotcomps)
        plotcomps = s.complist;
    end;
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    if isempty(frqlim)
      frqlim = [s.freqs(1) s.freqs(end)];
    end;
    fr = find(s.freqs >= frqlim(1) & s.freqs <= frqlim(end));
    if strcmp(largeone, 'on')% any number of IMs, can be templates or meanvar==1
        figure;
        row = round(sqrt(length(factors)*2));
        col = ceil(sqrt(length(factors)*2));
        if ~iseven(col)
            col = col+1;
        end;
        if row > 6 % don't plot too many per page
            row = 6;col = 6;
        end;
        pl= 1;
        minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
        maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
        for im = 1:length(factors)            
            for cp = 1:length(s.complist) % find largest variance IC
                tpact = activations(factors(im),length(s.freqs)*(cp-1)+1:length(s.freqs)*cp);
                mxval(1,cp) = max(abs(tpact(:)));
            end;
            [val plotic] = max(mxval); % index of largest IC
            if pl > row * col
                figure; pl = 1;
            end;
            % plot scalp map first:
            sbplot(row,col,pl);pl = pl+1;
            topoplot(EEG.icawinv(:,s.complist(plotic)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); 
            title(['IC ',int2str(s.complist(plotic))]);
            if meanvar==1% plot back projs
                lnwdth = 1; bkcol = [.6 0 1];
                sbplot(row,col,pl);pl = pl+1;
                backproj = winv(:,factors(im))*activations(factors(im),:) ; clear plotprj
                onebkprj = backproj(:,length(s.freqs)*(plotic-1)+1:length(s.freqs)*plotic); 
                plotprj = onebkprj + repmat(s.meanpwr(plotic,:),[size(onebkprj,1) 1]);
                if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                    ph = quadplot(s.freqs(fr),plotprj(:,fr),lnwdth,bkcol); hold on;
                    %set(ph,'color',cols(cmp,:));
                    ph = quadplot(s.freqs(fr),s.meanpwr(plotic,fr),lnwdth*2.5,'k');
                elseif strcmp(s.freqscale,'log') % log spaced
                    ph = logplot(s.freqs(fr),plotprj(:,fr),lnwdth,bkcol); hold on;
                    %set(ph,'color',cols(cmp,:));
                    ph = logplot(s.freqs(fr),s.meanpwr(plotic,fr),lnwdth*2.5,'k'); 
                else % otherwise linear
                    ph = plot(s.freqs(fr),plotprj(:,fr),'k-','linewidth',lnwdth,'color',bkcol); hold on;
                    %set(ph,'color',cols(cmp,:));
                    ph = plot(s.freqs(fr),s.meanpwr(plotic,fr),'k-','linewidth',lnwdth*2); 
                end;           
                set(gca,'xgrid','on');
                title(['IM ',int2str(factors(im))]);
                
            else % just plot template
                tpact = activations(factors(im),length(s.freqs)*(plotic-1)+1:length(s.freqs)*plotic);
                sbplot(row,col,pl);pl = pl+1;
                if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                    ph = quadplot(s.freqs(fr),tpact(:,fr),lnwdth); hold on;
                elseif strcmp(s.freqscale,'log') % log spaced
                    ph = logplot(s.freqs(fr),tpact(:,fr),lnwdth); hold on;
                else % otherwise linear
                    ph = logplot(s.freqs(fr),tpact(:,fr),lnwdth,cols(rcp,:)); 
                    hold on;  
                end;
                set(gca,'ylim',[minl maxl]); set(gca,'xgrid','on');  
                set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                plot([get(gca,'xlim')],[0 0],'r-');
                title(['IM ',int2str(factors(im))]);
            end;
        end;
    else% 'largeone','off'------------------------------------------------------
        if isempty(justone)% plot more than one IM
            minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
            maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
            if strcmp(superpos,'ics') % superimpose ic templates
                figure;row = round(sqrt(length(factors))); col = ceil(sqrt(length(factors)));
                cols = lines(length(plotcomps)); pl = 1;
                for tpp = 1:length(factors)
                    lnfac = (lnwdth-2)/length(plotcomps);
                    tp = factors(tpp);               
                    sbplot(row,col,pl)
                    for cp = 1:length(plotcomps)
                        rcp = find(plotcomps(cp) == s.complist); 
                        if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                            ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(cp,:)); hold on;
                        elseif strcmp(s.freqscale,'log') % log spaced
                            ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'b'); hold on;
                            set(ph,'color',cols(cp,:)); 
                        else % otherwise linear
                            ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                            hold on; 
                            set(ph,'color',cols(cp,:)); 
                        end;
                        lnwdth = lnwdth - lnfac;
                        if strcmp(scalelog,'log')
                            set(gca,'xscale','log');%%%%%%%%% **** added to put freqs on log scale
                        end;
                    end;
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'ylim',[minl maxl]); title(['IM ',int2str(factors(tpp))]);
                    set(gca,'box','off');
                    set(gca,'xgrid','on');
                    if pl == (row-1)*col+1
                        xlabel('Frequency (Hz)'); ylabel('Relative Power');
                    elseif pl > (row-1)*col+1
                        xlabel('Frequency (Hz)');
                    end;  
                    if pl <= col*(row-1)
                        set(gca,'xticklabel',[]);
                    end;pl = pl+1;
                end; textsc(['Superimposed IC Templates for single IMs: ',datpath],'title');
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                set(gcf,'color','w');
                axcopy
            elseif strcmp(superpos,'ims')  % superimpose IM templates for each IC
                figure;
                row = round(sqrt(length(plotcomps)*2)); 
                col = ceil(sqrt(length(plotcomps)*2));
                if ~iseven(col)
                    col = col+1; row = row-1;
                end;
                lnwdth = 2;
                cols = jet(length(factors)); 
                pl = 1;
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == s.complist); 
                    sbplot(row,col,pl); pl = pl+1;
                    topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                    set(gca,'fontsize',12);  title(['IC ',int2str(plotcomps(cp))]);
                    sbplot(row,col,pl); 
                    for tpp = 1:length(factors)
                        lnfac = (lnwdth-2)/length(plotcomps);
                        tp = factors(tpp);               
                        if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                            ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(tpp,:)); hold on;
                        elseif strcmp(s.freqscale,'log') % log spaced
                            ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'b'); hold on;
                            set(ph,'color',cols(tp,:)); 
                        else % otherwise linear
                            ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth);  hold on; 
                            set(ph,'color',cols(tpp,:)); 
                        end;
                        %lnwdth = lnwdth - lnfac;
                        if strcmp(scalelog,'log')
                            set(gca,'xscale','log');%%%%%%%%% **** added to put freqs on log scale
                        end;
                    end;
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'ylim',[minl maxl]); title(['IM templates']);
                    set(gca,'xgrid','on');
                    if pl == (row-1)*col+2
                        xlabel('Frequency (Hz)'); ylabel('Relative Power');
                    elseif pl > (row-1)*col+1
                        xlabel('Frequency (Hz)');
                    end;  
                    if pl <= col*(row-1)
                        set(gca,'xticklabel',[]);
                    end;pl = pl+1;
                end; 
                textsc(['Superimposed IM templates for single ICs: ',datpath(end-6:end),'-',savedat],'title');
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                set(gcf,'color','w');
                axcopy
               
            else  % do not super impose IMs or ICs
                figure;row = length(factors)+1; 
                if row > mxims + 1
                    row = round(row/2);
                    if row > mxims + 1
                        row = mxims + 1;
                    end;             
                end;            
                col = length(plotcomps)+1;
                pl = 2;           
                for cp = 1:length(plotcomps)
                    sbplot(row,col,pl)
                    topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
                    set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
                end;
                for tpp = 1:length(factors)
                    tp = factors(tpp);
                    if pl == row*col+1
                        set(gcf,'Position',[100 300 1400 900]);
                        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                        axcopy
                        ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                        set(ph,'fontsize',14);
                        figure; pl = 2;
                        for cp = 1:length(plotcomps)
                            sbplot(row,col,pl)
                            topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
                            set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
                        end;
                    end;
                    sbplot(row,col,pl)
                    hist(winv(:,tp),75);pl = pl+1;hold on;
                    set(gca,'fontsize',7);%set(gca,'xlim',[-2 6]);
                    plot([0 0],[get(gca,'ylim')],'r-');
                    set(gca,'yticklabel',[]);   %set(gca,'xticklabel',[]); 
                    title(['IM ',int2str(tp)]);
                    for cp = 1:length(plotcomps)
                        rcp = find(plotcomps(cp) == s.complist);
                        sbplot(row,col,pl)
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        if strcmp(s.freqscale,'quad') % quadratic spacing
                            quadplot(s.freqs(fr),tpact(:,fr),lnwdth,'b');
                            pl = pl+1;hold on;
                        elseif strcmp(s.freqscale,'log')  % log spacing
                            logplot(s.freqs(fr),tpact(:,fr),lnwdth,'b'); pl = pl+1;hold on;
                        else % otherwise linear
                            plot(s.freqs(fr),tpact(:,fr),'linewidth',lnwdth); pl = pl+1;hold on;
                            set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                        end;                    
                        set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                        set(gca,'ylim',[minl maxl]); 
                        set(gca,'xgrid','on');
                        set(gca,'fontsize',7);set(gca,'box','off');
                        set(gca,'ticklength',[.03 .03]);
                        plot([get(gca,'xlim')],[0 0],'r-');
                        if strcmp(scalelog,'log')
                            set(gca,'xscale','log');%%%%%%%%% **** added to put freqs on log scale
                        end;
                        if pl <= (row-1)*col+1 & tpp ~= length(factors)
                            set(gca,'xticklabel',[]);
                        end;                    
                        if pl <= (row-1)*col+2 | pl > (row-1)*col+3 & tpp ~= length(factors)
                            set(gca,'yticklabel',[]);
                        end;                    
                    end;
                end;
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                set(ph,'fontsize',14);
                axcopy
            end;
        else % if plotting 'justone' (in grid form instead of one row)
            cols = jet(length(plotcomps));
            minl = min(min(activations(justone,:)))-abs(min(min(activations(justone,:))))*.01;
            maxl = max(max(activations(justone,:)))+abs(max(max(activations(justone,:))))*.01;                
            figure;row = round(sqrt(length(plotcomps)*2))+1; col = round(sqrt(length(plotcomps)*2))+1;pl = 1;
            if ~iseven(col)
                col = col+1; row = row-1;
            end;       
            for cp = 1:length(plotcomps)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+2;
                set(gca,'fontsize',10);
                title(int2str(plotcomps(cp)));
            end;
            pl = 2;
            if meanvar == 1 % plot back projs
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %z = load('/data/common4/emotion/KmeansClustCoMods.mat' );
                            %s = load([datpath,'SpecCoModStuff.mat']);
                cols = jet(length(plotcomps));
                backproj = winv(:,justone)*activations(justone,:) ; clear plotprj
                for cmp = 1:length(plotcomps)
                    rcp =  find(plotcomps(cmp) == s.complist);
                    sbplot(row,col,pl)
                    onebkprj = backproj(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                    plotprj = onebkprj + repmat(s.meanpwr(rcp,:),[size(onebkprj,1) 1]);
                    if isempty(nlim)
                        nlim = min(plotprj(:))-(min(plotprj(:))*.1); 
                        mlim = max(plotprj(:))+(max(plotprj(:))*.1); 
                    end;
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs(fr),plotprj(:,fr),lnwdth); hold on;
                        set(ph,'color',cols(cmp,:));
                        ph = quadplot(s.freqs(fr),s.meanpwr(rcp,fr));set(ph,'color','k'); 
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = logplot(s.freqs(fr),plotprj(:,fr),lnwdth); hold on;
                        set(ph,'color',cols(cmp,:));
                        ph = logplot(s.freqs(fr),s.meanpwr(rcp,fr)); set(ph,'color','k'); 
                    else % otherwise linear
                        ph = plot(s.freqs(fr),plotprj(:,fr),'k-','linewidth',lnwdth); hold on;
                        set(ph,'color',cols(cmp,:));
                        ph = plot(s.freqs(fr),s.meanpwr(rcp,fr),'k-'); 
                    end;
                    %ph = plot(s.freqs,plotprj,'k-'); hold on;
                    %set(ph,'color',cols(cmp,:));
                    %ph = plot(s.freqs,s.meanpwr(rcp,:),'k-'); 
                    %set(gca,'ylim',[nlim mlim]);
                    fprintf('.');pl = pl+2; 
                    if cmp ~= length(plotcomps)
                        set(gca,'xticklabel',[]);
                    end;            
                    title(['IC ',int2str(plotcomps(cmp))]);
                end;
                if pl < (row * col)+1
                    sbplot(row,col,pl)
                    hist(winv(:,justone),75);hold on; pl = pl+1;
                    set(gca,'fontsize',10);%set(gca,'xlim',[-1.5 1.5]);
                    plot([0 0],[get(gca,'ylim')],'r-'); 
                    %set(gca,'yticklabel',[]); set(gca,'xticklabel',[]);
                    title('IM Weights');
                end;
                if pl < (row * col)+1
                    cols = jet(length(plotcomps));      
                    sbplot(row,col,pl)
                    for cp = 1:length(plotcomps)
                        rcp =  find(plotcomps(cp) == s.complist);
                        ph = plot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                        hold on; set(ph,'color',cols(cp,:));set(ph,'linewidth',lnwdth);
                        set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                    end;title(['IM ',int2str(justone)]); 
                    set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);
                    set(gca,'ylim',[min(activations(justone,:))-.1 max(activations(justone,:))+.1]);
                end;
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                ph=textsc(['Subj ',datpath,': ',savedat,'; IM ',int2str(justone)],'title');
                set(ph,'fontsize',14);
                axcopy
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            else  % plot templates (meanvar == 0) instead of back proj (meanvar==1)
                tp = justone;
                ymin = min(min(activations(tp,:)))-.05*min(min(activations(tp,:)));
                ymax = max(max(activations(tp,:)))+.05*max(max(activations(tp,:)));
                for cp = 1:length(plotcomps)
                    rcp =  find(plotcomps(cp) == s.complist);
                    sbplot(row,col,pl)
                    tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs(fr),tpact(:,fr),lnwdth); hold on;
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = logplot(s.freqs(fr),tpact(:,fr),lnwdth); hold on;
                    else % otherwise linear
                        ph = logplot(s.freqs(fr),tpact(:,fr),lnwdth,cols(rcp,:)); 
                        hold on;  
                    end;
                    set(gca,'ylim',[ymin ymax]);   
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'xgrid','on');
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'ticklength',[.03 .03]);
                    plot([get(gca,'xlim')],[0 0],'r-');
                    pl = pl+2;hold on;
                end;
                if pl < (row * col)+1 % plot IM weights
                    sbplot(row,col,pl)
                    hist(winv(:,tp),75);hold on; pl = pl+2;
                    plot([0 0],[get(gca,'ylim')],'r-'); 
                    title('IM Weights');
                end;
                if pl < (row * col)+ 2% plot all superimposed
                    cols = jet(length(plotcomps));      
                    sbplot(row,col,[pl pl+1])
                    for cp = 1:length(plotcomps)
                        rcp =  find(plotcomps(cp) == s.complist);
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                            ph = quadplot(s.freqs(fr),tpact(:,fr),lnwdth,cols(rcp,:)); hold on;
                        elseif strcmp(s.freqscale,'log') % log spaced
                            ph = logplot(s.freqs(fr),tpact(:,fr),lnwdth,cols(rcp,:)); hold on;
                        else % otherwise linear
                            ph = plot(s.freqs(fr),tpact(:,fr),'linewidth',lnwdth); 
                            hold on;  
                            set(ph,'color',cols(rcp,:));
                        end;                   
                    end;
                    title(['IM ',int2str(tp)]); 
                    set(gca,'ylim',[ymin ymax]);   
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'xgrid','on');
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'ticklength',[.03 .03]);
                    plot([get(gca,'xlim')],[0 0],'r-');
                end;
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                ph=textsc(['Spectral Templates Factor ',int2str(justone),' for ',datpath,': ',savedat],'title');
                set(ph,'fontsize',14);
                axcopy
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end;
        end;        
    end;
textsc([datpath(end-5:end),savedat],'title');
