% plots results from SpecCoModAnal.m
%
% SpecCoModPlot(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone);
%
% INPUTS:
% datset: dataset with EEG.icawinv to plot scalp maps
% datpath: full directory path where dataset can be found
% plotcomps: list of components to plot (must be in original 'complist')(default: all)
% factors: list of factor numbers to plot (default: all)
% savedat: filename that float files were saved as (same as input to SpecCoModAnal.m)
% frqlim: [minfreq maxfreq] to plot
% superpos -- ['y' or 'n'] 'y' will plot a single axis with all comps superimposed.
%             (can only be used for plotting all comps (justone needs to be []) 
% meanvar -- [0 | 1] if 1, then will plot mean spectra and trial backprojections, else will plot templates
% justone: (integer) factor number to plot alone. [] to plot all. Plots mean with spectral variations on top
% auxpath -- [string] if .mat is saved in a different folder from .set
%             this string should be the one pointing to the .mat
%

function SpecCoModTWPlot(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone,auxpath);
    

    scalelog = []; % if 'log' will scale x-axis in log
    lnwdth = 2.5;
    nlim = []; mlim = [];
    EEG = pop_loadset(datset, datpath);
    fprintf('\nLoading subject spectral decomposition data...\n');
    if exist('auxpath')
        datpath = auxpath;
    end;
    s = load([datpath,savedat,'.mat']);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph=floatread([datpath,savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    icamatall = floatread([datpath,savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws icamatall
    %%% for TW
    wv = activations'; 
    activations = winv'; winv = wv; clear wv
    %%%%%%%%%
    if isempty(plotcomps)
        plotcomps = s.complist;
    end;
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    if isempty(justone)
        minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
        maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
        if superpos == 'y'
            figure;row = round(sqrt(length(factors))); col = ceil(sqrt(length(factors)));
            cols = jet(length(plotcomps)); pl = 1;
            for tpp = 1:length(factors)
                lnfac = (lnwdth-2)/length(plotcomps);
                tp = factors(tpp);               
                sbplot(row,col,pl)
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == s.complist); 
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(cp,:)); hold on;
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); hold on;
                        set(gca,'xscale','log');
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
                end;pl = pl+1;
                set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                set(gca,'ylim',[minl maxl]); title(['Fac ',int2str(factors(tpp))]);
                set(gca,'box','off');
                set(gca,'xgrid','on');
                if tpp == length(factors)
                    xlabel('Frequency (Hz)');
                    ylabel('Relative Power');
                end;                
            end; textsc(['Superimposed Factor Spectra: ',datpath],'title');
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            set(gcf,'color','w');
            axcopy
        else  
            figure;row = length(factors)+1; 
            if row > 16
                row = round(row/2);
                if row > 16
                    row = 16;
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
                % this is for TW decomp to make uniform only within IM
                minl = min(min(activations(tp,:)))-abs(min(min(activations(tp,:))))*.01;
                maxl = max(max(activations(tp,:)))+abs(max(max(activations(tp,:))))*.01;
                if pl == row*col+1
                    set(gcf,'Position',[100 300 1400 900]);
                    set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
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
                set(gca,'fontsize',7);
                plot([0 0],[get(gca,'ylim')],'r-');
                set(gca,'yticklabel',[]);   set(gca,'xticklabel',[]); 
                title(['IM ',int2str(tp)]);
                
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == s.complist);
                    sbplot(row,col,pl)
                    if strcmp(s.freqscale,'quad') % quadratic spacing
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        quadplot(s.freqs(fr),tpact(:,fr),1.75,'b');
                        pl = pl+1;hold on;
                    elseif strcmp(s.freqscale,'log')  % log spacing
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        logplot(s.freqs(fr),tpact(:,fr),1.75,'b'); pl = pl+1;hold on;
                    else % otherwise linear
                        plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',1.75); pl = pl+1;hold on;
                        set(gca,'xtick',[10:10:frqlim(2)]);
                        set(gca,'xticklabel',{10 [] 30 [] []});
                        set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                    end;                    
                    %set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'ylim',[minl maxl]); 
                    set(gca,'xgrid','on');
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'ticklength',[.03 .03]);
                    plot([get(gca,'xlim')],[0 0],'r-');
                    if strcmp(scalelog,'log')
                        set(gca,'xscale','log');%%%%%%%%% **** added to put freqs on log scale
                    end;
                    if pl <= (row-1)*col+1
                        set(gca,'xticklabel',[]);
                        set(gca,'yticklabel',[]);
                    end;                    
                end;
            end;
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
        end;
    else
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
        if meanvar == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %z = load('/data/common4/emotion/KmeansClustCoMods.mat' );
            %s = load([datpath,'SpecCoModStuff.mat']);
            cols = jet(length(plotcomps));
            x=[1:size(activations,1)];x(justone) = [];                
            acts = activations; acts(x,:) = 0;
            backproj = winv*acts ; clear plotprj
            for cmp = 1:length(plotcomps)
                rcp =  find(plotcomps(cmp) == s.complist);
                sbplot(row,col,pl)
                onebkprj = backproj(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                for trls = 1:size(onebkprj,1)
                    plotprj(trls,:) = onebkprj(trls,:) + s.meanpwr(rcp,:);
                end;
                if isempty(nlim)
                    nlim = min(plotprj(:))+(min(plotprj(:))*.1); 
                    mlim = max(plotprj(:))+(max(plotprj(:))*.1); 
                end;
                ph = plot(freqs,plotprj,'k-'); hold on;
                set(ph,'color',cols(cmp,:));
                ph = plot(s.freqs,s.meanpwr(rcp,:),'k-'); 
                set(ph,'linewidth',2);set(gca,'xlim',[3 50]);
                set(gca,'ylim',[nlim mlim]);
                    if logyes == 1
                        set(gca,'xscale','log');set(gca,'xtick',[4 6 10 20 40]);
                        set(gca,'xtick',[5,10,30:20:frqlim(2)]);  
                        set(gca,'xticklabel',{5 10 30 [] []}); 
                    else
                        set(gca,'xtick',[10:10:frqlim(2)]);
                        set(gca,'xticklabel',{10 [] 30 [] []});
                    end;
                fprintf('.');pl = pl+2; 
                if cmp ~= length(plotcomps)
                    set(gca,'xticklabel',[]);
                end;            
                title(['Fc ',int2str(justone),';Cp ',int2str(plotcomps(cmp))]);
            end;
            if pl < (row * col)+1
                sbplot(row,col,pl)
                hist(winv(:,justone),75);hold on; pl = pl+1;
                set(gca,'fontsize',10);%set(gca,'xlim',[-1.5 1.5]);
                plot([0 0],[get(gca,'ylim')],'r-'); 
                %set(gca,'yticklabel',[]); set(gca,'xticklabel',[]);
                title('Factor Weights');
            end;
            if pl < (row * col)+1
                cols = jet(length(plotcomps));      
                sbplot(row,col,pl)
                for cp = 1:length(plotcomps)
                    rcp =  find(plotcomps(cp) == s.complist);
                    ph = plot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                    hold on; set(ph,'color',cols(cp,:));
                    set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                end;title(['Fac ',int2str(justone)]); 
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);
                set(gca,'ylim',[min(activations(justone,:))-.1 max(activations(justone,:))+.1]);
            end;
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates Factor ',int2str(justone),' for ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        else            
            tp = justone;
            ymin = min(min(activations(tp,:)))-.05*min(min(activations(tp,:)));
            ymax = max(max(activations(tp,:)))+.05*max(max(activations(tp,:)));
            for cp = 1:length(plotcomps)
                rcp =  find(plotcomps(cp) == s.complist);
                sbplot(row,col,pl)
                if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                    ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth); hold on;
                elseif strcmp(s.freqscale,'log') % log spaced
                    ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); hold on;
                else % otherwise linear
                    ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                    hold on;  
                end;
                pl = pl+2;hold on;
                plot([get(gca,'xlim')],[0 0],'k-'); 
                set(gca,'ylim',[ymin ymax]);   
                set(gca,'xgrid','on');
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
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(rcp,:)); hold on;
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); hold on;
                        set(ph,'color',cols(rcp,:));
                    else % otherwise linear
                        ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                        hold on;  
                        set(ph,'color',cols(rcp,:));
                    end;                   
                end;title(['IM ',int2str(tp)]); 
                set(gca,'ylim',[ymin ymax]); 
            end;
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates Factor ',int2str(justone),' for ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end;
    end;
    
