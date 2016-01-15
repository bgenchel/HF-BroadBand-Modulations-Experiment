% plots results from DecompMat.m on cross-correlation spectral data
%
% PlotCoModDecomp(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone);
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
% this differs from SpecCoModPlot only in that it has PCA matrix saved separately

function [winv activations] = PlotXcorrDecomp(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,auxpath);
    
    pageorient = 'landscape';

    scalelog = ''; % if 'log' will scale x-axis in log
    lnwdth = 2.5;
    nlim = []; mlim = [];
    fprintf('\nLoading subject spectral decomposition data...\n');
    if exist('auxpath')
        if ~isempty(auxpath)
            datpath = auxpath;
        end;
    end;
    s = load([datpath,savedat,'.mat']);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph=floatread([datpath,savedat,'.sph'],[s.pcmat s.pcmat],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcmat s.pcmat],[],0); 
    icamatall = floatread([datpath,savedat,'.fdt'],[s.pcmat s.numframes],[],0);    
    ws = wts*sph;      winv = pinv(ws); activations = ws*icamatall;  % templates 
    speceig = floatread([datpath,savedat,'MATEIG.fdt'],[length(s.rowmeans) s.pcmat],[],0);
    specwts = speceig*winv;   
    winv = specwts;    
    
    clear wts sph ws icamatall
    
    if isempty(plotcomps)
        plotcomps = s.complist;
    end;
    
    if size(plotcomps,2) > 21
        xxx=1;
        for xx = 0:20:size(plotcomps,2)
            if xx+20 > size(plotcomps,2)
                splitplotcomps{xxx} = plotcomps(:,[xx+1:end]);xxx=xxx+1;
            else
                splitplotcomps{xxx} = plotcomps(:,[xx+1:xx+20]);xxx=xxx+1;
            end;
        end;
    else
        splitplotcomps{1} = plotcomps;
    end;
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    if length(factors) == 1 % if only one IM requested
        justone = factors;  
    else
        justone = [];
    end;
    
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
    maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
    
    if isempty(justone) % plot multiple IMs
        for sp = 1:length(splitplotcomps)
            plotcomps = splitplotcomps{sp};
            if superpos == 'y'
                figure;row = round(sqrt(length(factors))); col = ceil(sqrt(length(factors)));
                cols = jet(size(plotcomps,2)); pl = 1;
                for tpp = 1:length(factors)
                    lnfac = (lnwdth-2)/length(plotcomps);
                    tp = factors(tpp);               
                    sbplot(row,col,pl)
                    for cp = 1:size(plotcomps,2)
                        rcp = intersect(find(plotcomps(1,cp) == s.complist(1,:)), find(plotcomps(2,cp) == s.complist(2,:)));
                        if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                            ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(cp,:)); hold on;
                        elseif strcmp(s.freqscale,'log') % log spaced
                            ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),1.75,'b'); hold on;
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
                    pl = pl+1;
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
                set(gcf,'PaperOrientation',pageorient);  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                set(gcf,'color','w');
                axcopy
            else  
                figure;row = length(factors); 
                if row > 16
                    row = round(row/2);
                    if row > 16
                        row = 16;
                    end;             
                end;            
                col = size(plotcomps,2)+1;
                pl = 1;           
                for tpp = 1:length(factors)
                    tp = factors(tpp);
                    if pl == row*col+1
                        set(gcf,'Position',[100 300 1400 900]);
                        set(gcf,'PaperOrientation',pageorient);  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                        axcopy
                        ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                        set(ph,'fontsize',14);
                        figure; pl = 1;
                    end;
                    sbplot(row,col,pl)
                    hist(winv(:,tp),75);pl = pl+1;hold on;
                    set(gca,'fontsize',7);%set(gca,'xlim',[-2 2]);
                    plot([0 0],[get(gca,'ylim')],'r-');
                    set(gca,'yticklabel',[]);   set(gca,'xticklabel',[]); 
                    title(['IM ',int2str(tp)]);
                    for cp = 1:size(plotcomps,2)
                        rcp = intersect(find(plotcomps(1,cp) == s.complist(1,:)), find(plotcomps(2,cp) == s.complist(2,:)));
                        sbplot(row,col,pl)
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        if strcmp(s.freqscale,'quad') % quadratic spacing
                            ph = quadplot(s.freqs(fr),tpact(:,fr),1.75,'b'); hold on;                   
                        elseif strcmp(s.freqscale,'log')  % log spacing
                            ph = logplot(s.freqs(fr),tpact(:,fr),1.75,'b');hold on;
                        else % otherwise linear
                            plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',1.75); hold on;
                            set(gca,'xtick',[10:10:frqlim(2)]);
                            set(gca,'xticklabel',{10 [] 30 [] []});
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
                        if pl <= col
                            title([int2str(plotcomps(1,cp)),'-',int2str(plotcomps(2,cp))]);
                        end;
                        if pl <= (row-1)*col+1 & tpp ~= length(factors)
                            set(gca,'xticklabel',[]);
                        end;                    
                        if pl <= (row-1)*col+2 | pl > (row-1)*col+3 & tpp ~= length(factors)
                            set(gca,'yticklabel',[]);
                        end;                    
                        pl = pl+1;
                    end;
                end;
                set(gcf,'Position',[100 300 1400 900]);
                set(gcf,'PaperOrientation',pageorient);  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                set(ph,'fontsize',14);
            end;
        end;
    else % plot only one IM with scalp maps.
        
        EEG = pop_loadset(datset, datpath); % for scalp maps
        npairs = size(s.complist,2); % number of IC pairs
        sbtot = 3*npairs; 
        eql = round(sqrt(sbtot));
        if eql < 5
            col = 3; row = ceil(sbtot/col);
        elseif eql >= 5 & eql < 8
            col = 6; row = ceil(sbtot/col);
        elseif eql >= 8 & eql < 11
            col = 9; row = ceil(sbtot/col);
        elseif eql >= 11 
            col = 12; row = ceil(sbtot/col);
        end;
        figure; pl = 1;
        for p = 1:npairs
            sbplot(row,col,pl); pl = pl+1; % first IC             
            topoplot(EEG.icawinv(:,s.complist(1,p)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); 
            set(gca,'fontsize',9);  title(int2str(s.complist(1,p)));
            
            
            sbplot(row,col,pl); pl = pl+1; % associated cross spectrum
            rcp = intersect(find(s.complist(1,p) == s.complist(1,:)), find(s.complist(2,p) == s.complist(2,:)));                    
            tpact = activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
            if strcmp(s.freqscale,'quad') % quadratic spacing
                ph = quadplot(s.freqs(fr),tpact(:,fr),1.75,'b'); hold on;                   
            elseif strcmp(s.freqscale,'log')  % log spacing
                ph = logplot(s.freqs(fr),tpact(:,fr),1.75,'b');hold on;
            else % otherwise linear
                plot(s.freqs,tpact(:,fr),'linewidth',1.75); hold on;
                set(gca,'xtick',[10:10:frqlim(2)]);
                set(gca,'xticklabel',{10 [] 30 [] []});
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
            end;                    
            set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
            set(gca,'ylim',[minl maxl]); 
            set(gca,'xgrid','on');
            set(gca,'fontsize',7);set(gca,'box','off');
            set(gca,'ticklength',[.03 .03]);
            plot([get(gca,'xlim')],[0 0],'r-');
            
            
            sbplot(row,col,pl); pl = pl+1; % second IC             
            topoplot(EEG.icawinv(:,s.complist(2,p)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); 
            set(gca,'fontsize',9);  title(int2str(s.complist(2,p)));
        end;             
        set(gcf,'Position',[100 300 1400 900]);
        set(gcf,'PaperOrientation',pageorient);  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        ph=textsc(['Cross Spectral Template ',int2str(justone),'; all IC pairs from ',datpath(end-4:end-1),': ',savedat,'.mat'],'title');
        set(ph,'fontsize',14);
        
    end;
