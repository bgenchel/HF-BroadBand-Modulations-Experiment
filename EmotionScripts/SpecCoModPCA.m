% plots the PCA templates that are used for subsequent
% decomposition in the WT fashion
%
%
%
%
%
%

function SpecCoModPCA(datset,datpath,savedat,plotcomps,factors);
    
    
    
    EEG = pop_loadset(datset, datpath);
    s = load([datpath,savedat,'.mat']);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = floatread([datpath,savedat,'.fdt'],[s.pcs s.numframes],[],0);    
   
    fr = find(s.freqs);
    if isempty(plotcomps)
        plotcomps = s.complist;
    end;
    if isempty(factors)
        factors = [1:size(data,1)];
    end;

    minl = min(min(data(factors,:)))-abs(min(min(data(factors,:))))*.01;
    maxl = max(max(data(factors,:)))+abs(max(max(data(factors,:))))*.01;
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
        sbplot(row,col,pl); ph = text(-.4,.5,['PC ',int2str(tp)]); 
        set(ph,'fontsize',14); axis('off');pl = pl+1;
        for cp = 1:length(plotcomps)
            rcp = find(plotcomps(cp) == s.complist);
            sbplot(row,col,pl)
            tpact = data(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
            if strcmp(s.freqscale,'quad') % quadratic spacing
                quadplot(s.freqs(fr),tpact(:,fr),1.75,'b');
                pl = pl+1;hold on;
            elseif strcmp(s.freqscale,'log')  % log spacing
                logplot(s.freqs(fr),tpact(:,fr),1.75,'b'); pl = pl+1;hold on;
            else % otherwise linear
                plot(s.freqs(fr),tpact(:,fr),'linewidth',1.75); pl = pl+1;hold on;
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
            end;                    
            %set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
            set(gca,'ylim',[minl maxl]); 
            set(gca,'xgrid','on');
            set(gca,'fontsize',7);set(gca,'box','off');
            set(gca,'ticklength',[.03 .03]);
            plot([get(gca,'xlim')],[0 0],'r-');
            if pl <= (row-1)*col+1 & tpp ~= length(factors)
                set(gca,'xticklabel',[]);
            end;                    
            if pl <= (row-1)*col+2 | pl > (row-1)*col+3 & tpp ~= length(factors)
                set(gca,'yticklabel',[]);
            end;                    
        end;
    end;
    set(gcf,'Position',[100 300 1400 900]);
    set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
    set(ph,'fontsize',14);
