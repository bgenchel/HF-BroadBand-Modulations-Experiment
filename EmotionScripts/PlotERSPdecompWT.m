% plots output of ERSPdecomp()
%
%
%
%
%

function PlotERSPdecomp(datset,name,fullpath,plotcomps,plotfacs);
    
    
    EEG = pop_loadset(datset, fullpath);
    s = load([fullpath,name,'.mat']);     
    for ff = 1:length(s.freqs)-1
        fspace(1,ff) = s.freqs(ff+1) - s.freqs(ff);
    end;
    if fspace(4) > fspace(1) & fspace(10) > fspace(1)
        logyes = 1;
    else
        logyes = 0;
    end;
    sph=floatread([fullpath,name,'.sph'],[s.numrows s.numrows],[],0); 
    wts=floatread([fullpath,name,'.wts'],[s.pcs s.numrows],[],0); 
    icamatall = floatread([fullpath,name,'.fdt'],[s.numrows s.numframes],[],0);    
    ws = wts*sph;    activations = ws*icamatall;    
    winv = pinv(ws); clear wts sph ws
    frqlim = [0 s.freqs(end)];
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    figure;row =length(plotfacs)+3; 
    if row > 18
        row = round(row/2);
    end;            
    col = length(plotcomps)+1;
    pl = 2;           
    for cp = 1:length(plotcomps)
        sbplot(row,col,pl)
        topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
        set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
    end;
    % plot mean baseline spectrum
    pl = pl+1;
    for cp = 1:size(s.meanpwr,1)
        sbplot(row,col,pl)
        blersp = s.meanpwr(cp,:);
        blersp = reshape(blersp,length(s.freqs),length(s.times));
        ph = plot(s.freqs,mean(blersp,2),'b-','linewidth',1.5);hold on;
        plot([10 10],[get(gca,'ylim')],'r-');
        set(gca,'xlim',[s.freqs(1) s.freqs(end)]);
        set(gca,'yticklabel',[]);
        pl = pl+1;
    end;
    
    
    % plot average ERSP
    pl = pl+1;
    linfreqs = linspace(s.freqs(1),s.freqs(end),length(s.freqs));
    [mesht meshf] = meshgrid(s.times,linfreqs);
    for cp = 1:size(s.meanpwr,1)
        sbplot(row,col,pl)
        blersp = s.meanpwr(cp,:);
        blersp = reshape(blersp,length(s.freqs),length(s.times));
        for t = 1:size(blersp,2)
            blersp(:,t) = blersp(:,t) - mean(blersp,2);
        end;        
        datout = griddata(s.times,s.freqs,blersp,mesht,meshf);
        lim = max(max(abs(datout)));
        ph = imagesc(s.times,linfreqs,datout,[-lim lim]);hold on;
        %ph = imagesc(s.times,s.freqs,blersp,[-lim lim]);hold on;
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        plot([0 0],[get(gca,'ylim')],'k-');
        plot([get(gca,'xlim')],[10 10],'k-');
        plot([get(gca,'xlim')],[30 30],'k-');
        set(gca,'ydir','norm');
        pl = pl+1;
    end;
    
           
    lim = max(max(abs(activations(plotfacs,:))));
    for tpp = 1:length(plotfacs)
        tp = plotfacs(tpp);
        if pl > row*col
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['ERSP Templates for  ',fullpath(end-4:end),': ',name],'title');
            set(ph,'fontsize',14); axcopy
            figure;
            pl = 2;           
            for cp = 1:length(plotcomps)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
                set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
            end;
        end;
        sbplot(row,col,pl)
        hist(winv(:,tp),75);pl = pl+1;hold on;
        set(gca,'fontsize',7);set(gca,'xlim',[-2 2]);
        plot([0 0],[get(gca,'ylim')],'r-');
        set(gca,'yticklabel',[]);                
        set(gca,'xticklabel',[]); 
        title(['IM ',int2str(tp)]);
        for cp = 1:length(plotcomps)
            rcp = find(plotcomps(cp) == s.comps);
            sbplot(row,col,pl)
            oneersp = activations(tp,(length(s.freqs)*length(s.times))*(rcp-1)+1:(length(s.freqs)*length(s.times))*rcp);
            oneersp = reshape(oneersp,length(s.freqs),length(s.times));
            datout = griddata(s.times,s.freqs,oneersp,mesht,meshf);
            imagesc(s.times,s.freqs,datout,[-lim lim]);                    
            pl = pl+1;hold on;
            set(gca,'ydir','norm');
            %if logyes == 1
                %set(gca,'yscale','log');
                %set(gca,'ytick',[4 6 10 20 40]);
            %end;
            set(gca,'fontsize',7);
            set(gca,'ticklength',[.03 .03]);
            if tpp ~= length(plotfacs) %| tpp ~= row-1
                set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
            end; 
            plot([0 0],[get(gca,'ylim')],'k-');
            plot([get(gca,'xlim')],[10 10],'k-');
            plot([get(gca,'xlim')],[30 30],'k-');
        end;
    end;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    ph=textsc(['ERSP Templates for ',fullpath(end-4:end),': ',name],'title');
    set(ph,'fontsize',14);
    axcopy
