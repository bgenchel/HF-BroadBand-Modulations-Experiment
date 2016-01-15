% same as PlotERSPdecomp() but no scalp maps


function PlotMatdecomp(name,fullpath,plotcomps,plotfacs);

   
    s = load([fullpath,name,'.mat']);     
    %for ff = 1:length(s.freqs)-1
    %    fspace(1,ff) = s.freqs(ff+1) - s.freqs(ff);
    %end;
%    if fspace(4) > fspace(1) & fspace(10) > fspace(1)
%        logyes = 1;
%    else
%        logyes = 0;      % 
%    end;
     logyes = 1;
     sph=floatread([fullpath,name,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpath,name,'.wts'],[s.pcs s.numtrials],[],0); 
    icamatall = floatread([fullpath,name,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*icamatall;    
    winv = pinv(ws); clear wts sph ws
    lim = max(max(abs(activations(plotfacs,:))));
    frqlim = [0 s.freqs(end)];
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    figure;row =length(plotfacs)+1; 
    if row > 16
        row = round(row/2);
    end;            
    col = length(plotcomps)+1;
    pl = 1;           
           
    for tpp = 1:length(plotfacs)
        tp = plotfacs(tpp);
        if pl > row*col
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['ERSP Templates for  ',fullpath(end-4:end),': ',name],'title');
            set(ph,'fontsize',14); axcopy
            figure;
            pl = 1;           
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
            imagesc(s.times,s.freqs,oneersp,[-lim lim]);                    
            pl = pl+1;hold on;
            set(gca,'ydir','norm');
            if logyes == 1
                set(gca,'yscale','log');
                set(gca,'ytick',[4 6 10 20 40]);
            end;
            set(gca,'fontsize',7);
            set(gca,'ticklength',[.03 .03]);
            if tpp == 1
                title(int2str(rcp));
            end;                    
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
