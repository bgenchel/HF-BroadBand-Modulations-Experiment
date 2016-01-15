% plots output of ERSPdecomp()
%
%
%
%
%

function PlotERSPdecompContext(datset,name,fullpath,plotcomps,plotfacs,addmat);
    
    ms = 3;   % markersize 
    nlett = 6; nfb = 21; nresp = 6;% numbers of questions in the 3 categories
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
    sph=floatread([fullpath,name,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpath,name,'.wts'],[s.pcs s.numtrials],[],0); 
    icamatall = floatread([fullpath,name,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;     acts = ws*icamatall;  
    activations = acts(:,1:length(s.comps)*length(s.freqs)*length(s.times));
    winv = pinv(ws); clear wts sph ws
    frqlim = [0 s.freqs(end)];
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    linfreqs = linspace(s.freqs(1),s.freqs(end),length(s.freqs));
    [mesht meshf] = meshgrid(s.times,linfreqs);
    % create mean context vectors for each dimension
    snum = 500; % shuffle 300 times 
    fprintf('Calculating %s shuffle permutations of trial weights',int2str(snum));
    for im = 1:size(winv,2)
        for trial = 1:size(addmat,2)
            wtqs(:,trial) = addmat(:,trial)*winv(trial,im);
        end;
        meanqs(:,im) = mean(wtqs,2);        
        for rep =1:snum        
            randwts = shuffle(winv(:,im));
            for trial = 1:size(addmat,2)
                randctx(:,rep,trial) = addmat(:,trial)*randwts(trial);
            end;
        end;
        randmean(:,:,im) = mean(randctx,3);
    end;
        
    % make signif masks for context trial plotting
    for im = 1:size(randmean,3)
        for q = 1:size(randmean,1)
            currvec = randmean(q,:,im);currvec = sort(currvec);
            minmask(q,im) = currvec(round(length(currvec)*.01));
            maxmask(q,im) = currvec(round(length(currvec) - (length(currvec)*.01)));
        end;    
    end;
    figure;row =length(plotfacs)+1; 
    if row > 18
        row = round(row/2);
    end;            
    col = length(plotcomps)+9;
    pl = 2;           
    for cp = 1:length(plotcomps)
        sbplot(row,col,pl)
        topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
        set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
    end;pl = col+1;
           
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
            end; pl = col+1;
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
        sbplot(row,col,[pl pl+7])
        ph = plot([1:nlett],meanqs(1:nlett,tp),'ro-','linewidth',2);hold on;
        set(ph,'markersize',ms);
        ph = plot([nlett+1:nlett+nfb],meanqs(nlett+1:nlett+nfb,tp),'bs-','linewidth',2); hold on;
        set(ph,'markersize',ms);
        ph = plot([nlett+nfb+1:nlett+nfb+nresp],meanqs(nlett+nfb+1:nlett+nfb+nresp,tp),'g^-','linewidth',2);
        set(ph,'markersize',ms+ms/2); hold on;
        set(gca,'xlim',[0 nlett+nfb+nresp+1]);%set(gca,'ylim',[-ctlim ctlim]);
        ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
        yl = get(gca,'ylim');starlev = (yl(2) - yl(1))*.1;
        for q = 1:size(minmask,1)
            if meanqs(q,tp)>minmask(q,tp)& meanqs(q,tp)<maxmask(q,tp)
                ph = plot(q,meanqs(q,tp),'r.');set(ph,'color',[.8 .8 .8]);
                set(ph,'markersize',ms*5);
            else
                ph = plot([q q],[yl(2)-starlev yl(2)-starlev],'r*');
                set(ph,'markersize',6);set(ph,'color',[1 .5 0]);
            end; 
        end;
        
        if pl > (row-1)*col | (pl> (row/2)*col& pl<(row/2+1)*col)
            CnxtLabelRewFB(max(yl),.3,7); % fontsize=9
        end;
        pl = pl+8;
    end;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    ph=textsc(['ERSP Templates for ',fullpath(end-4:end),': ',name],'title');
    set(ph,'fontsize',14);
    axcopy
    
