% Plots for TW spectral decomp a mean ERSP with +/- modulations for specified IM
% and specified time slices as over-plotted spectra
%
%
%
%
%
%
%



function PlotTWerspMean(filename,datset,fullpath,ic,im,labels,tmpnts);
    
    fs = 10;
    EEG = pop_loadset(datset, fullpath);
    name = [filename,int2str(ic)];  
    s = load([fullpath,name,'.mat']);   
    if s.pcs==0
        s.pcs= s.channo;    
    end;
    wts = floatread([fullpath,name,'.wts'],[s.pcs s.channo],[],0);
    sph = floatread([fullpath,name,'.sph'],[s.channo s.channo],[],0);     
    dat = floatread([fullpath,name,'.fdt'],[s.channo s.numframes],[],0);   
    ws = wts*sph;winv = pinv(ws);acts = ws*dat;
    mntf = reshape(s.icmean,length(s.freqs),length(s.times));  
    
    if isfield(s,'pcersp')
        if s.pcersp > 0               
            erspeig = floatread([fullpath,name,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcersp],[],0);
            erspdat = erspeig*dat(1:s.pcersp,:);
            ersptmpl = erspeig*winv(1:s.pcersp,im);  
            ersptmpl = reshape(ersptmpl,length(s.freqs),length(s.times)); 
        end;
    end;
    if isfield(s,'pcctx')
        if s.pcctx > 0
            cxteig = floatread([fullpath,name,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
            cxtdat = cxteig*dat(s.pcersp+1:end,im);  
            cxttmpl = cxteig*winv(s.pcersp+1:end,im);  
        end;
    end;
    figure;
    col = 6;
    row = 2 + ceil(length(tmpnts)/3);
    
    sbplot(row*2,col,1)
    topoplot(EEG.icawinv(:,ic),EEG.chanlocs(EEG.icachansind),'electrodes','off');
    title(['IC ',int2str(ic)]);
    sbplot(row*2,col,col+1)
    hist(acts(im,:),100); hold on;
    plot([0 0],[get(gca,'ylim')],'r-');
    set(gca,'xlim',[-2.9 2.9]);
        
    
    sbplot(row,col,[2 3])
    set(gca,'fontsize',fs);
    elim = max(abs(ersptmpl(:)));
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.times, s.freqs, ersptmpl,[-elim elim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.times, s.freqs, ersptmpl,[-elim elim]);hold on;
    else
        imagesc(s.times, s.freqs, ersptmpl,[-elim elim]);hold on;               
        set(gca,'ytick',[10:10:s.freqs(end)]);
    end;
    set(gca,'ydir','norm'); set(gca,'ticklength',[.04 .04]);
    hold on; plot([0 0],[get(gca,'ylim')],'k-');
    if isfield(s,'medwarpevs')
        wcols = lines(length(s.medwarpevs));
        for wp = 1:length(s.medwarpevs)
            ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
            set(ph,'color',wcols(wp,:));
        end;
    end;  
    title(['ERSP modulation template']);
    
    sbplot(row,col,[4 6])
    set(gca,'fontsize',fs);
    cols = hsv(length(labels)); cols(3,:) = [.8 .8 0];
    for ee = 1:length(cxttmpl)
        ph=plot(ee,cxttmpl(ee),'.'); hold on;
        set(ph,'markersize',15);                        
        set(ph,'color',cols(ee,:));
        ph=text(ee,cxttmpl(ee),labels{ee});;
        set(ph,'color',cols(ee,:));set(ph,'fontsize',14);
        set(ph,'rotation',90);
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.75 .75 .75]);                
    plot([1:length(cxttmpl)], cxttmpl','k-','linewidth',1);hold on;             
    set(gca,'xlim',[0 length(cxttmpl)+1]);
    if ~isempty(tmpnts)
    for tm = 1:length(tmpnts)
        ph = plot([tmpnts(tm) tmpnts(tm)],[get(gca,'ylim')],'k--','linewidth',3);
    end;
    end;
    set(gca,'box','off');
    title(['Context template IM',int2str(im)]);
    %%%%%----------------------------------------------------------
    % pull out max and min trials:
    
    [x y] = sort(acts(im,:),2); % y(1) most neg, y(end) most pos
    posvar = winv(1:s.pcersp,im) * acts(im,y(end));
    posvar = erspeig * posvar;
    %posvar = erspdat(:,y(end)); % from ersp data matrix instead
    posvar = reshape(posvar,length(s.freqs),length(s.times)); 
    posvar = mntf + posvar;
    
    negvar = winv(1:s.pcersp,im) * acts(im,y(1));
    negvar = erspeig * negvar;
    %negvar = erspdat(:,y(1)); % from ersp data matrix instead
    negvar = reshape(negvar,length(s.freqs),length(s.times)); 
    negvar = mntf + negvar;

    posersp = posvar- repmat(mean(mntf,2),[1 size(mntf,2)]);
    negersp = negvar- repmat(mean(mntf,2),[1 size(mntf,2)]);
    mnersp = mntf - repmat(mean(mntf,2),[1 size(mntf,2)]);
    
    lim1 = max(abs(posersp(:)));lim2 = max(abs(negersp(:)));lim = max([lim1 lim2]);

    %%%%  Plot the mean ERSP -----------------------------------------------
    sbplot(row,col,[col+3 col+4])
    set(gca,'fontsize',fs);
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.times, s.freqs, mnersp,[-lim lim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.times, s.freqs, mnersp,[-lim lim]);hold on;
    else
        imagesc(s.times, s.freqs, mnersp,[-lim lim]);hold on;               
        set(gca,'ytick',[10:10:s.freqs(end)]);
    end;
    set(gca,'ydir','norm'); set(gca,'ticklength',[.04 .04]);
    hold on; plot([0 0],[get(gca,'ylim')],'k-');
    if isfield(s,'medwarpevs')
        wcols = lines(length(s.medwarpevs));
        for wp = 1:length(s.medwarpevs)
            ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
            set(ph,'color',wcols(wp,:));
        end;
    end; 
    if ~isempty(tmpnts)
    for tm = 1:length(tmpnts)
        ph = plot([tmpnts(tm) tmpnts(tm)],[get(gca,'ylim')],'k--','linewidth',3);
    end;
    end;
    title(['Mean ERSP']);
    %%%% -----------------------------------------------
    % Plot the positive variation for specified IM:
    
    sbplot(row,col,[col+5 col+6])
    set(gca,'fontsize',fs);
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.times, s.freqs, posersp,[-lim lim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.times, s.freqs, posersp,[-lim lim]);hold on;
    else
        imagesc(s.times, s.freqs, posersp,[-lim lim]);hold on;               
        set(gca,'ytick',[10:10:s.freqs(end)]);
    end;
    set(gca,'ydir','norm'); set(gca,'ticklength',[.04 .04]);
    hold on; plot([0 0],[get(gca,'ylim')],'k-');
    if isfield(s,'medwarpevs')
        wcols = lines(length(s.medwarpevs));
        for wp = 1:length(s.medwarpevs)
            ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
            set(ph,'color',wcols(wp,:));
        end;
    end;  
    if ~isempty(tmpnts)
    for tm = 1:length(tmpnts)
        ph = plot([tmpnts(tm) tmpnts(tm)],[get(gca,'ylim')],'k--','linewidth',3);
    end;
    end;
    title(['Positively weighted trial']);
    cbar;
    %%%% -----------------------------------------------
    % Plot the negative variation for specified IM:
    [x y] = sort(acts(im,:),2); % y(1) most neg, y(end) most pos
    
    
    sbplot(row,col,[col+1 col+2])
    set(gca,'fontsize',fs);
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.times, s.freqs, negersp,[-lim lim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.times, s.freqs, negersp,[-lim lim]);hold on;
    else
        imagesc(s.times, s.freqs, negersp,[-lim lim]);hold on;               
        set(gca,'ytick',[10:10:s.freqs(end)]);
    end;
    set(gca,'ydir','norm'); set(gca,'ticklength',[.04 .04]);
    hold on; plot([0 0],[get(gca,'ylim')],'k-');
    if isfield(s,'medwarpevs')
        wcols = lines(length(s.medwarpevs));
        for wp = 1:length(s.medwarpevs)
            ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
            set(ph,'color',wcols(wp,:));
        end;
    end;
    if ~isempty(tmpnts)
    for tm = 1:length(tmpnts)
        ph = plot([tmpnts(tm) tmpnts(tm)],[get(gca,'ylim')],'k--','linewidth',3);
    end;
    end;
    ylabel('Frequency (Hz)');
    xlabel('Latency (ms)');
    title(['Negatively weighted trial']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Plot power spectra from specified time points  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(tmpnts)
    pl = col*2+1;
    for tm = 1:length(tmpnts)
        if pl == col*3 + 1
            pl = pl+1;
        end;
        tmidx = find(s.times >= tmpnts(tm));
        sbplot(row,col,[pl pl + 1]); pl = pl+2;
        set(gca,'fontsize',fs);
        %quadplot(s.freqs,mntf(:,tmidx(1)),2,'k'); hold on;
        quadplot(s.freqs,posvar(:,tmidx(1)) - mntf(:,tmidx(1)),1.5,'r');hold on;
        quadplot(s.freqs,negvar(:,tmidx(1)) - mntf(:,tmidx(1)),1.5,'b');
        ph = plot([get(gca,'xlim')],[0 0],'k-'); 
        title(['Freq. Modulation at ',int2str(tmpnts(tm)),' ms']);
        if tm == 1
            xlabel('Frequency (Hz)');
            ylabel('Power (db)');
        end;
    end;
    
        
    end;

