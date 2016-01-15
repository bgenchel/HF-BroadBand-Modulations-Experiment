% plots output of ConcatTWdecomp()
%
%
%
%
%
% addmat [matrix] (context questions x epochs) only not [] when context included
% labels -- if cxt(2) is 1, then include context labels for plotting (strings)
% plottype -- ['winv' or 'wtdmean'] winv is straight context templates. 

function [cxtout] = PlotConcatTW(filename,datset,fullpath,addmat,whichdims,labels,alpha,plottype,savettl,maxcol);
    
    
    
    if ~exist('maxcol')
        maxcol = 6;
    elseif isempty(maxcol)
        maxcol = 6;
    end;
    shuffnum = 500;
    context = 'false'; % initialize and then change later if true
    EEG = pop_loadset(datset, fullpath);
    s = load([fullpath,filename,'.mat']);     
    wts = floatread([fullpath,filename,'.wts'],[s.pcs s.channo],[],0);
    sph = floatread([fullpath,filename,'.sph'],[s.channo s.channo],[],0);     
    dat = floatread([fullpath,filename,'.fdt'],[s.channo s.numframes],[],0); 
    ws = wts*sph;winv = inv(ws);    acts = ws*dat;
    clear sph wts ws
    if isempty(whichdims)
        whichdims = [1:size(winv,2)];
    end;
    if isfield(s,'pcersp')
        if ~isempty(s.pcersp) & s.pcersp > 0
            erspeig = floatread([fullpath,filename,'EIGVEC.fdt'],[length(s.freqs)*length(s.times)*length(s.complist) s.pcersp],[],0);
            erspdat = erspeig*dat(1:s.pcersp,:);
            ersptmpl = erspeig*winv(1:s.pcersp,:);
        end;
    end;
    if isfield(s,'pcctx')
        if isfield(s,'useeps')
            addmat = addmat(:,s.useeps);
        end;
        if ~isempty(s.pcctx) & s.pcctx > 0
            delqs = []; % delete zero'd out questions
            for q = 1:size(addmat,1)
                if isempty(find(addmat(q,:)))
                    delqs = [delqs q];
                end;
            end;
            if ~isempty(delqs)
                addmat(delqs,:) = [];
            end;
            deleps = [];
            for ep = 1:size(addmat,2)
                x = find(addmat(:,ep) == 0);
                if ~isempty(x)
                    deleps = [deleps ep];
                end;
            end;
            if ~isempty(deleps)
                addmat(:,deleps)=[];
            end;
            cxteig = floatread([fullpath,filename,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
            cxtdat = cxteig*dat(s.pcersp+1:end,:);  
            cxttmpl = cxteig*winv(s.pcersp+1:end,:);  
            cols = hsv(size(addmat,1));
        end;
        context = 'true';
    end;
    
    clear dat
    if strcmp(plottype,'wtdmean')
        fprintf('\nAccumulating bootstrap distribution ...\n');
        %%%%%%%  calculate weighted ERSPs for each dim:
        wtdersp = (erspdat*acts')/size(erspdat,2);
        % multiply each trial element,sum then divide by ntrials 
        bootwts = zeros(size(wtdersp,1),size(wtdersp,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (erspdat*shuffle(acts,2)')/size(erspdat,2);
        end;    
        bootwts = sort(bootwts,3);
        limersp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limersp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        %%%%%%%  calculate context vectors for each dim:
        wtdctx = (cxtdat*acts')/size(erspdat,2); % use back-proj addmat
        %wtdctx = (addmat*acts')/size(erspdat,2);  % use original addmat
        % multiply each trial element,sum then divide by ntrials 
        bootwts = zeros(size(wtdctx,1),size(wtdctx,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (cxtdat*shuffle(acts,2)')/size(erspdat,2);
        end;
        bootwts = sort(bootwts,3);
        limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        mnctx = mean(bootwts,3);
        stdctx = std(bootwts,0,3);
    end;
    clear acts bootwts 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pg = 1;
    if strcmp(context,'true')
        row = length(s.complist)+1; col = length(whichdims)+1;
        if col > maxcol
            col = maxcol;
        end;
        for doned = 0:col-1:length(whichdims)
            diffd = length(whichdims) - doned; 
            if doned+1 <= length(whichdims)                
                set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                textsc([filename,'-',savettl],'title');
                figure;pl = 2; pg = pg+1;
            end;
            for dd = doned+1:doned+(col-1)
                if dd > length(whichdims)
                    pl = pl+(col-1) - diffd;break;
                end;
                d = whichdims(dd);
                sbplot(row,col,pl);
                ph = plot([0 size(cxttmpl,1)+1],[0 0],'k-'); hold on;
                set(ph,'color',[.75 .75 .75]);
                if strcmp(plottype,'winv')
                    for ee = 1:size(cxttmpl,1)
                        ph=plot(ee,cxttmpl(ee,d),'.'); hold on;
                        set(ph,'markersize',15); 
                        set(ph,'color',cols(ee,:));
                    end;
                    for ee = 1:size(cxttmpl,1)
                        ph=text(ee,cxttmpl(ee,d),labels{ee});
                        if length(labels) > 10
                        set(ph,'color',cols(ee,:));set(ph,'fontsize',5);
                        else
                        set(ph,'color',cols(ee,:));set(ph,'fontsize',9);
                        end;
                        set(ph,'rotation',90); 
                    end;
                    cxtout(:,d) = cxttmpl(:,d);
                elseif strcmp(plottype,'wtdmean')                        
                    tmpcomp = wtdctx(:,d);
                    for ee = 1:size(cxttmpl,1)
                        ph = plot([ee ee],[limmat(ee,d,1) limmat(ee,d,2)],'k-','linewidth',5);
                        set(ph,'color',[.8 .8 .8]);hold on;
                        if tmpcomp(ee)>limmat(ee,d,1)&tmpcomp<limmat(ee,d,2)
                            ph=plot(ee,tmpcomp(ee),'.');set(ph,'markersize',14); 
                        else
                            ph=plot(ee,tmpcomp(ee),'*');set(ph,'markersize',9); 
                        end;
                        
                        set(ph,'color',cols(ee,:));
                    end;
                    for ee = 1:size(cxttmpl,1)
                        ph=text(ee,tmpcomp(ee),labels{ee});
                        if length(labels) > 10
                            set(ph,'color',cols(ee,:));set(ph,'fontsize',5);
                        else
                            set(ph,'color',cols(ee,:));set(ph,'fontsize',9);
                        end;                        
                        set(ph,'rotation',90); 
                    end;
                    % mask the wtd mean by the bootstrap values(for output):
                    tmpcomp(find(tmpcomp>limmat(:,d,1)&tmpcomp<limmat(:,d,2))) = 0;
                    cxtout(:,d) = tmpcomp;
                end;
                set(gca,'xlim',[0 size(cxttmpl,1)+1]);
                title(['Dim ',int2str(d)]);
                pl = pl+1;
            end;
            for ic = 1:length(s.complist)
                if doned+1 <= length(whichdims)
                    sbplot(row,col,pl)
                    topoplot(EEG.icawinv(:,s.complist(ic)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                    title(['IC ',int2str(s.complist(ic))]); pl = pl+1;
                end;
                for dd = doned+1:doned+col-1
                    if dd > length(whichdims)
                        pl = pl+(col-1) - diffd;break;
                    end;
                    d = whichdims(dd);
                    lim = max(abs(erspdat(:,d)));
                    if strcmp(plottype,'winv')
                        onetmp = erspdat(length(s.freqs)*length(s.times)*(ic-1)+1:length(s.freqs)*length(s.times)*ic,d);
                    elseif strcmp(plottype,'wtdmean')      
                        onetmp = wtdersp(length(s.freqs)*length(s.times)*(ic-1)+1:length(s.freqs)*length(s.times)*ic,d);
                        onetmp(find(onetmp>limersp(length(s.freqs)*length(s.times)*(ic-1)+1:length(s.freqs)*length(s.times)*ic,d,1)&onetmp<limersp(length(s.freqs)*length(s.times)*(ic-1)+1:length(s.freqs)*length(s.times)*ic,d,2))) = 0;
                    end;
                    
                    lim = max(abs(onetmp));
                    if lim == 0
                        lim = 1;% if no sig points, make lim a valid value
                    end;
                    onetmp = reshape(onetmp,length(s.freqs),length(s.times));
                    sbplot(row,col,pl);
                    if strcmp(s.freqscale,'quad')
                        [realy] = quadimagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;                    
                    elseif  strcmp(s.freqscale,'log')                
                        mylogimagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;
                    else
                        imagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;               
                        set(gca,'ytick',[10:10:s.freqs(end)]);
                    end;
                    plot([0 0],[get(gca,'ylim')],'k-');
                    if isfield(s,'medwarpevs')
                        wcols = lines(length(s.medwarpevs));
                        for wp = 1:length(s.medwarpevs)
                            ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                            set(ph,'color',wcols(wp,:));
                        end;
                    end; 
                    pl = pl+1;
                end;
            end;
            if ~isempty(savettl)
                textsc([filename,'-',savettl],'title');
                str = ['print ',fullpath,savettl,int2str(pg),'.jpg -djpeg']; eval(str);
            end;
        end;
    else 
        figure; row = size(winv,2);col = length(s.complist);pl = 1;
        if row > maxcol
            row = maxcol;
        end;
        for ic = 1:length(s.complist)
            sbplot(row,col,pl)
            topoplot(EEG.icawinv(:,s.complist(ic)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
            title(['IC ',int2str(s.complist(ic))]); pl = pl+1;
        end;
        for dd = 1:length(whichdims)
            d = whichdims(dd);
            lim = max(winv(:,d));
            for ic = 1:length(s.complist)
                onetmp = winv(length(s.freqs)*length(s.times)*(ic-1)+1:length(s.freqs)*length(s.times)*ic,d);
                onetmp = reshape(onetmp,length(s.freqs),length(s.times));
                if pl > row*col
                    textsc(filename,'title');
                    figure; pl = 1;
                    for ic = 1:length(s.complist)
                        sbplot(row,col,pl)
                        topoplot(EEG.icawinv(:,s.complist(ic)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                        title(['IC ',int2str(s.complist(ic))]); pl = pl+1;
                    end;
                end;
                sbplot(row,col,pl)
                if strcmp(s.freqscale,'quad')
                    quadimagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;
                elseif  strcmp(s.freqscale,'log')                
                    mylogimagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;
                else
                    imagesc(s.times, s.freqs, onetmp,[-lim lim]);hold on;               
                    set(gca,'ytick',[10:10:s.freqs(end)]);
                end;
                if isfield(s,'medwarpevs')
                    wcols = lines(length(s.medwarpevs));
                    for wp = 1:length(s.medwarpevs)
                        ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                        set(ph,'color',wcols(wp,:));
                    end;
                end; 
                pl = pl+1;
            end;
        end;        
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    textsc([filename,'-',savettl],'title');
