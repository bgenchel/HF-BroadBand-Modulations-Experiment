% plots output of ConcatTWdecomp()
%
%
%
%
%
% addmat [matrix] (context questions x epochs) only not [] when context included
% labels -- if cxt(2) is 1, then include context labels for plotting (strings)
% plottype -- ['winv' or 'wtdmean'] winv is straight context templates. 
% comporder -- [cell array] the same size as gdcomps (and must match 'complist' in each
%              decomposition .mat file

function  PlotClustConcatTW(filename,datset,fullpaths,clusttempls,clustdims,plottype,labels,comporder,savettl);
    
    shuffnum = 500;
    alpha = .01;
    maxics = [];
    for dim = 1:size(clustdims,1)
        s = load([fullpaths{clustdims(dim,1)},filename,'.mat']);     
        maxics = [maxics length(s.complist)];
    end;
    maxics = max(maxics);
    figure; pg = 1; 
    row = size(clustdims,1)*2-1; 
    if row > maxics-2
        row = maxics-3;
        if iseven(row)
            row = row-1;
        end;
    end;
    col = maxics;
    sbplot(row,col,[1 col])
    cols = hsv(size(clusttempls,2));
    ph = plot([0 size(clusttempls,2)+1],[0 0],'k-','linewidth',1.5);hold on;
    ph = plot(mean(clusttempls,1),'k-','linewidth',1.5);            
    hold on;
    for c = 1:size(clusttempls,2)
        ph = plot(c,clusttempls(:,c)','k.');hold on;
        set(ph,'color',cols(c,:));set(ph,'markersize',15);         
        ph=text(c,mean(clusttempls(:,c)'),labels{c});;
        set(ph,'color',cols(c,:));set(ph,'fontsize',16);
        set(ph,'rotation',45);
        set(gca,'xlim',[0 size(clusttempls,2)+1]);set(gca,'xticklabel',[]);
        title('Cluster Context Templates');
    end;
    pl = col+1; pll = pl;
    context = 'false'; % initialize and then change later if true
    for dim = 1:size(clustdims,1)       
        if pl > row * col
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            textsc(filename,'title');
            if ~isempty(savettl)
                str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.jpg -djpeg']; eval(str);
                str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.eps -depsc -adobe -painters']; eval(str);pg=  pg+1;
            end;
            figure; pl = 1; pg = pg+1;
            sbplot(row,col,[1 col])
            cols = hsv(size(clusttempls,2));
            ph = plot([0 size(clusttempls,2)+1],[0 0],'k-','linewidth',1.5);hold on;
            ph = plot(mean(clusttempls,1),'k-','linewidth',1.5);            
            for c = 1:size(clusttempls,2)
                ph = plot(c,clusttempls(:,c)','k.');hold on;
                set(ph,'color',cols(c,:));set(ph,'markersize',15);         
                ph=text(c,mean(clusttempls(:,c)'),labels{c});;
                set(ph,'color',cols(c,:));set(ph,'fontsize',16);
                set(ph,'rotation',45);
                set(gca,'xlim',[0 size(clusttempls,2)+1]);set(gca,'xticklabel',[]);
                title(['Cluster Context Templates-- Page ',int2str(pg)]);
            end;
            pl = col+1; pll = pl;
        end;
        EEG = pop_loadset(datset, fullpaths{clustdims(dim,1)});
        s = load([fullpaths{clustdims(dim,1)},filename,'.mat']);     
        wts = floatread([fullpaths{clustdims(dim,1)},filename,'.wts'],[s.pcs s.channo],[],0);
        sph = floatread([fullpaths{clustdims(dim,1)},filename,'.sph'],[s.channo s.channo],[],0);     
        dat = floatread([fullpaths{clustdims(dim,1)},filename,'.fdt'],[s.channo s.numframes],[],0); 
        ws = wts*sph;winv = inv(ws); acts = ws*dat;  
        clear sph wts ws
        if ~isempty(comporder)
            s.complist = s.complist(comporder{clustdims(dim,1)});
        end;
        
        if isfield(s,'pcersp')
            if ~isempty(s.pcersp) & s.pcersp > 0
                erspeig = floatread([fullpaths{clustdims(dim,1)},filename,'EIGVEC.fdt'],[length(s.freqs)*length(s.times)*length(s.complist) s.pcersp],[],0);
                erspdat = erspeig*dat(1:s.pcersp,:);
                ersptmpl = erspeig*winv(1:s.pcersp,:);
            end;
        end;
        if isfield(s,'pcctx')
            cxteig = floatread([fullpaths{clustdims(dim,1)},filename,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
            cxtdat = cxteig*dat(s.pcersp+1:end,:);  
            cxttmpl = cxteig*winv(s.pcersp+1:end,:);  
            cols = hsv(size(cxtdat,1));
            context = 'true';
        end;
        clear dat limersp
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
        end;
        clear acts bootwts 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(context,'true');
            for ic = 1:length(s.complist)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,s.complist(ic)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                title(['IC ',int2str(s.complist(ic))]); pl = pl+1;
            end;
            pl = pll+col; pll  = pl;
            d = abs(clustdims(dim,2));
            if clustdims(dim,2) < 0 % if there was a flip
                multfac = -1; % then multiply ERSP by -1
            else
                multfac = 1;
            end;
            for ic = 1:length(s.complist)
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
                if clustdims(dim,2) < 0
                    onetmp = onetmp * -1;
                end;
                sbplot(row,col,pl);pl = pl+1;
                if strcmp(s.freqscale,'quad')
                    [realy] = quadimagesc(s.times, s.freqs, onetmp*multfac,[-lim lim]);hold on;                    
                elseif  strcmp(s.freqscale,'log')                
                    mylogimagesc(s.times, s.freqs, onetmp*multfac,[-lim lim]);hold on;
                else
                    imagesc(s.times, s.freqs, onetmp*multfac,[-lim lim]);hold on; 
                    set(gca,'ytick',[10:10:s.freqs(end)]);
                end;
                plot([0 0],[get(gca,'ylim')],'k-');
                set(gca,'xticklabel',[]);
                if isfield(s,'medwarpevs')
                    wcols = lines(length(s.medwarpevs));
                    for wp = 1:length(s.medwarpevs)
                        ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                        set(ph,'color',wcols(wp,:));
                    end;
                end;
            end; 
            pl = pll+col; pll = pl;
        else % if no cotext?
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
                        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                        if ~isempty(savettl)
                            str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.jpg -djpeg']; eval(str);pg=  pg+1;
                        end;
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
    end;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    textsc(filename,'title');
    if ~isempty(savettl)
        str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.jpg -djpeg']; eval(str);
        str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.eps -depsc -adobe -painters']; eval(str);pg=  pg+1;
    end;
    
