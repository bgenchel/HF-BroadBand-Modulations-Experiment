%  PlotClustTW(filename,datset,fullpaths,clusttempls,clustdims,labels,savettl);
% 
%
%
%
%
% plotinv -- if 1, will plot the inverse of context and ersp templates

function  PlotClustTW(filename,datset,fullpaths,clusttempls,clustdims,plottype,auxclust,labels,savettl,plotinv);
    
    shuffnum = 250; alpha = .05;
    if isempty(auxclust)
        auxclust = cell(1,max(clustdims(:,1)));
        for nx = 1:length(auxclust)
            if ~isempty(find(clustdims(:,1) == nx))
                auxclust{nx} = clustdims(find(clustdims(:,1) == nx),2)';
            end;
        end;
    end;
    
    frqlim = [3 128];
    col = 6; row = 6;
    figure; pg = 1; 
    sbplot(row,col,[1 round(col/2)])
    cols = hsv(size(clusttempls,2));
    ph = plot([0 size(clusttempls,2)+1],[0 0],'k-','linewidth',1.5);
    hold on;
    if plotinv == 1
        clusttempls = clusttempls * -1;
    end;
    medtempls = mean(clusttempls,1);
    ph = plot(medtempls,'k-'); set(ph,'color',[.5 .5 .5]);
    for c = 1:size(clusttempls,2)
        ph = plot(c,clusttempls(:,c)','k.');hold on;
        set(ph,'color',cols(c,:));set(ph,'markersize',15);         
        ph=text(c,mean(clusttempls(:,c)'),labels{c});
        set(ph,'color',cols(c,:));set(ph,'fontsize',12);
        set(ph,'rotation',90); set(gca,'box','off');
        set(gca,'xlim',[0 size(clusttempls,2)+1]);set(gca,'xticklabel',[]);
        set(gca,'ylim',[min(clusttempls(:))+.2*min(clusttempls(:)) max(clusttempls(:))+.5*max(clusttempls(:))]);
        title('Cluster Context Templates');
    end;
    pl = col+1; pll = pl; 
    context = 'false'; % initialize and then change later if true
 
    for dim = 1:size(clustdims,1)
        if ~isempty(auxclust{clustdims(dim,1)})&find(ismember(auxclust{clustdims(dim,1)},clustdims(dim,2)))
            if pl > row * col
                set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                if ~isempty(savettl)
                    str = ['print ',fullpaths{1}(1:end-9),savettl,int2str(pg),'.eps -depsc -painters -adobe']; eval(str);pg=  pg+1;
                    figure;
                else
                    figure;
                end;
                pl = 1; pg = pg+1;
                sbplot(row,col,[1 round(col/2)])
                cols = hsv(size(clusttempls,2));
                ph = plot([0 size(clusttempls,2)+1],[0 0],'k-','linewidth',1.5);
                hold on;
                for c = 1:size(clusttempls,2)
                    ph = plot(c,clusttempls(:,c)','k.');hold on;
                    set(ph,'color',cols(c,:));set(ph,'markersize',15);         
                    ph=text(c,mean(clusttempls(:,c)'),labels{c});
                    set(ph,'color',cols(c,:));set(ph,'fontsize',12);
                    set(ph,'rotation',45);
                    set(gca,'xlim',[0 size(clusttempls,2)+1]);set(gca,'xticklabel',[]);
                    title(['Cluster Context Templates-- Page ',int2str(pg)]);
                end;
                pl = col+1; pll = pl;
            end;
            EEG = pop_loadset(datset, fullpaths{clustdims(dim,1)});
            s = load([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'.mat']);     

            wts = floatread([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'.wts'],[s.numchans s.numchans],[],0);
            sph = floatread([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'.sph'],[s.numchans s.numchans],[],0);     
            dat = floatread([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'.fdt'],[s.numchans s.numframes],[],0);             
            ws = wts*sph;winv = pinv(ws);   acts = ws*dat;
            clear sph wts ws
            frs = find(s.freqs>frqlim(1)&s.freqs<frqlim(2));
            if isfield(s,'pcersp')
                if ~isempty(s.pcersp) & s.pcersp > 0
                    erspeig = floatread([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcersp],[],0);
                    erspdat = erspeig*dat(1:s.pcersp,:);
                    ersptmpl = erspeig*winv(1:s.pcersp,:);
                end;
            end;
            if isfield(s,'pcctx')
                cxteig = floatread([fullpaths{clustdims(dim,1)},filename,int2str(clustdims(dim,2)),'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
                cxtdat = cxteig*dat(s.pcersp+1:end,:);  
                cxttmpl = cxteig*winv(s.pcersp+1:end,:);  
                cols = hsv(size(cxtdat,1));
                context = 'true';
            end;
            if strcmp(plottype,'wtdmean')
                fprintf('\nAccumulating bootstrap distribution ...\n');
                clear limersp limmat 
                %%%%%%%  calculate weighted ERSPs for each dim:
                wtdersp = (erspdat*acts')/size(erspdat,2);% multiply each trial element,sum then divide by ntrials 
                bootwts = zeros(size(wtdersp,1),size(wtdersp,2),shuffnum);
                for b= 1:shuffnum
                    bootwts(:,:,b) = (erspdat*shuffle(acts,2)')/size(erspdat,2);
                end;
                bootwts = sort(bootwts,3);
                limersp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
                limersp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
                 clear bootwts
                %%%%%%%  calculate context vectors for each dim:
                wtdctx = (cxtdat*acts')/size(erspdat,2); 
                %wtdctx = (addmat*acts')/size(erspdat,2); 
                bootwts = zeros(size(wtdctx,1),size(wtdctx,2),shuffnum);
                for b= 1:shuffnum
                    bootwts(:,:,b) = (cxtdat*shuffle(acts,2)')/size(erspdat,2);
                end;
                bootwts = sort(bootwts,3);
                limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
                limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
            end;
            clear sph wts ws dat acts bootwts 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(context,'true');
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,clustdims(dim,2)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                title(['IC ',int2str(clustdims(dim,2))]); pl = pl+1;
                d = abs(clustdims(dim,3));
                if clustdims(dim,3) < 0 % if there was a flip
                    multfac = -1; % then multiply ERSP by -1
                else
                    multfac = 1;
                end;

                if strcmp(plottype,'winv')
                    tmpcomp =  ersptmpl(:,d);
                elseif strcmp(plottype,'wtdmean')
                    tmpcomp = wtdersp(:,d);
                    % mask the wtd mean by the bootstrap values:
                    tmpcomp(find(tmpcomp>limersp(:,d,1)&tmpcomp<limersp(:,d,2))) = 0;
                end; 
               tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));
                lim = max(abs(tmpcomp(:)));
                if lim == 0
                    lim = 1;
                end;            
                
                if clustdims(dim,3) < 0 | plotinv == 1 % if template was flipped
                    tmpcomp = tmpcomp * -1;
                end;
                sbplot(row,col,pl);pl = pl+1;
                if strcmp(s.freqscale,'quad')
                    [realy] = quadimagesc(s.times, s.freqs(frs), tmpcomp(frs,:)*multfac,[-lim lim]);hold on;
                elseif  strcmp(s.freqscale,'log')                
                    mylogimagesc(s.times, s.freqs(frs), tmpcomp(frs,:)*multfac,[-lim lim]);hold on;
                else
                    imagesc(s.times, s.freqs(frs), tmpcomp(frs,:)*multfac,[-lim lim]);hold on; 
                    set(gca,'ytick',[10:10:s.freqs(end)]);
                end;
                plot([0 0],[get(gca,'ylim')],'k-');
                set(gca,'xticklabel',[]); set(gca,'ticklength',[.04 .04]);
                title(['Subject ',int2str(clustdims(dim,1))]);
                if isfield(s,'medwarpevs')
                    wcols = lines(length(s.medwarpevs));
                    for wp = 1:length(s.medwarpevs)
                        ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                        set(ph,'color',wcols(wp,:));
                    end;
                end;
                %cbar;
            else 
                for dd = 1:length(whichdims)
                    if pl > row*col
                        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                        if ~isempty(savettl)
                            str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.jpg -djpeg']; eval(str);pg=  pg+1;
                            figure;
                        else
                            figure;
                        end;
                        pl = 1;
                        sbplot(row,col,pl)
                    end;                
                    d = whichdims(dd);
                    sbplot(row,col,pl); pl = pl+1;
                    topoplot(EEG.icawinv(:,clustdims(dim,2)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                    title(['IC ',int2str(clustdims(dim,2))]);
                    lim = max(winv(:,d));
                    onetmp = winv(length(s.freqs)*length(s.times),d);
                    onetmp = reshape(onetmp,length(s.freqs),length(s.times));
                    if clustdims(dim,2) < 0 | plotinv == 1
                        onetmp = onetmp * -1;
                    end;

                    sbplot(row,col,[pl pl+1])
                    if strcmp(s.freqscale,'quad')
                        quadimagesc(s.times, s.freqs(frs), onetmp(frs,:),[-lim lim]);hold on;
                    elseif  strcmp(s.freqscale,'log')                
                        mylogimagesc(s.times, s.freqs(frs), onetmp(frs,:),[-lim lim]);hold on;
                    else
                        imagesc(s.times, s.freqs(frs), onetmp(frs,:),[-lim lim]);hold on;               
                        set(gca,'ytick',[10:10:s.freqs(frs(end))]);
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
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    if ~isempty(savettl)
        str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.eps -depsc -painters -adobe']; eval(str);pg=  pg+1;
        str = ['print ',fullpaths{1}(1:end-5),savettl,int2str(pg),'.jpg -djpeg']; eval(str);pg=  pg+1;
    end;
    
