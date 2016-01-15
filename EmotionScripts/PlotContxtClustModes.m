function [keepersp1,keepersp2,keepidx1,keepidx2,times,freqs]= PlotContxtClustModes(stem,icidx,fullpaths,vertlines,clabels,nbins)
    
    keepidx1 = []; keepidx2 = [];
    keepersp1= []; keepersp2 = [];
    figure; row = 5; col = 5; pl = 1;
    d = 1; keepidx = []; keepersp=[];
    for im = 1:size(icidx,1)
        if pl > row*col
            figure; pl = 1;
        end;
        dim = icidx(im,3);
        nx = icidx(im,1);
        ic = icidx(im,2);
        wtsphname=[stem,int2str(icidx(im,2))];
        s = load ([fullpaths{nx},wtsphname,'.mat']);
        sph=floatread([fullpaths{nx},wtsphname,'.sph'],[s.numchans s.numchans],[],0); 
        wts=floatread([fullpaths{nx},wtsphname,'.wts'],[s.numchans s.numchans],[],0); 
        icamat=floatread([fullpaths{nx},wtsphname,'.fdt'],[s.numchans s.numframes],[],0);    
        ws = wts*sph;winv = pinv(ws);    acts = ws*icamat;  
        
        if ~isempty(s.pcersp)
            erspeig = floatread([fullpaths{nx},wtsphname,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcersp],[],0);
            erspdat = erspeig*icamat(1:s.pcersp,:);
        end;    
        if ~isempty(s.pcctx)
            idxeig = floatread([fullpaths{nx},wtsphname,'ADDEIG.fdt'],[length(s.cxtmean) s.pcctx],[],0);
            idxdat = idxeig*icamat(s.pcersp+1:end,:);
        end;
        
        if dim < 0
            mfac = -1;dim=abs(dim);
        else
            mfac=1;
        end;
        % find low pass of weights distribution hist:
        [n,x] = hist(acts(dim,:)*mfac,nbins); % low pass version of distribution to find peaks
        distrs(d,:) = n;% save distribution pattern for output
        distrtms(d,:) = x;
        % put zeros on ends to find peaks at edges:
        clear nn xx
        nn(1) = 0; nn = [nn n];
        nn(end+1) = 0; n = nn;
        xx(1) = x(1); xx = [xx x]; xx(end+1) = x(end); x = xx;
        peakval = []; % find peaks:
        for y = 2:length(n)-1
            if n(y) > n(y-1) & n(y) > n(y+1)
                peakval = [peakval x(y)]; % return the actual wt val
            end;
        end;        
        %if n(1) > n(2)
        %    peakval = [x(1) peakval];
        %end;
        %if n(end) > n(end-1)
        %    peakval = [peakval x(end)];
        %end;
        peaksurr = (max(x) - min(x))/12; % a 12th of total wt value span
                                         % pull out each trial set individually

        sbplot(row,col,pl); pl = pl+1;
        hist(acts(dim,:),nbins); hold on;
        ph = plot([0 0],[get(gca,'ylim')],'r-');
        title(['Dim ',int2str(dim),'-All Wts ']);
        for p = 1:length(peakval)
            if pl > row*col
                figure; pl = 1;
            end;
            trset = find(acts(dim,:)*mfac > peakval(p)-peaksurr & acts(dim,:)*mfac < peakval(p)+peaksurr);
            plotdat = mean(erspdat(:,trset),2)*mfac;
            plotidx = mean(idxdat(:,trset),2)*mfac; 
            keepe(p,:) = mean(erspdat(:,trset),2)'*mfac;
            keepi(p,:) = mean(idxdat(:,trset),2)'*mfac;
            lim = max(abs(plotdat));
                if lim == 0
                lim = 1;
            end;
            sbplot(row,col,pl); pl = pl+1;
            if strcmp(s.freqscale,'quad')
                quadimagesc(s.times,s.freqs,reshape(plotdat,length(s.freqs),length(s.times)),[-lim lim]);hold on;
            elseif strcmp(s.freqscale,'log')
                mylogimagesc(s.times,s.freqs,reshape(plotdat,length(s.freqs),length(s.times)),[-lim lim]);hold on;
                set(gca,'ydir','norm'); 
            end;
            for v = 1:length(vertlines)
                plot([vertlines(v) vertlines(v)],[get(gca,'ylim')],'k-');
            end;   
            title(['Dim ',int2str(dim),'-Mode ',int2str(p)]);
            cbar;
            clim = max(abs(plotidx)) + max(abs(plotidx))*.05;
             
            if clim == 0
                clim = 1;
            end;
            
            if pl > row*col
                figure; pl = 1;
            end;
            sbplot(row,col,pl); pl = pl+1;
            ph = plot([1:size(idxdat,1)],plotidx,'k-','linewidth',1);
            set(ph,'color',[.5 .5 .5]);
            set(gca,'xlim',[0 size(idxdat,1)+1]); hold on;
            ph = plot([get(gca,'xlim')],[0 0],'r-');
            for q = 1:size(idxdat,1)
                ph = plot(q,plotidx(q),'bo');
            end;
            set(gca,'ylim',[-clim clim]);
            for q = 1:size(idxdat,1)
                ph = text(q,-clim+abs(-clim*.05),clabels{q}); 
                set(ph,'rotation',90);
            end;
            title(['Peak Val ',num2str(peakval(p))]);
        end;
        keepersp1 = [keepersp1; keepe(1,:)];
        keepersp2 = [keepersp2; keepe(end,:)];
        keepidx1 = [keepidx1; keepi(1,:)];
        keepidx2 = [keepidx2; keepi(end,:)];
    end;
    textsc([wtsphname,'-- Selected Mode Trials'],'title');
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    times = s.times;
    freqs = s.freqs;
    
