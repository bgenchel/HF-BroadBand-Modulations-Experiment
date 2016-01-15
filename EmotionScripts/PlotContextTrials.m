% Plots single-trial context decompositions with hi and low averaged trials and % probability of context.
%
% PlotContextTrials(newttl, fullpath,sigtemps,erspmask,ctxsigs,ctxstd,shuffnum,tlim,elim,ctlim,clim,fontsz,ttl,savename,erspctx);
%
% newttl -- [string] data file title (will add 'stuff.mat' onto end to open file
% fullpath -- [sting] data path where data file can be found
% sigtemps -- vector of template indices to plot
% erspmask -- [templates X time/freq pnts] 1 when bootstrap distribution was sig, 0 otherwise.
% ctxsigs -- [templates x questions] matrix with 1's where context template is significant (by bootstrap)
% ctxstd -- if not empty, matrix of template X question with standard deviation of bootstrap values.
% shuffnum -- [integer] number of shuffle repetitions to perform in context probability
% tlim -- [number] color limit of ersp templates (default: 10)
% elim -- [number] color limit of average ersp plots (default: 5)
% ctlim -- [number] y limit of context template plots (default: .7)
% clim -- [number] y limit of % probability plots (default: .8)
% fontsz -- [number] fontsize for context labels
% ttl -- [string] title stem to add to figure top
% savename -- [string] name stem to save produced figures as. Will add .jpg extension. [] will not save.
%             Will save always as [savename,int2str(fignum),'.jpg'] to accommodate multiple figures
% erspctx -- [0 | 1 | 2] 0 if ERSP and context; 1 if ersp-only decomposition (only returns context 
%                        probability, not context templates); 2 if context only

function PlotContextTrials(newttl, fullpath,sigtemps,erspmask,ctxsigs,ctxstd,shuffnum,tlim,elim,ctlim,clim,fontsz,ttl,savename,erspctx);
    
    percplot = .1; % percent of trials to take for ersps
    gry = .4; % grey level of context labels
    fignum = 0;
    if isempty(fontsz)
        tlim = 9;
    end;
    if isempty(ttl)
        ttl = 'Component';
    end;
    
    if erspctx == 0
    %%%%%%%%%%%%%%%%%Plot hi/low weighted trials  %%%%%%%%%%%
    a = load([fullpath,'FBContextTrialMatrix.mat']);
    s = load([fullpath,newttl,'stuff.mat']);
    data= floatread([fullpath,newttl,'.fdt'],[s.channo s.frameno],[],0);
    sph = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.sph'],[s.channo s.channo],[],0);
    wts = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.wts'],[s.pcs s.channo],[],0);
    erspeig = floatread([fullpath,newttl,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcerspnum],[],0);
    erspdata = floatread([fullpath,newttl,'DAT.fdt'],[length(s.freqs)*length(s.times) inf],[],0);

    ws = wts*sph; winv = pinv(ws); activations = ws*data;
    nlett = 6; nfb = 21; nresp = 6;   % numbers of questions in the 3 categories
    ms = 3; starsize = 10;  
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%
    if isempty(tlim)
        tlim = max(max(abs(winv(1:s.pcerspnum,:)'*erspeig')))- max(max(abs(winv(1:s.pcerspnum,:)'*erspeig')))*.2;
    end;
    if isempty(elim)
        elim = max(max(abs(erspdata)))-max(max(abs(erspdata)))*.6;
    end;
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ctlim) % context template
        ctlim = max(max(abs(winv(s.pcerspnum+1:end,:))));
    end;
    repeat = 1;
    col = 6;  
    %%%%%%%%%%%%%%%%%
    % make signif masks for ersp trial plotting
    randvec = [1:size(erspdata,2)]; 
    for rep = 1:shuffnum
        picksome = shuffle(randvec); picksome = picksome(1:round(size(erspdata,2)*percplot));
        randersps(:,rep) = mean(erspdata(:,picksome),2);
    end;
    for tf = 1:size(randersps,1)
        currvec = randersps(tf,:);currvec = sort(currvec);
        minmask(tf,:) = currvec(round(shuffnum*.01));
        maxmask(tf,:) = currvec(length(currvec) - round(shuffnum*.01));
    end;
    %%%%%%%%%%%%%%%%%
    while repeat == 1
        figure; row = 6;pl = 1;
        if length(sigtemps) < 3
            row = 4;
        end;    
        for tmp = 1:row/2
            if tmp <= length(sigtemps)
                template = sigtemps(tmp);
                sortact = sort(activations(template,:));
                cutvalhi = sortact(round(length(sortact)*(1-percplot)));
                hivals = find(activations(template,:) > cutvalhi);
                cutvallo = sortact(round(length(sortact)*percplot));
                lovals = find(activations(template,:) < cutvallo);
                %plotersphi = mean(erspdata(1:length(s.freqs)*length(s.times),hivals),2);
                %plotersphi(find(plotersphi>minmask&plotersphi<maxmask)) = 0;
                %plotersphi = reshape(plotersphi,length(s.freqs),length(s.times));
                %plotersplo = mean(erspdata(1:length(s.freqs)*length(s.times),lovals),2);
                %plotersplo(find(plotersplo>minmask&plotersplo<maxmask)) = 0;
                %plotersplo = reshape(plotersplo,length(s.freqs),length(s.times));
                contxthi = mean(a.addmat(:,hivals),2);
                contxtlo = mean(a.addmat(:,lovals),2);
                [tsampperc,tbootperc,tallboots,hH] = CntxtBoot(a.addmat,a.addmat(:,hivals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperchi(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                [tsampperc,tbootperc,tallboots,lH] = CntxtBoot(a.addmat,a.addmat(:,lovals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperclo(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
               if isempty(clim)
                      clim = max([max(max(abs(tplotperchi))) max(max(abs(tplotperclo)))]);
                    starlev = clim*.25;
                else
                    starlev = clim*.25;   
                end;
                %%%%%  Plot Positive  ERSP Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                oneersp = winv (1:s.pcerspnum,template)'*erspeig';
                if ~isempty(erspmask)
                    oneersp(find(erspmask(template,:) == 0)) = 0;            
                end;
                imagesc(s.times,s.freqs,reshape(oneersp,length(s.freqs),length(s.times)),[-tlim tlim]); 
                hold on; set(gca,'ydir','norm');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]);     set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');title(['Template ',int2str(template)]);
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'xticklabel',[]);
                %%%%%  Plot Positive  ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                % multiply act row by 'data' matrix, then pass through eigvec
                epc = data(1:s.pcerspnum,:);
                for p = 1:size(epc,1)
                    epc(p,:) = epc(p,:).*activations(tmp,:);
                end;
                oneersp = erspeig*epc;
                oneersp = mean(oneersp,2);
                oneersp = reshape(oneersp,length(s.freqs),length(s.times))
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,oneersp,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title('hi 10% trials');
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        

                %%%%%  Plot Positive Context Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                allctxvals = winv(s.pcerspnum+1:end,template);
                ph = plot([1:nlett],winv(s.pcerspnum+1:s.pcerspnum+nlett,template),'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],winv(s.pcerspnum+nlett+1:s.pcerspnum+nlett+nfb,template),'bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(s.pcerspnum+nlett+nfb+1:s.pcerspnum+nlett+nfb+nresp,template),'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+2;hold on;
                set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                if ~isempty(ctxstd)
                    if template <= size(ctxstd,1)
                    for qq = 1:size(ctxstd,2)
                        ph = plot([qq qq],[allctxvals(qq,1)-ctxstd(template,qq) allctxvals(qq,1)+ctxstd(template,qq)],'k-');
                        set(ph,'color',[1 .5 0]);
                    end;
                    end;
                end; 
                if ~isempty(ctxsigs)
                    sigpnts = find(ctxsigs(template,:));
                    if ~isempty(sigpnts)
                        for hh = 1:length(sigpnts)
                            ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                            set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                        end;   
                    end;
                end;
                if  tmp == 2
                    CnxtLabelRewFB(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                %set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                if tmp == 1
                    title(['Context Template']);
                end;
                %%%%%  Plot Positive Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperchi(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperchi(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperchi(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+2; hold on; set(ph,'linewidth',2);
                if tmp == 1
                    title(['Context Probability']);
                end;
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                sigpnts = find(hH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;   
                end;
                 if  tmp == 2 
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
               %set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                
                %%%%%  Plot Negative ERSP Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,reshape(oneersp,length(s.freqs),length(s.times))*-1,[-tlim tlim]); 
                hold on; set(gca,'ydir','norm');    set(gca,'yscale','log'); pl = pl+1;    
                set(gca,'ticklength',[.02 .02]);    set(gca,'ytick',[5:5:s.freqs(end)]); 
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');title(['Template -',int2str(template)]);
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        
                
                %%%%%  Plot Negative ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,plotersplo,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title('low 10% trials');
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        
                
                %%%%%  Plot Negative Context Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                ph = plot([1:nlett],winv(s.pcerspnum+1:s.pcerspnum+nlett,template)*-1,'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],winv(s.pcerspnum+nlett+1:s.pcerspnum+nlett+nfb,template)*-1,'bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(s.pcerspnum+nlett+nfb+1:s.pcerspnum+nlett+nfb+nresp,template)*-1,'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+2;hold on;
                set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                set(ph,'markersize',ms);

                if  tmp == row/2 | tmp == row/4
                    CnxtLabelRewFB(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                if ~isempty(ctxsigs)
                    sigpnts = find(ctxsigs(template,:));
                    if ~isempty(sigpnts)
                        for hh = 1:length(sigpnts)
                            ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                            set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                        end;   
                    end;
                end;
                if ~isempty(ctxstd)
                    if template <= size(ctxstd,1)
                    allctxvals = winv(s.pcerspnum+1:end,template)*-1;
                    for qq = 1:size(ctxstd,2)
                        ph = plot([qq qq],[allctxvals(qq,1)-ctxstd(template,qq) allctxvals(qq,1)+ctxstd(template,qq)],'k-');
                        set(ph,'color',[1 .5 0]);
                    end;
                    end;
                end;            
                if tmp ~=row/2
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end;        

                %%%%%  Plot Negative Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperclo(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperclo(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperclo(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+2; hold on; set(ph,'linewidth',2);
                
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                if  tmp == row/2 | tmp == row/4
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                sigpnts = find(lH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;
                end;
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end;        
        end;% to if tmp <=length(sigtemps)
    end; % to for tmp loop
    set(gcf,'PaperOrientation','landscape');        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    textsc([ttl,' Decomp of ERSP/Context; Hi and low weighted trials; ERSPs and percent (probability) of context questions'],'title');  
    if tmp <= length(sigtemps)
        sigtemps(1:tmp) = [];
    else
        sigtemps = [];
    end;    
    if isempty(sigtemps)
        repeat = 0; 
    end;    
    fignum = fignum+1;
    if ~isempty(savename)
        makesave = ['print ',fullpath,savename,int2str(fignum),'.eps -depsc'];    eval(makesave);   
    end;
    end; % to while loop
    
    elseif erspctx == 1  % ersp only decomp
    %%%%%%%%%%%%%%%%%Plot hi/low weighted trials  %%%%%%%%%%%
    a = load([fullpath,'FBContextTrialMatrix.mat']);
    s = load([fullpath,newttl,'stuff.mat']);
    data= floatread([fullpath,newttl,'.fdt'],[s.channo s.frameno],[],0);
    sph = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.sph'],[s.channo s.channo],[],0);
    wts = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.wts'],[s.pcs s.channo],[],0);
    erspeig = floatread([fullpath,newttl,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcerspnum],[],0);
    erspdata = floatread([fullpath,newttl,'DAT.fdt'],[length(s.freqs)*length(s.times) inf],[],0);
    ws = wts*sph; winv = pinv(ws); activations = ws*data;
    nlett = 6; nfb = 21; nresp = 6;   % numbers of questions in the 3 categories
    ms = 3; starsize = 10;  
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%
    if isempty(tlim)
        tlim = max(max(abs(winv'*erspeig')))-max(max(abs(winv'*erspeig')))*.5;
    end;
    if isempty(elim)
        elim = max(max(abs(erspdata)))-max(max(abs(erspdata)))*.9;
    end;
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ctlim) % context template
        ctlim = max(max(abs(winv)));
    end;
    repeat = 1;
    col = 5;  
    %%%%%%%%%%%%%%%%%
    % make signif masks for ersp trial plotting
    randvec = [1:size(erspdata,2)]; 
    for rep = 1:shuffnum
        picksome = shuffle(randvec); picksome = picksome(1:round(size(erspdata,2)*percplot));
        randersps(:,rep) = mean(erspdata(:,picksome),2);
    end;
    for tf = 1:size(randersps,1)
        currvec = randersps(tf,:);currvec = sort(currvec);
        minmask(tf,:) = currvec(round(shuffnum*.01));
        maxmask(tf,:) = currvec(length(currvec) - round(shuffnum*.01));
    end;
    %%%%%%%%%%%%%%%%%
    while repeat == 1
        figure; row = 6;pl = 1;
        if length(sigtemps) < 3
            row = 4;
        end;    
        for tmp = 1:row/2
            if tmp <= length(sigtemps)
                template = sigtemps(tmp);
                sortact = sort(activations(template,:));
                cutvalhi = sortact(round(length(sortact)*(1-percplot)));
                hivals = find(activations(template,:) > cutvalhi);
                cutvallo = sortact(round(length(sortact)*percplot));
                lovals = find(activations(template,:) < cutvallo);
                plotersphi = mean(erspdata(1:length(s.freqs)*length(s.times),hivals),2);
                plotersphi(find(plotersphi>minmask&plotersphi<maxmask)) = 0;
                plotersphi = reshape(plotersphi,length(s.freqs),length(s.times));
                plotersplo = mean(erspdata(1:length(s.freqs)*length(s.times),lovals),2);
                plotersplo(find(plotersplo>minmask&plotersplo<maxmask)) = 0;
                plotersplo = reshape(plotersplo,length(s.freqs),length(s.times));
                contxthi = mean(a.addmat(:,hivals),2);
                contxtlo = mean(a.addmat(:,lovals),2);
                [tsampperc,tbootperc,tallboots,hH] = CntxtBoot(a.addmat,a.addmat(:,hivals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperchi(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                [tsampperc,tbootperc,tallboots,lH] = CntxtBoot(a.addmat,a.addmat(:,lovals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperclo(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                if isempty(clim)
                     clim = max([max(max(abs(tplotperchi))) max(max(abs(tplotperclo)))]);
                    starlev = clim*.25;
                else
                    starlev = clim*.25;
                end;
                %%%%%  Plot Positive  ERSP Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                oneersp = winv (:,template)'*erspeig';
                if ~isempty(erspmask)
                    oneersp(find(erspmask(template,:) == 0)) = 0;            
                end;
                imagesc(s.times,s.freqs,reshape(oneersp,length(s.freqs),length(s.times)),[-tlim tlim]); 
                hold on; set(gca,'ydir','norm');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]);     set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');title(['Template ',int2str(template)]);
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'xticklabel',[]);
                %%%%%  Plot Positive  ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,plotersphi,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title('hi 10% trials');
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        

                %%%%%  Plot Positive Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+2])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperchi(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperchi(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperchi(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+3; hold on; set(ph,'linewidth',2);
                if tmp == 1
                    title(['Context Probability']);
                end;
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                sigpnts = find(hH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;   
                end;
                if  tmp == 2 
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                %set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                
                %%%%%  Plot Negative ERSP Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,reshape(oneersp,length(s.freqs),length(s.times))*-1,[-tlim tlim]); 
                hold on; set(gca,'ydir','norm');    set(gca,'yscale','log'); pl = pl+1;    
                set(gca,'ticklength',[.02 .02]);    set(gca,'ytick',[5:5:s.freqs(end)]); 
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');title(['Template -',int2str(template)]);
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        
                
                %%%%%  Plot Negative ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,plotersplo,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title('low 10% trials');
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        
                
 
                %%%%%  Plot Negative Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+2])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperclo(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperclo(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperclo(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+3; hold on; set(ph,'linewidth',2);
                
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                if  tmp == row/2 | tmp == row/4
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                sigpnts = find(lH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;
                end;
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end;        
        end;% to if tmp <=length(sigtemps)
        end; % to for tmp loop
        set(gcf,'color','w');
    set(gcf,'PaperOrientation','landscape');        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    textsc([ttl,' Decomp of ERSP ONLY; Hi and low weighted trials; ERSPs and percent (probability) of context questions'],'title');  
    if tmp <= length(sigtemps)
        sigtemps(1:tmp) = [];
    else
        sigtemps = [];
    end;    
    if isempty(sigtemps)
        repeat = 0; 
    end;    
    fignum = fignum+1;
    if ~isempty(savename)
        makesave = ['print ',fullpath,savename,int2str(fignum),'.eps -depsc'];    eval(makesave);   
    end;
    end; % to while loop
    
    elseif erspctx == 2   % context only
        
        %%%%%%%%%%%%%%%%%Plot hi/low weighted trials  %%%%%%%%%%%
    a = load([fullpath,'FBContextTrialMatrix.mat']);
    s = load([fullpath,newttl,'stuff.mat']);
    data= a.addmat;
    sph = floatread([fullpath,newttl,'SPH.sph'],[s.channo s.channo],[],0);
    wts = floatread([fullpath,newttl,'WTS.wts'],[s.pcs s.channo],[],0);
    erspdata = floatread([fullpath,newttl,'DAT.fdt'],[length(s.freqs)*length(s.times) inf],[],0);
    ws = wts*sph; winv = pinv(ws); activations = ws*data;
    nlett = 6; nfb = 21; nresp = 6;   % numbers of questions in the 3 categories
    ms = 3; starsize = 10;  
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%
    if isempty(elim)
        elim = max(max(abs(erspdata)))-max(max(abs(erspdata)))*.9;
    end;
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ctlim) % context template
        ctlim = max(max(abs(winv)));
    end;
    repeat = 1;
    col = 5; 
    %%%%%%%%%%%%%%%%%
    % make signif masks for ersp trial plotting
    randvec = [1:size(erspdata,2)]; 
    for rep = 1:shuffnum
        picksome = shuffle(randvec); picksome = picksome(1:round(size(erspdata,2)*percplot));
        randersps(:,rep) = mean(erspdata(:,picksome),2);
    end;
    for tf = 1:size(randersps,1)
        currvec = randersps(tf,:);currvec = sort(currvec);
        minmask(tf,:) = currvec(round(shuffnum*.01));
        maxmask(tf,:) = currvec(length(currvec) - round(shuffnum*.01));
    end;
    
    while repeat == 1
        figure; row = 6;pl = 1;
        if length(sigtemps) < 3
            row = 4;
        end;    
        for tmp = 1:row/2
            if tmp <= length(sigtemps)
                template = sigtemps(tmp);
                sortact = sort(activations(template,:));
                cutvalhi = sortact(round(length(sortact)*(1-percplot)));
                hivals = find(activations(template,:) > cutvalhi);
                cutvallo = sortact(round(length(sortact)*percplot));
                lovals = find(activations(template,:) < cutvallo);
                plotersphi = mean(erspdata(1:length(s.freqs)*length(s.times),hivals),2);
                plotersphi(find(plotersphi>minmask&plotersphi<maxmask)) = 0;
                plotersphi = reshape(plotersphi,length(s.freqs),length(s.times));
                plotersplo = mean(erspdata(1:length(s.freqs)*length(s.times),lovals),2);
                plotersplo(find(plotersplo>minmask&plotersplo<maxmask)) = 0;
                plotersplo = reshape(plotersplo,length(s.freqs),length(s.times));
                contxthi = mean(a.addmat(:,hivals),2);
                contxtlo = mean(a.addmat(:,lovals),2);
                [tsampperc,tbootperc,tallboots,hH] = CntxtBoot(a.addmat,a.addmat(:,hivals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperchi(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                [tsampperc,tbootperc,tallboots,lH] = CntxtBoot(a.addmat,a.addmat(:,lovals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperclo(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                if isempty(clim)
                    clim = max([max(max(abs(tplotperchi))) max(max(abs(tplotperclo)))]);
                    starlev = clim*.25;
                else
                    starlev = clim*.25;   
                end;
                %%%%%  Plot Positive  ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,plotersphi,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title(['Fac ',int2str(template),': hi 10% trials']);
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        

                %%%%%  Plot Positive Context Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                allctxvals = winv(:,template);
                ph = plot([1:nlett],winv(1:nlett,template),'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],winv(nlett+1:nlett+nfb,template),'bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(nlett+nfb+1:nlett+nfb+nresp,template),'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+2;hold on;
                set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                if ~isempty(ctxstd)
                    for qq = 1:size(ctxstd,2)
                        ph = plot([qq qq],[allctxvals(qq,1)-ctxstd(template,qq) allctxvals(qq,1)+ctxstd(template,qq)],'k-');
                        set(ph,'color',[1 .5 0]);
                    end;
                end; 
                if ~isempty(ctxsigs)
                    sigpnts = find(ctxsigs(template,:));
                    if ~isempty(sigpnts)
                        for hh = 1:length(sigpnts)
                            ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                            set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                        end;   
                    end;
                end;
                %set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                 if  tmp == 2 
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                if tmp == 1
                    title(['Context Template']);
                end;
                %%%%%  Plot Positive Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperchi(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperchi(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperchi(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+2; hold on; set(ph,'linewidth',2);
                if tmp == 1
                    title(['Context Probability']);
                end;
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                sigpnts = find(hH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;   
                end;
                 if  tmp == 2 
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
               %set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                
                
                %%%%%  Plot Negative ERSP Average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,pl)
                set(gca,'fontsize',fontsz);
                imagesc(s.times,s.freqs,plotersplo,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
                plot([-617 -617],[get(gca,'ylim')],'k-.');
                plot([882 882],[get(gca,'ylim')],'k-.');
                set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
                set(gca,'ytick',[5:5:s.freqs(end)]); title(['Fac ',int2str(template),': low 10% trials']);
                set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                end;        
                
                %%%%%  Plot Negative Context Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                ph = plot([1:nlett],winv(1:nlett,template)*-1,'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],winv(nlett+1:nlett+nfb,template)*-1,'bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(nlett+nfb+1:nlett+nfb+nresp,template)*-1,'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+2;hold on;
                set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                set(ph,'markersize',ms);

                if  tmp == row/2
                    CnxtLabelRewFB(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                if ~isempty(ctxsigs)
                    sigpnts = find(ctxsigs(template,:));
                    if ~isempty(sigpnts)
                        for hh = 1:length(sigpnts)
                            ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                            set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                        end;   
                    end;
                end;
                if ~isempty(ctxstd)
                    allctxvals = winv(:,template)*-1;
                    for qq = 1:size(ctxstd,2)
                        ph = plot([qq qq],[allctxvals(qq,1)-ctxstd(template,qq) allctxvals(qq,1)+ctxstd(template,qq)],'k-');
                        set(ph,'color',[1 .5 0]);
                    end;
                end;            
                if tmp ~=row/2
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end;        

                %%%%%  Plot Negative Context Probability  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+1])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperclo(1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperclo(nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperclo(nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+2; hold on; set(ph,'linewidth',2);
                
                set(gca,'ylim',[-clim clim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                if  tmp == row/2
                    CnxtLabelRewFB(clim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                sigpnts = find(lH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[clim-starlev clim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;
                end;
                if tmp ~= row/2
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end;        
        end;% to if tmp <=length(sigtemps)
    end; % to for tmp loop
    set(gcf,'PaperOrientation','landscape');        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    textsc([ttl,' Hi and low weighted trials; Ctx templates and percent (probability)'],'title');  
    if tmp <= length(sigtemps)
        sigtemps(1:tmp) = [];
    else
        sigtemps = [];
    end;    
    if isempty(sigtemps)
        repeat = 0; 
    end;    
    fignum = fignum+1;
    if ~isempty(savename)
        makesave = ['print ',fullpath,savename,int2str(fignum),'.eps -depsc'];    eval(makesave);   
    end;
    end; % to while loop

    end; % to if erspctx 0,1,2
