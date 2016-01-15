% plots summaries of component factors clustered by correlation between context vectors 
% 
% PlotShapeCtxSum(mapset,fullpath,gdcomps,clustidx,erspmask,shuffnum,fontsz,sigmask,ttl,posneg,scale,tmplorprob)
%
% mapset -- [string] full dataset name to load for scalp map image (requires eeglab open)
% fullpath -- [string] full path of directory in which 'mapset' can be found
% gdcomps -- [vector] indices of components included in original decomposition
% clustidx -- [? X 2 matrix] first column indicates component index, second column = factor number
% addmat -- [Q's X trials] matrix of context indicators used in original decomposition
% filename -- [string] data filename stem to call in raw data for plotting
% sigmask -- [0 | 1] 1 to mask ersp by bootstrap
% ttl -- [string] title to be written at the top of the figure
% scale -- ['n' | 'y'] 0 will NOT scale context probability by number of positive events
% tmplorprob -- ['prob' or 'tmpl'] will plot probability of context or template
% erspmask -- [cell array] containing masking info (1's and 0's) for all templates
% ctxsigs -- [cell array] same as erspmask but for context 
% ctxstd -- [cell array] contains std from bootstrap to plot on context templates.  

function PlotShapeCtxSum(mapset,fullpath,gdcomps,clustidx,shuffnum,fontsz,sigmask,erspmask,ctxsigs,ctxstd,ttl,scale,tmplorprob)
    
    ctlim = [];
    percplot = .1; % percent of trials to take for ersps
    fignum = 0;
    if isempty(fontsz)
        fontsz = 9;
    end;
    if isempty(ttl)
        ttl = 'Cluster Summary';
    end;
    EEG = pop_loadset( mapset,fullpath);    
    a = load([fullpath,'FBContextTrialMatrix.mat']);
    % Calculate mult factor for context vector
    for ad = 1:size(a.addmat,1)
        percqs(1,ad) = 1-(length(find(a.addmat(ad,:) == 1))/length(a.addmat(ad,:)));
        if scale == 'n'
            percqs(1,ad) = 1;
        end;        
    end;
    %%%%%%%%%%%%%%%%%Plot hi/low weighted trials  %%%%%%%%%%%
    nlett = 9; nfb = 19; nresp = 6;   % numbers of questions in the 3 categories
    ms = 5; starsize = 10; 
    col = 7;  gry = .25;
    figure; pl = 1;
    row = size(clustidx,1);
    %%%%%%%%%%%%%%%%%
        for cmp = 1:size(clustidx,1)
            tlim = [];elim = []; 
            fprintf('.');
            comp = clustidx(cmp,1);        
            template = clustidx(cmp,2); 
            newttl = ['Cp',int2str(gdcomps(comp)),'CXTAllDat'];
            s = load([fullpath,newttl,'stuff.mat']);
            data= floatread([fullpath,newttl,'.fdt'],[s.channo s.frameno],[],0);
            sph = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.sph'],[s.channo s.channo],[],0);
            wts = floatread([fullpath,newttl,int2str(s.channo),'pc',int2str(s.pcs),'.wts'],[s.pcs s.channo],[],0);
            erspeig = floatread([fullpath,newttl,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcerspnum],[],0);
            %erspdata = floatread([fullpath,newttl,'DAT.fdt'],[length(s.freqs)*length(s.times) inf],[],0);
            ws = wts*sph; winv = pinv(ws); activations = ws*data;
            %backproj = winv(1:33,template)*activations(template,:);
            %backproj = erspeig*backproj;
            %%%%%%%%%%%%%%%%% %%%%%%%%%%%
            % run pca to reduce dims and backproj
            %[pc,eigvec,sv] = runpca(data,s.pcs);
            %newdata = eigvec*pc;
            %newdata = erspeig*newdata(1:s.pcerspnum,:);
            %%%%%%%%%%%%%%%%% %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% plot component scalp map
            sbplot(row,col,pl)
            set(gca,'fontsize',fontsz);
            topoplot(EEG.icawinv(:,gdcomps(comp)),EEG.chanlocs);pl = pl+1;
            title(['Comp ',int2str(gdcomps(comp))]);
            
            %% Plot sorted weights %%%%%%%%%%%%%%% 
            sbplot(row,col,pl) 
            set(gca,'fontsize',fontsz);
            hist(activations(template,:),100);hold on;pl = pl +1;
            %plot(sort(activations(template,:)),'linewidth',2);pl = pl +1;
            %set(gca,'xlim',[1 size(activations,2)]);hold on;
            plot([0 0],[get(gca,'ylim')],'g-');
            if cmp == 1
                title('Weights Distribution');
            end;
            %%%%%%%%%%%%%%%%%
            sortact = sort(activations(template,:));
            cutvalhi = sortact(round(length(sortact)*(1-percplot)));
            %hivals = find(activations(template,:) > cutvalhi);
            %plotersphi = mean(newdata(1:length(s.freqs)*length(s.times),hivals),2);
            %plotersphi = mean(erspdata(1:length(s.freqs)*length(s.times),hivals),2);
            %plotersphi = mean(backproj(1:length(s.freqs)*length(s.times),hivals),2);
            %if sigmask == 1
            %    plotersphi(find(plotersphi>minmask&plotersphi<maxmask)) = 0;
            %end;
            %plotersphi = reshape(plotersphi,length(s.freqs),length(s.times));
           
            %%%%%%%%%%%%%%%%%%%% factor ersp template
            sbplot(row,col,pl)
            set(gca,'fontsize',fontsz);
            oneersp = erspeig*winv (1:s.pcerspnum,template);
            if ~isempty(erspmask)
                if ~isempty(erspmask{comp})
                    oneersp(find(erspmask{comp}(template,:) == 0)) = 0;
                end;
            end;
            tlim = max(abs(oneersp));
            if tlim == 0
                tlim = .01;
            end;                
            imagesc(s.times,s.freqs,reshape(oneersp,length(s.freqs),length(s.times)),[-tlim tlim]); 
            hold on; set(gca,'ydir','norm');
            set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
            set(gca,'ytick',[5:5:s.freqs(end)]);     set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            hold on; plot([0 0],[get(gca,'ylim')],'k-');title(['Fac ',int2str(template)]);        
            if cmp == length(clustidx)
                cbar;
            end;        
            %%%%%%%%%%%%%%%%%%%% mean ersp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make signif masks for ersp trial plotting
            epc = data(1:s.pcerspnum,:);
            randacts = activations(template,:);
            randersps = zeros(shuffnum,size(erspeig,1));
            for rep = 1:shuffnum
                picksome = shuffle(randacts); 
                picksome = repmat(picksome,[size(epc,1) 1]);
                epc2 = epc.*picksome;
                %for p = 1:size(epc,1)
                    %epc2(p,:) = epc(p,:).*picksome;
                %end;
                oneersp = erspeig*epc2;
                oneersp = mean(oneersp,2);
                randersps(rep,:) = oneersp';
            end;
            for tf = 1:size(randersps,2)
                currvec = randersps(:,tf);currvec = sort(currvec);
                minmask(tf,:) = currvec(round(shuffnum*.01));
                maxmask(tf,:) = currvec(length(currvec) - round(shuffnum*.01));
            end;
            %%%%%%%%%%%%%%%
            sbplot(row,col,pl)
            set(gca,'fontsize',fontsz);
            % multiply act row by 'data' matrix, then pass through eigvec
            for p = 1:size(epc,1)
                epc(p,:) = epc(p,:).*activations(template,:);
            end;
            oneersp = erspeig*epc;
            oneersp = mean(oneersp,2);
            if sigmask == 1
                oneersp(find(oneersp>minmask&oneersp<maxmask)) = 0;
            end;
            oneersp = reshape(oneersp,length(s.freqs),length(s.times));
            set(gca,'fontsize',fontsz);
            elim = max(max(abs(oneersp)));
            if elim == 0
                elim = .01;
            end;            
            imagesc(s.times,s.freqs,oneersp,[-elim elim]); hold on; set(gca,'ydir','norm');hold on;
            plot([-617 -617],[get(gca,'ylim')],'k-.');        plot([882 882],[get(gca,'ylim')],'k-.');
            set(gca,'yscale','log'); pl = pl+1;    set(gca,'ticklength',[.02 .02]);
            set(gca,'ytick',[5:5:s.freqs(end)]); 
            if cmp == 1
                title(['Weighted Mean ERSP']);
            end;        
            set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            if cmp == length(clustidx)
                cbar;
            end;        
            
            %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if tmplorprob == 'prob'
                contxthi = mean(a.addmat(:,hivals),2);
                [tsampperc,tbootperc,tallboots,hH] = CntxtBoot(a.addmat,a.addmat(:,hivals),'t',shuffnum,.01); 
                for qs = 1:length(tsampperc)
                    tplotperchi(1,qs) = tsampperc(1,qs) - mean(tallboots(qs,:));
                end;
                if isempty(ctlim)
                    ctlim = max(max(abs(oneersp)));
                    starlev = ctlim*.25;
                else
                    starlev = ctlim*.25;   
                end;
                %%%%%%%%%%%%%%%%%%%% context probability 
                %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+2])
                set(gca,'fontsize',fontsz);
                ph = plot([1:nlett],tplotperchi(1:nlett).*percqs(1,1:nlett),'ro-'); hold on; set(gca,'fontsize',fontsz);set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],tplotperchi(nlett+1:nlett+nfb).*percqs(1,nlett+1:nlett+nfb),'bs-'); hold on; set(ph,'linewidth',2);
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],tplotperchi(nlett+nfb+1:end).*percqs(1,nlett+nfb+1:end),'g^-');
                set(ph,'markersize',ms+ms/2);
                pl=pl+3; hold on; set(ph,'linewidth',2);
                if cmp == 1
                    title(['Context Probability ']);
                end;
                if  cmp == row | cmp == row/2
                    CnxtLabelShape(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                set(gca,'ylim',[-ctlim ctlim]);set(gca,'xlim',[0 size(a.addmat,1)+1]);plot([get(gca,'xlim')],[0 0],'k-');
                sigpnts = find(hH);
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;   
                end;
                set(gca,'yticklabel',[]);
            elseif  tmplorprob == 'wtdm'
                %%%%%  Plot Weighted mean of all trials  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                
                epc = data(s.pcerspnum+1:end,:);  
                randacts = activations(template,:);
                randersps = zeros(shuffnum,size(epc,1));
                for rep = 1:shuffnum
                    picksome = shuffle(randacts); 
                    picksome = repmat(picksome,[size(epc,1) 1]);
                    epc2 = epc.*picksome;
                    oneersp = mean(epc2,2);
                    randersps(rep,:) = oneersp';
                end;
                siglev = .01;
                for tf = 1:size(randersps,2)
                    currvec = randersps(:,tf);currvec = sort(currvec);
                    cminmask(tf,:) = currvec(round(shuffnum*siglev));
                    cmaxmask(tf,:) = currvec(length(currvec) - round(shuffnum*siglev));
                end;
                %%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+2])
                set(gca,'fontsize',fontsz);
                % multiply act row by 'data' matrix, then pass through eigvec
                multmat = repmat(activations(template,:),[size(epc,1) 1]);
                epc = multmat.*epc;
                mnreps = mean(randersps,1);
                stdreps = std(randersps,1);
                twoersp = mean(epc,2)'; % for sig calc
                oneersp = mean(epc,2)'; % for zscore
                for tf = 1:length(oneersp) % don't include real weighting in zscore
                    oneersp(1,tf) = (oneersp(1,tf) - mnreps(tf))/stdreps(tf);% zscore
                end;                
                
                if isempty(ctlim)
                    ctlim = max(abs(oneersp));
                    starlev = ctlim*.25;
                else
                    starlev = ctlim*.25;   
                end;
                if sigmask == 1
                    sigpnts = find(twoersp < cminmask' | twoersp> cmaxmask');
                    grypnts = find(twoersp > cminmask' & twoersp < cmaxmask');
                end;
                ph = plot([1:nlett],oneersp(1,1:nlett),'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],oneersp(1,nlett+1:nlett+nfb)','bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],oneersp(1,nlett+nfb+1:end),'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+3;hold on;

                set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                if ~isempty(sigpnts)
                    for hh = 1:length(sigpnts)
                        ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                        set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                    end;   
                end;
                for hh = 1:length(grypnts)
                    ph = plot(grypnts(hh),oneersp(1,grypnts(hh)),'r.');
                    set(ph,'markersize',22);set(ph,'color',[.7 .7 .7]);
                end;   
                if   cmp == row | cmp == row/2
                    CnxtLabelShape(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);
                end;
                set(gca,'ylim',[-ctlim ctlim]);
                if cmp == 1
                    title(['Weighted Mean of Context (zscore); sig at p<',num2str(siglev)]);
                end;
            elseif  tmplorprob == 'tmpl'
                %%%%%  Plot Positive Context Template  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %%%%%%%%%%%%%%%
                sbplot(row,col,[pl pl+2])
                set(gca,'fontsize',fontsz);
                % multiply act row by 'data' matrix, then pass through eigvec
                allctxvals = winv(s.pcerspnum+1:end,template);
                if isempty(ctlim)
                    ctlim = max(max(abs(allctxvals)));
                    starlev = ctlim*.25;
                else
                    starlev = ctlim*.25;   
                end;
               ph = plot([1:nlett],allctxvals(1:nlett)'.*percqs(1,1:nlett),'ro-','linewidth',2);hold on;
                set(ph,'markersize',ms);
                ph = plot([nlett+1:nlett+nfb],allctxvals(nlett+1:nlett+nfb)'.*percqs(1,nlett+1:nlett+nfb),'bs-','linewidth',2); hold on;
                set(ph,'markersize',ms);                
                ph = plot([nlett+nfb+1:nlett+nfb+nresp],allctxvals(nlett+nfb+1:end)'.*percqs(1,nlett+nfb+1:nlett+nfb+nresp+1:end),'g^-','linewidth',2);
                set(ph,'markersize',ms+ms/2); pl = pl+3;hold on;
                
                 set(gca,'xlim',[0 nlett+nfb+nresp+1]);set(gca,'ylim',[-ctlim ctlim]);
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
                if ~isempty(ctxstd)
                    if template <= size(ctxstd{comp},1)
                    for qq = 1:size(ctxstd{comp},2)
                        ph = plot([qq qq],[allctxvals(qq,1)*percqs(1,qq)-ctxstd{comp}(template,qq) allctxvals(qq,1)*percqs(1,qq)+ctxstd{comp}(template,qq)],'k-');
                        set(ph,'color',[1 .5 0]);
                    end;
                    end;
                end; 
                if ~isempty(ctxsigs)
                    sigpnts = find(ctxsigs{comp}(template,:));
                    grypnts = [1:length(allctxvals)]; grypnts(sigpnts) = [];
                    if ~isempty(sigpnts)
                        for hh = 1:length(sigpnts)
                            ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-starlev ctlim-starlev],'r*');
                            set(ph,'markersize',starsize);set(ph,'color',[1 .5 0]);
                        end;   
                        for hh = 1:length(grypnts)
                            ph = plot(grypnts(hh),allctxvals(grypnts(hh),1),'r.');
                            set(ph,'markersize',22);set(ph,'color',[.7 .7 .7]);
                        end;   
                    end;
                end;
                if   cmp == row | cmp == row/2
                    CnxtLabelShape(ctlim,gry,fontsz);
                    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
                end;
                %set(gca,'xticklabel',[]);
                set(gca,'ylim',[-ctlim ctlim]);
                if cmp == 1
                    title(['Context Templates (ICA rep for signif)']);
                end;
            
            end;            
        end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
    textsc(ttl,'title');
    
    fprintf('\ndone. done.\n');
