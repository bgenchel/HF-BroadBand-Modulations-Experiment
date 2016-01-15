% Plot results from SnglTrialERSPDecomp.m with standard context conditions (Reward Twoback)
%
% PlotShapeCxtResults(pathname,filename,wtsmat,sphmat,savename,freqs,times,cxtlim,pcersp,epcanum,pccxt,cpcanum,pcstot,channo,ctxsigs,ctxstd,whichtmps);
%
% pathname: [string] full path where wts, sph can be found and savename will be saved
% datatitle -- [string] name used to save all float files
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% cxtlim: [decimal] Y-axis limit on context factor plots
% epcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for ERSPs
%          [] if no PCA
% cpcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for context matrix
%          [] if no PCA
% pcstot: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% scalefc -- [decimal] factor to scale ERSP data (post PCA) by in formula: 
%             alltrials = alltrials/max(max(alltrials))*scalefac;
% erspmask -- [templates X time/freq pnts] 1 when bootstrap distribution was sig, 0 otherwise.
% ctxsigs -- [templates x questions] matrix with 1's where context template is significant (by bootstrap)
% ctxstd -- if not empty, matrix of template X question with standard deviation of bootstrap values.
% whichtmps -- vector of indexes to templates of interest to plot (ie, [1:pcstot] for all)

function PlotShapeCxtResults(pathname, datatitle, savename,ctlim,erspmask,ctxsigs,ctxstd,whichtmps);
    
    s = load([pathname,datatitle,'stuff.mat']); 
    
    fs = 11;  % fontsize
    ms = 3;   % markersize 
    nlett = 9; nfb = 19; nresp = 6;   % numbers of questions in the 3 categories
    ntot = 34; % number of questions total
    %%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph = floatread([pathname,datatitle,int2str(s.channo),'pc',int2str(s.pcs),'.sph'],[s.channo s.channo],[],0);
    wts = floatread([pathname,datatitle,int2str(s.channo),'pc',int2str(s.pcs),'.wts'],[s.pcs s.channo],[],0);
    
    ws = wts*sph;winv = pinv(ws);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on ERSP NOT Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(s.pcerspnum) & isempty(s.pccntxnum) % ersp is pca's, but context IS NOT   
        erspeig = floatread([pathname,datatitle,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcerspnum],[],0);
        
        row = length(whichtmps); col = 6;
        figure; pl = 1;
        lim = max(max(winv(1:s.pcerspnum,whichtmps)'*erspeig'))-.5*max(max(winv(1:s.pcerspnum,whichtmps)'*erspeig'));
        for tmp = 1:length(whichtmps)
            template = whichtmps(tmp);
            tmpcomp =  winv (1:s.pcerspnum,template)'*erspeig';
            if ~isempty(erspmask)
            tmpcomp(find(erspmask(template,:) == 0)) = 0;            
            end;
            tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
            %lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            subplot(row,col,pl);
            imagesc(s.times,s.freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);set(gca,'yscale','log');
            set(gca,'ytick',[5:5:s.freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            hold on; plot([-617 -617],[get(gca,'ylim')],'m-.');
            hold on; plot([882 882],[get(gca,'ylim')],'m-.');
            set(gca,'fontsize',fs);  %title(['Cluster: ',int2str(template)]); %colorbar;
            if tmp ~= length(whichtmps)
                set(gca,'xticklabel',[]);
            end;
            
            subplot(row,col,[pl:pl+4])
            allctxvals = winv(s.pcerspnum+1:end,template);
            ph = plot([1:nlett],winv(s.pcerspnum+1:s.pcerspnum+nlett,template),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],winv(s.pcerspnum+nlett+1:s.pcerspnum+nlett+nfb,template),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(s.pcerspnum+nlett+nfb+1:s.pcerspnum+nlett+nfb+nresp,template),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);            
            set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
            if ~isempty(ctxstd)
                for qq = 1:size(ctxstd,2)
                    ph = plot([qq qq],[allctxvals(qq,1)-ctxstd(template,qq) allctxvals(qq,1)+ctxstd(template,qq)],'k-');
                    set(ph,'color',[1 .5 0]);
                end;
            end;            
            if tmp == length(whichtmps) | tmp == length(whichtmps)/2
                CnxtLabelShape(ctlim,.3,fs); end;  
            if ~isempty(ctxsigs)    
            sigpnts = find(ctxsigs(template,:));
            if ~isempty(sigpnts)
                for hh = 1:length(sigpnts)
                    ph = plot([sigpnts(hh) sigpnts(hh)],[ctlim-.15 ctlim-.15],'r*');
                    set(ph,'markersize',15);set(ph,'color',[1 .5 0]);
                end;   
            end;
            end;
        end;
        ph = textsc(['ERSP/Context decomposition: PCA of ERSP to',int2str(s.pcerspnum),' dims; '],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = ['print ',pathname,savename,'.eps -depsc'];    eval(makesave);   
        end;       
        axcopy
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on Context NOT ERSP %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isempty(pcersp) & ~isempty(pccxt) % ersp not pca'd, but context is   
        addeig = floatread([pathname,datatitle,'ADDEIG.fdt'],[(nlett+nfb+nresp) s.pccntxnum]);        

        row = pcstot; col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:s.pcs
            tmpcomp = winv (1:length(s.freqs)*length(s.times),template)';  %makes a freqs*times X 1
            tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
            
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            sbplot(row,col,pl);
            imagesc(s.times,s.freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
            set(gca,'ytick',[5:5:s.freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;
            clear newwinv
            newwinv = winv(s.pcerspnum+1:end,template)'*addeig;
            sbplot(row,col,pl:pl+4)
            ph = plot([1:nlett],newwinv(s.pcerspnum+1:epcanum+nlett),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],newwinv(s.pcerspnum+nlett+1:s.pcerspnum+nlett+nfb),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],newwinv(s.pcerspnum+nlett+nfb+1:s.pcerspnum+nlett+nfb+nresp),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template == 10
            CnxtLabelShape(ctlim,.3,fs); end;
        end;
        ph = textsc(['ERSP/Context decomposition: PCA of context NOT ERSP matrix before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = ['print ',pathname,savename,'.eps -depsc'];    eval(makesave);   
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on  ERSP AND Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ~isempty(pcersp) & ~isempty(pccxt)% Both ersp and context are pca'd
        addeig = floatread([pathname,datatitle,'ADDEIG.fdt'],[(nlett+nfb+nresp) s.pccntxnum]);        
        erspeig = floatread([pathname,datatitle,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcerspnum],[],0);
        
        row = pcstot; col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:s.pcs
            tmpcomp =  winv (1:s.pcerspnum,template)'*erspeig';
            tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            sbplot(row,col,pl);
            imagesc(s.times,s.freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
            set(gca,'ytick',[5:5:s.freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;

            clear newwinv
            for fac = 1:size(pccxt,2)
                newwinv = winv(s.pcerspnum+1:end,template)'*addeig';
            end;
            sbplot(row,col,pl:pl+4)
            ph = plot([1:nlett],newwinv(s.pcerspnum+1:s.pcerspnum+nlett),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],newwinv(s.pcerspnum+nlett+1:s.pcerspnum+nlett+nfb),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],newwinv(s.pcerspnum+nlett+nfb+1:s.pcerspnum+nlett+nfb+nresp),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template == 10
            CnxtLabelShape(ctlim,.3,fs); end;

        end;
        ph = textsc(['ERSP/Context decomposition: PCA of ERSP AND context matrices before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = ['print ',pathname,savename,'.eps -depsc'];    eval(makesave);   
        end;

    end;
