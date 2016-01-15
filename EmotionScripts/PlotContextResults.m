% Plot results from SnglTrialERSPDecomp.m with standard context conditions (Reward Twoback)
%
% PlotContextResults(pathname,wtsmat,sphmat,savename,freqs,times,cxtlim,pcersp,epcanum,pccxt,cpcanum,pcstot,channo);
%
% pathname: [string] full path where wts, sph can be found and savename will be saved
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% cxtlim: [decimal] Y-axis limit on context factor plots
% pcersp: [string] name of float file containing eigenvectors from PCA decomposition of ERSPs. Use same
%         string as used for SnglTrialERSPDecomp.         [] if no PCA on ERSPs. 
% epcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for ERSPs
%          [] if no PCA
% pccxt: [string] name of float file containing eigenvectors from PCA decomposition of context matrix.
%         Use same string as used for SnglTrialERSPDecomp.     [] if no PCA on context matrix
% cpcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for context matrix
%          [] if no PCA
% pcstot: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% scalefc -- [decimal] factor to scale ERSP data (post PCA) by in formula: 
%             alltrials = alltrials/max(max(alltrials))*scalefac;

function PlotContextResults(pathname,wtsmat,sphmat,savename,freqs,times,ctlim,pcersp,epcanum,pccxt,cpcanum,pcstot,channo,scalefac);
    fs = 11; ms = 5; nlett = 9; nfb = 19; nresp = 4;  ntot = 32;
    %%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph = floatread([pathname,sphmat],[channo channo],[],0);
    wts = floatread([pathname,wtsmat],[pcstot channo],[],0);
    
    ws = wts*sph;winv = pinv(ws);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Plot Reg context factors onto one page %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(pcersp) & isempty(pccxt)
        row = pcstot; col = 6;
        figure;pl = 1; 
        for template = 1:pcstot
            tmpcomp = winv (1:length(freqs)*length(times),template)';  %makes a freqs*times X 1
            lim = max(abs(tmpcomp))+max(abs(tmpcomp))*.5;
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            sbplot(row,col,pl); pl = pl+1;
            imagesc(times,freqs,tmpcomp,[-lim lim]); set(gca,'yscale','log');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;
            sbplot(row,col,pl:pl+4)
            ph = plot([1:nlett],winv(epcanum+1:epcanum+nlett,template),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],winv(epcanum+nlett+1:epcanum+nlett+nfb,template),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(epcanum+nlett+nfb+1:epcanum+nlett+nfb+nresp,template),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template == 10
            CnxtLabelRew(ctlim,.3,fs); end;
        end;
        ph = textsc(['ERSP/Context decomposition: NO PCA of either ERSP or context matrices before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on ERSP NOT Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ~isempty(pcersp) & isempty(pccxt) % ersp is pca's, but context IS NOT   
        erspeig = floatread([pathname,pcersp,'.fdt'],[length(freqs)*length(times) epcanum],[],0);
        
        row = pcstot; col = 6;
        figure; pl = 1;
        lim = max(max(winv(1:epcanum,:)'*erspeig'))-.5*max(max(winv(1:epcanum,:)'*erspeig'));
        for template = 1:pcstot
            tmpcomp =  winv (1:epcanum,template)'*erspeig';
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            %lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            sbplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);set(gca,'yscale','log');
            set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;
            set(gca,'xticklabel',[]);
            subplot(row,col,[pl:pl+4])
            ph = plot([1:nlett],winv(epcanum+1:epcanum+nlett,template),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],winv(epcanum+nlett+1:epcanum+nlett+nfb,template),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],winv(epcanum+nlett+nfb+1:epcanum+nlett+nfb+nresp,template),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
            if template == 10
            CnxtLabelRew(ctlim,.3,fs); end;
         end;
        ph = textsc(['ERSP/Context decomposition: PCA of ERSP (',int2str(epcanum),' dims) NOT context matrix before ICA; scalefac ',num2str(scalefac)],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = ['print ',pathname,savename,'.eps -depsc'];    eval(makesave);   
        end;       
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on Context NOT ERSP %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isempty(pcersp) & ~isempty(pccxt) % ersp not pca'd, but context is   
        addeig = floatread([pathname,pccxt,'.fdt'],[45 cpcanum]);        

        row = pcstot; col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:pcstot
            tmpcomp = winv (1:length(freqs)*length(times),template)';  %makes a freqs*times X 1
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            sbplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
            set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;
            clear newwinv
            newwinv = winv(epcanum+1:end,template)'*addeig;
            sbplot(row,col,pl:pl+4)
            ph = plot([1:nlett],newwinv(epcanum+1:epcanum+nlett),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],newwinv(epcanum+nlett+1:epcanum+nlett+nfb),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],newwinv(epcanum+nlett+nfb+1:epcanum+nlett+nfb+nresp),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template == 10
            CnxtLabelRew(ctlim,.3,fs); end;
        end;
        ph = textsc(['ERSP/Context decomposition: PCA of context NOT ERSP matrix before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on  ERSP AND Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ~isempty(pcersp) & ~isempty(pccxt)% Both ersp and context are pca'd
        addeig = floatread([pathname,pccxt,'.fdt'],[45 cpcanum]);        
        erspeig = floatread([pathname,pcersp,'.fdt'],[length(freqs)*length(times) epcanum],[],0);
        
        row = pcstot; col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:pcstot
            tmpcomp =  winv (1:epcanum,template)'*erspeig';
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            sbplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
            set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
            set(gca,'fontsize',fs);  title(['Cluster: ',int2str(template)]); colorbar;

            clear newwinv
            for fac = 1:size(pccxt,2)
                newwinv = winv(epcanum+1:end,template)'*addeig';
            end;
            sbplot(row,col,pl:pl+4)
            ph = plot([1:nlett],newwinv(epcanum+1:epcanum+nlett),'ro-','linewidth',2);hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+1:nlett+nfb],newwinv(epcanum+nlett+1:epcanum+nlett+nfb),'bs-','linewidth',2); hold on;
            set(ph,'markersize',ms);
            ph = plot([nlett+nfb+1:nlett+nfb+nresp],newwinv(epcanum+nlett+nfb+1:epcanum+nlett+nfb+nresp),'g^-','linewidth',2);
            set(ph,'markersize',ms+ms/2); pl = pl+5;hold on;
            set(gca,'xlim',[0 ntot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template == 10
            CnxtLabelRew(ctlim,.3,fs); end;

        end;
        ph = textsc(['ERSP/Context decomposition: PCA of ERSP AND context matrices before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;

    end;
