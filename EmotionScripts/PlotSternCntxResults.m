% Plot results from SnglTrialERSPDecomp.m with standard context conditions (Sternberg)
%
% PlotContextResults(pathname,wtsmat,sphmat,savename,freqs,times,pcersp,epcanum,pccxt,cpcanum,pcstot,channo,color);
%
% pathname: [string] full path where wts, sph can be found and savename will be saved
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% pcersp: [string] name of float file containing eigenvectors from PCA decomposition of ERSPs. Use same
%         string as used for SnglTrialERSPDecomp.         [] if no PCA on ERSPs. 
% epcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for ERSPs
%          [] if no PCA
% pccxt: [string] name of float file containing eigenvectors from PCA decomposition of context matrix.
%         Use same string as used for SnglTrialERSPDecomp.     [] if no PCA on context matrix
% cpcanum: [integer] Number of PCA dimensions requested from SnglTrialERSPDecomp.m for context matrix
%          [] if no PCA
% pcstot: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% color: [string] or vector of RGB values for context plot

function PlotContextResults(pathname,wtsmat,sphmat,savename,freqs,times,pcersp,epcanum,pccxt,cpcanum,pcstot,channo,color);
    fs = 7;  cttot = 19; ttl = wtsmat;
   %%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sph = floatread([pathname,sphmat],[channo channo]);
    wts = floatread([pathname,wtsmat],[pcstot channo]);
    
    ws = wts*sph;winv = pinv(ws);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Plot Reg context factors onto one page %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(pcersp) & isempty(pccxt)
        if pcstot > 16
            row = 8;
        else
            row = pcstot/2; 
        end;
        col = 6;
        figure;pl = 1; ms = 4;
        for template = 1:pcstot
            if pl == row * col + 1
                ph = textsc(['ERSP/Context decomposition: NO PCA of either ERSP or context matrices before ICA'],'title'); set(ph,'fontsize',14);
                figure; pl=1;
            end;            
            cnt = 0;
            tmpcomp = winv (1:length(freqs)*length(times),template)';  %makes a freqs*times X 1
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            subplot(row,col,pl); pl = pl+1;
            imagesc(times,freqs,tmpcomp,[-lim lim]); 
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Template: ',int2str(template)]); colorbar;
            if template == 1
                ctlim = max(max(abs(winv(length(freqs)*length(times)+1:length(freqs)*length(times)+cttot,:))));
            end;
            subplot(row,col,pl:pl+1)
            ph = plot(winv(length(freqs)*length(times)+1:length(freqs)*length(times)+cttot,template),'ro-','linewidth',2); 
            set(ph,'markersize',ms); set(ph,'color',color); pl = pl+2;hold on;
            set(gca,'xlim',[0 cttot+8]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]); xid = 1;
            if template==  row | template ==pcstot 
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'rotation',90);xid = xid+1;
            end;
        end;
        ph = textsc([ttl,'; NO PCA of either ERSP or context matrices before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on ERSP NOT Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ~isempty(pcersp) & isempty(pccxt) % ersp is pca's, but context IS NOT   
        erspwinv = floatread([pathname,pcersp,'.fdt'],[length(freqs)*length(times) epcanum]);
        if pcstot > 16
            row = 8;  
        else
            row = pcstot/2; 
        end;
        col = 6;
        figure; pl = 1; ms = 4;  mns = 0;
        for template = 1:pcstot
            if pl == row * col + 1
                ph = textsc(['ERSP/Context decomp: PCA of ERSP (',int2str(epcanum),' dims) NOT context (',int2str(shq),' Qs) before ICA'],'title'); set(ph,'fontsize',14);
                figure; pl=1; mns = template-1;
            end;            
            tmpcomp =  winv (1:epcanum,template)'*erspwinv';
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            lim = max(max(abs(tmpcomp)));
            subplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; set(gca,'ydir','norm'); 
            set(gca,'yscale','log'); set(gca,'ytick',[10:10:40]);
            set(gca,'yticklabel',[10:10:40]);    set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Template: ',int2str(template)]); %colorbar;
            if template-mns < row
                set(gca,'xticklabel',[]);
            end;    
            if template == 1
                ctlim = max(max(abs(winv(epcanum+1:epcanum+cttot,:))));
            end;
            subplot(row,col,pl:pl+1)
            ph = plot(winv(epcanum+1:epcanum+cttot,template),'ro-','linewidth',2); 
            set(ph,'markersize',ms); set(ph,'color',color); pl = pl+2;hold on;
            set(gca,'xlim',[0 cttot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);xid=1;
            if template == row*2 | template == row*2-1 |template == pcstot |template == pcstot-1 
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'rotation',90);xid = xid+1;
            end;
        end;
        
        ph = textsc([ttl,'; PCA of ERSP (',int2str(epcanum),' dims) NOT context (',int2str(cttot),' Qs) before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on Context NOT ERSP %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isempty(pcersp) & ~isempty(pccxt) % ersp not pca'd, but context is   
        addwinv = floatread([pathname,pccxt,'.fdt'],[45 cpcanum]);        
        if pcstot > 16
            row = 8;
        else
            row = pcstot/2; 
        end;
        col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:pcstot
            if pl == row * col + 1
                ph = textsc(['ERSP/Context decomposition: PCA of context NOT ERSP matrix before ICA'],'title'); set(ph,'fontsize',14);
                figure; pl=1;
            end;            
            tmpcomp = winv (1:length(freqs)*length(times),template)';  %makes a freqs*times X 1
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            subplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Template: ',int2str(template)]); colorbar;
            clear newwinv
            newwinv = winv(epcanum+1:end,template)'*addwinv';
            if template == 1
                ctlim = max(abs(newwinv));
            end;
            subplot(row,col,pl:pl+1)
            ph = plot(newwinv(1:cttot),'ro-','linewidth',2); 
            set(ph,'markersize',ms);  set(ph,'color',color);    pl = pl+2;hold on;
            set(gca,'xlim',[0 cttot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template==  row | template == pcstot
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'rotation',90);xid = xid+1;
            end;
        end;
        ph = textsc([ttl,'; PCA of context NOT ERSP matrix before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA on  ERSP AND Context %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ~isempty(pcersp) & ~isempty(pccxt)% Both ersp and context are pca'd
        addwinv = floatread([pathname,pccxt,'.fdt'],[46 cpcanum]);        
        erspwinv = floatread([pathname,pcersp,'.fdt'],[length(freqs)*length(times) epcanum]);
        if pcstot > 16
            row = 8;
        else
            row = pcstot/2; 
        end;
        col = 6;
        figure; pl = 1; ms = 4;
        for template = 1:pcstot
            if pl == row * col + 1
                ph = textsc(['ERSP/Context decomposition: PCA of ERSP AND context matrices before ICA'],'title'); set(ph,'fontsize',14);
                figure; pl=1;
            end;            
            tmpcomp =  winv (1:epcanum,template)'*erspwinv';
            tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
            lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
            subplot(row,col,pl);
            imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+1; %colormap('gray');
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['Template: ',int2str(template)]); colorbar;

            clear newwinv
            for fac = 1:size(pccxt,2)
                newwinv = winv(epcanum+1:end,template)'*addwinv';
            end;
            if template == 1
                ctlim = max(abs(newwinv));
            end;
            subplot(row,col,pl:pl+1)
            ph = plot(newwinv(1:cttot),'ro-','linewidth',2); 
            set(ph,'markersize',ms);  set(ph,'color',color);    pl = pl+2;hold on;
            set(gca,'xlim',[0 cttot+1]);set(gca,'ylim',[-ctlim ctlim]);
            ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
            if template==  row | template == pcstot
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'rotation',90);xid = xid+1;
            end;
        end;
        ph = textsc([ttl,'; PCA of ERSP AND context matrices before ICA'],'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        if ~isempty(savename)
            makesave = [pathname,savename];    eval(makesave);   
        end;
    end;
