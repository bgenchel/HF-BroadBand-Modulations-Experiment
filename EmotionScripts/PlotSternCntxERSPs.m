% Plot results from SnglTrialERSPDecomp.m with template, averaged highly weighted trials and context
%
% PlotContextResults(pathname,ttl,wtsmat,sphmat,freqs,times,pcersp,epcanum,pcstot,channo,fcs,color,percut);
%
% fcs = [integer vector]: list of factors to plot
% percut = Percent of highest and lowest weighted trials to average and plot
% color = [string] or vector of RGB values for context plots
% pathname: [string] full path where wts, sph can be found and savename will be saved
% ttl: [string] title stem that the data matrix was saved with (same as used for SnglTrialERSPDecomp.m)
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

function PlotContextResults(pathname,ttl,wtsmat,sphmat,freqs,times,pcersp,epcanum,pcstot,channo,fcs,color,percut);
    
    multfac = .6; % for plotting context template
    fs = 6;
    sph = floatread([pathname,sphmat],[channo channo]);
    wts = floatread([pathname,wtsmat],[pcstot channo]);    
    ws = wts*sph;winv = pinv(ws);
    data= floatread([pathname,ttl,'.fdt'],[channo inf]);activations = ws*data;
    erspwinv = floatread([pathname,pcersp,'.fdt'],[length(freqs)*length(times) epcanum]);
    projersp = erspwinv*data(1:epcanum,:);     
    
    cttot = 19;
    
    formula = 'mean(arg1,3);';
    figure; row = length(fcs); col =9; pl = 1; 
    for fac = 1:row
        fac = fcs(fac);
        subplot(row,col,pl);  pl = pl+1; 
        hist(activations(fac,:),100);hold on; plot([0 0],[get(gca,'ylim')],'r-');       

        tpact= sort(activations(fac,:));cutoff = tpact(round(length(tpact)*(1-percut/100))); % top 10%
        hitrials = find(activations(fac,:) > cutoff); 
        selcont = data(epcanum+1:end,hitrials); %selcont = mean(selcont,2);
        for q = 1:cttot
            totcount(1,q) = length(find(data(epcanum+q,:) == 1));
            newcont(1,q) = length(find(selcont(q,:) == 1));
        end;
        newcont = newcont./totcount;
        selersp = projersp(:,hitrials);
        selersp = reshape(selersp,length(freqs),length(times),size(selersp,2));
        [rsignif,rboot] = bootstat( selersp, formula, 'boottype', 'shuffle', 'bootside', 'both','alpha',.01, 'shufflemode','regular',  'basevect', [1:size(selersp,2)], 'dimaccu', 2);
        mnmask = rsignif(:,1); mnmask = repmat(mnmask,[1 size(selersp,2)]);
        mxmask = rsignif(:,2);    mxmask = repmat(mxmask,[1 size(selersp,2)]);
        selersp = mean(selersp,3); selersp(find(selersp>mnmask&selersp<mxmask)) = 0;
        
        subplot(row,col,pl)
        tmpcomp =  winv (1:1:epcanum,fac)'*erspwinv';
        tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
        lim = max(max(abs(tmpcomp))); 
        imagesc(times,freqs,tmpcomp,[-lim lim]); set(gca,'ydir','norm'); 
        set(gca,'yscale','log'); set(gca,'ytick',[10:10:40]);
        set(gca,'ytick',[5:5:45]); set(gca,'yticklabel',{[] 10 [] 20 [] 30 [] 40 []});
        set(gca,'ticklength',[.02 .02]);
        hold on; 
        if pl < row*col-col            
            set(gca,'xticklabel',[]);
        end;
        
        plot([0 0],[get(gca,'ylim')],'k-'); pl = pl+1; 
        title(['Template ',int2str(fac)]);

        subplot(row,col,pl)
        %if pl == 3
        lim = max(max(abs(selersp)))-5; 
        %end;
        imagesc(times,freqs,selersp,[-lim lim]); set(gca,'ydir','norm'); 
        set(gca,'yscale','log'); set(gca,'ytick',[10:10:40]);
        set(gca,'ytick',[5:5:45]); set(gca,'yticklabel',{[] 10 [] 20 [] 30 [] 40 []});
        set(gca,'ticklength',[.02 .02]); colorbar;
        hold on;
        if pl < row*col-col            
            set(gca,'xticklabel',[]);
        end;
        plot([0 0],[get(gca,'ylim')],'k-'); pl = pl+1; 
        title(['Mean ',int2str(fac),' Pos']);
        
        subplot(row,col,pl:pl+1)
        ph=plot(winv(epcanum+1:end,fac)*multfac,'y.-'); hold on;set(ph,'color',[1 .5 0]);
        ph=plot([1:cttot],newcont-(percut/100),'ro-'); hold on; set(ph,'color',color); 
        set(ph,'linewidth',2);set(ph,'markersize',5);
        set(gca,'xlim',[0 cttot+1]);        plot([get(gca,'xlim')],[0 0],'k-');
        if pl == 4
            ctlim = max(newcont)+.1;
        end;
        set(gca,'ylim',[-ctlim ctlim]);
        title(['Pos Trials (',int2str(length(hitrials)),')']);
            if pl == row*col - 5
                xid = 1; 
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
            end;
        pl = pl+2;
        
        tpact= sort(activations(fac,:));cutoff = tpact(round(length(tpact)*(percut/100))); % bottom 10%
        lotrials = find(activations(fac,:) < cutoff); 
        selcont = data(epcanum+1:end,lotrials); 
        for q = 1:cttot
            totcount(1,q) = length(find(data(epcanum+q,:) == 1));
            newcont(1,q) = length(find(selcont(q,:) == 1));
        end;
        newcont = newcont./totcount;
        selersp = projersp(:,lotrials);
        selersp = reshape(selersp,length(freqs),length(times),size(selersp,2));
        [rsignif,rboot] = bootstat( selersp, formula, 'boottype', 'shuffle', 'bootside', 'both','alpha',.01, 'shufflemode','regular',  'basevect', [1:size(selersp,2)], 'dimaccu', 2);
        mnmask = rsignif(:,1); mnmask = repmat(mnmask,[1 size(selersp,2)]);
        mxmask = rsignif(:,2);    mxmask = repmat(mxmask,[1 size(selersp,2)]);
        selersp = mean(selersp,3); selersp(find(selersp>mnmask&selersp<mxmask)) = 0;
        
        subplot(row,col,pl)
        lim = max(max(abs(tmpcomp))); 
        tmpcomp =  winv (1:epcanum,fac)'*erspwinv';
        tmpcomp = reshape(tmpcomp*-1,length(freqs),length(times));  %Good
        imagesc(times,freqs,tmpcomp,[-lim lim]); set(gca,'ydir','norm'); 
        set(gca,'yscale','log'); set(gca,'ytick',[10:10:40]);
        set(gca,'ytick',[5:5:45]); set(gca,'yticklabel',{[] 10 [] 20 [] 30 [] 40 []});
        set(gca,'ticklength',[.02 .02]);
        hold on; 
        if pl< row*col-col            
            set(gca,'xticklabel',[]);
        end;
        plot([0 0],[get(gca,'ylim')],'k-'); pl = pl+1; 
        title(['Template ',int2str(fac)]);
        subplot(row,col,pl)
        lim = max(max(abs(selersp)))-5; 
        imagesc(times,freqs,selersp,[-lim lim]); set(gca,'ydir','norm'); 
        set(gca,'yscale','log'); set(gca,'ytick',[10:10:40]);
        set(gca,'ytick',[5:5:45]); set(gca,'yticklabel',{[] 10 [] 20 [] 30 [] 40 []});
        set(gca,'ticklength',[.02 .02]); colorbar;
        hold on; 
        if pl < row*col-col            
            set(gca,'xticklabel',[]);
        end;
        plot([0 0],[get(gca,'ylim')],'k-'); pl = pl+1; 
        title(['Mean ',int2str(fac),' Neg']);
        
        subplot(row,col,pl:pl+1)
        ph=plot(winv(epcanum+1:end,fac)*multfac*-1,'y.-'); hold on;set(ph,'color',[1 .5 0]);
        ph=plot([1:cttot],newcont-(percut/100),'ro-'); hold on; set(ph,'color',color);  
        set(ph,'linewidth',2);set(ph,'markersize',5);
        set(gca,'xlim',[0 cttot+1]); set(gca,'ylim',[-ctlim ctlim]);
        plot([get(gca,'xlim')],[0 0],'k-');title(['Neg Trials (',int2str(length(lotrials)),')']);
        pl = pl+2;
            if pl >  row*col 
                xid = 1;
                ph = text(xid,-ctlim,' 0 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' 0 letter =  mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-1 letter =  mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter =  none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-2 letter = mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = none'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = ignore'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,'-3 letter = mem'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 0 '); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 1'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 2'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 3'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 4'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 5'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 6'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
                ph = text(xid,-ctlim,' load = 7'); set(ph,'fontsize',fs);
                set(ph,'rotation',90);xid = xid+1;
            end;

    end;
    textsc([wtsmat,' Factor Templates and Positive and Neg (',int2str(percut),'%) mean trials to verify factors'],'title');
    set(gcf,'PaperOrientation','landscape'); set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

    
    

