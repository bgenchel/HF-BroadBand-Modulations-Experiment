% plots results from SpecCoModAnal.m
%  PlotSpecCoModCtx(datpath,complist,savedat,frqlim,pcs,plotpcs,numtrials,numframes,freqs,emos,superpos);
%
% INPUTS:
% datpath: full directory path where dataset can be found
% savedat: filename that float files were saved as (same as input to SpecCoModAnal.m)
% frqlim: [minfreq maxfreq] to plot
% plotpcs: number of final PCA/ICA dimensions to plot
% superpos -- ['y' or 'n'] 'y' will plot a single axis with all comps superimposed.
%             (can only be used for plotting all comps (justone needs to be []) 

function PlotSpecCoModCtx(datpath,savedat,frqlim,plotpcs,emos,superpos);
    
    s = load([datpath,savedat,'Stuff.mat']); 

    sph=floatread([datpath,savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([datpath,savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
    pcaeig = floatread([datpath,savedat,'EIGVEC.fdt'],[length(s.freqs) s.pcnum],[],0);    
    ws = wts*sph;    
    if size(ws,1) == size(ws,2)
        winv = inv(ws); 
    else
        winv = pinv(ws);
    end;
    clear wts sph ws
    specnum = size(pcaeig,2);
    spectemps = pcaeig * winv(1:specnum,:);spectemps = spectemps';
    emowts = winv(specnum+1:end,:);
    
    if superpos == 'y'
        figure;cols = jet(plotpcs); 
        for tp = 1:plotpcs
            ph = plot(s.freqs,spectemps(tp,:)); 
            hold on; set(ph,'color',cols(tp,:));            
        end; textsc(['Superimposed Factor Spectra: ',savedat],'title');
        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        axcopy
    else  
        figure;row = 5; 
        if row > 16
            row = round(row/2);
        end;            
        col = 4;

        pl = 1;           
        hlim = max(spectemps(:));llim = min(spectemps(:));
        cols = jet(length(emos));
        for tp = 1:plotpcs
            if pl == row*col+1
                set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                set(ph,'fontsize',14);
                figure; pl = 1;
            end;
            sbplot(row,col,pl)
            ph = plot(s.freqs,spectemps(tp,:),'k-','linewidth',2); pl = pl+1;hold on;
            set(gca,'fontsize',7);set(gca,'box','off');
            set(gca,'xlim',frqlim); set(gca,'ylim',[llim hlim]);   
            set(gca,'xtick',[5:5:frqlim(2)]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40 [] []}); 
            set(gca,'xgrid','on'); set(gca,'yticklabel',[]);
            plot([get(gca,'xlim')],[0 0],'g'); title(['Factor ',int2str(tp)]);
            sbplot(row,col,pl)
            for e = 1:size(winv(specnum+1:end,:),1)
                ph = bar(e,emowts(e,tp)); hold on;
                set(gca,'ylim',[min(emowts(:,tp)) max(emowts(:,tp))]);
                set(ph,'facecolor',cols(e,:));
                yl = get(gca,'ylim');
                ph = text(e,yl(1),emos{e});
                set(ph,'rotation',90);%set(ph,'color',cols(e,:));
            end;pl = pl+1;
            set(gca,'xlim',[0 e+1]); 
        end;
        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        ph=textsc(['Spectral Templates and associated emotion indexes for ',datpath,': ',savedat],'title');
        set(ph,'fontsize',14);
        axcopy
    end;
    
