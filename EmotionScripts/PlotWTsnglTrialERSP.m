% Plot results from SnglTrialICAClustering.m using the 'wt' option and NO context
%
% PlotWTsnglTrialERSP(pathname,savedat,savename,complist,toplot,whichfacs);
%
% savedat: [string] title stem used in SnglTrialERSPDecomp.m
% pathname: [string] full path where wts, sph can be found and savename will be saved
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% pcs: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% toplot: [integer] Max number of factors to plot on one page.
% plotwts: if not [], will plot weights over the course of the session
%          0 will plot straight weights, any other number will plot weights smoothed by this factor
% whichfacs: vector of integers corresponding to factors to plot. [] for all

function PlotWTsnglTrialERSP(pathname,savedat,savename,toplot,plotwts,whichfacs);

    fs = 9; % fontsize
%%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = load([pathname,savedat,'.mat']);     
    wts = floatread([pathname,savedat,'.wts'],[s.pcmat s.numrows],[],0);
    sph = floatread([pathname,savedat,'.sph'],[s.numrows s.numrows],[],0);  
    data = floatread([pathname,savedat,'.fdt'],[s.numrows inf],[],0);  
    ws = wts*sph;winv = pinv(ws);activations = ws*data;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%% Plot Reg context factors onto one page %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(whichfacs)
        whichfacs = [1:s.pcmat];
    end;
    EEG = pop_loadset( 'sources.set',pathname);
    %figure;pl = 1; ms = 4;
    figure;pl = 2;
    lim = 5; row = toplot; col = length(s.complist)+1;
    for mp = 1:length(s.complist)
        sbplot(row,col,pl)
        topoplot(EEG.icawinv(:,s.complist(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
        title(int2str(s.complist(mp)));
    end;    
    for tp = 1:length(whichfacs)%pcs
        template = whichfacs(tp);
        if pl == row * col + 1
            ph = textsc([savedat,'; TW ERSP decomposition (co-modulation)'],'title'); set(ph,'fontsize',14);
            figure; pl=2;
            for mp = 1:length(s.complist)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,s.complist(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
                title(int2str(s.complist(mp)));
            end;    
        end;            
        sbplot(row,col,pl)
        hist(winv(:,template),50);pl = pl+1;hold on;set(gca,'xlim',[-2.5 2.5]);
        plot([0 0],[get(gca,'ylim')],'r-');set(gca,'fontsize',6);
            set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
                title(['IM ',int2str(template)]);
        cnt = 0;
        lim = max(max(activations(template,:)))-.1*max(max(activations(template,:)));
        for cp = 1:length(s.complist)
            clear onecomp
            onecomp = activations(template,(cp-1)*(length(s.times)*length(s.freqs))+1:cp*(length(s.times)*length(s.freqs)));
            onecomp = reshape(onecomp,length(s.freqs),length(s.times));
            sbplot(row,col,pl)
            imagesc(s.times, s.freqs, onecomp,[-lim lim]);pl = pl+1;hold on;
            set(gca,'ydir','norm');plot([get(gca,'xlim')],[10 10],'k:');
            set(gca,'ticklength',[.02 .02]);set(gca,'fontsize',fs);
            plot([0 0],[get(gca,'ylim')],'k-');set(gca,'ytick',[5:5:45]);
            set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
            if template == row - 1 | template == s.pcs 
                set(gca,'xticklabel',[0:500:1000]);set(gca,'ytick',[10:10:s.freqs(end)]);
            end;
        end;
    end;
    ph = textsc([savedat,'; WT ERSP decomposition (co-modulation)'],'title'); set(ph,'fontsize',14);
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    if ~isempty(savename)
        makesave = [pathname,savename];    eval(makesave);   
    end;

    
    if ~isempty(plotwts)
        clear wts sph ws
        figure;     pl=1;
        row = round(sqrt(length(whichfacs)));col = ceil(sqrt(length(whichfacs)));
        for im = 1:length(whichfacs)
            dim = whichfacs(im);
            sbplot(row,col,pl); pl = pl+1;
            if plotwts ~= 0
                newdata = zeros(0,size(winv,2));
                for em = 1:length(s.numtrials)
                    xwidth = round(length(sum(s.numtrials(1:em-1))+1:sum(s.numtrials(1:em)))/40*plotwts);
                    [outdata,outx] = movav(winv(sum(s.numtrials(1:em-1))+1:sum(s.numtrials(1:em)),dim),0,xwidth,0);
                    newdata(end+1:end+size(outdata,2),dim) = outdata;
                    intvls(em,:) = [size(outdata,2) length(sum(s.numtrials(1:em-1))+1:sum(s.numtrials(1:em))) - size(outdata,2)];
                end;
                for em = 1:size(intvls,1)
                    datlen = intvls(em,1);
                    datdiff = intvls(em,2);
                    
                    ph = plot([sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen],newdata(sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen,dim),'k-','linewidth',1); hold on;
                end;
                set(gca,'xlim',[0 sum(intvls(:,1))]);
                yl = get(gca,'ylim');
                for e = 1:length(s.numtrials)
                    ph = plot([sum(intvls(1:(e-1),1)) sum(intvls(1:(e-1),1))],[yl],'k-');
                end; 
            else            
                ph = plot(winv(:,dim));hold on;
                set(gca,'xlim',[0 sum(s.numtrials)]);
                yl = get(gca,'ylim');
                for e = 1:length(s.numtrials)
                    if labels == 1
                        ph = text(sum(s.numtrials(1:e-1))+1,2.5,e);
                        set(ph,'rotation',45);  
                    end;
                    ph = plot([sum(s.numtrials(1:(e-1))) sum(s.numtrials(1:(e-1)))],[yl],'k-');
                end; 
            end;
            ph = plot([get(gca,'xlim')],[0 0],'k--','linewidth',2);
            set(ph,'color','r');
            title([pathname(end-4:end-1),' IM ',int2str(dim)]);
        end;
    end;
    
