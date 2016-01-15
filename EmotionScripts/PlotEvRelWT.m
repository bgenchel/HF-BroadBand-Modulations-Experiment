% Plot results from EvRelSnglTrialWT.m 
%
% PlotEvRelWT(plotset,fullpath,name,toplot,plotwts,whichfacs);
%
% plotset -- [string] name of any dataset with matching ICA weights for scalp maps plotting
% fullpath: [string] full path where wts, sph can be found and savename will be saved
% name: [string] title stem used in SnglTrialERSPDecomp.m
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% pcs: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% toplot: [integer] Max number of factors to plot on one page.
% plotwts: if not [], will plot, on a new figure, weights over the trials as they were given to ICA
%          0 will plot straight weights, any other number will plot weights smoothed by this factor
% whichfacs: vector of integers corresponding to factors to plot ([] to plot all)
% multfac -- vector corresponding to the IMs in 'whichfacs' will scale each corresponding IM. Typically want 1 | -1. Default is all ones. 
function PlotEvRelWT(plotset,fullpath,name,toplot,plotwts,whichfacs,multfac);

    
    if exist('multfac')
       if isempty(multfac)
           multfac = ones(1,length(whichfacs));
       end;
    else
        multfac = ones(1,length(whichfacs));   
    end;
    
    fs = 7; % fontsize for time/freq plots
    
    %pop_editoptions( 'option_storedisk', 1, 'option_savematlab', 1, 'option_computeica', 1, 'option_rememberfolder', 1);

    %%%  load in data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = load([fullpath,name,'.mat']);     
    wts = floatread([fullpath,name,'.wts'],[s.pcs s.numtrials],[],0);
    sph = floatread([fullpath,name,'.sph'],[s.numtrials s.numtrials],[],0);  
    data = floatread([fullpath,name,'.fdt'],[s.numtrials inf],[],0);  
    ws = wts*sph;winv = pinv(ws);
    activations = ws*data; % create template matrix
    clear wts sph ws data  % save memory space
    
    %%%%%%  Check for log frequency spacing  %%%%%%%%%%%%%%%%
    tst1 = sqrt(s.freqs);
    tst2 = log(s.freqs);
    lfr = linspace(sqrt(s.freqs(1)),sqrt(s.freqs(end)),length(s.freqs));

    if tst1 == lfr
        frspace = 1; % quadratic spacing
    elseif tst2 == lfr 
        frspace = 2; % log spacing
    else
        frspace = 0; % linear spacing
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%% Plot Time/freq templates %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(whichfacs)
        whichfacs = [1:s.pcs]; % default is all IMs
    end;   
        
    EEG = pop_loadset( plotset,fullpath);
    %figure;pl = 1; ms = 4;
    figure;pl = 2;
    lim = 5; row = toplot+1; col = length(s.complist)+1;
    for mp = 1:length(s.complist)
        sbplot(row,col,pl)
        topoplot(EEG.icawinv(:,s.complist(mp)),EEG.chanlocs,'plotrad',.7);pl = pl+1;
        title(['IC',int2str(s.complist(mp))]);
    end;    
    lim = max(activations(:))-.1*max(activations(:));
    for tp = 1:length(whichfacs)%pcs
        template = whichfacs(tp);
        if pl == row * col + 1
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            ph = textsc([fullpath(end-4:end),': ',name,'; WT ERSP decomposition (co-modulation)'],'title'); set(ph,'fontsize',14);
            figure; pl=2;
            for mp = 1:length(s.complist)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,s.complist(mp)),EEG.chanlocs,'plotrad',.7);pl = pl+1;
                title(['IC',int2str(s.complist(mp))]);
            end;    
        end;            
        sbplot(row,col,pl)
        hist(winv(:,template),50);pl = pl+1;hold on;set(gca,'xlim',[-2.5 2.5]);
        plot([0 0],[get(gca,'ylim')],'r-');set(gca,'fontsize',6);
            set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
                title(['IM ',int2str(template)]);
        cnt = 0;
        for cp = 1:length(s.complist)
            clear onecomp
            onecomp = activations(template,(cp-1)*(length(s.times)*length(s.freqs))+1:cp*(length(s.times)*length(s.freqs)));
            onecomp = reshape(onecomp,length(s.freqs),length(s.times));
            sbplot(row,col,pl)
            if frspace == 1
                [yticks] = quadimagesc(s.times, s.freqs, onecomp*multfac(tp),[-lim lim]);pl = pl+1;hold on;
            elseif frspace == 2                
                mylogimagesc(s.times, s.freqs, onecomp*multfac(tp),[-lim lim]);pl = pl+1;hold on;
            else
                 imagesc(s.times, s.freqs, onecomp*multfac(tp),[-lim lim]);pl = pl+1;hold on;               
                 set(gca,'ytick',[10:10:s.freqs(end)]);
            end;
            set(gca,'ydir','norm');
            set(gca,'ticklength',[.04 .04]);set(gca,'fontsize',fs);
            plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'xtick',[500:500:s.times(end)]);
            set(gca,'xticklabel',{[] 1000 [] 2000});            
            if pl < (row-1)*col
                set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
            end;
            if isfield(s,'medwarpevs')
                cols = lines(length(s.medwarpevs));
                for wp = 1:length(s.medwarpevs)
                    ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                    set(ph,'color',cols(wp,:));
                end;
            end;              
            ph = plot([get(gca,'xlim')],[yticks(2) yticks(2)],'k:');
            ph = plot([get(gca,'xlim')],[yticks(4) yticks(4)],'k:');
        end;
    end;
    ph = textsc([fullpath(end-4:end),': ',name,'; WT ERSP decomposition (co-modulation)'],'title'); set(ph,'fontsize',14);
    set(gcf,'Position',[100 300 1400 900]);
    set(gcf,'PaperOrientation','landscape');    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

    if ~isempty(plotwts)
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
            title([fullpath(end-4:end-1),' IM ',int2str(dim)]);
        end;
    end;
    
