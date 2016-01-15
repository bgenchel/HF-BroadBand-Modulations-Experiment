% Plot results from SnglTrialICAClustering.m using the 'wt' option and NO context
%
% PlotTWsnglTrialERSP(mappath,path,filename,complist,freqs,times,channo,toplot,whichfacs,ttl);
%
% mappath -- cell array of [stings]: paths where 'sources.set' can be found for scalp maps
% filename: cell array of [strings] title stems data is saved as in SnglTrialERSPDecomp.m
% path -- [string]: full path where wts, sph can be found 
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% complist:list of components run in ERSP decomposition
% faclist --  (cell array of integers) indices of factors to plot for each subject
% toplot: [integer] Max number of factors to plot on one page.
% whichfacs: cell array of vectors: indices of factors to plot for each subject

function PlotContextResults(mappath,path,filename,complist,freqs,times,channo,toplot,whichfacs,ttl);

    if ~exist('ttl')
        ttl = 'WT ERSP decomposition (co-modulation)';
    end;    
    cd (path);
    fs = 9; % fontsize
    for nx = 1:length(complist)
        forcols(1,nx) = length(complist{nx});
    end;
    col = max(forcols)+1;
    row = toplot;  
    for rw = 1:row
    rowstarts(1,rw) = col*rw+1;
    end;
    lim = 5;
    figure;pl = 1; rw = 1;
    for nx = 1:length(complist)
        if ~isempty(whichfacs{nx})
        %%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s=load([filename{nx},'Stuff.mat']);
        data= floatread([path,filename{nx},'.fdt'],[s.channo inf]);
        sph = floatread([path,filename{nx},int2str(s.channo),'pc',int2str(s.pcs),'.sph'],[s.channo s.channo]);
        wts = floatread([path,filename{nx},int2str(s.channo),'pc',int2str(s.pcs),'.wts'],[s.pcs s.channo]);   
        ws = wts*sph;winv = pinv(ws);activations = ws*data;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%%%% Plot Reg context factors onto one page %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EEG = pop_loadset( 'sources.set',mappath{nx}); 
        if pl == row * col + 1
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            ph = textsc('WT ERSP decomposition (co-modulation)','title'); set(ph,'fontsize',14);
            figure; pl=1;rw = 1;
       end; 
        pl = pl+1;           
        for mp = 1:length(complist{nx})
            subplot(row,col,pl)
            topoplot(EEG.icawinv(:,complist{nx}(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
            title(int2str(complist{nx}(mp)));
        end;     
        pl = rowstarts(rw);rw = rw+1; 
        for tp = 1:length(whichfacs{nx})%pcs
            template = whichfacs{nx}(tp);
            if pl == row * col + 1
                set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                ph = textsc(ttl,'title'); set(ph,'fontsize',14);
                figure; pl=2;rw = 1;
                for mp = 1:length(complist{nx})
                    subplot(row,col,pl)
                    topoplot(EEG.icawinv(:,complist{nx}(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
                    title(int2str(complist{nx}(mp)));
                end;  
                pl = rowstarts(rw); rw = rw+1; 
            end;            
            subplot(row,col,pl)
            hist(winv(:,template),50);pl = pl+1;hold on;set(gca,'xlim',[-2.5 2.5]);
            plot([0 0],[get(gca,'ylim')],'r-');set(gca,'fontsize',6);
            set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
            title(['Subj: ',int2str(nx)]);
            cnt = 0;
            lim = max(max(activations(template,:)))-.1*max(max(activations(template,:)));
            for cp = 1:length(complist{nx})
                clear onecomp
                onecomp = activations(template,(cp-1)*(length(s.times)*length(s.freqs))+1:cp*(length(s.times)*length(s.freqs)));
                onecomp = reshape(onecomp,length(s.freqs),length(s.times));
                subplot(row,col,pl)
                imagesc(times, freqs, onecomp,[-lim lim]);pl = pl+1;hold on;
                set(gca,'ydir','norm');plot([get(gca,'xlim')],[10 10],'k:');
                set(gca,'yscale','log');set(gca,'ticklength',[.02 .02]);set(gca,'fontsize',fs);
                plot([0 0],[get(gca,'ylim')],'k-');set(gca,'ytick',[5:5:45]);
                set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
                if rw == length(rowstarts)
                    set(gca,'xticklabel',[0:500:1000]);set(gca,'yticklabel',{5 10 [] 20 [] [] [] [] 45});
                end;
            end;
            pl = rowstarts(rw);rw = rw+1;
        end;
        ph = textsc(ttl,'title'); set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        end;
    end;
    
    

