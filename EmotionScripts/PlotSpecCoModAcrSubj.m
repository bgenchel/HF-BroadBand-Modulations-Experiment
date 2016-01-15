% Plots spectral co-modulation templates, specified by cluster analysis, from multiple subjects
% 
%   [] = PlotSpecCoModAcrSubj(fullpaths,complist,facvec,freqs,plotrows,title);
%
%  mapset -- [string] full name of EEGLAB dataset that will be loaded for ICA scalp map plotting.
%  fullpaths --  cell array of strings where spectral data and dipole info dataset 'source.set' can be found
% savedat -- [string] file name stem to use opening 'Stuff' and float files
%  complist -- cell array of integers of components used in analysis for all subjects
%  facvec -- cell array of integers referring to all the clustered factors for each subject
%  freqs -- output of SpecCoModAnal.m
%  plotrows -- max number of rows to plot on each page
%  ttl -- [string] desired text for title of each page

function [] = PlotSpecCoModAcrSubj(mapset,fullpaths,savedat,complist,facvec,freqs,plotrows,ttl);

    manyfacs = 0;
    for nx = 1:length(complist)
        if ~isempty(facvec{nx})
            forcols(1,nx) = length(complist{nx});
            manyfacs = manyfacs + length(facvec{nx});
        end;
    end;
    col = max(forcols)+1;
    row = plotrows; 
    for rw = 1:row-1
    rowstarts(1,rw) = col*rw+1;
    end;
    fprintf('\nCollecting spectral co-modulation data to plot...\n');
    figure;pl = 1; rw = 1;
    
    for nx = 1:length(complist)
        if ~isempty(facvec{nx})
            s = load([fullpaths{nx},savedat,'.mat']); 
            sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
            %sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            %wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
            ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws
            EEG = pop_loadset(mapset, fullpaths{nx});
            if rw >= length(rowstarts) - 1
                set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                ph = textsc(ttl,'title'); set(ph,'fontsize',14);
                figure; pl=1;rw = 1;
            end; 
            pl = pl+1;           
            for mp = 1:length(complist{nx})
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,complist{nx}(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
                title(int2str(complist{nx}(mp)));
            end;     
            pl = rowstarts(rw);rw = rw+1; 
            for tp = 1:length(facvec{nx})%pcs
                template = facvec{nx}(tp);
                if pl > max(rowstarts)
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    ph = textsc(ttl,'title'); set(ph,'fontsize',14);
                    figure; pl=2;rw = 1;
                    for mp = 1:length(complist{nx})
                        sbplot(row,col,pl)
                        topoplot(EEG.icawinv(:,complist{nx}(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
                        title(int2str(complist{nx}(mp)));
                    end;  
                    pl = rowstarts(rw); rw = rw+1; 
                end; 
                sbplot(row,col,pl)
                hist(winv(:,facvec{nx}(tp)),75);pl = pl+1;hold on;
                set(gca,'fontsize',7);%set(gca,'xlim',[-1.5 1.5]);
                plot([0 0],[get(gca,'ylim')],'r-');
                set(gca,'yticklabel',[]); title([int2str(nx),'-',int2str(facvec{nx}(tp))]);
                for rcp = 1:length(complist{nx})
                    sbplot(row,col,pl)
                    plot(s.freqs,activations(facvec{nx}(tp),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                    pl = pl+1;hold on;
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'xlim',[s.freqs(1) s.freqs(end)]); set(gca,'ylim',[-2 10]);   
                    set(gca,'xtick',[5:5:s.freqs(end)]);  set(gca,'xticklabel',{[] 10 [] [] [] 30 [] [] [] []}); 
                    set(gca,'xgrid','on');
                end;
                if rw > length(rowstarts)
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    ph = textsc(ttl,'title'); set(ph,'fontsize',14);
                    figure; pl=2;rw = 1;
                    for mp = 1:length(complist{nx})
                        sbplot(row,col,pl)
                        topoplot(EEG.icawinv(:,complist{nx}(mp)),EEG.chanlocs,'plotrad',.5);pl = pl+1;
                        title(int2str(complist{nx}(mp)));
                    end;  
                    pl = rowstarts(rw); rw = rw+1; 
                else
                     pl = rowstarts(rw);rw = rw+1;                
                end; 
            end;       
            ph = textsc(ttl,'title'); set(ph,'fontsize',14);
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        end;
    end;
    
    
    
    
    
    
