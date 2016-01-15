% takes a cell array of clusters with factor numbers across subjects
% plots superimposed spectral modulation templates into the specified axis
% 
% PlotSpecCoModClusts(datpath,savedat,complist,facvec,figinfo);
%
% INPUTS:
% datpath -- full directory path where dataset can be found
% savedat -- filename that float files were saved as (same as input to SpecCoModAnal.m)
% complist -- full list of components in orginal decomposition
% facvec -- cell array of integers referring to all the clustered factors for each subject
% figinfo -- [nrows,ncols,place] number of rows, number of columns in figure and the subplot 
%            in which to start plotting.
%
%
%
%

function PlotSpecCoModClusts(fullpaths,savedat,complist,facvec,figinfo);
        
    
    pl = figinfo(3);
    fprintf('\nCollecting cluster spectral co-modulation data to plot...');
    allmins = zeros(1,0);
    allmaxs = zeros(1,0);
    allplots = zeros(1,0);
    for nx = 1:length(facvec)
        if ~isempty(facvec{nx})
            s = load([fullpaths{nx},savedat,'.mat']); 
            sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
            %sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            %wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
            ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws
            cols = hsv(length(complist{nx}));      
            for cp = 1:length(facvec{nx})
                allplots(end+1) = sbplot(figinfo(1),figinfo(2),pl);
                plot([s.freqs(1) s.freqs(end)],[0 0],'k--'); hold on;
                for rcp = 1:length(complist{nx})
                    ph = plot(s.freqs,activations(facvec{nx}(cp),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                    hold on; set(ph,'color',cols(rcp,:));
                    set(gca,'xlim',[s.freqs(1) s.freqs(end)]);
                end;title(['Sb ',int2str(nx),';Fac ',int2str(facvec{nx}(cp))]); 
                allmins(end+1) = min(activations(facvec{nx}(cp),:));
                allmaxs(end+1) = max(activations(facvec{nx}(cp),:));
                %set(gca,'ylim',[min(activations(facvec{nx}(cp),:))-.1 max(activations(facvec{nx}(cp),:))+.1]);
                pl = pl+1;
            end;
            fprintf('.');
        end;
    end;
    for sb = 1:length(allplots)
        set(allplots(sb),'ylim',[min(allmins) max(allmaxs)]);
    end;
    axcopy
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    fprintf('done. \n');
