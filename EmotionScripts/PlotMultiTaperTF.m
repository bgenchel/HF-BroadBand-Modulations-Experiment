% plots single subject results from 'runtimef.m'
%
% PlotRuntimefOUT(datafile,mapset,pathname,complist,origlist,mask,clim,frqlim,tmlim,ttls,sigdiff);
%
% datafile: full path and file name with .mat extension
% sigdiff: if 1 then assumes two conditions were run, 0 if single condition timef
% mask: if 1 then mask for significance from erspboot
% mapset: dataset name with .set extension from which scalp map data can be extracted
% path: full path where mapset can be found
% clim: color limit to use for all plots. will center the given number around 0 so green is 0
% frqlim: [min max] Frequency limits for plotting
% tmlim: [min max] Time limits to plot (ms)
% ttls: {string(s)}; Cell array of strings for single or two conditions timef plots (1 string per cond)
% complist: components to plot ersp for
% invdiff -- if 1, will plot the inverse (*-1) of the 'difference' image (sigdiff must be 1 as well)
% **  if time warping was used, will look for a vector of latencies in the data variables called 'warpevs'
% requires a field in datafile called 'complist' with all original ICs decomposed.
% plotopt -- [row col place] number of rows, columns in the figure. place is sbplot to plot into
%                            if [], will pop a figure and decide row and col
%
function [pl] = PlotMultiTaperTF(datafile,mapset,pathname,plotcomps,clim,frqlim,tmlim,plotopt,ttl);
    
            
if ~exist('plotopt')
  plotopt = [];
end;
if isempty(plotopt)
  figure;
  if ~isempty(mapset)
    plotopt = [round(sqrt(3*length(plotcomps))) ceil(sqrt(3*length(plotcomps))) 1]; 
    if plotopt(2) > 9
      plotopt(2) = 9; 
    end;
    if ~ismember(plotopt(2),[3 6 9])
      plotopt(2) = plotopt(2)+1;  
      if ~ismember(plotopt(2),[3 6 9])
        plotopt(2) = plotopt(2)+1;plotopt(1) = plotopt(1) - 1;
      end;
    end;
  else    
    plotopt = [round(2*sqrt(length(plotcomps))) ceil(sqrt(2*length(plotcomps))) 1]; 
    if ~ismember(plotopt(2),[2 4 6 8 10])
      plotopt(2) = plotopt(2)+1;
    end;
  end
end;
row = plotopt(1); col = plotopt(2); pl = plotopt(3);

% load dataset for scalp map if present
if ~isempty(mapset)
  EEG = pop_loadset(mapset,pathname); 
end;

% load TF data

s = load ([pathname,datafile]);
if isempty(frqlim)
  frqlim = [s.freqs(1) s.freqs(end)];
end;
if isempty(tmlim)
  tmlim = [s.times(1) s.times(end)];
end;

fr=find(s.freqs >= frqlim(1) & s.freqs<=frqlim(2));  
tms=find(s.times > tmlim(1) & s.times<tmlim(2));


for ic = 1:length(plotcomps)
  if ~isempty(mapset)
    sbplot(row,col,pl);pl=pl+1;
    topoplot(EEG.icawinv(:,plotcomps(ic)),EEG.chanlocs,'electrodes','off','plotrad',.65);   
    title([int2str(plotcomps(ic))]);
  end;
  data = 10*log10(s.icpwr{find(plotcomps(ic)==s.complist)}); baseline = mean(mean(data,3),2);
  data = data - repmat(baseline,[1 size(data,2) size(data,3)]);
  if isempty(clim)
    newlim = mean(data(fr,tms,:),3); newlim = max(abs(newlim(:)));
  else
    newlim = clim;
  end;
  sbplot(row,col,[pl pl+1]);  pl=pl+2;          
  imagesc(s.times(tms),s.freqs(fr),mean(data(fr,tms,:),3),[-newlim newlim]); hold on; set(gca,'ydir','norm'); 
  plot([0 0],[get(gca,'ylim')],'k-');title([pathname(end-4:end-1),'; IC ',int2str(plotcomps(ic))]);cbar;
end;

ph = textsc(ttl,'title'); %ph = set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    