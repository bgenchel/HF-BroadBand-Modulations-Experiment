% Uses WT spectral decomp and reconstructs an ERSP by
%
% [outersp,keeptrack,plotpnts] = PlotERSPfromWT(savedat,datset,datpath,plotdims,mask,allbigs,plotopt);
%
% plotdims -- [vector] dimension indices to plot, each with it's own figure.
% condvec -- [vector] vector of numbers the same size as the keeptrack,1 to
% indicate groupings of trials to average and plot separately as distinct
% conditions. [] will use first column of keeptrack to group conditions.
% frqlim -- [minfrq maxfrq] frequency limits for plots
% mask -- ['mask' or 'off] mask reconstructed ERSPs by permutation
% plotopt -- [place in current figure or 'off'] plot the results(integer) or just return output ('off')
%            Many conditions will default to place = 1; one condition will plot starting at place with
%            n columns equal to length(s.complist) and row = col-1;
%            Enter [] for default plotting parameters and to create a new figure;

function [meanersp,differsps] = PlotMeanERSPfromWT(savedat,datpaths,condvec,mask,clustcps,plotopt);

plttype = 'ersp'; % 'time' or 'ersp'
plotdiff = 'on';

if isempty(plotopt)
  plotopt = 'off';
end;
if ~isstr(plotopt)  % plot into existing figure with:
  row = plotopt(1); col = plotopt(2); place = plotopt(3);
elseif strcmp(plotopt,'on')
  ncomps = 0;
  for x = 1:length(clustcps)
    ncomps = ncomps + length(clustcps{x});
  end;  
  col = ceil(sqrt(ncomps+1));
  row = ceil(sqrt(ncomps));
  if col > 6
    col = 6; row = 6;
  end;
elseif strcmp(plotopt,'off')
  figure;place = 1; 
  row = 100; col = 100; % don't pop another figure;
end;
clim = [];
s = load([datpaths{2},savedat,'.mat']);
for cond = 1:length(condvec)
  meanersp{cond} = zeros(length(s.freqs),length(s.tmpoints),0);
  if strcmp(plotopt,'on')
    figure;place = 1; 
  end;
  for nx = 1:length(clustcps)
    if ~isempty(clustcps{nx})
      s = load([datpaths{nx},savedat,'.mat']);
      rawdat = floatread([datpaths{nx},savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);
      condidx = s.keeptrack(:,1)';   
      for ic = 1:length(clustcps{nx})
        ric = find(ismember(s.complist,clustcps{nx}(ic)));
        cpidx = [length(s.freqs)*(ric-1)+1:length(s.freqs)*ric];% one ic
        if strcmp(plotopt,'on')
          if place > row*col
            figure; place = 1;
          end;
        end;
        pl = 1;      
        pnts = find(condidx == condvec(cond));
        plotpnts = unique(s.keeptrack(pnts,3)); % time points to plot
        if strcmp(plttype,'time')
          %??
        elseif strcmp(plttype,'ersp')
          reersp = rawdat(pnts,cpidx);
          reersp = reshape(reersp,[size(reersp,1)/length(plotpnts) length(plotpnts) size(reersp,2) ]);
          if strcmp(mask,'mask') % create sig mask
            fprintf('\nCollecting bootstrap statistics for condition %s...\n',int2str(cond));
            [reersp] = GroupSig(reersp,[],.01,'permt');
          else
            reersp = mean(reersp,1); % take the mean!
            reersp = squeeze(reersp)';
          end;
        end;
        meanersp{cond}(:,:,end+1) = reersp;
        if isempty(clim)
          clim = max(max(abs(reersp)));
          if clim == 0
            clim = 1;
          end;
        end;
        if strcmp(plotopt,'on')
          if strcmp(plttype,'ersp')
            sbplot(row,col,place);place = place+1;
            if strcmp(s.freqscale,'quad')
              [realy,labely] = quadimagesc(plotpnts,s.freqs,reersp,[-clim clim]); hold on;
            elseif strcmp(s.freqscale,'log')
              [realy,labely] = mylogimagesc(plotpnts',s.freqs,reersp,[-clim clim]); hold on;
            else % linear
              imagesc(plotpnts',s.freqs,reersp,[-clim clim]); hold on;
              set(gca,'ydir','norm');
            end;
            ph = plot([0 0],[get(gca,'ylim')],'k-');
          elseif strcmp(plttype,'time')% plot over-plotted IM time courses
              %??
          end;
          title(['S ',int2str(nx),'; IC ',int2str(clustcps{nx}(ic))]);
        end;
      end;
    end;
  end;
  if strcmp(plotopt,'on')
    sbplot(row,col,place);
    if strcmp(s.freqscale,'quad')
      [realy,labely] = quadimagesc(plotpnts,s.freqs,mean(meanersp{cond},3),[-clim/1.5 clim/1.5]); hold on;
    elseif strcmp(s.freqscale,'log')
      [realy,labely] = mylogimagesc(plotpnts',s.freqs,mean(meanersp{cond},3),[-clim/1.5 clim/1.5]); hold on;
    else % linear
      imagesc(plotpnts',s.freqs,mean(meanersp{cond},3),[-clim/1.5 clim/1.5]); hold on;
      set(gca,'ydir','norm');
    end;
    ph = plot([0 0],[get(gca,'ylim')],'k-');
    title(['Average Clust ERSP']);
    textsc(['Condition ',s.datsets{cond}],'title');
  end;
end;
%plot difference if requested
if strcmp(plotdiff,'on') & length(meanersp) >1
  differsps = zeros(length(s.freqs),length(s.tmpoints),0);
  if strcmp(plotopt,'on')
  figure; place = 1;
  end
  for ic = 1:size(meanersp{1},3)
    differsps(:,:,end+1) = meanersp{1}(:,:,ic) - meanersp{2}(:,:,ic);
    if strcmp(plotopt,'on')
      sbplot(row,col,place);place = place+1;
      if strcmp(s.freqscale,'quad')
        [realy,labely] = quadimagesc(plotpnts,s.freqs,differsps(:,:,end),[-clim/2 clim/2]); hold on;
      elseif strcmp(s.freqscale,'log')
        [realy,labely] = mylogimagesc(plotpnts',s.freqs,differsps(:,:,end),[-clim/2 clim/2]); hold on;
      else % linear
        imagesc(plotpnts',s.freqs,differsps(:,:,end),[-clim/2 clim/2]); hold on;
        set(gca,'ydir','norm');
      end;
      ph = plot([0 0],[get(gca,'ylim')],'k-');
    end;
  end;
  if strcmp(plotopt,'on')
    sbplot(row,col,place);
    if strcmp(s.freqscale,'quad')
      [realy,labely] = quadimagesc(plotpnts,s.freqs,mean(differsps,3),[-clim/3 clim/3]); hold on;
    elseif strcmp(s.freqscale,'log')
      [realy,labely] = mylogimagesc(plotpnts',s.freqs,mean(differsps,3),[-clim/3 clim/3]); hold on;
    else % linear
      imagesc(plotpnts',s.freqs,mean(differsps,3),[-clim/3 clim/3]); hold on;
      set(gca,'ydir','norm');
    end;
    ph = plot([0 0],[get(gca,'ylim')],'k-');
    title(['Average Clust Difference']);
    textsc([s.datsets{1},' minus ',s.datsets{2}],'title');
  end;
end;
