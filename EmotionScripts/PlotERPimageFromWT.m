% Uses WT spectral decomp and reconstructs an ERP image by
%
% [outersp,keeptrack,plotpnts] = PlotERPimageFromWT(savedat,datset,datpath,plotdims,mask,allbigs,plotopt);
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

function PlotERPimageFromWT(savedat,datpath,plotdims,condvec,plotevents);

smoothfac = 10;
multfac = 1000; % sec to ms

rmpath('/data/MATLAB/binica/');rmpath('/data/MATLAB/binica_SCCN/')
s = load([datpath,savedat,'.mat']);
sph=floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0);
wts=floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);
icamatall = floatread([datpath,savedat,'.fdt'],[s.pcs s.numframes],[],0);
% pcamat = floatread([datpath,savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);
ws = wts*sph;    acts = ws*icamatall;    winv = inv(ws);
clear wts sph
speceig = floatread([datpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);
specwts = speceig*winv;
winv = specwts; clear speceig specwts  icamatall
%ws = pinv(winv);

if isempty(plotdims)
   plotdims = [1:size(winv,2)];
end;
plotpnts = unique(s.keeptrack(:,3)); % time points to plot
condidx = s.keeptrack(:,1)'; 
nepochs = length(condidx)/length(plotpnts);
epids = reshape(condidx,[length(condidx)/length(plotpnts) length(plotpnts)]);
nowcond = 1;
newevents = cell(1,max(condidx));
for ep = 1:nepochs
  if epids(ep,1) == nowcond
    newevents{nowcond}(end+1) = s.tackevents(ep);
  else
    nowcond = nowcond+1;
    newevents{nowcond}(end+1) = s.tackevents(ep); 
  end;
end;

row = round(sqrt(length(plotdims)));
col = ceil(sqrt(length(plotdims)));
if length(condvec) == 1 % one condition
  if row>5
    row = 4;col = 5;
  end;
elseif length(condvec) > 1 % more than one condition
  if ~iseven(col)
    col = col-1; row = col+1;
    if col>6
      row = 5;col = 6;
    end;
  end;
end;

figure;place = 1;
for dm = 1:length(plotdims)
  dim = plotdims(dm);
  if place > row*col
    figure; place = 1;
  end;
  for c = 1:length(condvec)
    cond = condvec(c);
    
    pnts = find(condidx == cond);
    imgpnts = winv(pnts,dim);
    imgpnts = reshape(imgpnts,[size(imgpnts,1)/length(plotpnts) length(plotpnts) ]);
%    sbplot(row,col,place);place=place+1;imagesc(s.tmpoints,[1:size(imgpnts,1)],imgpnts); hold on;plot([0 0],[get(gca,'ylim')],'k-');
    
    fnames=fieldnames(newevents{cond}(1)); % field names are constant
    for ev = 1:length(plotevents)
      evs(1,ev) = find(ismember(fnames,plotevents{ev}));
    end;
    clear evlats
    for ep = 1:length(newevents{cond})
      for ev = 1:length(plotevents)
        tmpev = getfield(newevents{cond}(ep),plotevents{ev});
        if isstr(tmpev)
          tmpev = str2num(tmpev);
        end;    
        evlats(ev,ep) = tmpev;
      end;
    end;  
    clim = iqr(imgpnts(:)) ;
    sbplot(row,col,place); place=place+1;    erpimage(imgpnts',evlats(1,:)*multfac, plotpnts, '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','caxis',[-clim clim],'auxvar',[multfac*evlats']);
    %sbplot(row,col,place); place=place+1;    erpimage(imgpnts',evlats(1,:)*multfac, plotpnts, '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','caxis',[-clim clim]);
    title(['Dim ',int2str(dim),'; ',s.datsets{cond}]);
  end;
end;

     