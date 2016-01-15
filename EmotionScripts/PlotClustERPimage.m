% Plot all activations from all ICs in a given cluster as a single ERPimage
%
% [outdata,outvar] = PlotClustERPimage(plotcps,filenames,fullpaths,tmlims,frqlims,smoothfac,sortby,plotevents,multfac);
%
% smoothfac -- [number] divide the number of trials by this number to 
%              determine smoothing factor for each panel. (default = 15)
% sortby -- ['var' or 'time'] 'var' will sort ERP image by the latency
%           of the variable listed first in 'plotevents' argument. 
%           'time' does not sort ERP images and leaves as time-on-trial
%           and blocked by subjects
% plotevents -- [cell array] list of strings giving event type(s) to plot
%               If more than one string, first is used as sortvar (if sortby
%               is set to 'var' instead of 'time'
% rev = [1 or -1] if -1, will flip all activations to opposite polarity
 

function [outdata,outvar,outtimes,keepidxs] = PlotClustERPimage(plotcps,filenames,fullpaths,tmlims,frqlims,smoothfac,sortby,plotevents,multfac,rev,priorev);

clim = [];
if smoothfac == 1 % no smoothing
  nosmooth = 1;
else
  nosmooth = 0;
end;

if isempty(priorev)
  priorev = 1;
end;
if isempty(multfac)
  multfac = 1;
end;
if isempty(tmlims)
  tmlims = [-1000 3000];
end;
if isempty(smoothfac)
  smoothfac = 15;
end;
if iscell(plotcps) % multiple groups, plot vs each other
  row = length(filenames); % conds as rows
  col = length(plotcps); % groups as columns
  if length(filenames)>2
    twofigs = 1;row = 2;
  else
    twofigs = 0; 
  end;  
  figure; pl = 1;
  for f = 1:length(filenames)
    if twofigs == 1 && f > 2
      figure; pl = 1;
    end;
    clear keepacts keepgrids keeplats keepidxs
    for g = 1:length(plotcps)
      allacts = [];clustlats = [];grids = zeros(67*67,0);subj=1;subjidx = [];
      for nx = 1:length(plotcps{g})
        if ~isempty(plotcps{g}{nx})
          ALLEEG=[];EEG=[];
          EEG = pop_loadset( 'filename', filenames{f}, 'filepath', fullpaths{nx});
          tms = find(EEG.times>tmlims(1)&EEG.times<tmlims(2));
          outtimes = EEG.times(tms);
          fnames=fieldnames(EEG.epoch(1)); % field names are constant
          for ic = 1:length(plotcps{g}{nx})
            [h grid]= topoplot(EEG.icawinv(:,plotcps{g}{nx}(ic)),EEG.chanlocs,'electrodes','off','noplot','on');
            grids(:,end+1) = reshape(grid,[67*67 1]);       
            if ~isempty(plotevents)
              for ev = 1:length(plotevents)
                evs(1,ev) = find(ismember(fnames,['event',plotevents{ev}]));
              end;
              clear evlats
              for ep = 1:length(EEG.epoch)
                zlat = find(ismember(cell2mat(EEG.epoch(ep).eventlatency),0));
                for ev = 1:length(plotevents)
                  tmpev = getfield(EEG.epoch(ep),['event',plotevents{ev}]);
                  tp = tmpev{zlat};
                  if isstr(tp)
                    tp = str2num(tp);
                  end;    
                  evlats(ev,ep) = tp*priorev; % -1 = latency before time-lock event
                end;
              end;  
              clustlats = [clustlats,evlats];
            else
              clustlats =  [];
            end;
            allacts = [allacts,squeeze(EEG.icaact(plotcps{g}{nx}(ic),tms,:))];% tms x trials
            subjidx = [subjidx,ones(1,size(EEG.icaact,3))*subj];subj = subj+1;
          end;
        end;
      end;
      keepacts{g} = allacts;
      keepgrids{g} = grids;
      keeplats{g} = clustlats;
      keepidxs{g} = subjidx;
    end;
    % Bring ALL (all groups) ICs into same polarity:--------
    [compin, polarities] = std_comppol(cat(2,keepgrids{:}));
    pos = 0;
    polarities = polarities*rev; % if rev is -1, will flip polarity
    for g = 1:length(keepgrids)
      if ~isempty(keepacts{g})
        polfac = [];
        for s = 1:max(keepidxs{g})
          polfac = [polfac,repmat(polarities(pos+s),[1 length(find(keepidxs{g}==s))])];
        end; pos = pos+s;
        polfac = repmat(polfac,[size(keepacts{g},1) 1]);
        clustlats = keeplats{g}; 
        allacts = keepacts{g}.*polfac; clear polfac 
        if nosmooth == 1% no smoothing
          smoothfac = size(allacts,2);nosmooth = 1;smoothfac
        end;
        % Now plot the activations or power:--------
        sbplot(row,col,pl);pl=pl+1;
        if isempty(frqlims) % activations
          if isempty(clim)
            clim = mean(iqr(mean(allacts,3)))-(.3*mean(iqr(mean(allacts,3))));
            if size(allacts,2) > 3500
              clim = clim - clim*.10; % large # trials
            end;
          end;
          if strcmp(sortby,'time')
            if ~isempty(clustlats)
              [outdata{g},outvar{g},outtrials] = erpimage(allacts,[], EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','auxvar',[multfac*clustlats'],'caxis',[-clim clim]); 
            else % don't plot aux
              [outdata{g},outvar{g},outtrials] = erpimage(allacts,[], EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','caxis',[-clim clim]); 
            end;
          else % sort by first event type
            [outdata{g},outvar{g},outtrials] = erpimage(allacts,clustlats(1,:)*multfac, EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','auxvar',[multfac*clustlats(2:end,:)'],'caxis',[-clim clim]);
          end;
        else % power
          if strcmp(sortby,'time')
            if ~isempty(clustlats)
              [outdata{g},outvar{g},outtrials] = erpimage(allacts,[], EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps','auxvar',[multfac*clustlats']);
            else % don't plot aux
              [outdata{g},outvar{g},outtrials] = erpimage(allacts,[], EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps');
            end;
          else % sort by first event type
            [outdata{g},outvar{g},outtrials] = erpimage(allacts,clustlats(1,:)*multfac, EEG.times(tms), '', size(allacts,2)/smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps','auxvar',[multfac*clustlats(2:end,:)']);
          end;
        end;
        title(['Group ',int2str(g),' (',int2str(max(keepidxs{g})),'); ',filenames{f}(1:end-4)]);
        if length(outtrials) == length(keepidxs{g}) % resort if not smoothed
          keepidxs{g} = keepidxs{g}(outtrials);
        end;
      else
        pl=pl+1; % keep empty for consistency
      end;
    end;
  end;  
else % all one group-------------------------------------------------------------
  row = 1; 
  col = length(filenames);
  if length(filenames)>3
    twofigs = 1;col = 3;
  else
    twofigs = 0;
  end;
  
  figure; pl = 1;
  for f = 1:length(filenames)
    if twofigs == 1 && f > 2
      figure; pl = 1;
    end;
    allacts = [];clustlats = [];grids = zeros(67*67,0);subj=1;subjidx = [];
    for nx = 1:length(plotcps)
      if ~isempty(plotcps{nx})
        ALLEEG=[];EEG=[];
        EEG = pop_loadset( 'filename', filenames{f}, 'filepath', fullpaths{nx});
        tms = find(EEG.times>tmlims(1)&EEG.times<tmlims(2));
        outtimes = EEG.times(tms);
       fnames=fieldnames(EEG.epoch(1)); % field names are constant
        for ic = 1:length(plotcps{nx})
          [h grid]= topoplot(EEG.icawinv(:,plotcps{nx}(ic)),EEG.chanlocs,'electrodes','off','noplot','on');
          grids(:,end+1) = reshape(grid,[67*67 1]);        
          for ev = 1:length(plotevents)
            evs(1,ev) = find(ismember(fnames,['event',plotevents{ev}]));
          end;
          clear evlats
          for ep = 1:length(EEG.epoch)
            zlat = find(ismember(cell2mat(EEG.epoch(ep).eventlatency),0));
            for ev = 1:length(plotevents)
              tmpev = getfield(EEG.epoch(ep),['event',plotevents{ev}]);
              evlats(ev,ep) = tmpev{zlat};
            end;
          end;  
          clustlats = [clustlats,evlats];
          allacts = [allacts,squeeze(EEG.icaact(plotcps{nx}(ic),tms,:))];% tms x trials
          subjidx = [subjidx,ones(1,size(EEG.icaact,3))*subj];subj = subj+1;
        end;
      end;
    end;
    % Bring all ICs into same polarity:--------
    [compin, polarities] = std_comppol(grids);
    polfac = [];
    for s = 1:max(subjidx)
      polfac = [polfac,repmat(polarities(s),[1 length(find(subjidx==s))])];
    end;
    polfac = repmat(polfac,[size(allacts,1) 1]);
    allacts = allacts.*polfac; clear polfac polarities
    if smoothfac == 1 | nosmooth == 1% no smoothing
      smoothfac = size(allacts,2);nosmooth = 1;smoothfac
    end;
    
    % Now plot the activations or power:--------
    sbplot(row,col,pl);pl=pl+1;
    if isempty(frqlims) % activations
      if isempty(clim)
        clim = mean(iqr(allacts'));
      end;
      if strcmp(sortby,'time')
        if ~isempty(clustlats)
          [outdata{f},outvar{f}] = erpimage(allacts,[], EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','auxvar',[clustlats'],'caxis',clim); 
        else % don't plot aux
          [outdata{f},outvar{f}] = erpimage(allacts,[], EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','caxis',clim); 
        end;
      else % sort by first event
        plottrials = find(clustlats(1,:) ~= 0);
        [outdata{f},outvar{f}] = erpimage(allacts(:,plottrials),clustlats(1,plottrials)*multfac, EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','auxvar',[multfac*clustlats(:,plottrials)'],'caxis',[-clim clim]);
      end;
    else % power
      if strcmp(sortby,'time')
        if ~isempty(clustlats)
          [outdata{f},outvar{f}] = erpimage(allacts,[], EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps','auxvar',[multfac*clustlats']);
        else % don't plot aux
          [outdata{f},outvar{f}] = erpimage(allacts,[], EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps');
        end;
      else % sort by first event
        plottrials = find(clustlats(1,:) ~= 0);
        [outdata{f},outvar{f}] = erpimage(allacts,clustlats(1,plottrials)*multfac, EEG.times(tms), '', smoothfac, 1 ,'yerplabel','','erp','on','cbar','on','coher',[frqlims .01],'plotamps','auxvar',[multfac*clustlats(:,plottrials)']);
      end;
    end;
    title(EEG.setname);
    if length(outtrials) == length(keepidxs{g}) % resort if not smoothed
      keepidxs{f} = subjidx(outtrials);
    end;
  end;
end;
%,'timewarp',{[1000*clustlats(1,:)',1000*clustlats(2,:)'],[2000 3500],{'c','m'}}