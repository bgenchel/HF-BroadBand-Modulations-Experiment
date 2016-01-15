% plot just the max activation vs latency for each trial
% in a dataset of ICs scaled by scalp RMS
% proj -- ['rms' or channel label]


function PlotMaxActVsLatERP(datset,fullpaths,complist,tmlims,bsln,rejs,proj);
allmax=[];alllat=[];
for nx = 1:length(complist)
  if ~isempty(complist{nx})
    EEG = pop_loadset('filename',datset,'filepath',fullpaths{nx});nchan=size(EEG.data,1);
    EEG = pop_rmbase( EEG,bsln );
    EEG = pop_select(EEG,'notrial',rejs{nx});% delete specified epochs
    tms = find(EEG.times>tmlims(1) & EEG.times < tmlims(2));
    clear dat
    for ic = 1:length(complist{nx})
      if strcmp(proj,'rms')
        fprintf('\nUsing full scalp RMS\n');
      rms = sqrt(mean(EEG.icawinv(:,complist{nx}(ic)).^2));
      else
        fprintf('\nProjecting to channel %s\n',proj);
        ch = find(ismember({EEG.chanlocs.labels},proj));
        rms = EEG.icawinv(ch,complist{nx}(ic));
      end;
      dat = EEG.icaact(complist{nx}(ic),tms,:)*rms;
      %tmpdat = EEG.icaact(complist(ic),tms,:)*rms;
      %dat(ic,:,:) = tmpdat - repmat(mean(tmpdat,2),[1 size(tmpdat,2) 1]);
      for ep = 1:size(dat,3)
        [tt maxlat(1,ep)] = max(abs(dat(1,:,ep)));
        trialmax(1,ep) = dat(1,maxlat(1,ep),ep); % real value (pos or neg)
      end;
      allmax = [allmax,trialmax];
      alllat = [alllat,maxlat];
    end;
  end;
end;
    
figure; 
  plot(alllat,allmax,'.','markersize',25);hold on;
    