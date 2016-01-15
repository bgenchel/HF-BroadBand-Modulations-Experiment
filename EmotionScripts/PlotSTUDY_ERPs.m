% Plots ERPs from STUDY
% overplots 'conds'
% if groups is not [], then overplots groups (currently disabled)
% 
function PlotSTUDY_ERPs(STUDY,ALLEEG,clust,groups,conds,chan,shuffnum,maskalpha,tlims,row,col,place)

for cnd = 1:length(conds)
  cond = conds(cnd);
  mnrts = [];clear acts chans
  for ic = 1:size(STUDY.cluster(clust).sets,2)    
    setidx = STUDY.cluster(clust).sets(cond,ic); 
    comp = STUDY.cluster(clust).comps(ic);
    CURRENTSTUDY = 0;[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, [1:length(STUDY.datasetinfo)] ,'retrieve',setidx,'study',1); 
     
    rts=[];
    for ep = 1:length(EEG.epoch)
        zpos = find(ismember(cell2mat(EEG.epoch(ep).eventlatency),0));
        
        if EEG.epoch(ep).eventtype{zpos-1} == 1|EEG.epoch(ep).eventtype{zpos-1} == 2
            rts = [rts EEG.epoch(ep).eventlatency{zpos-1}-EEG.epoch(ep).eventlatency{zpos-2}];
        else
            rts = [rts 0];% if no rt, make rt 0 ms
        end;
    end;   
    rts(find(rts==0))=[]; % delete no resp trials
    mnrts = [mnrts mean(rts)];
    rmsuv = sqrt(mean(EEG.icawinv(:,comp).^2));% RMS at scalp
    % collect winv, oriented by 'topopol' field in STUDY:
    %winvs{ic} = EEG.icawinv(:,comp)*STUDY.cluster(clust).topopol(ic);
    % collect activations, oriented by 'topopol' SAME field in STUDY:
    acts{ic} = EEG.icaact(comp,:)*rmsuv*STUDY.cluster(clust).topopol(ic);
    %chans{ic} = EEG.chanlocs; 
  end;  
  tms = EEG.times; % should be same for all
  condrts{cond} = mnrts;
  condacts{cond} = acts;
end;    
    
    
    