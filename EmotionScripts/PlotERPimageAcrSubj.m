% for P1 cluster analysis. Returns RTs and amplitudes
%
%      [allerps,allRTs,times] = PlotERPimageAcrSubj(paths,datset,datpaths, clustcps,tlims,flims,memtarg,elecs,ttl,ampst,phsest,vsort,rtsort,pwr,phs);
%
% paths -- [cell array of strings] fullpath for each subject where datset can be found
% datset -- [string] name of dataset to plot (must be same for all subjects)
% datpaths -- [cell array of strings] fullpath for each subject where baseline power data can be found
% tlims -- [mintime maxtime] in ms, time range of erp to plot and save
% flims -- {[minfreq maxfreq]} cell array of freq ranges to save data for using time domain time intervals
% memtarg -- [0|1] 1 to only choose memorize trials that are subsequent targets, else 0
% elecs -- cell array of strings corresponding to electrodes in which to search for highest variance for comp acts
% ttl -- cluster name to add to each plot
% ampst -- [latency in ms] if not empty, will plot an amplitude sorted plot in the flim range at specified time
% phsest -- [latency in ms] if not empty, will plot a phase-sorted plot in the flim range at specified time
% vsort -- [minms maxms] if not empty, will sort trials based on activtion value during specified time interval.
% rtsort -- [0|1] if 1, will plot an RT-sorted erpimage of all trials; requires a field in epoch-centered event called 'rt' which specifies reaction time for that trial
% pwr --  [0|1] if 1, will plot amplitudes of specified freq range, RMS activation amplitude otherwise
%               pwr is baselined using tasc_baseline.mat from subj directory.
% phs -- [integer between -180 and 180] Phase angle to start phase-sorted subplot
%

function  [allerps,allRTs,times] = PlotERPimageAcrSubj(paths,datset,datpaths,clustcps,tlims,flims,memtarg,elecs,ttl,ampst,phsest,vsort,rtsort,pwr,phs); 
    
    row = 2; col = 2;
    EEG = pop_loadset( datset,paths{1});        
    mntime = find(EEG.times>tlims(1) & EEG.times< tlims(2)); % time period to save mean erp for
    clear erpamps1 erpamps2 
    ALLEEG=[];EEG=[];allerps = zeros(length(mntime),0); allRTs = zeros(1,0);
    for nx = 1:length(clustcps)
        if ~isempty(clustcps{nx})
            ALLEEG=[];clear clusterps erps
            EEG = pop_loadset( datset,paths{nx});      ALLEEG=EEG; 
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,1);
            EEG = pop_rmbase( EEG, [tlims(1) 0]);
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,1);
            tm = find(EEG.times> tlims(1) & EEG.times < 0); % baseline for std 
            multfact = sqrt(mean(EEG.icawinv(:,clustcps{nx}).^2));
            clear acts
            w=1;clear picktrials RTs
            if memtarg == 1
                for ep = 1:length(EEG.epoch)
                    for ev = 1:length(EEG.epoch(ep).event)
                        if cell2mat(EEG.epoch(ep).eventinset(ev)) == 'y' & cell2mat(EEG.epoch(ep).eventlatency(ev)) == 0
                            RTs(1,w) = cell2mat(EEG.epoch(ep).eventrt(ev));
                            picktrials(1,w) = ep; w = w+1;
                        end;
                    end;
                end;
            else
                for ev = 1:length(EEG.event)
                    if EEG.event(ev).rt > 0
                        RTs(1,w) = EEG.event(ev).rt; w = w+1;
                    end;                 
                end;
                picktrials = [1:size(EEG.data,3)];
            end;     
            allRTs(1,end+1:end+length(RTs)) = RTs;
            % makes tm X epochs made into uV by rms of scalp map 
            % plot average erp with overlaying peak alpha freq sine wave
            oneerp =  squeeze(EEG.icaact(clustcps{nx},mntime,picktrials))*multfact;
            
            %%%%%%%%%%%%%%%%%%%%  find absolute orientation  %%%%%%%%%%%%%%%
            q=1;clear allelecs
            for d = 1:length(elecs)
                for c = 1:length(EEG.chanlocs)
                    if strcmp(EEG.chanlocs(c).labels,elecs{d})
                        allelecs(1,q) = c;q=q+1;
                    end;
                end;
            end;            
            topproj = EEG.icawinv(allelecs,clustcps{nx});
            [x y] = max(abs(topproj)); useelec = allelecs(y);
            compset = [1:size(EEG.data,1)]; compset(clustcps{nx}) = [];
            EEG = pop_subcomp( EEG,compset , 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'projected set');
            erps(1,:) = mean(EEG.data(useelec,find(EEG.times>0&EEG.times<600),:),3);
            erps(1,:) = erps(1,:)/std(erps(1,:));
            EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
            erps(2,:) = mean(EEG.icaact(clustcps{nx},find(EEG.times>0&EEG.times<600),:),3);
            erps(2,:) = erps(2,:)/std(erps(2,:));
            [corr,indx,indy,corrs] = matcorr(erps(1,:),erps(2,:));
            if corr < 0
                oneerp = oneerp*-1;
            end;  
            ALLEEG = pop_delset( ALLEEG, [2] );        clear erps
            allerps(:,end+1:end+size(oneerp,2)) = oneerp;
        end;
    end;
    %%%%%%%%%%%%%%%%%%%%
    if ~exist('phs')
        phs = 0;
    end;   
    divby = size(allerps,2);
    figure;
    if pwr == 1
        if ~isempty(ampst)
            subplot(row,col,1);        [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs, EEG.times(mntime),['Amp-sorted Activations (',int2str(flims(1)),'-',int2str(flims(2)),'Hz at ',int2str(ampst),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN  0 .9  NaN],'erp','noxlabel','ampsort',[ampst 0 flims(1) flims(2)],'coher',[flims(1) flims(2) .01],'plotamps','auxvar',allRTs*1000,'cbar');     
        end;
        if ~isempty(phsest)
            subplot(row,col,2);   [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs,  EEG.times(mntime),['Phase-sorted Activations (',int2str(flims(1)),'-',int2str(flims(2)),'Hz at ',int2str(ampst),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9 NaN],'erp','noxlabel','phasesort',[phsest 0 flims(1) flims(2) phs],'auxvar',allRTs*1000,'cbar','coher',[flims(1) flims(2) .01],'plotamps'); 
        end;
        if ~isempty(vsort)
            subplot(row,col,3);    [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs,  EEG.times(mntime),['Value-sorted Activations (',int2str(vsort(1)),'-',int2str(vsort(2)),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9 NaN],'erp','noxlabel','valsort' , [vsort(1) vsort(2) 1],'auxvar',allRTs*1000,'cbar','coher',[flims(1) flims(2) .01],'plotamps'); 
        end;
        if rtsort == 1
            subplot(row,col,4);    [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs*1000,EEG.times(mntime),['RT-sorted Activations '],size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9  NaN],'erp','noxlabel','cbar','coher',[flims(1) flims(2) .01],'plotamps'); 
        end;
    else
        if ~isempty(ampst)
            subplot(row,col,1);    [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs, EEG.times(mntime),['Amp-sorted Activations (',int2str(flims(1)),'-',int2str(flims(2)),'Hz at ',int2str(ampst),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN  0 .9  NaN],'erp','noxlabel','ampsort',[ampst 0 flims(1) flims(2)],'auxvar',allRTs*1000,'cbar');     
        end;
        if ~isempty(phsest)
            subplot(row,col,2);     [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs,  EEG.times(mntime),['Phase-sorted Activations (',int2str(flims(1)),'-',int2str(flims(2)),'Hz at ',int2str(ampst),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9 NaN],'erp','noxlabel','phasesort',[phsest 0 flims(1) flims(2) phs],'auxvar',allRTs*1000,'cbar'); 
        end;
        if ~isempty(vsort)
            subplot(row,col,3);      [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs,  EEG.times(mntime),['Value-sorted Activations (',int2str(vsort(1)),'-',int2str(vsort(2)),'ms) '] ,size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9 NaN],'erp','noxlabel','valsort' , [vsort(1) vsort(2) 1],'auxvar',allRTs*1000,'cbar'); 
        end;
        if rtsort == 1
            subplot(row,col,4);     [outdata1,outvar1,outtrials,limits,axhndls,erp1,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage( allerps,allRTs*1000,EEG.times(mntime),['RT-sorted Activations '],size(allerps,2)/divby,0,'limits',[tlims(1) tlims(2) NaN NaN  NaN NaN 0 .9  NaN],'erp','noxlabel','cbar'); 
        end;    
    end;    
    ph = textsc(ttl,'title');  set(ph,'fontsize',14);
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    
    times = EEG.times(mntime);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
