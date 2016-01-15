% Runs timef on all 'prebase' datasets to obtain spectral baseline for future newtimef computations

% remember to:
cd /data/common2/emotion/
mkdir ersps

eeglab
comment = 'Baseline power spectrum from "prebase.set" of 100 freqs, padratio, 2, wavelet,[3 .5]. ';
figure;
for nx = 35:35%length(gdcomps)    
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 

    EEG = pop_loadset( 'prebase.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    EEG.icaweights=wts;    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for cmp = 1:length(gdcomps{nx})
          [ersp,itc,powbase,times,freqs,erspboot,itcboot]=newtimef(EEG.icaact(gdcomps{nx}(cmp),:), EEG.pnts,[0 EEG.xmax*1000], 256,[3 .5],'plotitc','off', 'alpha',.01,'freqs',[1:.5:50],'baseline',EEG.xmax*1000, 'verbose','off');
        baseline(gdcomps{nx}(cmp),:) = powbase;
    fprintf('\nOne More Component Done: %i of %i\n',cmp,length(gdcomps{nx}));  clf  
    end;
    cd (['/data/common2/emotion/',paths{nx},'/ersps/']);
    save baseline.mat baseline freqs comment
    ALLEEG=[];EEG=[];clear baseline 
    fprintf('\nOne More Subj Done: %i\n',nx);    
end;


figure;    [ersp,itc,powbase,times,freqs,erspboot,itcboot]=newtimef(EEG.icaact(gdcomps{nx}(cmp),:), EEG.pnts,[0 EEG.xmax*1000], 256,[3 .5],'plotitc','off','freqs',[1:.5:50],'baseline',50, 'verbose','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get baseline from emotion period datasets:

emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'};

comment = 'Baseline power spectrum from randomly epoched from continuous data from all emotions (after initial button press),excluding baselines, to get a mean baseline of 100 freqs, wavelet,[3 .5]. ';
figure;
for nx =35:35%length(gdcomps) 
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
    allacts = zeros(length(gdcomps{nx}),1024,0);
    for k = 1:length(emos)        
        EEG = pop_loadset( [emos{k},'.set'], ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        ft = EEG.event(2).latency;
        for evtm = ft+512:1024:size(EEG.data,2)-512  % go up by 1 sec to create non-overlapping 1 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';        
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-2 2]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
        EEG = pop_rmbase( EEG,[-2000 2000]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        % load weights
        EEG.icaweights=wts;    EEG.icasphere=sph;EEG.icawinv=[];
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_eegthresh(EEG,0,gdcomps{nx} ,-400,400,EEG.xmin,EEG.xmax,0,1);
        EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
        EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
        EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
        EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
        % save activations of interest
        allacts(:,:,end+1:end+size(EEG.icaact,3)) = EEG.icaact(gdcomps{nx},:,:); 
        ALLEEG=[];EEG=[];    
    end;     
    for cmp = 1:length(gdcomps{nx})
        [ersp,itc,baseline(gdcomps{nx}(cmp),:),times,freqs,erspboot,itcboot]=newtimef(allacts(cmp,:,:), size(allacts,2),[-2000 2000], 256,[3 .5],'plotitc','off', 'alpha',.01,'freqs',[1:.5:50], 'verbose','off','baseline',2000,'winsize',256);
        fprintf('\nOne More Component Done: %i of %i\n',cmp,length(gdcomps{nx}));  clf  
    end;
    cd (['/data/common2/emotion/',paths{nx},'/ersps/']);
    save EmoBaseline.mat baseline freqs comment
    ALLEEG=[];EEG=[];clear baseline 
    fprintf('\nOne More Subj Done: %i\n',nx);    
end;


figure;    [ersp,itc,powbase,times,freqs,erspboot,itcboot]=newtimef(EEG.icaact(gdcomps{nx}(cmp),:), EEG.pnts,[0 EEG.xmax*1000], 256,[3 .5],'plotitc','off','freqs',[1:.5:50],'baseline',50, 'verbose','off');
,'baseline',500