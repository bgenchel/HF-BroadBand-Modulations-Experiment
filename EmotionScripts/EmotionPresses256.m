% this is a script to plot each button press trial on one page
%expt for jo74 on Jan 14th(?) new emotions, somewhat:
% awe;  frustration; joy;   anger; sad; happy; fear; love; jealousy; compassion; emabarrass;  content; grief; relief;
 
% import with event channel 265, then keep only channel 257 (button press)

EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [1:1099:1100] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1100');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [1101:1099:2200] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'emo 1101 2200');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [2201:1099:3300] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2201-end');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% merge and import log file
EEG = pop_mergeset( ALLEEG, [1  2  3], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole trial', 'overwrite', 'on');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',0);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'ButtonOnly.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

figure;plot(EEG.data); %to make sure you got the right thing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import emotion data in three chunks, need to merge log file each time by shifts.
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [1:1099:1100] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1100');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',0);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%  block 2   %%%%%%%%%
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [1077:876:2100] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'jo74 emo2 1077 2100');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-34);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%  block 3  %%%%%%%
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [2002:1137:3140] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 2002:3140');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-54);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to get last eyes open/closed
EEG = pop_readbdf('/data/common1/emotion/jo74/jo_emo.bdf', [2800:382:3183] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'last part');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-71);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/jo74/dummy-small-256.elp', 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 eeglab redraw

EEG = pop_saveset( EEG, 'emo-1-252.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

 % write into data info:
 jo74  - 29 yo female
data collected on Jan 14, 2004  using 256 (265 including analog) Biosemi acquisition system
First 1100 blocks of longer data file
merged with presentation log file (events aligned correctly and did not try to squeeze in all events that were not called in in these blocks)
filtered 1 50

Data  taken after emotion and Reward twoback expts. Used eyes_new.sce to get 2 min each of eyes open and closed. This dataset is only the open epochs. 


% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
emos1 = {'awe', 'frustration','joy','anger','sad'};
emos2 = {'happy','fear','love'};% ,'jealousy'}; had to be imported separately because spans gap
emos3 = {'compassion','emabarrass','content','grief','relief'};
for k = 1:length(emos1)
    EEG = pop_epoch( EEG, {  emos1{k}  }, [-.5  300], 'newname', emos1{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emos1{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:EEG.pnts); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,emos1{k} , '/data/common1/emotion/jo74/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 eeglab redraw
for k = 1:length(emos2)
    EEG = pop_epoch( EEG, {  emos2{k}  }, [-.5  300], 'newname', emos2{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emos2{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:EEG.pnts); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,emos2{k} , '/data/common1/emotion/jo74/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 eeglab redraw
for k = 1:length(emos3)
    EEG = pop_epoch( EEG, {  emos3{k}  }, [-.5  300], 'newname', emos3{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emos3{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:EEG.pnts); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,emos3{k} , '/data/common1/emotion/jo74/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 eeglab redraw
% for the last, you need to change length to 250 instead of 300
% select out data after stoprelax
emoset = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};
for k = 1:length(emoset)
EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    for ev = 1:length(EEG.event)
        if EEG.event(ev).Event_Type(1) == 'S'
            seltime = EEG.event(ev).latency;
        end;
    end;
    EEG = pop_select( EEG, 'point',[1 seltime] );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,emoset{k} , '/data/common1/emotion/jo74/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG=[];
end;
 eeglab redraw

% use events from buttonpressonly file to epoch real data 
eeglab
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bemabarrass','bcontent','bgrief','brelief'};
 EEG = pop_loadset( 'ButtonOnly.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

for k = 1:length(emos)
    if k == 14
        endt = 250;
    else
        endt=300;
    end;    
    EEG = pop_epoch( EEG, {  emos{k}  }, [-.5  endt], 'newname', emos{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emobuts{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:512); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 ALLEEG(1)=[];
%figure;plot(EEG.data); %to make sure you got the right thing
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    
    for ev = 1:length(EEG.event)
        if EEG.event(ev).Event_Type(1) == 'S'
            seltime = EEG.event(ev).latency;
        end;
    end;
    EEG = pop_select( EEG, 'point',[1 seltime] );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, k);
end;
%figure;plot(EEG.data); %to make sure you got the right thing
clear pressevents
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    maxdata = max(EEG.data);
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 1: size(EEG.data,2)
        if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press';
            EEG.event(end).Event_Type = 'Response';
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;


% to epoch on the button press only data:
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','emabarrassPress','contentPress','griefPress','reliefPress'};
for k = 1:length(emonames)   
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    EEG = pop_epoch( EEG,{'press'} , [-2  3], 'newname', emonames{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emonames{k},'overwrite','on');
    EEG = pop_rmbase( EEG,[0 100]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
end;
eeglab redraw

figure; 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    %diffpress = diff(diff(EEG.data));
    %diffpress = diff(EEG.data);
    diffpress = EEG.data;
    diffpress = squeeze(diffpress);
    subplot(4,4,k)
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 2500 NaN NaN NaN NaN NaN NaN],'caxis',[-max(max(diffpress(501:600,:)))-100 max(max(diffpress(501:600,:)))-100]);
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 550 NaN NaN  NaN NaN NaN NaN],'caxis',[-8500 8500]);
 erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 750  NaN NaN  NaN NaN NaN NaN]);
end;
figure;
%plset = [1,5,6,8,10,12,13,14];
plset = [2,3,4,7,9,11];
for k= 1:length(plset)
    subplot(2,3,k)
    k = plset(k);
     EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
     diffpress = diff(diff(EEG.data));
     %diffpress = diff(EEG.data);
     %diffpress = EEG.data;
     diffpress = squeeze(diffpress);
  erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 400  NaN NaN  NaN NaN NaN NaN],'caxis',[-15000 15000]);
end;
  erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 750 -1500000 1500000  NaN NaN NaN NaN]);

% Now create events at button presses looking at actual button presses
emoset = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','emabarrassPress','contentPress','griefPress','reliefPress'};
ALLEEG=[];
% This to epoch actual data based on the above events
for k = 1:length(emonames)   
    EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    % add press events from button only (above)
    EEG.event(end+1:end+length(pressevents{k})) = pressevents{k};
    EEG = pop_epoch( EEG,{'press'} , [-2  3], 'newname', emonames{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emonames{k});
    EEG = pop_rmbase( EEG, [-2000 3000]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
    EEG = pop_saveset( EEG,emonames{k} , '/data/common1/emotion/jo74/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG=[];
end;
eeglab redraw

    


