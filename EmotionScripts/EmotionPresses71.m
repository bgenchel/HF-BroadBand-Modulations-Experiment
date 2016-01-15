% finds press events from cont data and epochs on button presses

eeglab

emoset = {'awe.set', 'frust.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','content.set','grief.set','relief.set'};
datset = {'awepress.set','frustpress.set','joypress.set','angerpress.set','sadpress.set','happypress.set','fearpress.set','lovepress.set','jealousypress.set','compassionpress.set','contentpress.set','griefpress.set','reliefpress.set'};
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','contentPress','griefPress','reliefPress'};
paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/kl80/','/data/common1/emotion/rr80/'};
% use events from buttonpressonly file to epoch real data 
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief'};
for nx = 1:length(paths)
    for k = 1:length(emoset)
        EEG = pop_loadset( emoset{nx}, paths{nx});
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %EEG = pop_select( EEG, 'channel',72);
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');

        EEG.data = rmbase(EEG.data(72,:),EEG.pnts,1:512); 
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
      
        for ev = 1:length(EEG.event)
            if EEG.event(ev).stimulus(1) == 'R'
                seltime = EEG.event(ev).latency;
                break
            end;
        end;
        EEG = pop_select( EEG, 'nopoint',[seltime size(EEG.data,2)-1] );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG,CURRENTSET );

        %figure;plot(EEG.data); %to make sure you got the right thing
        maxdata = max(EEG.data);
        thresh = maxdata*.1;  % sets threshold at 10% of max value
        r=1; 
        for g = 1: size(EEG.data,2)
            if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
                EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
                EEG.event(end).type = 'press';
                EEG.event(end).stimulus = 'Response';
            end;    
        end;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
        pressevents{k} = EEG.event;
        % to epoch on the button press only data:
        EEG = pop_epoch( EEG,{'press'} , [-2  2], 'newname', emonames{k}, 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emonames{k},'overwrite','on');
        EEG = pop_rmbase( EEG,[0 100]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
        EEG = pop_saveset( EEG, datset{k}, paths{nx});
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        ALLEEG=[]; EEG=[];
    end;
end;

eeglab redraw











nx=2;
for k = 1:length(emoset)
EEG = pop_loadset( emoset{k}, paths{nx});
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG.data = rmbase(EEG.data,EEG.pnts,1:1024); 
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end;

clear pressevents
for k = 1:length(emoset)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    maxdata = max(EEG.data(72,:,:));
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 1: size(EEG.data,2)
        if EEG.data(72,g)>thresh & EEG.data(72,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press';
            EEG.event(end).Event_Type = 'Response';
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;
% to epoch on the button press only data:
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','contentPress','griefPress','reliefPress'};
for k = 1:length(emonames)   
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    EEG = pop_epoch( EEG,{'press'} , [-2  3], 'newname', emonames{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emonames{k},'overwrite','on');
    EEG = pop_rmbase( EEG,[0 100]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
end;
eeglab redraw


% done saving new button datasets. From now on call these in:
datset = {'awepress.set','frustpress.set','joypress.set','angerpress.set','sadpress.set','happypress.set','fearpress.set','lovepress.set','jealousypress.set','compassionpress.set','contentpress.set','griefpress.set','reliefpress.set'};
paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/rr80'};
nx=1;
for k = 1:length(datset)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
EEG = pop_saveset( EEG, datset{k}, paths{nx});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end;
 
% Call in button press datasets
eeglab
datset = {'awepress.set','frustpress.set','joypress.set','angerpress.set','sadpress.set','happypress.set','fearpress.set','lovepress.set','jealousypress.set','compassionpress.set','contentpress.set','griefpress.set','reliefpress.set'};
paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/rr80'};
nx=1;
for k = 1:length(datset)
    EEG = pop_loadset( datset{k}, paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;


emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
figure; 
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    diffpress = diff(diff(EEG.data(72,:,:)));
    %diffpress = diff(EEG.data);
    %diffpress = EEG.data(72,:,:);
    diffpress = squeeze(diffpress);
    subplot(4,4,k)
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 2500 NaN NaN NaN NaN NaN NaN],'caxis',[-max(max(diffpress(501:600,:)))-100 max(max(diffpress(501:600,:)))-100]);
 erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 400  NaN NaN NaN NaN NaN NaN]);
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 750  NaN NaN  NaN NaN NaN NaN]);
end;


figure;  k=9;
     EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
     diffpress = diff(diff(EEG.data(72,:,:)));
     %diffpress = diff(EEG.data);
     %diffpress = EEG.data;
     diffpress = squeeze(diffpress);
  erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 400  NaN NaN  NaN NaN NaN NaN],'caxis',[-90 90]);
%erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 400   NaN NaN  NaN NaN NaN NaN]);

