% finds avgerage button-press duration for each emotion
eeglab
path = '/data/common1/emotion/jo74/';
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bemabarrass','bcontent','bgrief','brelief'};

cd (path)
EEG = pop_loadset('ButtonOnly.set',path);
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
 eeglab redraw
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
 eeglab redraw
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
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    time  = EEG.times;
    tm = find(time>-100 & time<3000);
    time = time(tm);
    subplot(4,4,k)
    meg = mean(EEG.data,3);meg = meg(1,tm);
    plot(time,meg);
    top = find(meg == (max(meg)));
    stp = find(meg(1,top:end)<350000);
    stm = time(stp(1)+(top-1));
    set(gca,'xtick',[-100:100:3000]);
    set(gca,'xticklabel',{[] 0  1 2  3 4  5  6  7  8  9 10 [] [] [] [] 15 [] [] [] [] 20 [] [] [] []  25 [] [] [] [] 30});
    set(gca,'xgrid','on');
    set(gca,'xlim',[-100 stm+500]);
    set(gca,'fontsize',14);
    title(emos{k});    
    lenset(1,k) = EEG.times(1,stp(1)+(top-1)+(tm(1)-1)+50);
end;

% result:
lenset =

   1.0e+03 *

  Columns 1 through 7 

awe-fear:    2.4141    0.8086    0.7109    0.7656    2.1875    1.0469    1.0117

  Columns 8 through 14 

love-relief:    1.4219    0.8047    1.8242    0.9219    1.3516    1.5352    1.7070
