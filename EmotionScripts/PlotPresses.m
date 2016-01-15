% plot emotion presses superimposed (from data epoched on beg of press)
% eeglab
emos = {'awepress.set' 'awe2press.set' 'frustpress.set' 'frust2press.set' 'joypress.set' 'joy2press.set' 'angerpress.set' 'anger2press.set' 'sadpress.set' 'sad2press.set' 'surprisepress.set' 'surprise2press.set' 'happypress.set' 'happy2press.set' 'fearpress.set' 'fear2press.set' 'lovepress.set' ...
        'love2press.set' 'jealousypress.set' 'jealousy2press.set' 'compassionpress.set' 'compassion2press.set' 'contentpress.set' 'content2press.set' 'griefpress.set' 'grief2press.set' 'reliefpress.set' 'relief2press.set'}; %ap80
%emos = {'frustpress.set'  'joypress.set'  'angerpress.set'  'sadpress.set' 'surprisepress.set' 'happypress.set'   'lovepress.set'  'jealousypress.set' 'compassionpress.set'  'contentpress.set' 'griefpress.set' 'reliefpress.set' }; % rr80
emos = { 'frustpress.set'  'joypress.set'  'angerpress.set' 'sadpress.set'  'surprisepress.set' 'happypress.set'  'fearpress.set'  'love2press.set'  'jealousy2press.set'  'compassion2press.set'  'content2press.set'  'grief2press.set' }; %ap80
figure;pl = 1; 
for em = 1:length(emos)
    EEG = pop_loadset( emos{em},'/data/common1/emotion/ap80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    tm = find (EEG.times>-200 & EEG.times<1500);
%    if em == 5 
%        figure;pl = 1;    end;
%    if em == 9
%        figure;pl = 1;    end;
%    if em == 13 
%        figure;pl = 1;    end;
%    if em == 17 
%        figure;pl = 1;    end;
%    if em == 20
%        figure;pl = 1;    end;
%    if em == 17
%        figure;pl = 1;    end;
    for p = 3:20%:length(EEG.epoch)
        subplot(round(sqrt(length(emos)))+1,round(sqrt(length(emos))),pl)
        plot(EEG.times(1,tm(1:end-1)),diff(EEG.data(72,tm,p)));
        %plot(EEG.times(1,tm(1:end)),EEG.data(72,tm,p));
        hold on;       
    end;
    set(gca,'xlim',[-200 1000]);
y = get(gca,'ylim');
    set(gca,'ylim',[0 y(2)]);
    set(gca,'yticklabel',[]);
    set(gca,'xticklabel',[]);
    title(emos{em});
    pl = pl+1; 
    ALLEEG = pop_delset( ALLEEG, [1] );
end;
 

%%%%  for jo74 (256 channel data)
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
% now save avg presses for other scripts:

for k = 1:length(emos)
        EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
avpr = mean(EEG.data,3);
allavg(k,:) = avpr;
end;
