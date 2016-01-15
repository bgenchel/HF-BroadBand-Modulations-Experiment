% procedures to import emotion data on ps82
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
% ref'ed to 135,212.  (E11,G20)
% chan 265 is the button press
% he was instructed to press the button at emotion surges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common2/emotion/ps82/ps82emot.bdf';
spath = '/data/common2/emotion/ps82/';

EEG = pop_readbdf(rawdat, [1:999:1000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
pack
EEG = pop_readbdf(rawdat, [1001:999:2000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'emo 2');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 2);
pack
EEG = pop_readbdf(rawdat, [2001:999:3000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 3');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 3);
pack
EEG = pop_readbdf(rawdat, [3001:999:4000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,  'setname', 'emo 4');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 4);
pack
EEG = pop_readbdf(rawdat, [4001:999:5000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,  'setname', 'emo 5');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 5);
EEG = pop_readbdf(rawdat, [5001:999:6000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,  'setname', 'emo 6');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 6);

eeglab redraw
EEG = pop_mergeset( ALLEEG, [1:6 ], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole session', 'overwrite', 'off');
% remove boundary events and rename numbered events
for ev = length(EEG.event):-1:1   % done
    if EEG.event(ev).type(1) == 'b'
        EEG.event(ev) = [];
    end;
end; 
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
      %%%%%%%%%%%%%%%%%%%%%
EEG = pop_saveset( EEG, 'ButtonOnly.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear relevents
for ev = 2:length(EEG.event)   % need to take out boundaries first to make non-strings
    if EEG.event(ev).type == '30' & EEG.event(ev-1).type == '2'
        relevents(1,1) = ev;
    end;
    if EEG.event(ev).type == '30' & EEG.event(ev-1).type == '25'
        relevents(1,2) = ev;
    end;
    if EEG.event(ev).type == '9' 
        relevents(1,3) = ev;
    end;
    if EEG.event(ev).type == '10' 
        relevents(1,4) = ev;
    end;
    if EEG.event(ev).type == '11' 
        relevents(1,5) = ev;
    end;
    if EEG.event(ev).type == '12' 
        relevents(1,6) = ev;
    end;
    if EEG.event(ev).type == '13' 
        relevents(1,7) = ev;
    end;
    if EEG.event(ev).type == '14' 
        relevents(1,8) = ev;
    end;
    if EEG.event(ev).type == '15' 
        relevents(1,9) = ev;
    end;
    if EEG.event(ev).type == '16' 
        relevents(1,10) = ev;
    end;
    if EEG.event(ev).type == '17' 
        relevents(1,11) = ev;
    end;
    if EEG.event(ev).type == '18' 
        relevents(1,12) = ev;
    end;
    if EEG.event(ev).type == '19' 
        relevents(1,13) = ev;
    end;
    if EEG.event(ev).type == '20' 
        relevents(1,14) = ev;
    end;
    if EEG.event(ev).type == '21' 
        relevents(1,15) = ev;
    end;
    if EEG.event(ev).type == '22'
        relevents(1,16) = ev;
    end;
    if EEG.event(ev).type == '23' 
        relevents(1,17) = ev;
    end;
end;
% make sure all relevents are non-zero
evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'}; %this is emo track order 9-23: DO NOT CHANGE
for ev = 1:length(relevents)
    EEG.event(relevents(ev)).type = evtypes{ev};
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'ButtonOnly.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bhappy','bsad','blove' ,'bfear','bcompassion','bjealousy','bcontent','bgrief','brelief','bexcite','bdisgust'};
ALLEEG=[];EEG=[];
EEG = pop_loadset( 'ButtonOnly.set', '/data/common2/emotion/ps82/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
% need to estimate end for excite (EEG cut off)
for k = 1:length(emos)
        endt = (EEG.event(end).latency-EEG.event(end-4).latency)/256;
    EEG = pop_epoch( EEG, {  emos{k}  }, [-.5  endt], 'newname', emos{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emobuts{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:512); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 eeglab redraw
 ALLEEG(1)=[];
%figure;plot(EEG.data); %to make sure you got the right thing
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    
    for ev = 1:length(EEG.event)
        if str2num(EEG.event(ev).type) == 24
            seltime = EEG.event(ev).latency; break
        end;
    end;
    EEG = pop_select( EEG, 'point',[1 seltime] );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, k);
end;
 eeglab redraw
%figure;plot(EEG.data); %to make sure you got the right thing
clear pressevents
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    maxdata = max(EEG.data);
    thresh = maxdata*.2;  % sets threshold at 10% of max value
    r=1; 
    for g = 2: size(EEG.data,2)
        if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press';
            EEG.event(end).Event_Type = 'Response'; 
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  ***  Next use Button info to import and epoch real data  ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Import emotion data in three chunks, need to merge log file each time by shifts.
eeglab  % need to estimate end for excite (EEG cut off)
rawdat = '/data/common2/emotion/ps82/ps82emot.bdf';
spath = '/data/common2/emotion/ps82/';
elpfile = '/data/common2/emotion/ps82/ps82-256.elp';

%%%%  block 1   %%%%%%%%%
EEG = pop_readbdf(rawdat, [1:679:680] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 1 ');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{elpfile , 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-1-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[]; EEG=[]; clear ans
pack
%%%%  block 2   %%%%%%%%%
EEG = pop_readbdf(rawdat, [681:759:1440] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 2  ');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'emo-2-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 3  %%%%%%%
EEG = pop_readbdf(rawdat, [1441:989:2430] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 3');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-3-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 4  %%%%%%%
EEG = pop_readbdf(rawdat, [2431:889:3320] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 4');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-4-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 5  %%%%%%%
EEG = pop_readbdf(rawdat, [3321:679:4000] ,265,[135:77:212]); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 5');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-5-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
%%%%%%  block 6  %%%%%%%
EEG = pop_readbdf(rawdat, [4001:1004:5005] ,265,[135:77:212]); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 6');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-6-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack  
%%%%%%  block 7  %%%%%%%
EEG = pop_readbdf(rawdat, [5006:987:5993] ,265,[135:77:212]); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 6');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-7-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack  
exit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set','emo-6-254.set','emo-7-254.set'};
for r = 1:length(datset)
EEG = pop_loadset( datset{r},spath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, datset2{r} , spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % write into data info:
 ab75  - ~29 yo male
data collected on Oct 13, 2004  using 256 (265 including analog) Biosemi acquisition system
referenced to channels 135,212 (E11 and G20) 
file includes subject debrief and photos of electrodes. Data collected by jenny, julie
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually
Performed  Button Emotion task # 2

% Look at data and find bad channels
%remove = [191,224];
remove = [69,77,162];
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set','emo-6-254.set','emo-7-254.set'};
datset2 = {'emo-1-252.set','emo-2-252.set','emo-3-252.set','emo-4-252.set','emo-5-252.set','emo-6-252.set','emo-7-252.set'};
datset3 = {'emo-1-249.set','emo-2-249.set','emo-3-249.set','emo-4-249.set','emo-5-249.set','emo-6-249.set','emo-7-249.set'};
floatset = {'emo1.fdt','emo2.fdt','emo3.fdt','emo4.fdt','emo5.fdt','emo6.fdt','emo7.fdt'};
spath = '/data/common2/emotion/ps82/';
for r = 1:length(datset2)
EEG = pop_loadset( datset2{r},spath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG = pop_select( EEG, 'nochannel',remove );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG,datset{r} , spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for evtm = 256:512:size(EEG.data,2)-256  % go up by 3 sec to create 3 sec epochs
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 1000;        
    end;
    EEG = pop_epoch( EEG,{1000} , [-1 1]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1000 1000]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejtrend(EEG,1,[1:size(EEG.data,1)] ,750,65,0.4,2,1);
    EEG = pop_rejkurt(EEG,1,[1:size(EEG.data,1)] ,4,4,0,1);        
    EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,4,4,0,1);
    EEG.data = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
floatwrite(EEG.data, [spath,floatset{r}]);
numframes(1,r) = size(EEG.data,2);
ALLEEG=[]; EEG=[];
end;
numframes =
 94720   +   111104   +   155136  +    158720  +    111616  +    169472+  169984


sum = 970752  
2:7 = 876032

% now in linux

cat  emo2.fdt  emo3.fdt   emo4.fdt   emo5.fdt   emo6.fdt  emo7.fdt > emoall.fdt
%  
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/ps82/psemoICA.sc



eeglab;  chan = 249;  spath = '/data/common2/emotion/ps82/';  sname = 'ps82';
EEG = pop_loadset( ['emo-1-',int2str(chan),'.set'],spath );
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);


sph=floatread([spath, 'sph',int2str(chan),'pc100.sph'],[chan chan]); 
wts=floatread([spath, 'wts',int2str(chan),'pc100.wts'],[100 chan]); 
EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_topoplot(EEG,0, [1:36] , [sname,' Emotion Components; PCA to 100 '],[6 6] ,0, 'electrodes', 'off', 'masksurf', 'on');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
print -Pcoloring -dpsc2 [-painters]

subj = [1:3,5,7:11,13,14,16,17,19,20,23,27,28];
figure; pl=1;  row = ceil(length(subj)/3); col = 9;
for c = 1:length(subj)
    cc=subj(c);    subplot(row,col,pl);
    topoplot(EEG.icawinv(:,cc),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.5);
    title(int2str(cc));    subplot(row,col,pl+1:pl+2);    
    psd(EEG.icaact(cc,:),512,EEG.srate,512);hold on;
    title(int2str(cc));    set(gca,'xtick',[5:5:40]);
    set(gca,'xlim',[1 40]);    set(gca,'xgrid','on');    pl=pl+3;
end;
ph=textsc([sname,' Emotion Component Spectra (PCA to 100) from ',int2str(chan),' chan dataset 1'],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

print -Pcoloring -dpsc2 [-painters]




% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
spath = '/data/common2/emotion/ps82/';
datsets = {'emo-1-249.set','emo-2-249.set','emo-3-249.set','emo-4-249.set','emo-5-249.set','emo-6-249.set','emo-7-249.set'};

evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'}; %this is emo track order 9-23
ALLEEG=[];EEG=[];
% do prebase and post base separately (without button press)
%%**  get pressevents first to mark first button press
for ds = 1:length(datsets)
    EEG = pop_loadset( datsets{ds}, spath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    for ev = 1:length(EEG.event) 
        if ev > 1
            if EEG.event(ev).type == 30 & EEG.event(ev-1).type == 2
                relevents(1,1) = ev;
            end;
            if EEG.event(ev).type == 30 & EEG.event(ev-1).type == 25
                relevents(1,2) = ev;
            end;
        end;
        if EEG.event(ev).type == 9 
            relevents(1,3) = ev;
        end;
        if EEG.event(ev).type == 10 
            relevents(1,4) = ev;
        end;
        if EEG.event(ev).type == 11 
            relevents(1,5) = ev;
        end;
        if EEG.event(ev).type == 12 
            relevents(1,6) = ev;
        end;
        if EEG.event(ev).type == 13 
            relevents(1,7) = ev;
        end;
        if EEG.event(ev).type == 14 
            relevents(1,8) = ev;
        end;
        if EEG.event(ev).type == 15 
            relevents(1,9) = ev;
        end;
        if EEG.event(ev).type == 16 
            relevents(1,10) = ev;
        end;
        if EEG.event(ev).type == 17 
            relevents(1,11) = ev;
        end;
        if EEG.event(ev).type == 18 
            relevents(1,12) = ev;
        end;
        if EEG.event(ev).type == 19 
            relevents(1,13) = ev;
        end;
        if EEG.event(ev).type == 20 
            relevents(1,14) = ev;
        end;
        if EEG.event(ev).type == 21 
            relevents(1,15) = ev;
        end;
        if EEG.event(ev).type == 22
            relevents(1,16) = ev;
        end;
        if EEG.event(ev).type == 23 
            relevents(1,17) = ev;
        end;
    end;
    for ev = 1:length(relevents)
        if relevents(ev) ~= 0
            EEG.event(relevents(ev)).type = evtypes{ev};
        end;
    end;clear relevents
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,datsets{ds} , spath);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG=[];EEG=[];
end;

%%%%%%***********  get pressevents FIRST!!!********
ALLEEG=[];EEG=[];
% need to estimate end for excite (EEG cut off)
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'};
for ds = 1:length(datsets)
    EEG = pop_loadset( datsets{ds}, spath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  

    for emo = 1:length(emos)
        if ~isempty(find(ismember({EEG.event.type}, emos{emo})))
            stf = EEG.event(find(ismember({EEG.event.type}, emos{emo}))).latency;
            endf = EEG.event(find(ismember({EEG.event.type}, emos{emo}))+1).latency;
            endt = endf/256;  stt = stf/256;  eptime = endt-stt;
            EEG = pop_epoch( EEG, {  emos{emo}  }, [-.5  eptime], 'newname', emos{emo}, 'epochinfo', 'yes');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emos{emo});
            EEG.data = rmbase(EEG.data,EEG.pnts,1:EEG.pnts); 
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG.event(end+1).type = 'press';            
            EEG.event(end).latency = pressevents{emo}(1).latency;
            EEG.event(end).epoch = 1;            
           EEG = eeg_checkset(EEG, 'eventconsistency');            
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            EEG = pop_saveset( EEG,emos{emo} , spath);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            ALLEEG = pop_delset( ALLEEG, [2] );
        end;
        EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    end;
    ALLEEG=[];EEG=[];
end;

emos = {'prebase','postbase'};
for ds = 1:length(datsets)-1:length(datsets)
    EEG = pop_loadset( datsets{ds}, spath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  

    for emo = 1:length(emos)
        if ~isempty(find(ismember({EEG.event.type}, emos{emo})))
            stf = EEG.event(find(ismember({EEG.event.type}, emos{emo}))).latency;
            endf = EEG.event(find(ismember({EEG.event.type}, emos{emo}))+1).latency;
            endt = endf/256;  stt = stf/256;  eptime = endt-stt;
            EEG = pop_epoch( EEG, {  emos{emo}  }, [-.5  eptime], 'newname', emos{emo}, 'epochinfo', 'yes');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emos{emo});
            EEG.data = rmbase(EEG.data,EEG.pnts,1:EEG.pnts); 
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = pop_saveset( EEG,emos{emo} , spath);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            ALLEEG = pop_delset( ALLEEG, [2] );
        end;
        EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    end;
    ALLEEG=[];EEG=[];
end;
