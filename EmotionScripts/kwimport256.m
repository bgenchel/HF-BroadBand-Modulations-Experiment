% procedures to import emotion data on kw78
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
% ****  changed because G21 is bad:  ref'ed to 135,212.  (E11,G20)
% chan 265 is the button press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
EEG = pop_readbdf(rawdat, [100:100:200] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1000');
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common1/emotion/kw78/kw78.bdf';
spath = '/data/common1/emotion/kw78/';

EEG = pop_readbdf(rawdat, [1:999:1000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
pack
EEG = pop_readbdf(rawdat, [1001:999:2000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'emo 1001 2000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 2);
pack
EEG = pop_readbdf(rawdat, [2001:999:3000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2001-3000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 3);
pack
EEG = pop_readbdf(rawdat, [3001:999:4000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,  'setname', 'emo 3001-4000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 4);
pack
EEG = pop_readbdf(rawdat, [4001:954:4955] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,  'setname', 'emo 5');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 5);
pack
EEG = pop_mergeset( ALLEEG, [1  2  3  4  5], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole session', 'overwrite', 'on');
% remove boundary events and rename numbered events
for ev = length(EEG.event):-1:1   % done
    if EEG.event(ev).type(1) == 'b'
        EEG.event(ev) = [];
    end;
end; % will exceed matrix dimensions if it takes any out
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
      %%%%%%%%%%%%%%%%%%%%%
EEG = pop_saveset( EEG, 'ButtonOnly.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ev = 2:length(EEG.event)  
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
evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'}; %this is emo track order 9-23
for ev = 1:length(relevents)
    EEG.event(relevents(ev)).type = evtypes{ev};
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    
eeglab
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bhappy','bsad','blove' ,'bfear','bcompassion','bjealousy','bcontent','bgrief','brelief','bexcite','bdisgust'};

EEG = pop_loadset( 'ButtonOnly.set', '/data/common2/emotion/kw78/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

for k = 1:length(emos)
    if k == length(emos)
        endt = (EEG.event(end).latency-EEG.event(end-4).latency)/256;
    else
        endt=300;
    end;    
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
    EEG.event(end+1).latency = (EEG.xmax-180)*256;  % make event 2 min before end
    EEG.event(end).type = 'press1';
    EEG.event(end).Event_Type = 'Response';
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
    pressevents{k} = EEG.event(end);
end;
% used above instead because button press was so close to end of each emotion
clear pressevents
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    maxdata = max(EEG.data);
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 1: size(EEG.data,2)
        if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press1';
            EEG.event(end).Event_Type = 'Response';
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;

% to epoch on the button press only data:
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','contentPress','griefPress','reliefPress','disgustPress','excitePress'};
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
    %diffpress = diff(diff(EEG.data));
    %diffpress = diff(EEG.data);
    diffpress = EEG.data;
    diffpress = squeeze(diffpress);
    subplot(4,4,k)
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 2500 NaN NaN NaN NaN NaN NaN],'caxis',[-max(max(diffpress(501:600,:)))-100 max(max(diffpress(501:600,:)))-100]);
% erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 2500 -20000 20000 NaN NaN NaN NaN],'caxis',[-4500 4500]);
 erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 2500 -3000000 3000000 NaN NaN NaN NaN],'caxis',[-3000000 3000000]);
end;
figure;  k=9;
     EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
     diffpress = diff(diff(EEG.data));
     %diffpress = diff(EEG.data);
     %diffpress = EEG.data;
     diffpress = squeeze(diffpress);
  erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts-2),emos{k},0, 1 ,'erp','limits',[-50 650 -2000 2000 NaN NaN NaN NaN],'caxis',[-5000 5000]);
     erpimage( diffpress, ones(1,size(diffpress,2))*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k},0, 1 ,'erp','limits',[-50 650 -3000000 3000000  NaN NaN NaN NaN]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  ***  Next use Button info to import and epoch real data  ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Import emotion data in three chunks, need to merge log file each time by shifts.
eeglab
rawdat = '/data/common4/emotion/kw78/kw78.bdf';
spath = '/data/common4/emotion/kw78/';
elpfile = '/data/common4/emotion/kw78/kw78-256.elp';

savedat = 'Emo-HP';

%ImportEmoRaw('kw78.bdf','kw78-256.elp',spath,savedat,[],1 );% get presses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  block 1 (done!)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [1 1010] ,'ref' ,[131 215]);% E3 and G23
%EEG = pop_readbdf(rawdat, [1:1009:1010] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [spath(end-4:end-1),' Emo 1 ']);

EEG.button = EEG.data(257,:);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG=pop_chanedit(EEG,  'load',{elpfile , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_eegfilt( EEG, .5, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG = pop_saveset( EEG, [savedat,'-1.set'], spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[]; EEG=[]; clear ans
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  block 2 (done!)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [1011 2020] ,'ref' ,[131 215]);% E3 and G23
%EEG = pop_readbdf(rawdat, [1011:1009:2020] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [spath(end-4:end-1),' Emo 2 ']);

EEG.button = EEG.data(257,:);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
eeglab redraw

EEG=pop_chanedit(EEG,  'load',{elpfile , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_eegfilt( EEG, .5, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
EEG = pop_saveset( EEG, [savedat,'-2.set'], spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  block 3  %%%%%%% (done!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [2021 2935] ,'ref' ,[131 215]);% E3 and G23
%EEG = pop_readbdf(rawdat, [2021:914:2935] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [spath(end-4:end-1),' Emo 3 ']);

EEG.button = EEG.data(257,:);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG=pop_chanedit(EEG,  'load',{elpfile , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_eegfilt( EEG, .5, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
EEG = pop_saveset( EEG, [savedat,'-3.set'], spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  block 4  %%%%%%%(done!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [2935 3940] ,'ref' ,[131 215]);% E3 and G23
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [spath(end-4:end-1),' Emo 4 ']);

EEG.button = EEG.data(257,:);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG=pop_chanedit(EEG,  'load',{elpfile , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_eegfilt( EEG, .5, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG = pop_saveset( EEG, [savedat,'-4.set'], spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  block 5 (done!) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [3941 4955] ,'ref' ,[131 215]);% E3 and G23
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [spath(end-4:end-1),' Emo 5 ']);

EEG.button = EEG.data(257,:);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG=pop_chanedit(EEG,  'load',{elpfile , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_eegfilt( EEG, .5, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');

EEG = pop_saveset( EEG, [savedat,'-5.set'], spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans

exit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete low and bad chans:
datset = {[savedat,'.set']};
numchans = 200;
zfactor = -.07;
rmref = {'E3','G23'};
[EEG,delchans,nchan] = Cut2ChanX(datset,spath,numchans,zfactor,rmref,[],[],1,1);

datset = {[savedat,'-1-',int2str(nchan),'.set'],[savedat,'-2-',int2str(nchan),'.set']};

wtstem = 'HP';
[EEG,numframes,nchan] = MakeFloats(datset,spath,wtstem,[],[]);% outputs # of 'EEG' chans used
pcs = 120;ext=1;
LinuxICA(wtstem,spath,sum(numframes),nchan,pcs,ext);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new import info:
 kw78  - ~26 yo female
data collected on May 20, 2004  using 256 (265 including analog) Biosemi acquisition system
referenced to channels 131,214 (E3 and G23) 
channel G21(213) was very bad the whole time, couldn't get the noise down with any tricks.
file includes subject debrief and photos of electrodes. Data collected by julie, jenny
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered above .5 Hz only (no low pass)
no log file merge, need to rename events manually
Performed NO button Emotion expt number 3
Kept button response as EEG.button and ECG channels retained as 'type' 'ECG'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % write into data info (original import)
 kw78  - ~26 yo female
data collected on May 20, 2004  using 256 (265 including analog) Biosemi acquisition system
referenced to channels 135,212 (E11 and G20) 
channel G21(213) was very bad the whole time, couldn't get the noise down with any tricks.
file includes subject debrief and photos of electrodes. Data collected by julie, jenny
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually
Performed NO button Emotion expt number 3

% Look at data and find bad channels
remove = [74,75,213,245];
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set'};
datset2 = {'emo-1-250.set','emo-2-250.set','emo-3-250.set','emo-4-250.set','emo-5-250.set'};
floatset = {'/data/common1/emotion/kw78/emo1.fdt','/data/common1/emotion/kw78/emo2.fdt','/data/common1/emotion/kw78/emo3.fdt','/data/common1/emotion/kw78/emo4.fdt','/data/common1/emotion/kw78/emo5.fdt'};
spath = '/data/common1/emotion/kw78/';
for r = 1:length(datset)
  EEG = pop_loadset( datset{r},spath); 
  EEG.icachansind = [1:214];  
  EEG = pop_select( EEG, 'nochannel',[215:216] );
  [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
  EEG = pop_select( EEG, 'nochannel',remove );
  [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
  EEG = pop_saveset( EEG,datset{r} , spath);
  [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
  for evtm = 384:768:size(EEG.data,2)-384  % go up by 3 sec to create 3 sec epochs
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 1000;        
    end;
    EEG = pop_epoch( EEG,{1000} , [-1.5 1.5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1500 1500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,1,[1:size(EEG.data,1)] ,4,4,0,1);        
    EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,4,4,0,1);
    EEG.data = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
floatwrite(EEG.data, floatset{r});
numframes(1,r) = size(EEG.data,2);
ALLEEG=[]; EEG=[];
end;
numframes =

      141312  +    145920  +    134400  +   158208   +    158976


sum = 738816  


% now in linux

cat  emo1.fdt  emo2.fdt   emo4.fdt  emo5.fdt  > emoall.fdt
%1,2,4,5:604416  
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/kw78/kwemoICA250.sc


% in matlab, call all three in and cat there.


figure;
for c = 1:size(EEG.data,1)
psd(EEG.data(c,:),512,EEG.srate,512);hold on;
end;

[corr] = PlotWinvCorr(wv1,wv2,chanlocs1,chanlocs2,[]);


eeglab
EEG = pop_loadset( [savedat,'-',int2str(nchan),'.set'], spath);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,1);

sph=floatread([spath,'HP',int2str(nchan),'.sph'],[nchan nchan]); 
wts =floatread([spath,'HP',int2str(nchan),'.wts'],[pcs nchan]); 

wts=floatread('/data/common2/emotion/kw78/wts250pc110.wts',[110 250]); 
sph=floatread('/data/common2/emotion/kw78/sph250pc110.sph',[250 250]); 
wts=floatread('/data/common2/emotion/kw78/wts250pc110.wts',[110 250]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
EEG = eeg_checkset(EEG);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_topoplot(EEG,0, [37:55] , 'kw78 Emotion Components; PCA to 100 ',[] ,0, 'electrodes', 'off', 'masksurf', 'on');
comps = [3,4,5,7,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26];

, 'plotrad',0.5
subj = [1:5,7:10,11,12:14,17:19,20,22,24,25,30,37,43,45,52];
subj = [1:6,9:12,14:19,22:25,28,34];
figure; pl=1;  row = ceil(length(subj)/3); col = 9;
for c = 1:length(subj)
    cc=subj(c);
    subplot(row,col,pl);
    topoplot(EEG.icawinv(:,cc),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.65);
    title(int2str(cc));
    subplot(row,col,pl+1:pl+2);    
    psd(EEG.icaact(cc,:),512,EEG.srate,512);hold on;
    title(int2str(cc));
    set(gca,'xtick',[5:5:40]);
    set(gca,'xlim',[1 40]);
    set(gca,'xgrid','on');
    pl=pl+3;
end;
ph=textsc('kw78 Emotion Component Spectra(PCA to 110) from 250 chan data, set 1','title');
set(ph,'fontsize',14);




% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
spath = '/data/common2/emotion/kw78/';
datsets = {'emo-1-250.set','emo-2-250.set','emo-3-250.set','emo-4-250.set','emo-5-250.set'};

evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'}; %this is emo track order 9-23

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
            EEG.event(end+1).type = 'press1';            
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
