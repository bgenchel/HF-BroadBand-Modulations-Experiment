% procedures to import emotion data on sg75
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
% ****  changed because G21 is bad:  ref'ed to 135,212.  (E11,G20)
% chan 265 is the button press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common2/emotion/sg75/sg75.bdf';
spath = '/data/common2/emotion/sg75/';

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
EEG = pop_readbdf(rawdat, [4001:932:4933] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,  'setname', 'emo 5');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 5);
pack
EEG = pop_mergeset( ALLEEG, [1  2  3  4  5], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole session', 'overwrite', 'on');
% remove boundary events and rename numbered events
for ev = length(EEG.event):-1:1   %  need to do this before labelling events to make non-strings
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
evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'}; %this is emo track order 9-23
for ev = 1:length(relevents)
    EEG.event(relevents(ev)).type = evtypes{ev};
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

eeglab
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bhappy','bsad','blove' ,'bfear','bcompassion','bjealousy','bcontent','bgrief','brelief','bexcite','bdisgust'};

EEG = pop_loadset( 'ButtonOnly.set', '/data/common2/emotion/sg75/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

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
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 2: size(EEG.data,2)
        if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press1';
            EEG.event(end).Event_Type = 'Response'; break  % only logs first press
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  ***  Next use Button info to import and epoch real data  ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Import emotion data in three chunks, need to merge log file each time by shifts.
eeglab
rawdat = '/data/common2/emotion/sg75/sg75.bdf';
spath = '/data/common2/emotion/sg75/';
elpfile = '/data/common2/emotion/sg75/sg75-256.elp';

EEG = pop_readbdf(rawdat, [1:949:950] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 1 ');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{elpfile , 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-1-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[]; EEG=[]; clear ans
pack
%%%%  block 2   %%%%%%%%%
EEG = pop_readbdf(rawdat, [951:989:1940] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 2  ');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'emo-2-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 3  %%%%%%%
EEG = pop_readbdf(rawdat, [1941:929:2870] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 3');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-3-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 4  %%%%%%%
EEG = pop_readbdf(rawdat, [2871:998:3869] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 4');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-4-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 5  %%%%%%%
EEG = pop_readbdf(rawdat, [3869:1031:4900] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 5');
EEG = pop_select( EEG, 'nochannel',[255:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-5-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans

exit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % write into data info:
 sg75  - ~29 yo female
data collected on June 8,2004  using 256 (265 including analog) Biosemi acquisition system
referenced to channels 135,212 (E11 and G20) 
file includes subject debrief and photos of electrodes. Data collected by julie,jenny
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually
Performed No Button Emotion task # 1

% Look at data and find bad channels
remove = [13,28,30,92,102,141,185,207];
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set'};
datset2 = {'emo-1-246.set','emo-2-246.set','emo-3-246.set','emo-4-246.set','emo-5-246.set'};
floatset = {'/data/common2/emotion/sg75/emo1.fdt','/data/common2/emotion/sg75/emo2.fdt','/data/common2/emotion/sg75/emo3.fdt','/data/common2/emotion/sg75/emo4.fdt','/data/common2/emotion/sg75/emo5.fdt'};
spath = '/data/common2/emotion/sg75/';
for r = 1:length(datset)
    EEG = pop_loadset( datset{r},spath); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'nochannel',remove );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    EEG = pop_saveset( EEG,datset2{r} , spath);
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

      120576   +   144384  +    136704   +   156672   +   166656


sum = 724992  


% now in linux

cat  emo1.fdt  emo2.fdt emo3.fdt   emo4.fdt  emo5.fdt  > emoall.fdt
%1:5:  
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/sg75/sgemoICA.sc


% in matlab, call all three in and cat there.


eeglab
EEG = pop_loadset( 'emo-1-246.set', '/data/common2/emotion/sg75/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

sph=floatread('/data/common2/emotion/sg75/sph246pc100.sph',[246 246]); 
wts=floatread('/data/common2/emotion/sg75/wts246pc100.wts',[100 246]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_topoplot(EEG,0, [1:36] , 'sg75 Emotion Components; PCA to 100 ',[6 6] ,0, 'electrodes', 'off', 'masksurf', 'on');
, 'plotrad',0.5
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
print -Pcoloring -dpsc2 [-painters]

subj = [2:5,7:11,14:18,25];
figure; pl=1;  row = ceil(length(subj)/3); col = 9;
for c = 1:length(subj)
    cc=subj(c);
    subplot(row,col,pl);
    topoplot(EEG.icawinv(:,cc),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.5);
    title(int2str(cc));
    subplot(row,col,pl+1:pl+2);    
    psd(EEG.icaact(cc,:),512,EEG.srate,512);hold on;
    title(int2str(cc));
    set(gca,'xtick',[5:5:40]);
    set(gca,'xlim',[1 40]);
    set(gca,'xgrid','on');
    pl=pl+3;
end;
ph=textsc('sg75 Emotion Component Spectra (PCA to 100) from 246 chan data, set 1','title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

print -Pcoloring -dpsc2 [-painters]




% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
spath = '/data/common2/emotion/sg75/';
datsets = {'emo-1-246.set','emo-2-246.set','emo-3-246.set','emo-4-246.set','emo-5-246.set'};

evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'}; %this is emo track order 9-23

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
