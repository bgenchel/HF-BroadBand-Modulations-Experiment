% procedures to import emotion data on tv81
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
% ****  changed because G21 is bad:  ref'ed to 135,212.  (E11,G20)
% chan 265 is the button press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common1/emotion/tv81/tv81.bdf';
spath = '/data/common1/emotion/tv81/';

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
EEG = pop_readbdf(rawdat, [2001:1049:3050] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2001-3000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 3);
pack
EEG = pop_readbdf(rawdat, [3051:1028:4079] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,  'setname', 'emo 3001-4000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 4);
pack
EEG = pop_mergeset( ALLEEG, [1  2  3  4], 0);
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
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bhappy','bsad','blove' ,'bfear','bcompassion','bjealousy','bcontent','bgrief','brelief','bexcite','bdisgust'};

EEG = pop_loadset( 'ButtonOnly.set', '/data/common2/emotion/tv81/');
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
    maxdata = max(EEG.data);
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 1: size(EEG.data,2)
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
rawdat = '/data/common2/emotion/tv81/tv81.bdf';
spath = '/data/common2/emotion/tv81/';
elpfile = '/data/common2/emotion/tv81/tv81-256.elp';

EEG = pop_readbdf(rawdat, [1:899:900] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 1 ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{elpfile , 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
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
EEG = pop_readbdf(rawdat, [901:969:1870] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 2  ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
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
EEG = pop_readbdf(rawdat, [1871:879:2750] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 3');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
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
EEG = pop_readbdf(rawdat, [2751:749:3500] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 4');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
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
EEG = pop_readbdf(rawdat, [3501:578:4079] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 5');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
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
tv81  - ~23 yo male
data collected on May 25, 2004  using 256 (265 including analog) Biosemi acquisition system
referenced to channels 135,212 (E11 and G20) 
channel G21(213) is dead. flat lined; C12 started badly, but flattend out by the end
file includes subject debrief and photos of electrodes. Data collected by julie, jenny
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually
Performed NO button emotion task # 4

% Look at data and find bad channels
remove = [1,21,25,76,120,121,125,139,143,153,154,158,190,200,203,204,208,213];
below =  [15,65,66,74,75,86,87,89,90,91,92,93,94,95,96,128,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,193,194,195,196,197,198,199,  200,201,202,203,204,205,206,207,208,209,211,212,213,219,220,221,222,224]; % found by taking all locs below z = .04 (see below)
remove = [below remove];
remove = unique(remove);
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set'};
datset2 = {'emo-1-180.set','emo-2-180.set','emo-3-180.set','emo-4-180.set','emo-5-180.set'};
floatset = {'/data/common2/emotion/tv81/emo1.fdt','/data/common2/emotion/tv81/emo2.fdt','/data/common2/emotion/tv81/emo3.fdt','/data/common2/emotion/tv81/emo4.fdt','/data/common2/emotion/tv81/emo5.fdt'};
spath = '/data/common2/emotion/tv81/';
for r = 1:length(datset)
EEG = pop_loadset( datset2{r},spath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%alllocs = {EEG.chanlocs.Z};
%alllocs = cell2mat(alllocs);
%below = find(alllocs<-.04);    %% finds electrodes BELOW the 'center' (all Z's below zero)
%below(find(below == 249))=[];   %% get back the below the eye channels 
%below(find(below == 250))=[];   % only works the first time around (with all 256 channels)
%below(find(below == 251))=[];
%below(find(below == 252))=[];
%remove = [below remove];
%remove = unique(remove);
%EEG = pop_select( EEG, 'nochannel',remove );
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
%EEG = pop_saveset( EEG,datset2{r} , spath);
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for evtm = 384:768:size(EEG.data,2)-384  % go up by 3 sec to create 3 sec epochs
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 1000;        
    end;
    EEG = pop_epoch( EEG,{1000} , [-1.5 1.5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1500 1500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,1,[1:size(EEG.data,1)] ,3,3,0,1);        
    EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,3,3,0,1);
    EEG.data = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
floatwrite(EEG.data, floatset{r});
numframes(1,r) = size(EEG.data,2);
ALLEEG=[]; EEG=[];
end;
numframes =

      105216   +   105216   +   102144    +   92160   +    68352
sum = 473088



% now in linux

cat  emo1.fdt  emo2.fdt   emo3.fdt  emo4.fdt  emo5.fdt  > emoall.fdt
% use pca  to 130
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/tv81/tvemoICA252.sc


% in matlab, call all three in and cat there.


figure;
for c = 1:size(EEG.data,1)
psd(EEG.data(c,:),512,EEG.srate,512);hold on;
end;


eeglab
EEG = pop_loadset( 'emo-1-180.set', '/data/common2/emotion/tv81/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

sph=floatread('/data/common2/emotion/tv81/sph180-90.sph',[180 180]); 
wts=floatread('/data/common2/emotion/tv81/wts180-90.wts',[90 180]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_topoplot(EEG,0, [1:36] , 'tv81 Emotion Components; 180 chan; PCA to 90 ',[6 6] ,0, 'electrodes', 'off', 'masksurf', 'on');
, 'plotrad',0.5
subj = [1:6,8:11,13:15,17,19,20,22,23,25,27,29,30,32,35];
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
ph=textsc('tv81 Emotion Component Spectra from 180 chan data (pca to 90), set 1','title');
set(ph,'fontsize',14);




% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
spath = '/data/common2/emotion/tv81/';
datsets = {'emo-1-180.set','emo-2-180.set','emo-3-180.set','emo-4-180.set','emo-5-180.set'};

evtypes = {'prebase','postbase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'}; %this is emo track order 9-23

% do prebase and post base separately (without button press)
%%**  get pressevents first to mark first button press
ALLEEG=[];EEG=[];
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
emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'};
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

