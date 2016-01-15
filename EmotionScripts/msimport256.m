% procedures to import emotion data on ms82
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
%  ref'ed to 135,213.  (E11,G21)
% chan 265 is the button press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common1/emotion/ms82/ms82_2.bdf';
spath = '/data/common1/emotion/ms82/';

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
EEG = pop_readbdf(rawdat, [4001:891:4892] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,  'setname', 'emo 4001-4892');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 5);
% merge and import log file
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

eeglab
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief','bdisgust','bexcite'};

EEG = pop_loadset( 'ButtonOnly.set', '/data/common2/emotion/ms82/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
 eeglab redraw

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
rawdat = '/data/common1/emotion/ms82/ms82_2.bdf';
spath = '/data/common1/emotion/ms82/';
elpfile = '/data/common1/emotion/ms82/ms82-256.elp';

EEG = pop_readbdf(rawdat, [1:1107:1108] ,265,[135:77:212]);
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
EEG = pop_readbdf(rawdat, [1109:1012:2121] ,265,[135:77:212]);
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
EEG = pop_readbdf(rawdat, [2122:1013:3135] ,265,[135:77:212]);
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
EEG = pop_readbdf(rawdat, [3136:1058:4194] ,265,[135:77:212]);
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
EEG = pop_readbdf(rawdat, [4195:697:4892] ,265,[135:77:212]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 5');
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
EEG = pop_saveset( EEG, 'emo-5-254.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
 eeglab redraw

 exit
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % write into data info:
 tl81  - ~23 yo female
data collected on May 12, 2004  using 256 (265 including analog) Biosemi acquisition system
merged with presentation log file 
file includes subject debrief and photos of electrodes. Data collected by julie, jenny
file includes subject debrief and photos of electrodes. Data collected by julie,
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually

% Look at data and find bad channels
remove = [23,24,25,29,31,61,67,81,133,134,137,140,144,156,192,207];
%remove = [11,17,23,24,25,29,61,80,126,140,144,156,192];
spath = '/data/common2/emotion/ms82/';
datset = {'emo-1-254.set','emo-2-254.set','emo-3-254.set','emo-4-254.set','emo-5-254.set'};
datset2 = {'emo-1-238.set','emo-2-238.set','emo-3-238.set','emo-4-238.set','emo-5-238.set'};
floatset = {'/data/common1/emotion/ms82/emo1.fdt','/data/common1/emotion/ms82/emo2.fdt','/data/common1/emotion/ms82/emo3.fdt','/data/common1/emotion/ms82/emo4.fdt','/data/common1/emotion/ms82/emo5.fdt'};
for r = 1:length(datset)
EEG = pop_loadset( datset2{r},spath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
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
floatwrite(EEG.data, floatset{r});
numframes(1,r) = size(EEG.data,2);
ALLEEG=[]; EEG=[];
end;
numframes =
      205824   +   179712  +    183552   +   171264  +    135168

sum =  875520 

numframes = with autoreject at 3,3

      160512 +     122112  +    129024   +   121344  +    115968

sum =  648960

% now in linux

cat  emo1.fdt  emo2.fdt  emo3.fdt  emo4.fdt  emo5.fdt   > emoall.fdt
% sets 1,2,3,5: 704256 frames
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/ms82/msemoICA238.sc


eeglab
EEG = pop_loadset( 'emo-1-238.set', '/data/common2/emotion/ms82/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

%sph=floatread('/data/common2/emotion/ms82/sph238-145.sph',[238 238]); % tried and deleted, ugly
%wts=floatread('/data/common2/emotion/ms82/wts238-145.wts',[145 238]); 
%sph=floatread('/data/common2/emotion/ms82/sph238-135.sph',[238 238]); 
%wts=floatread('/data/common2/emotion/ms82/wts238-135.wts',[135 238]); 
sph=floatread('/data/common2/emotion/ms82/sph238-110.sph',[238 238]); 
wts=floatread('/data/common2/emotion/ms82/wts238-110.wts',[110 238]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

pop_topoplot(EEG,0, [1:20] , 'ms82 Emotion Components; 238 chan; PCA to 110 ',[6 6] ,0, 'electrodes', 'off', 'masksurf', 'on');
37:72   ,'plotrad',0.5
1:36

subj = [5,7,8,9,10,11,16,17,18,20,21,30,33,37,40,58];
subj = [5:7,9,10,11,12,13,14,22,27,32];
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
    set(gca,'xticklabel',{5 [] 15 [] 25 [] 35 []});
    set(gca,'xlim',[1 40]);
    set(gca,'xgrid','on');
    pl=pl+3;
end;
ph=textsc('ms82 Emotion Component Spectra (PCA to 110) from 238 chan data, set 1','title');
set(ph,'fontsize',14);


% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
spath = '/data/common2/emotion/ms82/';
datsets = {'emo-1-238.set','emo-2-238.set','emo-3-238.set','emo-4-238.set','emo-5-238.set'};

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
% epoch on emotion til stop event 
% long emo length = 235 sec
% need to rename events first
spath = '/data/common2/emotion/ms82/';
datsets = {'emo-1-238.set','emo-2-238.set','emo-3-238.set','emo-4-238.set','emo-5-238.set'};
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite'};
%%**  get pressevents first to mark first button press

    ALLEEG=[];EEG=[];
for ds = 1:length(datsets)
    EEG = pop_loadset( datsets{ds}, spath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    for emo = 1:length(emos)
        if find(ismember({EEG.event.type}, emos{emo}))
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
            for rr = 1:length(pressevents{emo})
                EEG.event(end+1).type = 'press';            
                EEG.event(end).latency = pressevents{emo}(rr).latency;
                EEG.event(end).epoch = 1; 
            end;      
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
    ALLEEG=[];EEG=[];
for ds = 1:length(datsets)-1:length(datsets)
    EEG = pop_loadset( datsets{ds}, spath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    for emo = 1:length(emos)
        if find(ismember({EEG.event.type}, emos{emo}))
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
