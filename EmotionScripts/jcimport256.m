% procedures to import emotion data on jc82
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
%  ref'ed to 135,213.
% chan 265 is the button press
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
EEG = pop_readbdf('/data/common1/emotion/jc82/jc82emo.bdf', [1:1099:1100] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1100');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_readbdf('/data/common1/emotion/jc82/jc82emo.bdf', [1101:1099:2200] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'emo 1101 2200');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_readbdf('/data/common1/emotion/jc82/jc82emo.bdf', [2201:798:3000] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2201-3300');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_readbdf('/data/common1/emotion/jc82/jc82emo.bdf', [3001:674:3676] ,265,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2201-3300');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% merge and import log file
EEG = pop_mergeset( ALLEEG, [1  2  3  4], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole session', 'overwrite', 'on');
EEG = pop_importpres(EEG,  '/data/common1/emotion/jc82/jc82-emotionimage.log', 'Code', 'Time',0);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'ButtonOnly.set', '/data/common1/emotion/jc82/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG(1)=[];
ALLEEG(1)=[];
ALLEEG(1)=[];
eeglab redraw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
emos = {'relief','grief','content','embarrass','compassion','jealousy','love','fear','happy','sad','anger','joy', 'frustration' ,'awe'};
emobuts = {'brelief','bgrief','bcontent','bembarrass','bcompassion','bjealousy','blove','bfear','bhappy','bsad','banger','bjoy', 'bfrustration' ,'bawe'};
 EEG = pop_loadset( 'ButtonOnly.set', '/data/common1/emotion/jc82/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
 eeglab redraw

for k = 12:length(emos)
    if k == 1
        endt = 150;
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
        if str2num(EEG.event(ev).type) == 250
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
emonames = {'reliefPress','griefPress','contentPress','embarrassPress','compassionPress','jealousyPress','lovePress','fearPress','happyPress','sadPress','angerPress','joyPress', 'frustrationPress' ,'awePress'};
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
 jc82  - ~22 yo female
data collected on Feb 4, 2004  using 256 (265 including analog) Biosemi acquisition system
merged with presentation log file 
file includes subject debrief and photos of electrodes. Data collected by julie
Had 60 Hz noise in data when ran through analog box. Sometimes alleviated by grounding, sometimes it made it worse.
filtered 1 50



















% Remove all analog channels:
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG(1), EEG, CURRENTSET, 'setname', 'emotion 1101 2200 - 256' , 'overwrite', 'on');
 eeglab redraw
% load in chanlocs
% EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/jo74/jo256.elp', 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/home/julie/dummy-small-256.elp', 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% filter the data
EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
 eeglab redraw
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
 eeglab redraw

% gets rid of low electrodes
alllocs = {EEG.chanlocs.Z};
alllocs = cell2mat(alllocs);

below = find(alllocs<-.035);    %% finds electrodes BELOW the 'center' (all Z's below zero)
%below = find(alllocs<-.033);    %% finds electrodes BELOW the 'center' (all Z's below zero)

below(find(below == 249))=[];   %% get back the below the eye channels 
below(find(below == 250))=[];   % only works the first time around (with all 256 channels)
below(find(below == 251))=[];
below(find(below == 252))=[];

%below(find(below == 253))=[];
%below(find(below == 254))=[];
%below(find(below == 255))=[];
%below(find(below == 256))=[];

EEG = pop_select( EEG, 'nochannel',below);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
% gives 189 channels
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);   % hi pass filter over 1 Hz
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'emotion 2201-3200', 'overwrite', 'on');

EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
% need to epoch with no overlap and run rejection because there are lots of noisy time periods
% first remove these channels:
%[140,141,177,207];
%EEG = pop_select( EEG, 'nochannel',[140,141,177,207] );
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
for evtm = 256:512:size(EEG.data,2)-256  % go up by 1.5 sec to create 3 sec epochs til 2 sec before button
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
    EEG.event(end).latency = evtm;
    EEG.event(end).type = 'fake';        
end;
EEG = pop_epoch( EEG,{'fake'} , [-1 1]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
EEG = pop_rmbase( EEG,[-1000 1000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_rejkurt(EEG,1,[1:248] ,5,4,0,1);
EEG = pop_jointprob(EEG,1,[1:248] ,5,4,0,1);
x=reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
floatwrite(x, '/data/common1/emotion/jc82/emo1-248.fdt');
% emo1:    emo2:  

cat emo1-248.fdt emo2-248.fdt  > emo248-1-2.fdt
/data/common/matlab/ica_linux2.4 < EmoEegIca.sc



EEG = pop_saveset( EEG, 'emo-1-248.set', '/data/common1/emotion/jc82/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

 floatwrite(EEG.data, '/data/common1/emotion/jo74/alldataeps.fdt')

 floatwrite(EEG.data, '/data/common1/RewTwoback/jo74/tb2eps.fdt')
% in linux terminal: 
 cat emo1epochs.fdt emo2epochs.fdt emo3epochs.fdt  > allemoeps.fdt
 cat tb1eps.fdt tb2eps.fdt   > alltbeps.fdt
 cat emo1-252.fdt emo2-252.fdt emo3-252.fdt   > emo252-1-2.fdt
  % run ica
/data/common/matlab/ica_linux2.4 < allemo3.sc
/data/common/matlab/ica_linux2.4 < alltb.sc
/data/common/matlab/ica_linux2.4 < emo256ica.sc
/data/common/matlab/ica_linux2.4 < emo252ica2.sc
/data/common/matlab/ica_linux2.4 < allimageica.sc

% in matlab, call all three in and cat there.
dat3=floatread('emo3.fdt',[189 inf]); 
data=floatread('allemoeps.fdt',[189 inf]) ;
alldata = cat(2,alldata,dat3);
% or in matlab
[wts,sph] =binica(alldata);
 % makes 194 chan total

EEG = pop_loadset( 'emo-1-252.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
for index = 1:length(EEG.event)
    if length(EEG.event(index).type)>2
   if EEG.event(index).type(1:3) == 'ecl'
      EEG.event(end+1) =  EEG.event(index);% appends events to the end
      EEG.event(end).latency = EEG.event(index).latency+3*EEG.srate;
      EEG.event(end).type = 'fix';
      EEG.event(end+1) =  EEG.event(index);% appends events to the end
      EEG.event(end).latency = EEG.event(index).latency+7*EEG.srate;
      EEG.event(end).type = 'fix';
   end;
   end;
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% Reorder events according to latency
EEG = pop_editeventvals(EEG, 'sort',{ 'latency',0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


EEG = pop_loadset( 'emo-3-252.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG = pop_epoch( EEG, {  'eclose'  }, [0 10], 'newname', 'eclose3', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'eclose3', 'overwrite', 'on');
EEG = pop_rmbase( EEG, [0 10000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


  % make fake events every sec throughout  nf expts  (continuous data)for ica
tmpts = size(EEG.data,2); pl=1;
for index = EEG.srate:EEG.srate*2:tmpts
      %EEG.event(end+1) =  EEG.event(end);% appends events to the end
      EEG.event(pl).latency = index;
      EEG.event(pl).type = 'fake'; pl = pl+1;
end;
EEG = eeg_checkset(EEG, 'eventconsistency');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% Reorder events according to latency
EEG = pop_editeventvals(EEG, 'sort',{ 'latency',0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


[outeeg,indices] = pop_epoch (EEG, {'fake'}, [-1 1], 'epochinfo','yes');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, outeeg,CURRENTSET);
EEG = pop_rmbase( EEG, [-1000 1000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw

figure;
for c = 1:size(EEG.data,1)
psd(EEG.data(c,:),512,EEG.srate,512);hold on;
end;


sph=floatread('sph155.sph',[155 155]); 
wts=floatread('wts155.wts',[155 155]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

sph=floatread('sph189-pc80.sph',[189 189]); 
wts=floatread('wts189-pc80.wts',[80 189]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

sph=floatread('sph252-160.sph',[252 252]); 
wts=floatread('wts252-160.wts',[160 252]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


taken out after scale down to 189
cut = [3     4     6    18  25  44    52    53    58    63    73    76    80    81 93    94    98   102 103 104 117  119   120   122   143   148   153   154   156 157  161 165   166   178]


    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
pop_prop( EEG, 0, 16);
    EEG = eeg_retrieve(ALLEEG, 2); CURRENTSET = 2;

pop_prop( EEG, 0, 8);


 [trmeans] = changrads(EEG.icawinv,'/data/common1/emotion/jo74/mylocs.loc');