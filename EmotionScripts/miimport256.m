% procedures to import emotion data on mi83
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
%  ref'ed to 135,213.  (E11,G21)
% chan 265 is the button press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  *** Import Button Only   ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % import with event channel 265, then keep only channel 257 (button press)
eeglab
rawdat = '/data/common4/emotion/mi83/mi83.bdf';
spath = '/data/common4/emotion/mi83/';

EEG = pop_biosig(rawdat, 'blockrange', [1 1000]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo 1 1000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);

EEG = pop_biosig(rawdat, 'blockrange', [1001 2000]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,  'setname', 'emo 1001 2000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 2);

EEG = pop_biosig(rawdat, 'blockrange', [2001 3000]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,  'setname', 'emo 2001-3000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 3);

EEG = pop_biosig(rawdat, 'blockrange', [3001 4000]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,  'setname', 'emo 3001-4000');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 4);

EEG = pop_biosig(rawdat, 'blockrange', [4001 5000]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,  'setname', 'emo 4001-4513');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 5);

EEG = pop_biosig(rawdat, 'blockrange', [5001 5497]);  % CHECK FOR LAST FRAME IN RAW DATA FILE!!!!
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,  'setname', 'emo 5001-5497');
EEG = pop_select( EEG, 'channel',257);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 6);
% merge and import log file
EEG = pop_mergeset( ALLEEG, [1  2  3  4  5 6], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'Button Press, whole session', 'overwrite', 'on');
% remove boundary events and rename numbered events
for ev = 1:length(EEG.event)   % done
    if EEG.event(ev).type(1) == 'b'
        EEG.event(ev) = [];
    end;
end; % will exceed matrix dimensions if it takes any out
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
      %%%%%%%%%%%%%%%%%%%%%
EEG = pop_saveset( EEG, 'ButtonOnly.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
%%5485.996%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff = {EEG.event.type};
for f = 1:length(ff)
    fff(1,f) = str2num(ff{f});
end;

chname = find(fff ==1);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'intro';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==2);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'prebase';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%chname = find(fff ==3);
%for f = 1:length(chname)
%EEG.event(chname(f)).type = 'instruct';
%end;
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% either or:
chname = find(fff ==4);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'instruct';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==5);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'instruct2';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==6);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'relax';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==7);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'instruct3';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==8);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'intoemo';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==9);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'awe';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==10);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'frustration';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==11);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'joy';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==12);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'anger';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==13);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'happy';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==14);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'sad';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==15);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'love';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==16);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'fear';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==17);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'compassion';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==18);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'jealousy';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==19);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'content';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==20);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'grief';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==21);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'relief';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==22);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'excite';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==23);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'disgust';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==26);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'thanks';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

chname = find(fff ==100);
for f = 1:length(chname)
EEG.event(chname(f)).type = 'endvoice';
end;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);



eeglab
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief','bdisgust','bexcite'};

EEG = pop_loadset( 'ButtonOnly.set', '/data/common1/emotion/mi83/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
 eeglab redraw

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
            EEG.event(end).type = 'press';
            EEG.event(end).Event_Type = 'Response';
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;
save /data/common1/emotion/mi83/pressevents.mat pressevents

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
rawdat = '/data/common4/emotion/mi83/mi83.bdf';
spath = '/data/common4/emotion/mi83/';
elpfile = '/data/common4/emotion/mi83/mi83-256.elp';

%EEG = pop_readbdf(rawdat, [1:1048:1049] ,265,[135:78:213]);
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [1 1049] ,'ref' ,[133 214]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 1 ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{elpfile , 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-1-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[]; EEG=[]; clear ans
pack
%%%%  block 2   %%%%%%%%%
%EEG = pop_readbdf(rawdat, [1050:865:1915] ,265,[135:78:213]);
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [1050 1915] ,'ref' ,[133 214]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 2  ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'emo-2-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 3  %%%%%%%
%EEG = pop_readbdf(rawdat, [1916:934:2850] ,265,[135:78:213]);
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [1916 2850] ,'ref' ,[133 214]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 3');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-3-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 4  %%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [2851 3830] ,'ref' ,[133 214]);
%EEG = pop_readbdf(rawdat, [2851:979:3830] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 4');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-4-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 5  %%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [3831 4640] ,'ref' ,[133 214]);
%EEG = pop_readbdf(rawdat, [3831:809:4640] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 5');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-5-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
%%%%%%  block 6  %%%%%%%
EEG = pop_biosig(rawdat,'channels' ,[1:265], 'blockrange' , [4641 5497] ,'ref' ,[133 214]);
%EEG = pop_readbdf(rawdat, [4641:856:5497] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 6');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B9'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[255:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-6-254-NF.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
 eeglab redraw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % write into data info:
 mi83  - ~ yo male
data collected on May 13, 2004  using 256 (265 including analog) Biosemi acquisition system
merged with presentation log file 
file includes subject debrief and photos of electrodes. Data collected by julie, jenny
Had some 60 Hz noise in data when ran through analog box. Partly alleviated by grounding, 
filtered 1 50
no log file merge, need to rename events manually

% Look at data and find bad channels
remove = [21,66,69,82,87,92];%[66 87 92];
datset = {'emo-1-254-NF.set','emo-2-254-NF.set','emo-3-254-NF.set','emo-4-254-NF.set','emo-5-254-NF.set','emo-6-254-NF.set'};
datset2 = {'emo-1-248-NF.set','emo-2-248-NF.set','emo-3-248-NF.set','emo-4-248-NF.set','emo-5-248-NF.set','emo-6-248-NF.set'};
floatset = {'/data/common4/emotion/mi83/emo1.fdt','/data/common4/emotion/mi83/emo2.fdt','/data/common4/emotion/mi83/emo3.fdt','/data/common4/emotion/mi83/emo4.fdt','/data/common4/emotion/mi83/emo5.fdt','/data/common4/emotion/mi83/emo6.fdt'};
for r = 2:length(datset)
EEG = pop_loadset( datset{r},'/data/common4/emotion/mi83/'); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG = pop_select( EEG, 'nochannel',remove );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_saveset( EEG,datset2{r} , '/data/common4/emotion/mi83/');
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
      122112   +   111360  +    136704   +   155136   +   119040  +    129792


sum =  774144


% now in linux

cat  emo1.fdt  emo3.fdt  emo4.fdt  emo5.fdt  emo6.fdt   > emoall.fdt
% sets 1,3:6 = 662784 frames
/data/common/matlab/ica_linux2.4 < /data/common4/emotion/mi83/miemoICA248.sc

emo3.fdt    emo5.fdt 

% in matlab, call all three in and cat there.


figure;
for c = 1:size(EEG.data,1)
psd(EEG.data(c,:),512,EEG.srate,512);hold on;
end;


eeglab
EEG = pop_loadset( 'emo-1-248-NF.set', '/data/common4/emotion/mi83/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,1);

sph=floatread('/data/common4/emotion/mi83/sph248-2.sph',[248 248]); 
wts=floatread('/data/common4/emotion/mi83/wts248-2.wts',[100 248]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_topoplot(EEG,0, [1:50] , 'mi83 Emotion Components 248 chan data (FINAL) ',[] ,0, 'electrodes', 'off', 'plotrad',0.65, 'masksurf', 'on');
EEG.artifact.Veyes = [18];
EEG.artifact.Heyes = [11];
EEG.artifact.muscle = [3,22,25,29,31,33,34,38,42,44,45,49,50];
EEG.artifact.heart = [16,19];

subj = [1,2,4:9,12:24,26,30,31,33,34,46,51];
figure; pl=1;  row = 12; col = 12;mxfrq = 50; 
for c = 1:length(subj)
    cc=subj(c);
    sbplot(row,col,pl);
    topoplot(EEG.icawinv(:,cc),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.65);
    title(int2str(cc));
    sbplot(row,col,[pl+1 pl+2]);     
    [Pxx freqs] = pwelch(EEG.icaact(cc,:),256,128,512,EEG.srate);hold on;        
    % log this data
    Pxx = 10*log10(Pxx);
    ph = plot(freqs(find(freqs>3&freqs<mxfrq)),Pxx(find(freqs>3&freqs<mxfrq)),'b-','linewidth',2);
    
    title(int2str(cc));
    set(gca,'xtick',[0:10:mxfrq]);
    set(gca,'xlim',[3 mxfrq]);    set(gca,'xgrid','on');    pl=pl+3;
    plot([10 10],[get(gca,'ylim')],'r-');
    
end;
ph=textsc('mi83 Emotion Component Spectra from 248 chan data, set 1','title');set(ph,'fontsize',14);


% epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
eeglab
spath = '/data/common1/emotion/mi83/';
datsets = {'emo-1-248.set','emo-2-248.set','emo-3-248.set','emo-4-248.set','emo-5-248.set','emo-6-248.set'};
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite'};

%%**  get pressevents first to mark first button press
    ALLEEG=[];EEG=[];
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

