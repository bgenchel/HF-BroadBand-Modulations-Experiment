% procedures to import emotion data on jo74
% 1600 blocks seems to be the limit on import, at least with 2GB
%  event channel  =265
%  ref'ed to 135,213.
% chan 257 is the button press
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  ***  Next use Button info to import and epoch real data  ***  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Import emotion data in three chunks, need to merge log file each time by shifts.
eeglab
rawdat = '/data/common1/emotion/jo74/jo_emo.bdf';
spath = '/data/common1/emotion/jo74/';
elpfile = '/data/common1/emotion/jo74/dummy-small-256.elp';

EEG = pop_readbdf(rawdat, [1:833:834] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 1 ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{elpfile , 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',0);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-1-252.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[]; EEG=[]; clear ans
pack
%%%%  block 2   %%%%%%%%%
EEG = pop_readbdf(rawdat, [835:995:1830] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 2  ');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-29);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG, 'emo-2-252.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 3  %%%%%%%
EEG = pop_readbdf(rawdat, [1833:830:2663] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 3');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-49);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-3-252.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
%%%%%%  block 4  %%%%%%%
EEG = pop_readbdf(rawdat, [2664:519:3183] ,265,[135:78:213]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', 'emo part 4');
EEG = pop_select( EEG, 'nochannel',[257:264] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ elpfile, 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'B12'}, 'convert',{ 'chancenter',[],0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',-69);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'nochannel',[253:256] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 EEG = pop_eegfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
EEG = pop_saveset( EEG, 'emo-4-252.set', spath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[]; clear ans
pack
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                                                                                    

 % in linux terminal: 
 cat emo1epochs.fdt emo2epochs.fdt emo3epochs.fdt  > allemoeps.fdt
 cat tb1eps.fdt tb2eps.fdt   > alltbeps.fdt
 cat emo1-252.fdt emo2-252.fdt emo3-252.fdt   > emo252-1-2-3.fdt
  % run ica
/data/common/matlab/ica_linux2.4 < allemo3.sc
/data/common/matlab/ica_linux2.4 < alltb.sc
/data/common/matlab/ica_linux2.4 < emo256ica.sc
/data/common/matlab/ica_linux2.4 < emo252ica2.sc
/data/common/matlab/ica_linux2.4 < allimageica.sc


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


EEG = pop_loadset( 'emo-1-252.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG = pop_epoch( EEG, {  'eclose'  }, [0 10], 'newname', 'eclose3', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', 'eclose3', 'overwrite', 'on');
EEG = pop_rmbase( EEG, [0 10000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);



figure;
for c = 1:size(EEG.data,1)
psd(EEG.data(c,:),512,EEG.srate,512);hold on;
end;



sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); % final**
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
EEG.icaweights=wts;
EEG.icasphere=sph;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

subj = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
subj = [1:30];
figure; pl=1;  row = 7; col = 8;
for c = 1:length(subj)
    cc=subj(c);
    subplot(row,col,pl);
    topoplot(EEG.icawinv(:,cc),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.5);
    title(int2str(cc));
    subplot(row,col,pl+1);    
    psd(EEG.icaact(cc,:),512,EEG.srate,512);hold on;
    title(int2str(cc));
    set(gca,'xtick',[5:5:40]);
    set(gca,'xlim',[1 40]);
    set(gca,'xgrid','on');
    pl=pl+2;
end;
ph=textsc('jo74 Emotion Component Spectra from 252 chan data, set 1','title');
set(ph,'fontsize',14);




 [trmeans] = changrads(EEG.icawinv,'/data/common1/emotion/jo74/mylocs.loc');
 

 
 % epoch on emotion til 300 seconds (or shorter, data permitting)
% long emo length = 235 sec
eeglab
spath = '/data/common1/emotion/jo74/';
datsets = {'emo-1-252.set','emo-2-252.set','emo-3-252.set','emo-4-252.set'};
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
%emos = {'prebase','postbase'};
% do prebase and post base separately (without button press)
%%**  get pressevents first to mark first button press
ALLEEG=[];EEG=[];
for ds = 3:length(datsets)
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

