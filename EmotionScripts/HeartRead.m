% imports emo data for heart beat and button press channels only.

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed

% Variables

datfiles = {'tl81.bdf','mi83.bdf','ms82-2.bdf','js75.bdf','kw78.bdf','jo82.bdf','kl80.bdf','ar81.bdf','eb79.bdf','dg75.bdf','an82.bdf','jw84.bdf','tv81.bdf','sr81.bdf','an70.bdf','sg75.bdf','mr72.bdf','dk74.bdf','dn86.bdf','mr71.bdf','md85.bdf','mr72-2.bdf','cj82-2.bdf','kc66.bdf','ts79.bdf','es76.bdf','mm78.bdf','ab75.bdf','hs83.bdf','ps82.bdf','as82.bdf','ef76.bdf','jl83.bdf' , 'rr83.bdf', 'jw77.bdf'};

 
           
%% (done and saved)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
easysubjs = [1:22,24:35]; 
for nxx = 1:length(easysubjs)
    nx=easysubjs(nxx);
    clear blocks
    x=openbdf([fullpaths{nx},datfiles{nx}]);
    numframes = x.Head.NRec;
    strt = 0;
    for bl = 1:ceil(numframes/1000)
        if bl ~= ceil(numframes/1000)
            blocks{bl} = [strt strt+999];
        else
            blocks{bl} = [strt numframes];
        end;
        strt = strt + 1000;
    end;    
    ALLEEG=[];EEG=[]; 
    for bl = 1:length(blocks)
        EEG = pop_biosig([fullpaths{nx},datfiles{nx}],'channels' ,[1:x.Head.NS], 'blockrange' , blocks{bl} );% E3 and G23
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, bl+4,  'setname', [fullpaths{nx}(end-4:end-1),' Heart Channels']);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, bl+4);
        EEG = pop_select( EEG, 'nochannel',[257:264] );
        EEG = pop_select( EEG, 'nochannel',[1:254] ); % take out heart because 'types' is ridiculous
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, bl+4,  'overwrite', 'on');
        delevs = [];
        for ev = 1:length(EEG.event) 
            if EEG.event(ev).type == 30 
                for evv = ev-1:-1:1
                    if strcmp(EEG.event(evv).type,'prebase_instruct')
                        EEG.event(ev).type = 'prebase'; break
                    elseif  strcmp(EEG.event(evv).type,'postbase_instruct')
                        EEG.event(ev).type = 'postbase'; break
                    end;
                end;
            end;       
            if EEG.event(ev).type == 1
                EEG.event(ev).type = 'instruct1';
            end;
            if EEG.event(ev).type == 2
                EEG.event(ev).type = 'prebase_instruct';
            end;
            if EEG.event(ev).type == 3
                EEG.event(ev).type = 'instruct2';
            end;
            if EEG.event(ev).type == 5
                EEG.event(ev).type = 'instruct3';
            end;
            if EEG.event(ev).type == 7
                EEG.event(ev).type = 'instruct4';
            end;
            if EEG.event(ev).type == 25
                EEG.event(ev).type = 'postbase_instruct';
            end;
            if EEG.event(ev).type == 6
                EEG.event(ev).type = 'relax';
            end;
            if EEG.event(ev).type == 8
                EEG.event(ev).type = 'enter';
            end;
            if EEG.event(ev).type == 24
                EEG.event(ev).type = 'exit';
            end;
            if EEG.event(ev).type == 9 
                EEG.event(ev).type = 'awe';
            end;
            if EEG.event(ev).type == 10 
                EEG.event(ev).type = 'frustration';
            end;
            if EEG.event(ev).type == 11 
                EEG.event(ev).type = 'joy';
            end;
            if EEG.event(ev).type == 12 
                EEG.event(ev).type = 'anger';
            end;
            if EEG.event(ev).type == 13 
                EEG.event(ev).type = 'happy';
            end;
            if EEG.event(ev).type == 14 
                EEG.event(ev).type = 'sad';
            end;
            if EEG.event(ev).type == 15 
                EEG.event(ev).type = 'love';
            end;
            if EEG.event(ev).type == 16 
                EEG.event(ev).type = 'fear';
            end;
            if EEG.event(ev).type == 17 
                EEG.event(ev).type = 'compassion';
            end;
            if EEG.event(ev).type == 18 
                EEG.event(ev).type = 'jealousy';
            end;
            if EEG.event(ev).type == 19 
                EEG.event(ev).type = 'content';
            end;
            if EEG.event(ev).type == 20 
                EEG.event(ev).type = 'grief';
            end;
            if EEG.event(ev).type == 21 
                EEG.event(ev).type = 'relief';            
            end;
            if EEG.event(ev).type == 23 % yes, this is the way the tracks were numbered
                EEG.event(ev).type = 'disgust';
            end;
            if EEG.event(ev).type == 22 
                EEG.event(ev).type = 'excite';
            end;
            if EEG.event(ev).type > 30000
                delevs = [delevs ev];
            end;    
        end;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, bl+4);
        EEG = pop_editeventvals(EEG, 'delete',delevs);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, bl+4);
    end;
    EEG = pop_mergeset( ALLEEG, [1:length(ALLEEG)], 0);
    EEG = pop_saveset( EEG, 'Heart2.set', newpaths{nx});  
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% events saved in Heart2 files in 'newpaths' which were taken from ButtonEvents.set ****

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % fullpaths is a cell array of strings of form {'/data/common4/emotion/tl81/'...}
% emos is a cell array of emotion event names, such as {'awe','sad','anger',...} 
% press1 is the event to indicate that they are feeling emotion
% press an event where subjects are trying to express the emotion through the button

% for emotion task
addpath('/home/julie/MatlabScripts/emotion') 
DataInfo    % this matlab file loads all subject info needed 
addpath('/data/common2/emotion')
fullpaths = newpaths;


datset = 'Heart2.set'; % same name for all subjects, directory defines subject

% to find initial button press (press1) and select data between that and 'exit' button press
for nx = 1:35  %length(fullpaths) % loop through all subjects    
    EEG = pop_loadset( datset,fullpaths{nx} );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
    for em = 1:length(emos) % loop through all emotions
        StartTime = [];
        StopTime= [];
        for ev = 1:length(EEG.event)
            if strcmp(emos{em},EEG.event(ev).type)
                for evv = ev+1:length(EEG.event) % look for subsequent events to emo start
                    if strcmp('press1',EEG.event(evv).type)
                        StartEvent = evv; % option 1: epoch on events
                        StartTime = EEG.event(evv).latency; % option 2: select times
                        for evvv = evv+1:length(EEG.event)% look for subsequent events to feeling
                            if strcmp('exit',EEG.event(evvv).type)
                                StopTime = EEG.event(evvv).latency;
                                break
                            end;
                        end; 
                        break;
                    end;                    
                end;
                break; % only one event should be found, so just quit when it finds it.
            end;
        end;
        if ~isempty(StopTime) % option 2
            EEG = pop_select(EEG,'point',[StartTime StopTime]);
            EEG = eeg_checkset( EEG );
            EEG = pop_rmbase( EEG, [], [1:EEG.pnts] );
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
            EEG = pop_saveset( EEG,[emos{em},'Heart.set'],fullpaths{nx}); % specify name
            ALLEEG = pop_delset( ALLEEG, [2] );
            [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
        end;       
        
        %if ~isempty(StartEvent) % else, for option 1:
        %    EEG = pop_epoch(EEG,{},[0 StopTime],'eventindices',StartEvent); % one epoch
        %    EEG = pop_rmbase( EEG,[0 StopTime]);
        %    EEG = eeg_checkset( EEG );
        %    EEG = pop_saveset( EEG,[emos{em},'Heart.set'],fullpaths{nx}); % specify name
        %    ALLEEG = pop_delset( ALLEEG, [2] );
        %    [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
        %end;
    end;    
end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find heart beats %%
for nx=1:35 %length(fullpaths) 
  chn = 1; chunksize = EEG.srate*5; %also try at chunksize of 10
                                    
  for em  = 1:length(emos)
    EEG = pop_loadset( [emos{em},'Heart.set'], fullpaths{nx});
    clear hrate
    m=1;
    index = [];
    for x = 1:EEG.srate*.25:size(EEG.data,2)-chunksize  %
      dat = EEG.data(chn,x:x+chunksize);
      dat = dat - mean(dat);
      keeppnt = [];
      for xx = 2:length(dat)-1
        val = dat(xx); 
        preval = dat(xx-1);
        postval = dat(xx+1);
        if val < preval & val < postval & val < mean(dat) - .5*max(abs(dat))
          keeppnt = [keeppnt x+xx];
          index = [index x+xx];
        end;      
      end;
      %keeppnt = unique(keeppnt);
      index = unique(index);
      hrate(1,m) = length(keeppnt); m=m+1; %hrate will be overall bpm
    end;
  str = ['save ', fullpaths{nx},[emos{em}, 'HeartInfo.mat hrate index']]; eval(str);
  end;
                                                  
end;

%fullpaths{9} incomplete -- ended at awe (12/15)
%skip fullpaths{3} and fullpaths{23}
%still need to save time delay between peaks (entries in keeppnts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate HRV %%
srate = 256;
maxfrq = 5;
for nx = 1:35  %length(fullpaths)
    for em = 1:length(emos)
        hrvBPM = []; PS = [];
        str = ['load ', fullpaths{nx},[emos{em}, 'HeartInfo.mat']];
        eval(str);
        for i = 2:length(index)
            hrvBPS = (index(i) - index(i-1)) / 256 * 60; %beats per sec
            hrvBPM = [hrvBPM hrvBPS];
        end;
        for sample = 1:max(hrvBPM)
            place =find(index > sample);place = place(1);
            hrdat(1,sample) = hrvBPM(place);
        end;        
        [Pxx,F] = pwelch(hrdat,srate*100,srate/4,srate*80,srate);
        Pxx = 10*log10(Pxx); % convert to log power
        figure; plot(F(find(F < maxfrq)),Pxx(find(F < maxfrq)));        
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for jo74 trial
EEG = pop_importpres(EEG,  '/data/common1/emotion/jo74/jo-1-14-emotionimage.log', 'Code', 'Time',0);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
siglevel = .01;wsize = 128; freqs = [.0001:.0001:.0033]
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief','bdisgust','bexcite'};
ALLEEG=[];
for nx = 1:length(paths)
    EEG = pop_loadset( 'HeartOnly.set',['/data/common2/emotion/',paths{nx}] );
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG.data = rmbase(EEG.data,EEG.pnts,1:size(EEG.data,2)); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % 
    for k = 1:length(emos)
        endt = (EEG.event(end).latency-EEG.event(end-4).latency)/256;
        EEG = pop_epoch( EEG, {  emos{k}  }, [-.5  endt], 'newname', emos{k}, 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emobuts{k});
        EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    end;
    ALLEEG(1)=[];
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
    for k = 1:15
              EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
  
    fprintf('\n\nEmotion: %s',emos{k});
    fprintf('\nMean: %s',num2str(mean(EEG.data(1,:))));
    fprintf('\nStd: %s',num2str(std(EEG.data(1,:))));
    end;
    
    figure; 
    for k = 1:length(emos)
        EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
        subplot(4,4,k)
%[Pxx freqs] = psd(EEG.data(ht,:),81920,256,5120);
%figure;plot(freqs(find(freqs<.0033)),Pxx(find(freqs<.0033)))
%figure;plot(freqs(find(freqs>.0033&freqs<.05)),Pxx(find(freqs>.0033&freqs<.05)))
%figure;plot(freqs(find(freqs>.05 &freqs<.15)),Pxx(find(freqs>.05 &freqs<.15)))
%figure;plot(freqs(find(freqs>.15 &freqs<.45)),Pxx(find(freqs>.15&freqs<.45)))
ht=1;
figure; [ersp,itc,powbase,times,freqs,erspboot,itcboot]=newtimef( EEG.data(ht,:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000],EEG.srate,0,'freqs',[0:1:4],'baseline',EEG.xmax*1000,'plotitc','off');
        subplot(4,4,k)
        imagesc(times,freqs,ersp,[-40 40]);      
        title(emos{k})
    end;
        comp_ersp(:,:,ht)=ersp;
        ersp_boot(:,:,ht)= erspboot';
            
ALLEEG=[];EEG=[];

        
       , 'alpha',siglevel ,'winsize',wsize
 , 'alpha',siglevel
 
 
 
 % to epoch on the button press only data:
emonames = {'awePress', 'frustrationPress','joyPress','angerPress','sadPress','happyPress','fearPress','lovePress','jealousyPress','compassionPress','emabarrassPress','contentPress','griefPress','reliefPress'};
for k = 1:length(emonames)   
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    EEG = pop_epoch( EEG,{'press'} , [-2  2], 'newname', emonames{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emonames{k},'overwrite','on');
    EEG = pop_rmbase( EEG,[-2000 2000]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
end;
eeglab redraw

%  Plot ERP images of button, and both heart channels
figure;  pl=1;
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
%    heart = EEG.data(1,:)-EEG.data(2,:);
%    subplot(4,4,pl)
 %erpimage( heart, ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k}, 0, 1 ,'erp','limits',[-1000 2000 NaN NaN NaN NaN NaN NaN],'phasesort',[100 0 1 2],'caxis',  [-200 200]);pl=pl+1;
    subplot(7,4,pl)
 erpimage( EEG.data(1,:), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k}, 0, 1 ,'erp','limits',[-1000 2000 NaN NaN NaN NaN NaN NaN],'phasesort',[100 0 1 2],'caxis',  [-200 200]);pl=pl+1;
    subplot(7,4,pl)
 erpimage( EEG.data(2,:), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts),emos{k}, 0, 1 ,'erp','limits',[-1000 2000 NaN NaN NaN NaN NaN NaN],'phasesort',[100 0 1 2],'caxis',  [-200 200]);pl=pl+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epoch a heart dataset on the low spike of each heart beat
ALLEEG=[];EEG=[];
EEG = pop_loadset( 'HeartOnly.set', ['/data/common2/emotion/',paths{nx}]);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
figure; 
m=1;pl=1;
for r = 250:250:size(EEG.data,2)-250  % takes a LONG time
    mn = mean(EEG.data(1,r:r+250));
    newEEG = EEG.data(1,:)-mn;
    [a b] =min(newEEG(1,r:r+250));
    plot(newEEG(1,r:r+250))
    EEG.event(end+1) =  EEG.event(end);% appends events to the end
    EEG.event(end).latency = (r-1)+b;
    EEG.event(end).type = 'fake'; 
end;
EEG = eeg_checkset(EEG, 'eventconsistency');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

[outeeg,indices] = pop_epoch (EEG, {'fake'}, [-.3 .5], 'epochinfo','yes');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, outeeg,CURRENTSET);
EEG = pop_rmbase( EEG, [-300 500]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

figure;
for k = 1:size(EEG.data,3)
    plot(EEG.times,EEG.data(1,:,k),'g');hold on;
end;
ph=plot(EEG.times,mean(EEG.data(1,:,:),3),'r');
set(ph,'linewidth',2);
kern = mean(EEG.data(1,:,:),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use a single heart beat kernel
ALLEEG=[];EEG=[];
EEG = pop_loadset( 'HeartOnly.set', ['/data/common2/emotion/',paths{nx}]);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
output = zeros(1,size(EEG.data,2)-length(kern));
for x = 1:size(EEG.data,2)-length(kern)
    output(1,x:(x-1)+length(kern)) = kern.*EEG.data(1,x:(x-1)+length(kern));
if x == size(EEG.data,2)/4| x == size(EEG.data,2)/3| x == size(EEG.data,2)/2
    fprintf('\nDone %s out of %s blocks',int2str(x), int2str(size(EEG.data,2)));
    end;
end;
figure; plot(output);

% import heart beat events to continuous emo1

EEG = pop_loadset( 'emo-1-252.set', '/data/common/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
for p = 1:1500
EEG.event(end+1) = hbeat{p};
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end;
[outeeg,indices] = pop_epoch (EEG, {'fake'}, [-.2 .5], 'epochinfo','yes');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, outeeg,CURRENTSET);
EEG = pop_rmbase( EEG, [-199 0]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);



