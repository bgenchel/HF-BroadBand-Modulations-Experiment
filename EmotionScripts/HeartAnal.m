% takes heart-only datasets from emotion expts for analysis

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
fullpaths = newpaths;

datset = 'Heart2.set'; % same name for all subjects, directory defines subject

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create datasets for each emotion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find initial button press (press1) and select data between that and 'exit' button press

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
%skip fullpaths{3} and fullpaths{23}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find large spike direction for each subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdsubjs=[1:2,4:22,24:31,33:35];
for nxx = 1:length(gdsubjs)  % skips 3 & 23,32
    nx = gdsubjs(nxx);
    EEG = pop_loadset( [emos{11},'Heart.set'], fullpaths{nx});
    chunk = EEG.data(1,256:256+1024*2);
    chunk = chunk - mean(chunk);
    exvals(1) = min(chunk); exvals(2) = max(chunk);
    if abs(exvals(1)) > exvals(2)
        ori(nx) = -1;% flip EEG data to analyze later
    end;
    chunk2 = EEG.data(2,256:256+1024*2);
    chunk2 = chunk2 - mean(chunk2);
    figure; sbplot(1,2,1);
    plot(chunk); title('Channel 1');
    sbplot(1,2,2); plot(chunk2); title('Channel 2');
    chns(nx) = input('which electrode is better? 1 or 2?'); 
    ori(nx) = input('and which orientation? 1 or -1?'); close
end;
str = ['save ', fullpaths{nx}(1:end-5),'HeartInfo.mat chns ori']; eval(str);
% three means that neither channel has a spike opposite from 
% long broad peak, thus will be confusing. So ignore these 
% or use special parameters on chan 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find large spike direction for each subject (using pulse IC, not heart channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdsubjs=[1:2,4:10,12:15,18:20,23:26,28:35];
for nxx = 5:length(gdsubjs)  % skips 3 & 23,32
    nx = gdsubjs(nxx);  ALLEEG=[];EEG=[];
    EEG = pop_loadset( 'sources.set', fullpaths{nx});
    heart = EEG.heart;
    EEG = pop_loadset( [emos{11},'.set'], fullpaths{nx});
    chunk = EEG.icaact(heart(1),256:256+1024*2);
    chunk = chunk - mean(chunk);
    exvals(1) = min(chunk); exvals(2) = max(chunk);
    if abs(exvals(1)) > exvals(2)
        ori(nx) = -1;% flip EEG data to analyze later
    end;
    figure; 
    plot(chunk); set(gca,'xlim',[1 length(chunk)]);title('IC 1');
    if length(heart) == 2
        chunk2 = EEG.icaact(heart(2),256:256+1024*2);
        chunk2 = chunk2 - mean(chunk2);
        figure; plot(chunk2); title('IC 2');
    end;
    chns(nx) = input('which IC is better? 1 or 2?'); 
    ori(nx) = input('and which orientation? 1 or -1?'); 
    close;close;
end;
str = ['save ', fullpaths{nx}(1:end-5),'PulseInfo.mat h chns ori']; eval(str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find heart rates (from heart channel data)%%
str = ['load ', fullpaths{nx}(1:end-5),'HeartInfo.mat chns ori']; eval(str);
nsecs = 6;% 6 sec chunks to find peaks (needs to be even)
chn = 1;   % channel of interest                
gdsubjs=[1:2,4:14,16:17,19:22,25:28,30:31,33:35];

for nxx = 2:length(gdsubjs)  % skips 3 & 23,32
    nx = gdsubjs(nxx);
    if chns(nx) < 3
    clear hrate movinghrate
    for em  = 1:length(emos)
        EEG = pop_loadset( [emos{em},'Heart.set'], fullpaths{nx});
        chunksize = EEG.srate*nsecs;
        allfrms = [];movinghb = [];lostframes = [];
        for x = chunksize/2:EEG.srate/2:size(EEG.data,2)-chunksize/2             
             try
                datchunk= EEG.data(chns(nx),x-chunksize/2+1:x+chunksize/2)*ori(nx);
                datchunk = datchunk - mean(datchunk); % subtract baseline
                b = regress(datchunk',[ones(1,length(datchunk))' [1:length(datchunk)]']); 
                y = b(2)*[1:length(datchunk)] + b(1);
                datchunk = datchunk-y; % take out linear trend
                [mag frm] = findpeaks(datchunk,'minpeakheight',max(abs(datchunk))/2);
                allfrms = [allfrms frm+(x-1)];
                movinghb = [movinghb length(frm)/nsecs*60];%div by nsec in chunk * by 60 sec/min
            catch
                lostframes = [lostframes x];
            end;
        end;
        allfrms = unique(allfrms);
        hrate(1,em) = length(allfrms)/(size(EEG.data,2)/EEG.srate)*60;%hrate will be overall bpm
        % remove spurious values
        delidx = find(movinghb< (mean(movinghb) - std(movinghb)*1.5)|movinghb> mean(movinghb) + std(movinghb)*1.5);
        movinghrate{em} = movinghb;
        hrate(2,em) = mean(movinghb); % another way to calculate it
    end;
    str = ['save ', fullpaths{nx},'HeartInfoJO.mat hrate movinghrate lostvec']; eval(str);
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find heart rates (from pulse ICs)%%
str = ['load ', fullpaths{nx}(1:end-5),'PulseInfo.mat h chns ori']; eval(str);

for nxx = 1:length(gdsubjs)  % skips 3 & 23,32
    nx = gdsubjs(nxx);
    clear hrate 
    for em  = 1:length(emos)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset( [emos{em},'.set'], fullpaths{nx});
        if ori(nx) < 0
             EEG.icaact(h{nx}(chns(nx)),:) =  EEG.icaact(h{nx}(chns(nx)),:)*-1;% invert
        end;
        run = [];beat = [];lastx = 19;
        for x = 20:size(EEG.icaact,2) - 20            
            if EEG.icaact(h{nx}(chns(nx)),x) > EEG.icaact(h{nx}(chns(nx)),x - 1)
                if x == lastx + 1 % continues a run
                    run = [run x];
                else
                    if ~isempty(run)
                        stepsize = EEG.icaact(h{nx}(chns(nx)),run(end)) - EEG.icaact(h{nx}(chns(nx)),run(1));
                        if stepsize > 1.5 % large enough to be a beat?
                            beat = [beat run(end)];
                        end;
                        run = [];
                    end;                    
                end;
                lastx = x;
            end;            
        end;
        hrate(1,em) = (length(beat)/(size(EEG.icaact,2)/EEG.srate))*60;%hrate will be overall bpm
    end;      
    str = ['save ', fullpaths{nx},'PulseInfoJO.mat hrate movinghrate']; eval(str);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantify hrate for each emotion across subjects

badnx = [];clear emrates emrates2
for nxx = 1:length(gdsubjs)  % skips 3 & 23
    nx = gdsubjs(nxx);
    if chns(nx) < 3
    str = ['load ', fullpaths{nx},'PulseInfoJO.mat hrate movinghrate']; eval(str);    
    %str = ['load ', fullpaths{nx},'HeartInfoJO.mat hrate movinghrate lostvec']; eval(str);    
    %if mean(hrate(2,:)) > 100
    %    badnx = [badnx nx];
    %end;
     for em = 1:length(emos)    
        emrates(em,nxx) = hrate(1,em);
        %emrates2(em,nxx) = hrate(2,em);
    end;
    end;
end;
figure; bar(emrates);
figure; bar(emrates2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot moving hr for each subject

gdsubjs=[1:2,4:14,16:17,19:22,25:28,30:31,33:35];
badnx = [];clear emrates emrates2
figure; row = 6; col = 1; pl = 1;
for nxx = 1:length(gdsubjs)  % skips 3 & 23
    nx = gdsubjs(nxx);
    str = ['load ', fullpaths{nx},'HeartInfoJO.mat hrate movinghrate lostvec']; eval(str);  
    hrates(nx,:) = hrate(2,:);
    if ~isempty(movinghrate{1})
        if pl > row*col
            figure; pl = 1;
        end;
        sbplot(row,col,pl); pl = pl+1;
        allhr = [];
        for em = 1:length(emos) 
            allhr = [allhr,movinghrate{em}]; % collect all hr for total var
            sbjvars(nxx,em) = var(movinghrate{em});
            ph = plot([1:length(movinghrate{em})],movinghrate{em}); hold on;
            set(ph,'color',cols(em,:));lens(1,em) = length(movinghrate{em});
        end;
        set(gca,'xlim',[0 max(lens)+1]);
        set(gca,'ylim',[50 150]);
        totalvar(1,nxx) = var(allhr);
    end;
end;
str = ['save ', fullpaths{nx}(1:end-5),'AllSubjHeartInfo.mat hrates']; eval(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate overall hr with emotion valence/arousal (no correlation found)


gdsubjs=[1:2,4:14,16:17,19:22,25:28,30:31,33:35];
w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
erating = w.emoval;rw = 1;  % w.emoval or w.emoactiv
%erating = w.emoactiv; rw = 2; % w.emoval or w.emoactiv
clear cr
useemos = [1:7,9:15];% take out compassion for correlation
            
for nxx = 1:length(gdsubjs)  % skips 3 & 23
    nx = gdsubjs(nxx);
    str = ['load ', fullpaths{nx},'HeartInfoJO.mat']; eval(str);  

    hrcr(1,nx) = corr2(hrate(1,useemos),w.emoval(useemos));
    hrcr(2,nx) = corr2(hrate(1,useemos),w.emoactiv(useemos));
end;
