% procedure to follow for Emotion data prep

% Read in raw data file
% change EEG.srate appropriately for merge with log
clear wrlat rtlat
nbevents = length(EEG.event);
for index = 1:nbevents
rtlat (1,index) =   EEG.event(index).latency;
end;rtlat
EEG.srate = 500.127;
% merge with log and check again.
% if correct, delete channels 3 and 74
% change 
EEG.srate = 500;
% resample to 250 Hz
% filter 0 50
% name and fill in About this dataset
% save

% add zeros to end so you can get the last event
EEG.data(:,end+1:end+50000)=zeros(size(EEG.data,1),50000);
EEG.pnts = size(EEG.data,2);
EEG.xmax = EEG.pnts/EEG.srate;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw

% epoch on each emotion and save

name = {'awe' 'frust'  'joy'  'anger'  'sad'   'happy' 'surprise' 'hatred' 'fear' 'love'    'compassion'   'jealousy'  'relief' }; % for jo74
filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set'   'happy.set' 'surprise.set' 'hatred.set' 'fear.set' 'love.set'    'compassion.set'   'jealousy.set'  'relief.set' }; % for jo74



name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'  'hatred'   'grief'   'relief' }; % for kl80
filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'hatred.set'   'grief.set'   'relief.set' }; % for kl80


% load cont data set to epoch in slot 1 (others empty)
paths = {'/data/common1/emotion/kl80/','/data/common1/emotion/ap80/','/data/common1/emotion/rr80/','/data/common1/emotion/jo74/'};
subj = 1;
for em = 12:length(name)  % start at 1 or 2 for first or second session
    [outeeg,indices] = pop_epoch (EEG, {name{em}}, [-1 275], 'epochinfo','yes','newname',name{em});  
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, outeeg);
    EEG.data(:,:) = rmbase(EEG.data(:,:),0,1:250); % remove baseline from first 200 samples
% Create 'emotion' field--------------------------------------------------
    placehold = ones(length(EEG.event),1);
    EEG = pop_editeventfield(EEG,'indices',[1:length(EEG.event)],'emotion',placehold);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% fill all emotion fields with an emotion name
    for e = 1:length(EEG.event)
        EEG.event(e).emotion = name{em};
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    EEG = pop_saveset( EEG,filnm{em} , paths{subj});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG = pop_delset( ALLEEG, [2] );
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;

%  replace all data after space bar press with nans (one epoch datasets)
for em = 1:length(filnm) 
    EEG = pop_loadset( filnm{em},paths{subj});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG.data(:,EEG.event(2).latency:size(EEG.data,2))=NaN;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG,filnm{em} , paths{subj});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG = pop_delset( ALLEEG, [1] );
end;
% to check:
figure; plot(EEG.data(72,1:size(EEG.data,2)));  % check for button threshold
%%%%%%%%%