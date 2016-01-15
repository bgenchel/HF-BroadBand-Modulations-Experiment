% Plots several emotion Time Freqs for a given component, plus the topo map and spectrum

%  goes through all emotions and does single trial timef on all of them
subj1 = [2,3,5,7,8,9,10,11,12,14,17,18,20,23,25,27,29,30,37,38,40];    %kl80
subj2 = [7:12,14,17:20,22,24,26];  %ap80
subj3 = [2,3,5:18,20,21,23:25,27,28,29,31,33,36,37,38,41,43,48,49,50];  % rr80
    gdcomps = {subj1 subj2 subj3 };% subj3 subj4};
             
paths = {'/data/common1/emotion/kl80/','/data/common1/emotion/ap80/','/data/common1/emotion/rr80/','/data/common1/emotion/jo74/'};


%  ap80 
name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'   'grief'   'relief' }; % for rr80

%filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'grief.set'   'relief.set' }; % for rr80
%filnm = { 'awepress.set' 'frustpress.set'  'joypress.set'  'angerpress.set'  'sadpress.set' 'surprisepress.set'   'happypress.set' 'fearpress.set' 'lovepress.set' 'jealousypress.set' 'compassionpress.set'  'griefpress.set'   'relief2press.set' }; % for ap80
filnm = { 'awe2press.set' 'frust2press.set'  'joy2press.set'  'anger2press.set'  'sad2press.set' 'surprise2press.set'   'happy2press.set' 'fear2press.set' 'love2press.set' 'jealousy2press.set' 'compassion2press.set'  'grief2press.set'   'relief2press.set' }; % for rr80

%  rr80 
name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'   'grief'   'relief' }; % for rr80

%filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'grief.set'   'relief.set' }; % for rr80
filnm = { 'frustpress.set'  'joypress.set'  'angerpress.set'  'sadpress.set' 'surprisepress.set'   'happypress.set' 'lovepress.set' 'jealousypress.set' 'compassionpress.set'  'griefpress.set'   'reliefpress.set' }; % for rr80

%  Plot each component, all emotions
% load in chanlocs and adjust   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
subj = 2;  p=3;
cd (paths{subj})
load wts 
load sph
EEG = pop_loadset( filnm{p},paths{subj});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
   
    % now recalculate ica.act
    EEG.icaweights = wts;
    EEG.icasphere = sph;
    EEG.icaact = []; EEG.icawinv=[];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
%EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/rr80/rr80.elp', 'filetype',''}, 'forcelocs',{0, 'X', 'Cz'}, 'convert',{ 'chancenter',[],1});
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/ap80/ap80.elp', 'filetype',''}, 'forcelocs',{0, 'X', 'Cz'}, 'convert',{ 'chancenter',[],1});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load AllpressTF.mat
%load Allpress2TF.mat
%load AllemoTF.mat
fr = find(freqs<30);

%eoi = {'Awe','Frustration','Anger','Sadness'};%'Sadness'};'Surprise','Happy'
eoi = {'Frustration','Joy','Anger','Sadness','Love','Jealousy','Relief'};%'Sadness'};'Surprise','Happy'
figure; pl = 2;  lim = 6.5;
%pp =  [1,2,3,4,7,8,11] ;    % rr80 % set of numbers of emotions to plot
pp =  [2,3,4,5,9,10,13] ;     % ap80% set of numbers of emotions to plot
cp = 1 ;     % component index to plot above button
col = length(pp)+1;
row = 2;
EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
subplot(row,col,1)
topoplot(EEG.icawinv(:,gdcomps{subj}(cp)),EEG.chanlocs,'electrodes','off');hold on
set(gca,'fontsize',[14]);
title(int2str(gdcomps{subj}(cp)));
load baseline.mat
frb = find(freqs<30);
subplot(row,col,col+1)
plot(freqs(frb),baseline(gdcomps{subj}(cp),frb));
hold on
set(gca,'xlim',[2 30]);
set(gca,'yticklabel',[]);
set(gca,'xgrid','on');
set(gca,'fontsize',[14]);
title('Spectrum')

for em = 1:length(pp)
    p = pp(em);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    load AllpressTF.mat
    tmpemo =  AllpressTF{1,p};  % 1 is ersp, 2 is itc
    tmpboot = Allpressboot{1,p};% 1 is ersp, 2 is itc
    %if em == length(pp)
    %    load Allpress2TF.mat
    %    tmpemo =  Allpress2TF{1,p};  % 1 is ersp, 2 is itc
    %    tmpboot = Allpress2boot{1,p};% 1 is ersp, 2 is itc
    %end;
    %load AllemoTF.mat
    fr = find(freqs<30);
    minboot = tmpboot(1,:,:);
    maxboot = tmpboot(2,:,:);
    minmask = repmat (minboot,[200,1,1]);
    maxmask = repmat (maxboot,[200,1,1]);
    tmpemo(find(tmpemo(:,:,gdcomps{subj}(cp))<minmask(:,:,gdcomps{subj}(cp))'&tmpemo(:,:,gdcomps{subj}(cp))>maxmask(:,:,gdcomps{subj}(cp))'))=0;
    subplot(row,col,pl)        
    imagesc(times,freqs(fr),tmpemo(fr,:,gdcomps{subj}(cp)),[-lim lim]); hold on
    set(gca,'xticklabel',[]);
    set(gca,'fontsize',[14]);
    title(eoi{em});
    if em == length(pp)
        colorbar; end;
        ALLEEG = pop_delset( ALLEEG, [2] );
        % load dataset for button press
        EEG = pop_loadset( filnm{p},paths{subj});
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        subplot(row,col,pl+col)
        plot(EEG.times,EEG.data(72,:,15)); hold on   
        set(gca,'fontsize',[14]);
        set(gca,'xlim',[-500 3000]);
        set(gca,'xticklabel',[0,1,2 3]);
        set(gca,'yticklabel',[]);
        title(eoi{em});
        pl = pl+1;
end;

