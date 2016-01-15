% runs crossf on jo74 emotion data
eeglab
path = '/data/common1/emotion/jo74/';
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
%emos = {'awePress.set', 'frustrationPress.set','joyPress.set','angerPress.set','sadPress.set','happyPress.set','fearPress.set','lovePress.set','jealousyPress.set','compassionPress.set','emabarrassPress.set','contentPress.set','griefPress.set','reliefPress.set'};

comment = 'jo74 All emos in order of presentation. Crossf on randomly epoched epochs: -1 1 crossf(EEG.icaact(subj1(index1),:),EEG.icaact(subj1(index2),:), EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000], EEG.srate, [3 .5], alpha, 0.01,padratio, 4,winsize,256); 63 freqs starting at 3 hz; epochs cut down to ~75 by, if over 75, removing 2/3 of diff from beginning and 1/3 from end'; 

sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
% use events from buttonpressonly file to epoch real data 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bemabarrass','bcontent','bgrief','brelief'};
 EEG = pop_loadset( 'ButtonOnly.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
 eeglab redraw

for k = 1:length(emos)
    if k == 14
        endt = 250;
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
        if EEG.event(ev).Event_Type(1) == 'S'
            seltime = EEG.event(ev).latency;
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

% run EmotionPresses256.m to get pressevents{}
% frust=71
ALLEEG=[]; EEG=[];
cd (path)
emos = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};
figure;
for e = 3:4%length(emos)
    EEG = pop_loadset( emos{e},path);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);   
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
  % make fake events every sec throughout  nf expts  (continuous data)for ica
tmpts = size(EEG.data,2);
for index = pressevents{e}(1).latency:EEG.srate*2:tmpts
      EEG.event(end+1) =  EEG.event(end);% appends events to the end
      EEG.event(end).latency = index;
      EEG.event(end).type = 'fake'; 
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
EEG = pop_rejkurt(EEG,0,subj1 ,5,4,0,1);
EEG = pop_jointprob(EEG,0,subj1 ,5,4,0,1);
%Before rejects:  109    71    64    78    79    90    86    98    62   106    66    74    85   87

%after rejects:    99    67    61    70    73    84    77    88    56    98    64    68    78   77
% mean and median number of epochs = 75
if length(EEG.epoch)>75
    num = length(EEG.epoch)-75;
    rembeg = round(num/3)*2;
    remend = round(num/3);
    EEG = pop_select( EEG, 'notrial',1:rembeg );
   [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
   endnum = length(EEG.epoch);   
    EEG = pop_select( EEG, 'notrial',endnum-(remend-1):endnum );
   [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end;
    
    coher=zeros(63,200,subj1(end));
    crsboot=zeros(63,2,subj1(end));
    crsangle=zeros(63,200,subj1(end));
    cohercell = cell (1,subj1(end));
    crsbootcell = cell (1,subj1(end));
    crsanglecell = cell (1,subj1(end));
    %-- always list lowest to highest comp numbers-------------
    %done =
    for index1=1:length(subj1)-1
        for index2=index1+1: length(subj1)
            [coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher, alltfX, alltfY] = crossf(EEG.icaact(subj1(index1),:),EEG.icaact(subj1(index2),:), EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000], EEG.srate, [3 .5], 'alpha', 0.01,'padratio', 4,'winsize',256);
            coher(:,:,subj1(index2)) = coh;
            cohercell{1,subj1(index1)} = coher;
            crsboot(:,:,subj1(index2)) = cohboot;
            crsbootcell{1,subj1(index1)} = crsboot;
            crsangle(:,:,subj1(index2)) = cohangles;
            crsanglecell{1,subj1(index1)}= crsangle;
            clf
        fprintf('\nShowing Comp: %i\n',subj1(index1));
        fprintf('\n vs     Comp: %i\n',subj1(index2));
        fprintf('\n in dataset: %s\n',emos{e});
        end
    end;
    emocoh{e} = cohercell;
    emoboot{e} = crsbootcell;
    emoangles{e} = crsanglecell;
    ALLEEG=[]; EEG=[];
end;
f =3;  s = 4;
ecoh{1} = emocoh{f};
ecoh{2} = emocoh{s};
eboo{1} = emoboot{f};
eboo{2} = emoboot{s};
eang{1} = emoangles{f};
eang{2} = emoangles{s};

load  /data/common1/emotion/jo74/RndEpochsXcoh.mat
load /data/common1/emotion/jo74/RndEpochsXang.mat  

emocoh{f} = ecoh{1};
emocoh{s} = ecoh{2};
emoboot{f} = eboo{1};
emoboot{s} = eboo{2};
emoangles{f} = eang{1};
emoangles{s} = eang{2};
    
save /data/common1/emotion/jo74/RndEpochsXcoh.mat emocoh freqsout timesout comment emoboot
save /data/common1/emotion/jo74/RndEpochsXang.mat emoangles freqsout timesout comment 


save /home/julie/EmoTF/jo74CrsEmos.mat emocoh emoboot  freqsout timesout comment
save /home/julie/EmoTF/jo74CrsAngEmos.mat  emoangles freqsout timesout comment


% Visualize Crosses
load  /data/common1/emotion/jo74/RndEpochsXcoh.mat
load /home/julie/EmoTF/jo74CrsEmos.mat 
%load /data/common1/emotion/jo74/ButtonLength.mat
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
%emos = {'jo74:  awe  (Button-Press Duration only)', 'jo74:  frustration  (Button-Press Duration only)','jo74:  joy  (Button-Press Duration only)','jo74:  anger  (Button-Press Duration only)','jo74:  sad  (Button-Press Duration only)','jo74:  happy  (Button-Press Duration only)','jo74:  fear  (Button-Press Duration only)','jo74:  love  (Button-Press Duration only)' ,'jo74:  jealousy  (Button-Press Duration only)','jo74:  compassion  (Button-Press Duration only)','jo74:  emabarrass  (Button-Press Duration only)','jo74:  content  (Button-Press Duration only)','jo74:  grief  (Button-Press Duration only)','jo74:  relief  (Button-Press Duration only)'};
emos = {'jo74:  AWE -- Randomly epoched during button presses', 'jo74:  FRUSTRATION -- Randomly epoched during button presses','jo74:  JOY -- Randomly epoched during button presses','jo74:  ANGER -- Randomly epoched during button presses','jo74:  SAD -- Randomly epoched during button presses','jo74:  HAPPY -- Randomly epoched during button presses','jo74:  FEAR -- Randomly epoched during button presses','jo74:  LOVE -- Randomly epoched during button presses' ,'jo74:  JEALOUSY -- Randomly epoched during button presses','jo74:  COMPASSION -- Randomly epoched during button presses','jo74:  EMBARRASS -- Randomly epoched during button presses','jo74:  CONTENT -- Randomly epoched during button presses','jo74:  GRIEF -- Randomly epoched during button presses','jo74:  RELIEF -- Randomly epoched during button presses'};

compsrun = subj1; ; 
lim = .7;
f = find (freqsout<30);
t = find (timesout>-2000 & timesout<3000);
for e = 9:10%length(emos)
    figure; clear forsum
    % first plot scalp maps in top row
    for mp = 2:length(compsrun)   
        subplot(length(compsrun),length(compsrun),mp)
        topoplot(EEG.icawinv(:,compsrun(mp)),EEG.chanlocs,'electrodes','off');
        title(int2str(compsrun(mp)));
    end;
    subplot(length(compsrun),length(compsrun),length(compsrun)+1)
    topoplot(EEG.icawinv(:,compsrun(1)),EEG.chanlocs,'electrodes','off');
    title(int2str(compsrun(1)));
    for kp = 2:length(compsrun)-1
        subplot(length(compsrun),length(compsrun),length(compsrun)*kp+kp)
        topoplot(EEG.icawinv(:,compsrun(kp)),EEG.chanlocs,'electrodes','off');
        title(int2str(compsrun(kp)));
    end;    
    p = ((length(compsrun)-1)/2)*length(compsrun);
    for w=1:length(compsrun)-1
        a = compsrun(w);
        b = compsrun(w+1:end);
        cohercell = emocoh{e};
        Diff = cohercell{1,a};
        crsbootcell = emoboot{e};
        bootmat1 = crsbootcell {1,a};   %chooses a 3D array from cell array
        minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
        minmask = repmat (minboot,[1,size(timesout,2),1]); %makes (63,200,71) matrix of minboot
        Diff(find(Diff <= minmask)) = 0;%finds and zeros cohmats less than minmask
        dd = Diff(f,t,:);
        forsum(:,:,w) = sum(dd,3);
        for m=w+1:length(b)+w
            % subplot(length(compsrun)-1,length(compsrun)-1,((length(compsrun)-1)*(w-1))+(m-1))
            subplot(length(compsrun),length(compsrun),((length(compsrun))*(w-1))+(m-1)+length(compsrun)+1)
            imagesc(timesout(t),freqsout(f),Diff(f,t,b(m-w)),[-lim lim]);
            hold on
            set(gca,'ticklength',[0.05 0.05]); 
            set(gca,'tickdir','out');      
            set(gca,'xtick',[-500:250:500]); 
            set(gca,'xticklabel',[]);    
            set(gca,'fontsize',5);    
            set(gca,'ytick',[5:5:30]);
            set(gca,'yticklabel',[]); 
            %set(gca,'xlim',[-lenset(1,e)/2 lenset(1,e)]);
            %plot([0 0],[0 30],'k-');
            title(int2str(b(m-w)))
        end;
    end;
   %set(gca,'xticklabel',[{-2 [] -1 [] 0 [] 1 [] 2}]);    
    set(gca,'yticklabel',[5:5:30]);    
    colorbar
    textsc(emos{e},'title');
    %axcopy
    newlim = max(max(sum(forsum)))-.4;
    sbplot(5,10,32:33)
    imagesc(timesout(t),freqsout(f),sum(forsum,3),[-newlim newlim]); colorbar
            set(gca,'ticklength',[0.03 0.03]); 
     title('Sum of Cross Cohs');
end;

%  Pull out features of crosses
%load  /data/common1/emotion/jo74/RndEpochsXcoh.mat
%load /home/julie/EmoTF/jo74CrsEmos.mat 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
% freqsout(1:20)
%    3.0000    3.7500    4.5000    5.2500    6.0000    6.7500    7.5000
%    8.2500    9.0000    9.7500   10.5000   11.2500   12.0000   12.7500
%   13.5000   14.2500   15.0000   15.7500   16.5000   17.2500

frqrange = find(freqsout > 14.75 & freqsout < 15.25);
tmrange = find(timesout > -1000 & timesout < 1000);
cutoff = .3;
for e = 1:length(emos)
    subj1a = [];subj1b=[];
    for cp1 = 1:length(subj1)-1
        for cp2 = 2:length(subj1)
            cohercell = emocoh{e};
            onecoh = cohercell{1,subj1(cp1)};
            crsbootcell = emoboot{e};
            bootmat1 = crsbootcell {1,subj1(cp1)};   %chooses a 3D array from cell array
            minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
            minmask = repmat (minboot,[1,size(timesout,2),1]); 
            onecoh(find(onecoh <= minmask)) = 0;%finds and zeros cohmats less than minmask
            onecoh = onecoh(frqrange,tmrange,subj1(cp2));
            oneval = mean(mean(onecoh));
            if oneval > cutoff
                subj1a(end+1) = subj1(cp1);
                subj1b(end+1) = subj1(cp2);
            end;
        end;
    end;
fprintf('\nif e ==');
fprintf('%i ',e);
fprintf('\n%% ');
fprintf('%s: ',emos{e});
fprintf('\nsubj1a = [%s];',int2str(subj1a));
fprintf('\nsubj1b = [%s];',int2str(subj1b));
fprintf('\nend;');
end;


%coh cutoff = .2; 8.25 Hz;  -1 1 sec around button press
if e ==1 
% awe: 
subj1a = [5];
subj1b = [13];
end;
if e ==2 
% frustration: 
subj1a = [5 10];
subj1b = [13  18];
end;
if e ==3 
% joy: 
subj1a = [5  5  5];
subj1b = [10  13  46];
end;
if e ==4 
% anger: 
subj1a = [10];
subj1b = [15];
end;
if e ==5 
% sad: 
subj1a = [];
subj1b = [];
end;
if e ==6 
% happy: 
subj1a = [];
subj1b = [];
end;
if e ==7 
% fear: 
subj1a = [6];
subj1b = [23];
end;
if e ==8 
% love: 
subj1a = [];
subj1b = [];
end;
if e ==9 
% jealousy: 
subj1a = [];
subj1b = [];
end;
if e ==10 
% compassion: 
subj1a = [];
subj1b = [];
end;
if e ==11 
% emabarrass: 
subj1a = [10 10];
subj1b = [15  18];
end;
% content: 
subj1a = [5  5  5];
subj1b = [12  15  18];
end;
% grief: 
subj1a = [5];
subj1b = [18];
end;
% relief: 
subj1a = [];
subj1b = [];
end;

%coh cutoff = .2; 9 Hz;  -1 1 sec around button press
if e ==1 
% awe: 
subj1a = [5   5  10  13];
subj1b = [13  46  46  18];
end;
if e ==2 
% frustration: 
subj1a = [5   7  10  10  10  10  11  12  13  13  15  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [13  10  10  13  18  46  10  10  10  18  10  18  46  10  10  10  10  10  10  10  10  10  10];
end;
if e ==3 
% joy: 
subj1a = [5   5   5  11  13  15  15];
subj1b = [10  13  46  36  46  18  46];
end;
if e ==4 
% anger: 
subj1a = [5  10  15  15];
subj1b = [13  15  24  46];
end;
if e ==5 
% sad: 
subj1a = [5];
subj1b = [46];
end;
if e ==6 
% happy: 
subj1a = [5  10  10  15];
subj1b = [46  15  18  24];
end;
if e ==7 
% fear: 
subj1a = [];
subj1b = [];
end;
if e ==8 
% love: 
subj1a = [5  5  5];
subj1b = [13  18  46];
end;
if e ==9 
% jealousy: 
subj1a = [5  10  10  10  10  15  15];
subj1b = [15  13  15  18  24  18  24];
end;
if e ==10 
% compassion: 
subj1a = [];
subj1b = [];
end;
if e ==11 
% emabarrass: 
subj1a = [5   5   7  10  10  10  11  12  12  13  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [12  15  10  10  15  18  10  10  15  10  10  18  10  10  10  10  10  10  10  10  10  10];
end;
if e ==12 
% content: 
subj1a = [5 10];
subj1b = [15  18];
end;
if e ==13 
% grief: 
subj1a = [5];
subj1b = [18];
end;
if e ==14 
% relief: 
subj1a = [5  12];
subj1b = [46  15];
end;


%cutoff = .2; 9.75 Hz;  -1 1 sec around button press
if e ==1 
%awe: 
subj1a = [5   5  10  11  13];
subj1b = [13  46  46  36  18];
end;
if e ==2 
%frustration: 
subj1a = [5   5   7   7  10  10  10  10  10  10  11  12  13  13  15  15  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [13  18  10  18  10  13  15  18  24  46  10  10  10  18  10  18  24  46  10  10  10  10  10  10 10  10  10  10];
end;
if e ==3 
%joy: 
subj1a = [5   5   5   5   5   6   7   7  10  10  10  11  11  12  12  12  13  13  15  15  15  15  16  18  18  19  22  23  24  25  28  36  40];
subj1b = [10  12  13  18  46   7   7  46   7  15  46   7  36   7  15  46   7  46   7  18  24  46   7   7  46   7   7   7   7   7   7   7   7];
end;
if e ==4 
%anger: 
subj1a = [5   5   5   7  10  10  10  10  11  12  13  15  15  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [13  15  46  10  10  15  24  46  10  10  10  10  18  24  46  10  10  10  10  10  10  10  10  10  10];
end;
if e ==5 
%sad: 
subj1a = [5   5  15];
subj1b = [13  46  18];
end;
if e ==6 
%happy: 
subj1a = [5   5   7  10  10  10  10  11  12  13  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [13  46  10  10  15  18  24  10  10  10  10  24  10  10  10  10  10  10  10  10  10  10];
end;
if e ==7 
%fear: 
subj1a = [5  5];
subj1b = [13  46];
end;
if e ==8 
%love: 
subj1a = [5  5 10];
subj1b = [13  46  15];
end;
if e ==9 
%jealousy: 
subj1a = [5   7   7  10  10  10  10  10  10  11  12  13  15  15  15  16  18  19  22  23  24  24  25  28  36  40];
subj1b = [15  10  18  10  13  15  18  24  46  10  10  10  10  18  24  10  10  10  10  10  10  46  10  10  10  10];
end;
if e ==10 
%compassion: 
subj1a = [5 10];
subj1b = [46  18];
end;
if e ==11 
%emabarrass: 
subj1a = [5   5   7  10  10  10  10  11  12  12  13  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [12  18  10  10  12  15  18  10  10  15  10  10  18  10  10  10  10  10  10  10  10  10  10];
end;
if e ==12 
%content: 
subj1a = [5 10];
subj1b = [46  18];
end;
if e ==13 
%grief: 
subj1a = [5];
subj1b = [46];
end;
if e ==14 
%relief: 
subj1a = [5  5];
subj1b = [10  46];
end;

%cutoff = .2;  10.5 Hz;  -1 1 sec around button press
if e ==1 
% awe: 
subj1a = [5   5   5  11];
subj1b = [13  18  46  36];
end;
if e ==2 
% frustration: 
subj1a = [5   5   5   7   7  10  10  10  10  11  11  12  13  13  15  15  15  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [12  13  18  10  18  10  15  18  46  10  36  10  10  18  10  18  24  46  10  10  10  10  10  10  10  10  10  10];
end;
if e ==3 
% joy: 
subj1a = [5   5   5   5   5   6   7   7   7  10  10  10  11  11  12  12  12  13  15  15  15  15  16  18  18  19  22  23  24  25  28  36  40];
subj1b = [10  12  13  18  46   7   7  18  46   7  15  46   7  36   7  15  46   7   7  18  24  46   7   7  46   7   7   7   7   7   7   7   7];
end;
if e ==4 
% anger: 
subj1a = [5   5   5   5   7   7  10  10  10  13  15  15  15];
subj1b = [13  15  18  46  18  46  15  24  46  18  18  24  46];
end;
if e ==5 
% sad: 
subj1a = [5   5  15];
subj1b = [13  46  18];
end;
if e ==6 
% happy: 
subj1a = [5   5   7  10  15  15];
subj1b = [13  46  18  15  24  46];
end;
if e ==7 
% fear: 
subj1a = [5   5   5   5   5  10  15  15  18];
subj1b = [10  12  13  18  46  46  18  46  46];
end;
if e ==8 
% love: 
subj1a = [5 10];
subj1b = [46  15];
end;
if e ==9 
% jealousy: 
subj1a = [5   5   5   7   7  10  10  10  10  10  10  11  12  13  15  15  15  16  18  19  22  23  24  24  25  25  28  28  36  36  40  40];
subj1b = [13  18  46  10  18  10  13  15  18  24  46  10  10  10  10  18  24  10  10  10  10  10  10  25  10  25  10  25  10  25  10  25];
end;
if e ==10 
% compassion: 
subj1a = [5  7];
subj1b = [46  46];
end;
if e ==11 
% emabarrass: 
subj1a = [7   7  10  10  10  11  12  13  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [10  18  10  15  18  10  10  10  10  10  10  10  10  10  10  10  10  10  10];
end;
if e ==12 
% content: 
subj1a = [5   5   7  10  13];
subj1b = [13  46  46  18  18];
end;
if e ==13 
% grief: 
subj1a = [5  7];
subj1b = [46  18];
end;
if e ==14 
% relief: 
subj1a = [5  5  7];
subj1b = [10  46  46];
end;
%cutoff = .2;  11.25 Hz;  -1 1 sec around button press
if e ==1 
% awe: 
subj1a = [5  5];
subj1b = [13  18];
end;
if e ==2 
% frustration: 
subj1a = [5   5   5  10  11  12  12  13  15  15  15];
subj1b = [13  18  46  46  36  15  46  18  18  24  46];
end;
if e ==3 
% joy: 
subj1a = [5   5   5   5   6   7   7   7  10  11  11  12  12  13  15  15  15  16  18  18  19  22  23  24  25  28  36  40];
subj1b = [10  13  18  46   7   7  18  46   7   7  36   7  46   7   7  18  46   7   7  46   7   7   7   7   7   7   7   7];
end;
if e ==4 
% anger: 
subj1a = [5   5   5   5   6   7   7  10  13  15  15  18];
subj1b = [13  15  18  46  25  18  46  46  18  18  46  46];
end;
if e ==5 
% sad: 
subj1a = [15];
subj1b = [18];
end;
if e ==6 
% happy: 
subj1a = [5   5   7   7  13  15  15];
subj1b = [13  46  18  46  18  18  24];
end;
if e ==7 
% fear: 
subj1a = [5   5   5   5   6   7  10  12  15  15  18];
subj1b = [10  13  18  46  10  46  46  46  18  46  46];
end;
if e ==8 
% love: 
subj1a = [5];
subj1b = [46];
end;
if e ==9 
% jealousy: 
subj1a = [5   5   5   7  10  10  15  24  25  28  36  40];
subj1b = [13  18  46  18  18  46  18  25  25  25  25  25];
end;
if e ==10 
% compassion: 
subj1a = [5  7];
subj1b = [46  46];
end;
if e ==11 
% emabarrass: 
subj1a = [5   5   7   7  10  10  10  11  12  13  13  15  16  18  19  22  23  24  25  28  36  40];
subj1b = [13  46  10  18  10  15  18  10  10  10  18  10  10  10  10  10  10  10  10  10  10  10];
end;
if e ==12 
% content: 
subj1a = [5   5   7   7  13];
subj1b = [13  46  18  46  18];
end;
if e ==13 
% grief: 
subj1a = [5];
subj1b = [46];
end;
if e ==14 
% relief: 
subj1a = [5  5  5];
subj1b = [10  13  46];
end;
%cutoff = .2;  12 Hz;  -1 1 sec around button press
if e ==1 
% awe: 
subj1a = [5];
subj1b = [13];
end;
if e ==2 
% frustration: 
subj1a = [5   5   5   7  10  13  15];
subj1b = [13  18  46  46  46  18  18];
end;
if e ==3 
% joy: 
subj1a = [5   5   5   5   7  11  15  15];
subj1b = [10  13  15  18  46  36  18  46];
end;
if e ==4 
% anger: 
subj1a = [5   5   5  10  13  15];
subj1b = [13  18  46  46  18  18];
end;
if e ==5 
% sad: 
subj1a = [5];
subj1b = [15];
end;
if e ==6 
% happy: 
subj1a = [5   5   5   6   7   7  13];
subj1b = [13  18  46  25  18  46  18];
end;
if e ==7 
% fear: 
subj1a = [5   5   5   5   5   6   6   7  10  10  11  12  13  13  13  15  16  18  18  19  22  23  24  25  28  36  40];
subj1b = [6  10  13  18  46   6  10   6   6  46   6   6   6  18  46   6   6   6  46   6   6   6   6   6   6   6   6];
end;
if e ==8 
% love: 
subj1a = [];
subj1b = [];
end;
if e ==9 
% jealousy: 
subj1a = [5   5   5  13  15];
subj1b = [13  18  46  18  18];
end;
if e ==10 
% compassion: 
subj1a = [];
subj1b = [];
end;
if e ==11 
% emabarrass: 
subj1a = [5   5   5  13];
subj1b = [13  18  46  18];
end;
if e ==12 
% content: 
subj1a = [5  7];
subj1b = [13  46];
end;
if e ==13 
% grief: 
subj1a = [5];
subj1b = [13];
end;
if e ==14 
% relief: 
subj1a = [5  5];
subj1b = [13  46];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tl =  {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
figure; params = 'cutoff = .2;   11.25 Hz;    -1 1 sec around button press';
for e = 1:14
% insert comp pairs here:
if e ==1 
% awe: 
subj1a = [5];
subj1b = [13];
end;
if e ==2 
% frustration: 
subj1a = [5   5   5   7  10  13  15];
subj1b = [13  18  46  46  46  18  18];
end;
if e ==3 
% joy: 
subj1a = [5   5   5   5   7  11  15  15];
subj1b = [10  13  15  18  46  36  18  46];
end;
if e ==4 
% anger: 
subj1a = [5   5   5  10  13  15];
subj1b = [13  18  46  46  18  18];
end;
if e ==5 
% sad: 
subj1a = [5];
subj1b = [15];
end;
if e ==6 
% happy: 
subj1a = [5   5   5   6   7   7  13];
subj1b = [13  18  46  25  18  46  18];
end;
if e ==7 
% fear: 
subj1a = [5   5   5   5   5   6   6   7  10  10  11  12  13  13  13  15  16  18  18  19  22  23  24  25  28  36  40];
subj1b = [6  10  13  18  46   6  10   6   6  46   6   6   6  18  46   6   6   6  46   6   6   6   6   6   6   6   6];
end;
if e ==8 
% love: 
subj1a = [];
subj1b = [];
end;
if e ==9 
% jealousy: 
subj1a = [5   5   5  13  15];
subj1b = [13  18  46  18  18];
end;
if e ==10 
% compassion: 
subj1a = [];
subj1b = [];
end;
if e ==11 
% emabarrass: 
subj1a = [5   5   5  13];
subj1b = [13  18  46  18];
end;
if e ==12 
% content: 
subj1a = [5  7];
subj1b = [13  46];
end;
if e ==13 
% grief: 
subj1a = [5];
subj1b = [13];
end;
if e ==14 
% relief: 
subj1a = [5  5];
subj1b = [13  46];
end;
%plot on mri head
    gdcompsa = {subj1a};
    gdcompsb = {subj1b};
    if ~isempty(subj1a)
        % find sources for gdcompsa set
        allsourcesa =  EEG.dipfit.model(gdcompsa{1}(1));
        for ss = 1:length(gdcompsa)
            if ~isempty(gdcompsa{ss})        
                EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
                dipsources= EEG.dipfit.model(gdcompsa{ss}(1));
                for w = 1:length(gdcompsa{ss})
                    dipsources(1,w)= EEG.dipfit.model(gdcompsa{ss}(w));
                end;           
                allsourcesa(end+1:end+size(dipsources,2)) = dipsources;
                dipsources=[];
            end; 
        end;
        allsourcesa(1) = [];
        % find sources for gdcompsb set
        allsourcesb =  EEG.dipfit.model(gdcompsa{1}(1));
        for ss = 1:length(gdcompsb)
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            if ~isempty(gdcompsb{ss})        
                dipsources= EEG.dipfit.model(gdcompsb{ss}(1));
                for w = 1:length(gdcompsb{ss})
                    dipsources(1,w)= EEG.dipfit.model(gdcompsb{ss}(w));
                end;           
                allsourcesb(end+1:end+size(dipsources,2)) = dipsources;
                dipsources=[];
            end; 
        end;
        allsourcesb(1) = [];
        %%%
        allab = allsourcesa;
        allab(end+1:end+size(allsourcesb,2)) = allsourcesb;
        subplot(4,4,e)
        dipplot(allab,'dipolelength',0,'gui','off','dipolesize',25,'view',[0 0 1],'color',{'r'},'image','mri');
        % idx 1:3 non sources, the rest are backward order(last:first of b, then last to first of a; continuous numbering) , last to first, last idx is empty
        pl=1;count=0; 
        htmp = get(gca,'children'); 
        for idx = 1:length(htmp)-1
            ctmp = get(htmp(idx),'userdata');
            if ~isstruct(ctmp)
                count = count+1;
            end;    
        end;
        pl=1;clear  totnum 
        for idx = count+1:length(htmp)-1
            ctmp = get(htmp(idx),'userdata');    
            totnum(idx-count,:) = ctmp.rv(3:4);
        end;
        totnum = unique(totnum,'rows');
        totsource = size(totnum,1);  %***

        pl=1;p=1;findnext = 1;clear newstruct
        for idx = count+1:2:length(htmp)-1
            ctmp = get(htmp(idx),'userdata');
            if  findnext==1
                if ctmp.posxyz(2,1) ~= 0 
                    findnext = 0;
                end;        
                if pl<=totsource/2 
                    ctmp.cnum = pl;
                else
                    ctmp.cnum = pl-(totsource/2);
                end; pl = pl+1;
            else
                if pl<=totsource/2
                    ctmp.cnum = pl-1;
                else
                    ctmp.cnum = pl-1-(totsource/2);
                end;  findnext = 1;     
            end;
            newstruct(p) = ctmp;p=p+1;
        end;


        ff={newstruct.cnum};
        ff=cell2mat(ff); clear allcoo 
        pl=1;bilines=[];
        for cc = 1:totsource/2
            clear coords
            ids = find(ff==cc);
            tt ={newstruct(ids).pos3d};
            tp = tt{1};
            coords(1,:) = tp(1,:);;
            tp = tt{2};
            coords(2,:) = tp(1,:);;
            if size(tt,2)>2
                tp = tt{2};
                tpbi(1,:) =  tp(1,:);
                tp = tt{3};
                tpbi(2,:) =  tp(1,:);
                bilines{pl} =  tpbi;pl=pl+1;
            end;
            allcoo{cc} = coords;
        end;

        for cc = 1:size(allcoo,2)
            tp = allcoo{cc};    
            line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
        end;
        if ~isempty(bilines)
            for cc = 1:size(bilines,2)
                tp = bilines{cc};    
                ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
                set(ph,'linestyle','--');
            end;
        end;
        
        ph=title(tl{e});
        set(ph,'color','r'); 
        set(gca,'fontsize',12);
    else
        subplot(4,4,e)
        set(gca,'xcolor','w');
        set(gca,'ycolor','w');
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        ph=title(tl{e});
        set(ph,'color','r'); 
    end;
end;
ph=textsc(params,'title');
set(ph,'color','r'); 
set(gcf,'color','w');
