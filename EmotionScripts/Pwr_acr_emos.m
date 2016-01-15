%% find power in comps across diff emotions

subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % jo74 emot
emos = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set' ,'jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};
%%%%%%%%%%%%%%  Channel Spectra   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(emos)
    EEG = pop_loadset( emos{k},'/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    if k == 1 | k == 5 | k == 9 | k == 13
        figure;pl = 1;
    end;
    subplot(2,2,pl)
    [spectra,freqs,speccomp,contrib,specstd] = spectopo(EEG.data, EEG.pnts, EEG.srate,'winsize',512,'nfft',512,'chanlocs',EEG.chanlocs, 'freq',[6 8.8 10],'freqrange',[2 16],'limits',[4000 NaN NaN NaN NaN NaN],'electrodes','off' );
    title( emos{k});
    pl = pl+1; ALLEEG=[];    
end;
%%%%%%%%%%%%%%  Component Spectra   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % jo74 emot
sph=floatread('sph252-160.sph',[252 252]); 
wts=floatread('wts252-160.wts',[160 252]); 
figure; 
for nx = 1:length(emos)
    EEG = pop_loadset( emos{nx},'/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for cmp = 1:length(subj1)
        %[Pxx,F] = psd(EEG.icaact(subj1(cmp),2560:end-1280),256,EEG.srate,256);% much smoother
        [Pxx,F] = psd(EEG.icaact(subj1(cmp),512:end),256,EEG.srate,256);% much smoother
        cspec(subj1(cmp),:) = log10(Pxx(:,1))';
    end;
    allemopwr{1,nx} = cspec;
    clf;ALLEEG=[];EEG=[];
end;
comment = ['jo74 psd power spectra in good comps (5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46) during each emotion (the different cell arrays) in this order: awe, frustration,joy,anger,sad,happy,fear,love ,jealousy,compassion,emabarrass,content,grief,relief; starts at 10 sec to end of emotin(space bar)';
freqs = F;
save /data/common1/emotion/jo74/allemopwr.mat allemopwr freqs comment
save /data/common1/emotion/jo74/PostEyesClosed.mat eyespec  freqs % after whole session (incl 2back) eyes closed
save /data/common1/emotion/jo74/PostEyesOpen.mat eyesopenspec  freqs % after whole session (incl 2back) eyes closed
save /data/common1/emotion/jo74/EmoEyes.mat emoeyes freqs % after whole session (incl 2back) eyes closed

%%%%%
load /data/common1/emotion/jo74/allemopwr.mat
emonames = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief','EyesClosed'};
emoorder = [1,3,6,8,10,12,14,2,4,7,9,11,5,13,15];
fr=find(freqs>=2&freqs<=30);  % need to load a dataset first to get this
xf = freqs(fr);
color = {[    1.0000         0         0],[    1.0000    0.5000         0],[.75    .75         0],[    0.5000    1.0000    0.5000],[         0    1.0000    1.0000],[         0    0.5000    1.0000],[         0         0    1.0000],[    1.0000         0         0],[    1.0000    0.5000         0],[.75 ...
                    .75         0],[    0.5000    1.0000    0.5000],[         0    1.0000    1.0000],[         0    0.5000    1.0000],[         0         0    1.0000],[.5 .5 .5]};
row = 2;
col = 4;
tl = 'jo74: Emotion Power';
figure; pl = 2; p = 1;clear meanspec norm
for g = 1:2%length(subj1)
if g == 3 |g == 5 |g == 7 |g == 9 |g == 11 |g == 13 |g == 15 |g == 17 |g==19
    figure;pl = 2;p = 1;
end;
    sbplot(row,col,p)
    topoplot(EEG.icawinv(:,subj1(g)),EEG.chanlocs,'electrodes','off');hold on;
    set(gca,'fontsize',14);
    title(int2str(subj1(g)));
    p = p+4;  
    subplot(row,col,pl:pl+2)
    for nx=1:length(emos)        
        emopwr = allemopwr{1,emoorder(nx)};
        %emopwr = allemopwr{1,nx};
        %norm = mean(emopwr(subj1(g),fr));
        %meanspec(pl,:) = postsess(subj1(g),fr)-norm;
        meanspec(pl,:) = emopwr(subj1(g),fr); 
        ph=plot (freqs(fr),meanspec(pl,:));hold on;
        set(ph,'color',color{nx});
        set(ph,'linewidth',1.75);
        if nx > 7
            set(ph,'linestyle','--');
        end;        
    end;
    ph = plot(eyesopenspec(subj1(g),fr));
        set(ph,'color',color{end});
        set(ph,'linewidth',2);
            set(ph,'linestyle',':');
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'xlim',[xf(1) 30]);
    %set(gca,'ylim',[-.4 .8]);
    set(gca,'yticklabel',{[] [] 0 [] [] [] .8});
    set(gca,'fontsize',14);
    title(int2str(subj1(g)))
    legend(emonames(emoorder));
    pl = pl+4;
end;
textsc(tl,'title');axcopy;

load /data/common1/emotion/jo74/EmoEyes.mat
load /data/common1/emotion/jo74/PostEyesClosed.mat
load /data/common1/emotion/jo74/PostEyesOpen.mat
figure; row = 5; col = 4;pl = 2;p=1;
for g = 1:length(subj1)
    if g == 11
        figure;pl = 2;p=1;
    end;    
    subplot(row,col,p)
    topoplot(EEG.icawinv(:,subj1(g)),EEG.chanlocs,'electrodes','off');hold on;
    set(gca,'fontsize',14);
    title(int2str(subj1(g)));
    p = p+2;  
    subplot(row,col,pl)
    ph = plot(emoeyes(subj1(g),fr),'r');hold on;
    set(ph,'linewidth',2);
    ph = plot(eyesopenspec(subj1(g),fr),'b');hold on;
    set(ph,'linewidth',2);
    ph = plot(eyespec(subj1(g),fr),'k');hold on;
    set(ph,'linewidth',2);
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'xlim',[xf(1) 30]);
    %set(gca,'ylim',[-.4 .8]);
    set(gca,'yticklabel',{[] [] 0 [] [] [] .8});
    set(gca,'fontsize',14);
    title(int2str(subj1(g)))
    pl = pl+2;
end;
legend('EmoEyes','Post-Open','Post-Close');

%%%%%%%%%%%%%%%%%%%%   OLD 71 Channel Analysis:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Emotion expt
emos = {'awepress.set' 'awe2press.set' 'frustpress.set' 'frust2press.set' 'joypress.set' 'joy2press.set' 'angerpress.set' 'anger2press.set' 'sadpress.set' 'sad2press.set' 'surprisepress.set' 'surprise2press.set' 'happypress.set' 'happy2press.set' 'fearpress.set' 'fear2press.set' 'lovepress.set' 'love2press.set' 'jealousypress.set' 'jealousy2press.set' 'compassionpress.set' 'compassion2press.set' 'contentpress.set' 'content2press.set' 'griefpress.set' 'grief2press.set' 'reliefpress.set' 'relief2press.set'};
cd /data/common1/emotion/ap80/
load wts
load sph
load aplocs
a = [7:12,14,17:20,22,24,26];
figure;
for e = 1:length(emos)
    EEG = pop_loadset( emos{e},'/data/common1/emotion/ap80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    EEG.chanlocs = chan;    
    % now recalculate ica.act
    EEG.icaweights = wts;
    EEG.icasphere = sph;
    EEG.icaact = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
    [specmem freqs compspec c] = pop_spectopo(EEG, 0, [-500  3000], 'EEG' , 'plot', 'off', 'plotchan',0,'percent', 100, 'icacomps', a, 'freqrange',[2 30],'freq',[10]);
    clf;ALLEEG=[];
    emocompspc (:,:,e) = compspec;
end;
comment = 'pop_spectopo on each emotion press set (1 and 2 separately) in the order of the first data file events (awe,frust... etc) 1,2 alternating.a and freqs included in .mat spec from -500 to 3000ms around each press';
% save comppwrspec.mat emocompspc a freqs comment
% got rid of last epoch of all sets because they are nans
for e = 4:2:length(emos)
    EEG = pop_loadset( emos{e},'/data/common1/emotion/ap80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    if isnan(EEG.data(1,end,end))==1
        EEG = pop_select( EEG, 'notrial',size(EEG.data,3));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
        if isnan(EEG.data(1,end,end))==1
            EEG = pop_select( EEG, 'notrial',size(EEG.data,3));
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
        end;
    EEG = pop_saveset( EEG,emos{e} , '/data/common1/emotion/ap80/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end; 
    ALLEEG=[];           
end;

% Plot Power across emotions
dats = {'awe' 'awe2' 'frust' 'frust2' 'joy' 'joy2' 'anger' 'anger2' 'sad' 'sad2' 'surprise' 'surprise2' 'happy' 'happy2' 'fear' 'fear2' 'love' 'love2' 'jealousy' 'jealousy2' 'compassion' 'compassion2' 'content' 'content2' 'grief' 'grief2' 'relief' 'relief2'};
nemos = length(dats);
for q = 1:length(a)
    figure;
    for p = 1:nemos
        subplot(round(sqrt(nemos))+1,round(sqrt(nemos)),p);
        plot (freqs(1:60,1),emocompspc(a(q),1:60,p),'b-');
        hold on; 
        set(gca,'xlim',[1 30]);
        set(gca,'ylim',[25 75]);
        title(dats{p});
    end;
    textsc(int2str(a(q)),'title');
end;
%  plot power from all emos on one plot
dats = {'awe' 'awe2' 'frust' 'frust2' 'joy' 'joy2' 'anger' 'anger2' 'sad' 'sad2' 'surprise' 'surprise2' 'happy' 'happy2' 'fear' 'fear2' 'love' 'love2' 'jealousy' 'jealousy2' 'compassion' 'compassion2' 'content' 'content2' 'grief' 'grief2' 'relief' 'relief2'};
nemos = length(dats);
col = {'b' 'r' 'g' 'k' 'm' 'c'}; %'b' 'r' 'g' 'k' 'm' 'c' 'b' 'r' 'g' 'k' 'm' 'c' 'b' 'r' 'g' 'k' 'm' 'c' 'b' 'r' 'g' 'k'};
for q = 1:length(a)
    figure; b=1;
    for p = 2:2:12%nemos
        plot (freqs(1:60,1),emocompspc(a(q),1:60,p),col{b}); b=b+1;
        hold on
        set(gco,'linewidth',3);
    end;
    set(gca,'xlim',[1 30]);
    set(gca,'ylim',[25 75]);
    textsc(int2str(a(q)),'title');
end;

% resample all 2 datasets
dats = { 'frust2.set'  'joy2.set' 'anger2.set'  'sad2.set' 'surprise2.set' 'happy2.set'  'fear2.set' 'love2.set' 'jealousy2.set' 'compassion2.set' 'content2.set' 'grief2.set' 'relief2.set'};
cd /data/common1/emotion/ap80/
for fix = 1:length(dats)
    EEG = pop_loadset( dats{fix},'/data/common1/emotion/ap80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG.srate  =500;
    EEG = pop_resample( EEG, 250);
    EEG = pop_saveset( EEG,dats{fix} , '/data/common1/emotion/ap80/');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG=[];           
end;

