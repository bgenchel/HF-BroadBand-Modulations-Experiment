% Top part for new 256 data
% run timef on all relevant epochs:
% find baseline power first:
comment = 'jo74 timef( dat, size(dat,2),[0 size(dat,2)/EEG.srate*1000], 256,[3 .5], winsize,256,padratio,4,timesout,200,alpha,.01,type, coher); baseine is compsXfreqs'; 
tms = EEG.event(1).latency:EEG.event(2).latency;
figure;
for n = 1:length()
    dat = EEG.icaact(rbcomps(n),tms);
    [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( dat, size(dat,2),[0 size(dat,2)/EEG.srate*1000], 256,[3 .5], 'winsize',256,'padratio',4,'timesout',200,'alpha',.01,'type', 'phasecoher');
    baseline(rbcomps(n),:) = powbase;clf
end;
save /data/common/projects/Gimagery/rb74/ersps/baseline.mat baseline freqs times comment


rbcomps = [3,4,5,6,7,9,12,16];
figure; p = 1; 
for k = 1:length(EEG.event)-1
    tms = EEG.event(k).latency:EEG.event(k+1).latency;
    if size(tms,2) > 4096  
        comp_ersp=zeros(63,200,rbcomps(end));
        comp_itc=zeros(63,200,rbcomps(end));
        ersp_boot = zeros(2,63,rbcomps(end));
        itc_boot= zeros(63,rbcomps(end));
        for n = 1:length(rbcomps)
        dat =  EEG.icaact(rbcomps(n),tms);
            [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef(dat, size(dat,2),[0 size(dat,2)/EEG.srate*1000], 256,[3 .5], 'winsize',256,'padratio',4,'plotitc','off','timesout',200,'powbase',baseline(rbcomps(n),:),'alpha',.01,'type', 'phasecoher');
            comp_ersp(:,:,rbcomps(n))=ersp;
            comp_itc(:,:,rbcomps(n))=itc;
            clf
        end;
    calcset(1,p) = k; p = p+1;
    erspcell{k} = comp_ersp;
    itccell{k} = comp_itc;
    end;
end;
comment = 'rb74 ersp and itc cells of comp_ersps with all good comps of a certain time period (single trials) cells are events. Calcset is a set of numbers representing the time stretches that were plotted (longer than 4096 pnts)';;
save /data/common/projects/Gimagery/rb74/ersps/AllEventErspsItcs.mat erspcell itccell freqs times comment calcset;
 
%%%%%%%%%  Plot the Results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
load /data/common/projects/Gimagery/rb74/ersps/baseline.mat
load /data/common/projects/Gimagery/rb74/ersps/AllEventErspsItcs.mat 
rbcomps = [3,4,5,6,7,9,12,16];
row = length(rbcomps);
col = 10;%size(calcset,2);
figure; lim = 30;
fr = find(freqs< 30);pl = 1;
for k = 1:length(rbcomps)    
    subplot(row,col,pl);
    topoplot(EEG.icawinv(:,rbcomps(k)),EEG.chanlocs,'electrodes','off');hold on;
    title(int2str(rbcomps(k))); pl=pl+1;
    subplot(row,col,pl);
    plot(freqs,baseline(rbcomps(k),:)');
    set(gca,'xlim',[3 20]); pl=pl+1;
    for n = 34:39            %10-17,18-25, 26:33,34-39
        subplot(row,col,pl);
        oneersp = erspcell{calcset(n)};
        imagesc(times,freqs(fr),oneersp(fr,:,rbcomps(k)),[-lim lim]);hold on;
        plot([0 0],[50 0],'k-');pl = pl+1;
        set(gca,'xticklabel',[]);
        if n < length(calcset)+1
            title(calcset(n));
        end;
    end; pl=pl+2;
end;
colorbar




%  goes through all emotions and does single trial timef on all of them
subj1 = [2,3,5,7,8,9,10,11,12,14,17,18,20,23,25,27,29,30,37,38,40];    %kl80
subj2 = [7:12,14,17:20,22,24,26];  %ap80
subj3 = [2,3,5:18,20,21,23:25,27,28,29,31,33,36,37,38,41,43,48,49,50];  % rr80
    gdcomps = {subj1 subj2 subj3 };% subj3 subj4};
             
paths = {'/data/common1/emotion/kl80/','/data/common1/emotion/ap80/','/data/common1/emotion/rr80/','/data/common1/emotion/jo74/'};

%  kl80
name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'  'hatred'   'grief'   'relief' }; % for kl80
filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'hatred.set'   'grief.set'   'relief.set' }; % for kl80
datafiles = {'awe.mat' 'frust.mat'  'joy.mat'  'anger.mat'  'sad.mat' 'surprise.mat'   'happy.mat' 'fear.mat' 'love.mat' 'jealousy.mat' 'compassion.mat'  'hatred.mat'   'grief.mat'   'relief.mat'

%  ap80 
name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'   'grief'   'relief' }; % for rr80

%filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'grief.set'   'relief.set' }; % for rr80
%filnm = { 'awepress.set' 'frustpress.set'  'joypress.set'  'angerpress.set'  'sadpress.set' 'surprisepress.set'   'happypress.set' 'fearpress.set' 'lovepress.set' 'jealousypress.set' 'compassionpress.set'  'griefpress.set'   'reliefpress.set' }; % for ap80
filnm = { 'awe2press.set' 'frust2press.set'  'joy2press.set'  'anger2press.set'  'sad2press.set' 'surprise2press.set'   'happy2press.set' 'fear2press.set' 'love2press.set' 'jealousy2press.set' 'compassion2press.set'  'grief2press.set'   'relief2press.set' }; % for rr80
%  rr80 
name = {'awe' 'frust'  'joy'  'anger'  'sad' 'surprise'   'happy' 'fear' 'love' 'jealousy' 'compassion'   'grief'   'relief' }; % for rr80

%filnm = {'awe.set' 'frust.set'  'joy.set'  'anger.set'  'sad.set' 'surprise.set'   'happy.set' 'fear.set' 'love.set' 'jealousy.set' 'compassion.set'  'grief.set'   'relief.set' }; % for rr80
filnm = { 'frustpress.set'  'joypress.set'  'angerpress.set'  'sadpress.set' 'surprisepress.set'   'happypress.set' 'lovepress.set' 'jealousypress.set' 'compassionpress.set'  'griefpress.set'   'reliefpress.set' }; % for rr80



subj = 2;
cd (paths{subj})
load wts
load sph
figure;
load baseline.mat
for p = 7:length(filnm)
    EEG = pop_loadset( filnm{p},paths{subj});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    % for continuous data
    %nonans = isnan(EEG.data(1,:));
    %EEG.data(:,nonans)=[];
    %EEG.pnts = size(EEG.data,2);
   
    % now recalculate ica.act
    EEG.icaweights = wts;
    EEG.icasphere = sph;
    EEG.icaact = []; EEG.icawinv=[];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
% calculate timef 
    comp_ersp=zeros(65,200,50);
    for n=1:length(gdcomps{subj})
        [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(gdcomps{subj}(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 250,[3 .5], 'winsize',256,'padratio',4,'plotitc','off','timesout',200,'powbase',baseline(gdcomps{subj}(n),:),'alpha',.01,'type', 'coher');
        comp_ersp(:,:,gdcomps{subj}(n))=ersp;
        comp_itc( :,:,gdcomps{subj}(n))=itc;
        ersp_boot(:,:,gdcomps{subj}(n)) =  erspboot;
        itc_boot(:,gdcomps{subj}(n)) = itcboot';
    clf    
    end;
    AllpressTF{1,p} = comp_ersp;
    Allpressboot{1,p} = ersp_boot;
    AllpressTF{2,p} = comp_itc;
    Allpressboot{2,p} = itc_boot;
    ALLEEG=[];
end;
    comment = 'ap80 Emotion pilot.  epoched on beginning of FIRST round of emotion press.  Baseline is from whole session.winsize=256,padratio=4;timesout=200; all emotions saved in one cell array in the following order: awe frust  joy  anger  sad surprise   happy fear love jealousy compassion   grief   relief; with first row of cells being ersp and secon row being itc';
    save AllpressTF.mat AllpressTF Allpressboot times freqs comment baseline;


%  Plot each component, all emotions
% load in chanlocs and adjust   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/rr80/rr80.elp', 'filetype',''}, 'forcelocs',{0, 'X', 'Cz'}, 'convert',{ 'chancenter',[],1});
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/ap80/ap80.elp', 'filetype',''}, 'forcelocs',{0, 'X', 'Cz'}, 'convert',{ 'chancenter',[],1});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row = 17%length(gdcomps{subj});
col = length(filnm)+2;
figure;pl=3;tp = 1;
for cp = 18:length(gdcomps{subj})
%load AllpressTF.mat
%load Allpress2TF.mat
%load AllemoTF.mat
fr = find(freqs<30);
    if cp == 18
        figure;pl=3;tp=1;
    end;
    
    for p = 1:length(filnm)
        tmpemo =  AllpressTF{1,p};  % 1 is ersp, 2 is itc
        tmpboot = Allpressboot{1,p};% 1 is ersp, 2 is itc
minboot = tmpboot(1,:,:);
maxboot = tmpboot(2,:,:);
minmask = repmat (minboot,[200,1,1]);
maxmask = repmat (maxboot,[200,1,1]);
        tmpemo(find(tmpemo(:,:,gdcomps{subj}(cp))<minmask(:,:,gdcomps{subj}(cp))'&tmpemo(:,:,gdcomps{subj}(cp))>maxmask(:,:,gdcomps{subj}(cp))'))=0;
        subplot(row,col,pl)        
        imagesc(times,freqs(fr),tmpemo(fr,:,gdcomps{subj}(cp)),[-10 10]);
        pl=pl+1;
        %set(gca,'xticklabel',[]);
    end;
    colorbar
    pl=pl+1;
    subplot(row,col,tp)
    topoplot(EEG.icawinv(:,gdcomps{subj}(cp)),EEG.chanlocs,'electrodes','off');
    title(int2str(gdcomps{subj}(cp)));tp=tp+1;
load baseline.mat
frb = find(freqs<30);
    subplot(row,col,tp)
plot(freqs(frb),baseline(gdcomps{subj}(cp),frb));
hold on
set(gca,'xlim',[2 30]);
set(gca,'yticklabel',[]);
title(int2str(gdcomps{subj}(cp)))
 tp=tp+col-1;   
    pl=pl+1;
end; 
textsc('awe frust joy anger sad surprise happy fear love jealousy compassion grief relief','title');

% avg "good" and "bad" emotions
%awe frust  joy  anger  sad surprise   happy fear love jealousy compassion   grief   relief
gdemo = [1,3,7,9,11,13];
bdemo = [2,4,8,10];
sdemo = [4,12];
row = 8;
col = 8;
figure;pl=1; 
%load /data/common1/emotion/rr80/AllemoTF.mat
fr = find(freqs<30);
for cp = 1:length(gdcomps{subj})
    subplot(row,col,pl)
    topoplot(EEG.icawinv(:,gdcomps{subj}(cp)),EEG.chanlocs,'electrodes','off');
    title(int2str(gdcomps{subj}(cp)));pl = pl+1;
    for p = 1:length(datafiles)
        x= AllemoTF{1,gdemo(1)};
        for y = 2:length(gdemo)
            xx= AllemoTF{1,gdemo(y)};           
            x = x+xx;
        end; x = x/length(gdemo);
    end;
    subplot(row,col,pl)
    imagesc(1:250,freqs(fr),x(fr,300:550,gdcomps{subj}(cp)),[-20 20]);
    pl=pl+1;
    set(gca,'xticklabel',[]); title('Good');
    % avg bad emotions
    for p = 1:length(datafiles)
        x= AllemoTF{1,bdemo(1)};
        for y = 2:length(bdemo)
            xx= AllemoTF{1,bdemo(y)};           
            x = x+xx;
        end; x = x/length(bdemo);
    end;
    subplot(row,col,pl)
    imagesc(1:250,freqs(fr),x(fr,300:550,gdcomps{subj}(cp)),[-20 20]);
    pl=pl+1;
    set(gca,'xticklabel',[]); title('Bad');
    % avg sad emotions
    for p = 1:length(datafiles)
        x= AllemoTF{1,sdemo(1)};
        for y = 2:length(sdemo)
            xx= AllemoTF{1,sdemo(y)};           
            x = x+xx;
        end; x = x/length(sdemo);
    end;
    subplot(row,col,pl)
    imagesc(1:250,freqs(fr),x(fr,300:550,gdcomps{subj}(cp)),[-20 20]);
    pl=pl+1;
    set(gca,'xticklabel',[]); title('Sad/Grief');
    colorbar
end;

% Runtimef on tmp4ica data for baseline  (whole trial, artifacts removed)

% make fake events every sec throughout    (continuous data)for ica
tmpts = size(EEG.data,2);
for index = EEG.srate:EEG.srate*4:tmpts
      EEG.event(end+1) =  EEG.event(end);% appends events to the end
      EEG.event(end).latency = index;
      EEG.event(end).type = 'fake';
end;
EEG = eeg_checkset(EEG, 'eventconsistency');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% Reorder events according to latency
EEG = pop_editeventvals(EEG, 'sort',{ 'latency',0});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


[outeeg,indices] = pop_epoch (EEG, {'fake'}, [-2 2], 'epochinfo','yes');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, outeeg,CURRENTSET);
EEG = pop_rmbase( EEG, [-2000 2000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
figure;   subj = 2;
for n=1:length(gdcomps{subj})
    [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(gdcomps{subj}(n),:), EEG.pnts,[-2000 2000], 250,[3 .5], 'winsize',256,'padratio',4,'plotitc','off','baseline',2000);
    baseline(gdcomps{subj}(n),:) = powbase;
    clf
end;
comment = 'baseline from artifact corrected data from whole trial (all emotions) using fake epochs from -2 2 throughout session, winsize=256,padratio=4,baseline=2000(-2 2)';
save baseline.mat baseline freqs times comment
 
figure;
load baseline.mat
fr = find(freqs<30);
for pl = 1:length(gdcomps{subj})
subplot(6,6,pl)
plot(freqs(fr),baseline(gdcomps{subj}(pl),fr));
hold on
set(gca,'xlim',[2 30]);
title(int2str(gdcomps{subj}(pl)))
end;
