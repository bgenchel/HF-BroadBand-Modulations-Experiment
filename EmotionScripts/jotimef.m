% runs timef on all emos from jo74
eeglab
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
EEG = pop_loadset( 'emo-3-252.set', '/data/common/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
for evtm = 384:768:size(EEG.data,2)-384  % go up by 1.5 sec to create 3 sec epochs til 2 sec before button
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
    EEG.event(end).latency = evtm;
    EEG.event(end).type = 'fake';        
end;
EEG = pop_epoch( EEG,{'fake'} , [-1.5 1.5]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
EEG = pop_rmbase( EEG,[-1500 1500]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_select( EEG, 'trial',[1:2:220] );  % 173 trials,176,177
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_rejkurt(EEG,0,subj1 ,5,4,0,1);
EEG = pop_jointprob(EEG,0,subj1 ,5,4,0,1);
%save as 3 tmp.sets
EEG = pop_saveset( EEG, 'tmp3.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
ALLEEG=[];EEG=[];clear sph wts

EEG = pop_loadset( 'tmp1.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%merge and save as 'BaselineAllEmos.set'

%now run timef to get a whole trial baseline:
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
comment = 'Whole-session baseline measurement. Random epochs (-1.5 1.5s) were selected from whole session for this dataset of 526 epochs. [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,padratio,4,timesout,200,baseline,1500,alpha,.01,type, coher; '; 
    comp_ersp=zeros(63,200,subj1(end));
    comp_itc=zeros(63,200,subj1(end));
    ersp_boot = zeros(2,63,subj1(end));
    itc_boot= zeros(63,subj1(end));
    baseline = zeros(subj1(end),63);
figure;
    for n = 1:length(subj1)
        [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], 'winsize',256,'padratio',4,'baseline',1500,'alpha',.01,'type', 'phasecoher');
        comp_ersp(:,:,subj1(n))=ersp;
        ersp_boot(:,:,subj1(n)) =  erspboot;
        baseline(subj1(n),:) = powbase;
        clf
    end;
save /data/common1/emotion/jo74/ersps/WholeTrialBaseline.mat comp_ersp  ersp_boot  times freqs comment baseline


        eeglab
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
emoset = {'awePress.set', 'frustrationPress.set','joyPress.set','angerPress.set','sadPress.set','happyPress.set','fearPress.set','lovePress.set','jealousyPress.set','compassionPress.set','emabarrassPress.set','contentPress.set','griefPress.set','reliefPress.set'};

%comment = 'Pre-session eyes closed segments to be used for baseline.[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(k),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,padratio,4,timesout,200,baseline,2000,alpha,.01,type, coher; '; 
%comment = 'eyes closed segments of eyes taken after emotion AND twoback on the same day. to be used for baseline comparison to emotion data because it is a longer period of time, maybe more stable.cut up into -1.5 1.5 sec epochs and cleaned.[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(k),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,padratio,4,timesout,200,baseline,1500,alpha,.01,type, coher; '; 
comment = 'All emos in order of presentation using baseline from timef on random epochs from whole session(WholeTrialBaseline.mat)  eyes closed epochs.[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,powbase,baseline(subj1(n),:),padratio,4,alpha,.01,type, coher'; 

sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
load /data/common1/emotion/jo74/ersps/WholeTrialBaseline.mat
%load /data/common1/emotion/jo74/ersps/pre-baseline.mat
%load /data/common1/emotion/jo74/ersps/EpochedErspLateEyes.mat
figure;
for k = 1:length(emoset)
    EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_rejkurt(EEG,0,subj1 ,4,4,0,1);        
        EEG = pop_jointprob(EEG,0,subj1 ,4,4,0,1);
    comp_ersp=zeros(63,200,subj1(end));
    comp_itc=zeros(63,200,subj1(end));
    ersp_boot = zeros(2,63,subj1(end));
    itc_boot= zeros(63,subj1(end));
    
    for n = 1:length(subj1)
        [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], 'winsize',256,'padratio',4,'powbase',baseline(subj1(n),:),'alpha',.01,'type', 'phasecoher');
        comp_ersp(:,:,subj1(n))=ersp;
        comp_itc(:,:,subj1(n))=itc;
        ersp_boot(:,:,subj1(n)) =  erspboot;
        itc_boot(:,subj1(n)) = itcboot';
        baseline(subj1(n),:) = powbase;
        clf
    end;
    erspcell{k} = comp_ersp;
    itccell{k} = comp_itc;
    ebootcell{k} = ersp_boot;
    ibootcell{k} = itc_boot;
    ALLEEG=[]; EEG=[];
    %[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], 'winsize',256,'padratio',4,'plotitc','off','baseline',1500,'alpha',.01,'type', 'phasecoher');
    %   baseline(subj1(n),:) = powbase;
    %   clf   
end;
save /data/common1/emotion/jo74/ersps/ERSPsWholeTrialBsln.mat erspcell itccell ebootcell ibootcell times freqs comment baseline
save /data/common1/emotion/jo74/ersps/EpochedErspLateEyes.mat erspcell itccell ebootcell ibootcell times freqs comment baseline
save /data/common1/emotion/jo74/ersps/EpochedErspPreBase.mat erspcell itccell ebootcell ibootcell times freqs comment baseline
save /data/common1/emotion/jo74/ersps/laterEyes-baseline.mat baseline comment freqs 
% Plot ersps with topo and baseline
load /data/common1/emotion/jo74/ersps/EpochedErspLateEyes.mat
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
EEG = pop_loadset( 'sources.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
tm = find(times>-500 & times< 2000);
fr = find(freqs< 20);
%subj1 = [6,10,12,13,15,18,19,22,23,24,25];
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
emoset = [1:14];
gdcomps = {subj1};
figure; k=1;pl=1; lim = 5;
for k = 1:length(subj1)    
    subplot(length(subj1),length(emoset)+2,pl);
    topoplot(EEG.icawinv(:,subj1(k)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
    title(int2str(subj1(k))); pl=pl+1;
    subplot(length(subj1),length(emoset)+2,pl);
    plot(freqs,baseline(subj1(k),:)');
    set(gca,'xlim',[3 30]); pl=pl+1;
    set(gca,'yticklabel',[]);
    set(gca,'xgrid','on');
    for n = 1:length(emoset)
        subplot(length(subj1),length(emoset)+2,pl);
        oneersp = erspcell{emoset(n)};
        plotersp = oneersp(:,:,subj1(k));
        oneboot = ebootcell{n};
        minmask = oneboot(1,:,:);
        maxmask = oneboot(2,:,:);
        minmask = repmat(minmask,[200 1 1]);
        maxmask = repmat(maxmask,[200 1 1]);
        onemin = minmask(:,:,subj1(k))';
        onemax = maxmask(:,:,subj1(k))';
        %plotersp(find(plotersp>onemin & plotersp < onemax))=0;
        imagesc(times(tm),freqs(fr),plotersp,[-lim lim]); hold on;
        plot([0 0],[50 0],'k-');pl = pl+1;
        set(gca,'ydir','norm');
        if k == 1
            ph =  title(emos{emoset(n)});
            set(ph,'rotation',45);
            set(ph,'fontsize',14);
        end;
        if k == length(subj1)
            set(gca,'xticklabel',{[0 1 2]});
        else
            set(gca,'xticklabel',[]);
        end;        
    end;
end;

    
save /data/common1/emotion/jo74/ersps/baseline.mat baseline freqs times comment

load /data/common1/emotion/jo74/ersps/pre-baseline.mat
figure;
for k = 1:length(subj1)
subplot(5,4,k)
ln = plot(freqs,baseline(subj1(k),:)','b');hold on;
set(ln,'linewidth',2);
end;
load /data/common1/emotion/jo74/ersps/post-baseline.mat
for k = 1:length(subj1)
subplot(5,4,k)
ln = plot(freqs,baseline(subj1(k),:)','r');hold on;
set(ln,'linewidth',2);
end;
load /data/common1/emotion/jo74/ersps/baseline.mat
for k = 1:length(subj1)
subplot(5,4,k)
ln = plot(freqs,baseline(subj1(k),:)','g');hold on;
set(ln,'linewidth',2);
set(gca,'xlim',[3 40]);
title(int2str(subj1(k)));
end;
%%%  plot pwr with collapsed time in a range   
colors = { [0         0    0.75], [ 0 0    1],     [0    0.25    1],     [0    0.5   1],     [0    0.75    1],     [0    1    1],    [0.25    1    0.75],    [0.5    1    0.5],    [0.75 1 0.25],    [1    1 0],[1 0.75 0],    [1    0.5 0],    [1 0.25 0],    [1 0 0]};
figure; p=1;
for comp = 1:length(subj1) 
    if comp==11
        figure;p=1;
    end;    
    subplot(5,4,p); 
    topoplot(EEG.icawinv(:,subj1(comp)),EEG.chanlocs, 'electrodes', 'off', 'shrink', 'skirt');
    title(int2str(subj1(comp)));
    subplot(5,4,p+1);
    for k = 1:length(emoset)   
        tmpersp = erspcell{k};
        oneersp = tmpersp(:,:,subj1(comp));
        oneersp = mean(oneersp,2);
        ph= plot(freqs,oneersp);hold on;
        set(ph,'color',colors{k});
    end;
    title(subj1(comp));
    p=p+2;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Run timef on whole trial of each emotion
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
emoset = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};

comment = 'jo74 All emos in order of presentation using baseline from timef on pre and post session eyes closed epochs.Single trial timef on each emotion period. [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,powbase,baseline(subj1(n),:),padratio,4,alpha,.01,type, coher,timesout,400'; 

figure;
load /data/common1/emotion/jo74/ersps/baseline.mat
for k = 1:length(emoset)
    EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    sph=floatread('sph252-160.sph',[252 252]); 
    wts=floatread('wts252-160.wts',[160 252]); 
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    comp_ersp=zeros(63,600,subj1(end));
    comp_itc=zeros(63,600,subj1(end));
    ersp_boot = zeros(2,63,subj1(end));
    itc_boot= zeros(63,subj1(end));
    
    for n = 1:length(subj1)
        fprintf('\nWorking on Comp %i\n',subj1(n));
        fprintf('\n...of dataset: %s\n',emoset{k});
     [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], 'winsize',256,'padratio',4,'type', 'phasecoher','timesout',600 ,'powbase',baseline(subj1(n),:));
     comp_ersp(:,:,subj1(n))=ersp;
         clf
    end;
    erspcell{k} = comp_ersp;
    ALLEEG=[]; EEG=[];
end;
save /home/julie/EmoTF/allemoerspitc.mat erspcell times freqs comment baseline
% Now image these timefs
load /home/julie/EmoTF/allemoerspitc.mat
fr = find(freqs<30);
% image as comps across emos
subj1 = [5,6,7,10,11,12,13,15,18,19,22,23,24,25];
figure;pl = 1;
for cmp = 1:length(subj1)
    for k = 1:length(emos)
        subplot(length(subj1),length(emos),pl)
        oneemo = erspcell{k};
        imagesc(times,freqs(fr),oneemo(fr,:,subj1(cmp)),[-10 10]); pl = pl+1;
    end;
end;


% image as single components across 14 emotions
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for comp = 1:length(subj1)   
    figure;   lim=20; smfac = 7;
     subplot(4,4,1)
    pop_topoplot(EEG,0, subj1(comp), '',[1 1] ,0, 'electrodes', 'off', 'shrink', 'skirt', 'masksurf', 'on');
    for k = 2:length(emoset)+1    
        subplot(4,4,k)
        tmpersp = erspcell{k-1};
        tmpersp(find(tmpersp(:,:,subj1(comp))>minmask' &tmpersp(:,:,subj1(comp))<maxmask'));
        oneersp = tmpersp(fr,:,subj1(comp));   clear smersp     
        for f = 1:size(oneersp,1)
        pl=1;
        for sm = 1:round(size(oneersp,2)/(smfac+1))
                smersp(f,sm) = mean(oneersp(f,pl:pl+smfac));pl = pl+smfac;% smooth by 'smfac' time frames
            end;
        end;        
        %imagesc(times,freqs(fr), oneersp,[-lim lim]);
        imagesc([1:size(smersp,2)],freqs(fr), smersp,[-lim lim]);
        title(emos{k-1})
    end;
    colorbar;
end;
% image as power in a freq range with all emotions on one plot
fr = find(freqs>5 & freqs<7);
tm = find(times>8000);
colors = { [0         0    0.75], [ 0 0    1],     [0    0.25    1],     [0    0.5   1],     [0    0.75    1],     [0    1    1],    [0.25    1    0.75],    [0.5    1    0.5],    [0.75 1 0.25],    [1    1 0],[1 0.75 0],    [1    0.5 0],    [1 0.25 0],    [1 0 0]};
figure; p=1;
for comp = 1:length(subj1)   
    subplot(7,4,p); pl=1;
    pop_topoplot(EEG,0, subj1(comp), '',[1 1] ,0, 'electrodes', 'off', 'shrink', 'skirt', 'masksurf', 'on');
    for k = 1:length(emoset)   
        tmpersp = erspcell{k};
        oneersp = tmpersp(fr,tm,subj1(comp));
        oneersp = mean(oneersp,1);
        oneersp = mean(oneersp,2);
        allpow(1,pl) = oneersp; pl=pl+1;
    end;
    subplot(7,4,p+1);
    bar(allpow);
    title(comp);
    colorbar;p=p+2;
end;
%%%  plot pwr with collapsed time in a range   
tm = find(times>10000);
figure; p=1;
for comp = 1:length(subj1) 
    if comp==11
        figure;p=1;
    end;    
    subplot(5,4,p); 
    topoplot(EEG.icawinv(:,subj1(comp)),EEG.chanlocs, 'electrodes', 'off', 'shrink', 'skirt');
    title(int2str(subj1(comp)));
    subplot(5,4,p+1);
    for k = 1:length(emoset)   
        tmpersp = erspcell{k};
        oneersp = tmpersp(:,tm,subj1(comp));
        oneersp = mean(oneersp,2);
        ph= plot(freqs,oneersp);hold on;
        set(ph,'color',colors{k});
    end;
    title(subj1(comp));
    p=p+2;
end;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runtimef on actual button presses:
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
comment = 'Timef on Button presses only. saved as freqxtimexemotion in order of appearance in expt. winsize,256,padratio,4,alpha,.01,type, phasecoher,baseline,2000';
figure; 
for k = 1:length(emos)
     EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
   [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.data(1,:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,3, 'winsize',256,'padratio',4,'alpha',.01,'type', 'phasecoher','baseline',3000, 'maxfreq',128);
    comp_ersp(:,:,k)=ersp;
   comp_itc(:,:,k)=itc;
   ersp_boot(:,:,k) =  erspboot;
   itc_boot(:,k) = itcboot';
   baseline(k,:) = powbase;
   clf
  
end;
save /data/common1/emotion/jo74/ersps/ButtonErspsItcs.mat comp_ersp ersp_boot comp_itc itc_boot times freqs comment baselineones(1,167)
x = ones(1,167);
x = 2*x;
diffplot = zeros(1,size(EEG.data,2)-2,size(EEG.data,3));
for r = 1:size(EEG.data,3)
    diffplot(1,:,r) = diff(diff(EEG.data(1,:,r)));
end;
diffplot = squeeze(diffplot);
diffplot = [zeros(1,size(diffplot,2));diffplot;zeros(1,size(diffplot,2))];
diffplot = reshape(diffplot,1,size(diffplot,1)*size(diffplot,2));
figure;   [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef(diffplot , EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,3, 'winsize',256,'padratio',4,'type', 'phasecoher','powbase',x, 'maxfreq',128);

% Plot ERSP
tm = find(times>-100 & times< 1000);
figure;
for k = 1:length(emos)
    subplot(4,4,k)
    imagesc(times(tm),freqs,comp_ersp(:,tm,k),[-25 15]);hold on;
    plot([0 0],[50 0],'k-');
    title(emos{k});
end;
colorbar

% Plot ITC
tm = find(times>-300 & times< 1000);
figure;
for k = 1:length(emos)
    subplot(4,4,k)
    minboot = itc_boot(:,k);
    minmask = repmat(minboot,[1 200]);
    plotitc = comp_itc(:,:,k);
    plotitc(find(plotitc<minmask))=0;
    imagesc(times(tm),freqs,plotitc(:,tm),[0 1]);hold on;
    plot([0 0],[50 0],'k-');
    title(emos{k});
end;
colorbar

%Plot actual trajectories
tm = find(EEG.times>-300 & EEG.times< 1000);
figure;
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    subplot(4,4,k)
    for n = 1:size(EEG.data,3)
        plot(EEG.times(tm),EEG.data(1,tm,n));hold on;
    end;
    plot([0 0],[-2000000000 20000000000],'k-');
    title(emos{k});
    set(gca,'xlim',[EEG.times(tm(1)) EEG.times(tm(end))]);
    set(gca,'ylim',[-100000 max(max(EEG.data(1,tm,:)))+1000000]);
end;
colorbar

%%%%%%%%%%%%%%%%% timef on 2nd derivative of button  %%%%%%%%%%%%%%%%%%%%%
% run part of EmotionPresses256.m to epoch ButtonOnly.set into individual emos
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
%bsln = ones(145,1);
%bsln = 2*bsln;
figure;
for k = 1:length(emos)
    
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    
    diffplot = zeros(1,size(EEG.data,2)-2,size(EEG.data,3));
    for r = 1:size(EEG.data,3)
        diffplot(1,:,r) = diff(diff(EEG.data(1,:,r)));
    end;
    diffplot = squeeze(diffplot);
    diffplot = [zeros(1,size(diffplot,2));diffplot;zeros(1,size(diffplot,2))];
    diffplot = reshape(diffplot,1,size(diffplot,1)*size(diffplot,2));
    
  [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( diffplot, EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,3.5,'type', 'phasecoher', 'maxfreq',128 ,'padratio',2,'winsize',512);%,'powbase',bsln
    comp_ersp(:,:,k)=ersp;
    comp_itc(:,:,k)=itc;
    baseline(k,:) = powbase;
    clf
end;
comment = 'timef on 2nd derivative of button press using wavelets: 3.5,type, phasecoher, maxfreq,128 ,padratio,2,winsize,512; Data from ButtonOnly.set, epoched into emotions and button presses. From jo74-256 channel data; baseline is automatically calculated from -2000 to 0; epochs from -2 4 sec;Apr 15,2004;';
save /data/common/emotion/jo74/ersps/AccErspsItcs.mat freqs times comp_ersp comp_itc baseline
    
exit
exit

% plot the results
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
fr = find(freqs < 80);
figure;
for k = 1:length(emos)   
subplot(4,4,k)
imagesc(times,freqs(fr),comp_ersp(fr,:,k),[-40 40]);hold on;
plot([0 0],[0 128]);
set(gca,'ytick',[10:10:80]);
set(gca,'xtick',[-1000:200:1000]);
set(gca,'yticklabel',{10 [] 30 [] 50 [] 70 []});
set(gca,'xticklabel',{[] -.8 [] -.4 [] 0 [] .4 [] .8 []});
set(gca,'fontsize',14);
title(emos{k});
end;;
colorbar;

    ,'mtaper',[2 2]
    
    whos times freqs
    times(1)
    times(end)
freqs(1) 
freqs(end) 

FFT on frustration 2nd derivative:
times: 200   -1750   1.7461e+03
freqs: 128   1 128
FFT with mtaper [1 1]:
times:200 -1500  1.4961e+03
freqs:128    1 128
FFT with mtaper [2 2]:
times:200  -1000 996.0938
freqs:256    .5 128
FFT with 2's baseline
wavelet 3:
times:200   -1750  1.7461e+03
freqs:41    6 126
  padratio:4 
  times:200   -1750  1.7461e+03
  freqs:82    6 127.5
  padratio:4  winsize:256
  times:200   -1500 1.4961e+03
  freqs:167   3 127.5
  padratio:4  winsize:512
  times:200   -1000 996.0938
  freqs:338   1.5 127.875
mtaper doesn't show a major shift in pattern from others, but power varies wildly (.7 in [1 1] to 4 in [2 2])
wavelet 3 is probably best for itc. ersp does not show the low freq shift, though. maybe freqs don't go low enogh. With larger winsize, starting to get hints of it. I think I believe the wavelets the most. But who knows what's really the best? if any...
            
            
           
