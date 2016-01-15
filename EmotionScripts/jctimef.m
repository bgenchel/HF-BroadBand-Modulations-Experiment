% runs timef on all emos from jo74

subj1 = [];
emoset = {'relief.set','grief.set','content.set','embarrass.set','compassion.set','jealousy.set','love.set','fear.set','happy.set','sad.set','anger.set','joy.set', 'frustration.set' ,'awe.set'};

%comment = 'Pre-session eyes closed segments used as baseline.[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(k),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,padratio,4,timesout,200,baseline,2000,alpha,.01,type, coher'; 
comment = 'All emos in order of presentation using baseline from timef on pre and post session eyes closed epochs.[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(k),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], winsize,256,powbase,baseline(subj1(n),:),padratio,4,alpha,.01,type, coher'; 

figure;
load /data/common1/emotion/jc82/ersps/baseline.mat
for k = 1:length(emoset)
    EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    sph=floatread('sph252-160.sph',[252 252]); 
    wts=floatread('wts252-160.wts',[160 252]); 
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
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
    %[ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(subj1(k),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,[3 .5], 'winsize',256,'padratio',4,'plotitc','off','baseline',2000,'alpha',.01,'type', 'phasecoher');
    %   baseline(subj1(k),:) = powbase;
    %   clf   
end;
save /data/common1/emotion/jc82/ersps/allemoerspitc.mat erspcell itccell ebootcell ibootcell times freqs comment baseline
% Plot ersps with topo and baseline
emos = {'relief','grief','content','embarrass','compassion','jealousy','love','fear','happy','sad','anger','joy', 'frustration' ,'awe'};
figure; k=1;pl=1; lim = 5;
EEG = pop_loadset( emoset{k}, '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
sph=floatread('sph252-160.sph',[252 252]); 
wts=floatread('wts252-160.wts',[160 252]); 
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
tm = find(times>-500 & times< 2500);
fr = find(freqs< 30);
for k = 1:length(subj1)    
    subplot(length(subj1),length(emoset)+2,pl);
    topoplot(EEG.icawinv(:,subj1(k)),EEG.chanlocs,'electrodes','off');hold on;
    title(int2str(subj1(k))); pl=pl+1;
    subplot(length(subj1),length(emoset)+2,pl);
    plot(freqs,baseline(subj1(k),:)');
    set(gca,'xlim',[3 30]); pl=pl+1;
    for n = 1:length(emoset)
        subplot(length(subj1),length(emoset)+2,pl);
        oneersp = erspcell{n};
        imagesc(times(tm),freqs(fr),oneersp(fr,tm,subj1(k)),[-lim lim]);
        plot([0 0],[50 0],'k-');pl = pl+1;
        if n < length(emoset)+1
            title(emos{n});
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runtimef on actual button presses:
emos = {'relief','grief','content','embarrass','compassion','jealousy','love','fear','happy','sad','anger','joy', 'frustration' ,'awe'};
comment = 'Timef on second derivative of Button presses only. saved as freqXtimeXemotion in order of appearance in expt. winsize,256,padratio,4,alpha,.01,type, phasecoher,powbase is 2s';
figure; 
x = ones(1,167);
x = 2*x;
for k = 1:length(emos)
     EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
     %dat = reshape(EEG.data,1,size(EEG.data,2)*size(EEG.data,3));
diffplot = zeros(1,size(EEG.data,2)-2,size(EEG.data,3));
for r = 1:size(EEG.data,3)
    diffplot(1,:,r) = diff(diff(EEG.data(1,:,r)));
end;
diffplot = squeeze(diffplot);
diffplot = [zeros(1,size(diffplot,2));diffplot;zeros(1,size(diffplot,2))];
diffplot = reshape(diffplot,1,size(diffplot,1)*size(diffplot,2));
figure;   [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( diffplot, EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,3, 'winsize',256,'padratio',4,'type', 'phasecoher', 'maxfreq',128 );%,'mtaper',[1 2]);%,'baseline',3000
    comp_ersp(:,:,k)=ersp;
   comp_itc(:,:,k)=itc;
   baseline(k,:) = powbase;
   clf
 ,'powbase',x
end;
save /data/common1/emotion/jc82/ersps/DDiffButtonErspsItcs.mat comp_ersp ersp_boot comp_itc itc_boot times freqs comment baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = ones(1,128);
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
tm = find(times>-300 & times< 1000);
fr = find(freqs<128);
figure;
for k = 1:length(emos)
    subplot(4,4,k)
    imagesc(times(tm),freqs(fr),comp_ersp(fr,tm,k),[-15  10]);hold on;
    plot([0 0],[128 0],'k-');
    title(emos{k});
end;
colorbar

% Plot ITC
tm = find(times>-300 & times< 1000);
figure;
for k = 1:length(emos)
    subplot(4,4,k)
    imagesc(times(tm),freqs,comp_itc(:,tm,k),[-.8 1]);hold on;
    plot([0 0],[128 0],'k-');
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



x = ones(1,167);
x = 2*x;
diffplot = zeros(1,size(EEG.data,2)-2,size(EEG.data,3));
for r = 1:size(EEG.data,3)
    diffplot(1,:,r) = diff(diff(EEG.data(1,:,r)));
end;
diffplot = squeeze(diffplot);
diffplot = [zeros(1,size(diffplot,2));diffplot;zeros(1,size(diffplot,2))];
diffplot = reshape(diffplot,1,size(diffplot,1)*size(diffplot,2));

figure;   [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( diffplot, EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], 256,0, 'winsize',256,'padratio',4,'type', 'phasecoher', 'maxfreq',128 ,'mtaper',[2 2]);%,'baseline',3000
