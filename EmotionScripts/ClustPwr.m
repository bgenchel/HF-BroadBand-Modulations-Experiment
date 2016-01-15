% clusters power from timef(same baseline) across different components using ICA
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
%subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];

load /data/common1/emotion/jo74/ersps/EpochedErspPreBase.mat
load /data/common1/emotion/jo74/ersps/allemoerspitc.mat
load /data/common1/emotion/jo74/ersps/ERSPsWholeTrialBsln.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = find(freqs<30);
% makes a comps X emo*t/f
% this is fine for more than one subject
allpwr  =zeros(length(subj1),length(fr)*length(emos));
for comp = 1:length(subj1)
    tmppwr  =zeros(1,0);
    for n = 1:length(emos)
        tmpersp = erspcell{n};
        oneersp = tmpersp(fr,:,subj1(comp));
        oneersp = mean(oneersp,2);
        mn = mean(oneersp,1);
        oneersp = oneersp - mn;
        tmppwr(1,end+1:end+size(oneersp,1)) = oneersp';  % 
    end;  
    allpwr(comp,:) = tmppwr;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform for comps X emos*t/f
icadata = zeros (length(emos),length(fr),size(weights,1));
for n = 1:size(weights,1)
    tmpcomp = activations (n,:);
    tmpcomp = reshape(tmpcomp,length(fr),length(emos));  %Good
    tmpcomp = tmpcomp';
    icadata(:,:,n) = tmpcomp;
end; 
%%%%  Image the results   %%%%%%%
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};

% image as color
figure;pl=1;
for comp = 1:size(weights,1)
    subplot(4,4,pl)
    imagesc(freqs(fr),[1:size(icadata,1)],icadata(:,:,comp),[-6 6]);hold on;
    title(int2str(comp)); 
    pl=pl+1;
end;
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Or use emos as input to cluster between emos
fr = find(freqs<50);
tm = find(times>-400 & times< 1500);
allpwr  =zeros(length(emos),length(fr)*length(tm)*length(subj1));
for n = 1:length(emos)
    tmppwr  =zeros(1,0);
    for comp = 1:length(subj1)
        tmpersp = erspcell{n};
        oneersp = tmpersp(fr,tm,subj1(comp));
        oneersp = reshape(oneersp,1,size(oneersp,1)*size(oneersp,2));
        tmppwr(1,end+1:end+size(oneersp,2)) = oneersp;  % 
    end;  
    allpwr(n,:) = tmppwr;
end;

% Run ICA

[weights,sphere,compvars,bias,signs,lrates,activations] = runica (allpwr);

%weights = icacomps(clusters) X # emotions
% winv = # emotions X # emo clusters
winv = weights*sphere;

%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% or for emos X comps*freqs*time
for k = 1:length(subj1)
    icadata = zeros (length(fr),length(tm),length(emos));
    for n = 1:length(emos)
        g = k-1;
        tmpcomp = activations (n,(length(fr)*length(tm))*g+1:(length(fr)*length(tm))*k);
        tmpcomp = reshape(tmpcomp,length(fr),length(tm));  %Good
        icadata(:,:,n) = tmpcomp;
    end; 
    alldat{k} = icadata;
end;

% image as color
EEG = pop_loadset( 'sources.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
for comp = 1:length(subj1)
    icadata = alldat{comp};
    figure;pl=1;
        subplot(4,4,pl)
        topoplot(EEG.icawinv(:,subj1(comp)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
        title(int2str(subj1(comp))); pl=pl+1;    
    for n = 1:length(emos) 
        subplot(4,4,pl)
        imagesc(times(tm),freqs(fr),icadata(:,:,n),[-6 6]);hold on;
        title(int2str(n)); 
        set(gca,'ydir','norm');
        pl=pl+1;
    end;   
    colorbar
    textsc(int2str(subj1(comp)),'title');
end;
% scatter plot:
figure; 
for n = 1:length(emos)
subplot(4,4,n)
plot(winv(:,n),'k.');
for k = 1:length(emos)
text(k,winv(k,n),emos{k});
end;
title(int2str(n));
end;

% visualize a single cluster (all comps)
for n = 1:length(emos) 
    figure;pl=1;
    for comp = 1:length(subj1)
    icadata = alldat{comp};
        subplot(4,4,pl)
        imagesc(times(tm),freqs(fr),icadata(:,:,n),[-6 6]);hold on;
        title(int2str(subj1(comp))); 
        pl=pl+1;
    end;   
    colorbar
    textsc(['Cluster: ',int2str(n)],'title');
end;

% plot as single clusters with scatter and scalp maps
EEG = pop_loadset( 'sources.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
for clust = 1:14
    figure; pl = 3; fac = 1; % -1 for inversion neg = [1,13
    subplot(6,6,1:2)
    plot(winv(:,clust)*fac,'k.'); hold on;
    for k = 1:length(emos)
        text(k,winv(k,clust)*fac,emos{k});
    end;
    set(gca,'fontsize',16);
    title(['Cluster: ',int2str(clust)]);
    plot([0 20],[0 0],'r');

    for comp = 1:length(subj1)
        subplot(6,6,pl)
        topoplot(EEG.icawinv(:,subj1(comp)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
        title(int2str(subj1(comp))); pl=pl+1;
        icadata = alldat{comp};
        subplot(6,6,pl)
        imagesc(times(tm),freqs(fr),icadata(:,:,clust)*fac,[-6 6]);hold on;
        title(int2str(subj1(comp))); 
        set(gca,'ydir','norm');
        set(gca,'ticklength',[.02 .02]);
        pl=pl+1;
    end;
end;

%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% run PCA on each activation cluster
for n = 1:size(weights,1)
onepca = activations(n,:);
onepca = reshape(onepca,length(emos),length(fr));
[pc,eigvec,sv] = runpca(onepca');% 
%%% pc = #comps(pc clusters) X #comps(input comps): gives you eigenvectors (inv(weights))
%%  rows are pc comps, columns are input components
%%% eigvec= totalength X #comps(pc clusters): gives you the ersp comp maps once re-transformed
pcs{n} = pc;
eigs{n} = eigvec;
end;
figure;
for n=1:length(pcs)
    subplot(4,4,n)
    ploteig = eigs{n};
    imagesc(freqs(fr),[1:14],ploteig',[-10 10]);
end;


figure;
for n=1:length(pcs)
    subplot(4,4,n)
    plotpc = pcs{n};
    bar(plotpc(1,:));
    title(n);
end;


%%%%%%%%%%%%%%%%%%%
winv = inv(weights*sphere);
figure;
for n = 1:size(weights,1)
    subplot(4,4,n)
    g=bar(winv(:,n));hold on;
    set(gca,'ylim',[-.5 .5]);
    set(gca,'xlim',[0 17]);
    set(gca,'xtick',[1:16]);
    set(gca,'xgrid','on');
    set(gca,'fontsize',7);
        %set(g,'color',cb);
end;

figure;
for n = 1:size(weights,1)
    subplot(4,4,n)
    g=bar(weights(n,:),'g');hold on;
    set(gca,'ylim',[-1 1]);
    set(gca,'xlim',[0 17]);
    set(gca,'xtick',[1:16]);
    set(gca,'xgrid','on');
    set(gca,'fontsize',7);
        %set(g,'color',cb);
    title(int2str(n)); 
end;

% plot 3 winvs versus each other
p1 = icadata(:,7:11,1);    % cluster 1   7.5-10.5 Hz
p1=mean(p1,2);

p2 = icadata(:,11:13,2);    % cluster 2   10.5 - 12 Hz
p2=mean(p2,2);

p3 = icadata(:,8:11,3);    % cluster 3    8.25-10.5  Hz
p3=mean(p3,2);

p4 = icadata(:,33:36,4);    % cluster 4    27-30 Hz
p4=mean(p4,2);

p5 = icadata(:,1:3,5);    % cluster 5     3-4.5 Hz
p5=mean(p5,2);

p6 = icadata(:,17:19,6);    % cluster 6    15-16.5 Hz
p6=mean(p6,2);

p7 = icadata(:,7:10,7);    % cluster 7    7.5-9.75  Hz
p7=mean(p7,2);

p8 = icadata(:,27:29,8);    % cluster 8   22.5 - 24 Hz
p8=mean(p8,2);

p9 = icadata(:,5:6,9);    % cluster 9      6-6.75 Hz
p9=mean(p9,2);

p10 = icadata(:,21:23,10);    % cluster 10   18-19.5  Hz
p10=mean(p10,2);

p11 = icadata(:,18:20,11);    % cluster 11    15.75-17.25  Hz
p11=mean(p11,2);

p12 = icadata(:,24:27,12);    % cluster 12    20.25-22.5 Hz
p12=mean(p12,2);

p13 = icadata(:,3:4,13);    % cluster 13    4.5-5.25  Hz
p13=mean(p13,2);

p14 = icadata(:,3:5,14);    % cluster 14    4.5-6 Hz
p14=mean(p14,2);

p15 = icadata(:,15:16,15);    % cluster 15    13.5-14.25 Hz
p15=mean(p15,2);  

p16 = icadata(:,31:32,16);    % cluster 16     25.5-26.25  Hz
p16=mean(p16,2);

% Plot in 3D
allclust = {p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16};
for kk = 1:length(allclust)-2
    if kk == 1| kk == 5| kk == 9| kk == 13
        figure;  pl=1;
    end;    
    subplot(2,2,pl)
    plot3(allclust{kk},allclust{kk+1},allclust{kk+2},'*');hold on;
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'zgrid','on');    
    for k=1:14
        text(allclust{kk}(k),allclust{kk+1}(k),allclust{kk+2}(k),emos(k));
    end; 
    title([int2str(kk),' vs ',int2str(kk+1),' vs ',int2str(kk+2)]); pl = pl+1;
end;

% Plot in 2D
allclust = {p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16};
for kk = 1:length(allclust)-4
    if kk == 1| kk == 5| kk == 9| kk == 13
        figure;  pl=1;
    end;    
    subplot(2,2,pl)
    plot(allclust{kk},allclust{kk+4},'*');hold on;
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    for k=1:14
        text(allclust{kk}(k),allclust{kk+4}(k),emos(k));
    end; 
    title([int2str(kk),' vs ',int2str(kk+4)]); pl = pl+1;
end;

%%%%%%%%%%%%%%%%%%%
cb = jet(14);
figure; 
for r = 1:size(weights,2)
g=plot3(weights(1,emoorder(r)),weights(2,emoorder(r)),weights(3,emoorder(r)),'.');hold on;
set(g,'color',cb(r,:));
set(g,'MarkerSize',30);
end;
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');
g=colorbar;
set(g,'ytick',[2:4.7:65]);
set(g,'yticklabel',{'sad','grief','frustration','anger','fear','jealousy','emabarrass','awe', 'joy','happy','love' ,'compassion','content','relief'});


figure;
for n = 1:13
    subplot(4,4,n)
    for r = 1:size(weights,2)
        g=plot(weights(n,emoorder(r)),weights(n+1,emoorder(r)),'.');hold on;
        set(g,'color',cb(r,:));
        set(g,'markersize',20);
    end;
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'zgrid','on');
    title(int2str(n));
end;
g=colorbar;
set(g,'ytick',[2:4.7:65]);
set(g,'yticklabel',{'sad','grief','frustration','anger','fear','jealousy','emabarrass','awe', 'joy','happy','love' ,'compassion','content','relief'});

figure;
for n = 6:14
    subplot(3,4,n-5)
    for r = 1:size(weights,2)
        g=plot(weights(5,emoorder(r)),weights(n,emoorder(r)),'.');hold on;
        set(g,'color',cb(r,:));
        set(g,'markersize',20);
    end;
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'zgrid','on');
    title(int2str(n));
end;
g=colorbar;
set(g,'ytick',[2:4.7:65]);
set(g,'yticklabel',{'sad','grief','frustration','anger','fear','jealousy','emabarrass','awe', 'joy','happy','love' ,'compassion','content','relief'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Instead cluster single trial timef's from each emotion 
eeglab  
 subj1 =  [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46]; % jo74
'num','on','dipolesize',30
gdcomps = {subj1};

paths = {'/data/common1/emotion/jo74/'};
baseset = 'pre-post-eyes.set';
datset = {'awePress.set', 'frustrationPress.set','joyPress.set','angerPress.set','sadPress.set','happyPress.set','fearPress.set','lovePress.set' ,'jealousyPress.set','compassionPress.set','emabarrassPress.set','contentPress.set','griefPress.set','reliefPress.set','pre-post-eyes.set'};
lotime = -1500;
hitime = 1500;
maxfreq = 40;
wsize = 256;
btlim = [-2000 2000];
tlim = [-2000 3000];
prat = 2;
wave = [3 .5];
pl = 1;
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
nx=1;n=1;
EEG = pop_loadset( baseset,paths{nx}); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);            
data = squeeze(EEG.icaact(gdcomps{nx}(n),:,1:10));
[tf, freqs, times, itcvals] = timefreq(data, EEG.srate,'wavelet',wave, 'tlimits',btlim,'winsize',wsize,'freqs',[0 maxfreq],'padratio',prat);
tm = find(times> lotime & times< hitime); 
%alltrials = zeros(length(freqs)*length(tm),0);
alltrials = zeros(length(freqs)+length(freqs),0);
 ALLEEG=[]; EEG=[];
 p = 1;
for nx = 1:length(paths)    
    if ~isempty(gdcomps{nx})
        EEG = pop_loadset( baseset,paths{nx}); 
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        for n = 1:length(gdcomps{nx})
            data = squeeze(EEG.icaact(gdcomps{nx}(n),:,:));
            [tf, freqs, times, itcvals] = timefreq(data, EEG.srate,'wavelet',wave, 'tlimits',btlim,'winsize',wsize,'freqs',[0 maxfreq],'padratio',prat);
            basetfdata{n} = abs(tf).^2;
        end;
        ALLEEG=[]; EEG=[];clear tf itcvals data 
        fprintf('\nDone Fixation for Subject %s\n',int2str(nx));
        for ds = 8:length(datset)
            EEG = pop_loadset( datset{ds},paths{nx}); 
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            EEG.icaweights=wts;
            EEG.icasphere=sph;EEG.icawinv=[];
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);            
            for n = 1:length(gdcomps{nx})                
                data = squeeze(EEG.icaact(gdcomps{nx}(n),:,:));
                [tf, freqs, times, itcvals] = timefreq(data, EEG.srate,'wavelet',wave, 'tlimits',tlim,'winsize',wsize,'freqs',[0 maxfreq],'padratio',prat);
                tfdata{n} = abs(tf).^2;
            end;
            ALLEEG=[]; EEG=[];clear tf itcvals tf1 tf2 normtf data tfnorm
            %if p>1
            %    load /tmp/temp.mat
            %end;
            tm = find(times> lotime & times< hitime); 
            for nn = 1:length(gdcomps{nx})
                clear normtf
                tf1 = tfdata{nn};  %  tf = freqs X times X epochs
                tf1 = tf1(:,tm,:);
                tf1 = 10*log10(tf1);
                tf2 = basetfdata{nn};
                tf2 = 10*log10(tf2);
                tf2 = tf2(:,tm,:);
                tf2 = mean(tf2,3);
                tf2 = mean(tf2,2);
                for w = 1:size(tf2,1)
                    normtf(w,:,:) = tf1(w,:,:) - tf2(w);   % gives dB diff of cond-fix  (w=freqs)
                end;
                sdofmean = std(normtf,1,2);
                mnormtf = mean(normtf,2);
                mnormtf = squeeze(mnormtf);
                sdofmean = squeeze(sdofmean);
                comb = [mnormtf; sdofmean];
                 %normtf = reshape(normtf,size(normtf,1)*size(normtf,2),size(normtf,3));
                % makes into freqs*times X trials
                alltrials(:,end+1:end+size(comb,2)) = comb;
            end;
            %save /tmp/temp.mat alltrials
            %clear alltrials
            clear  tf1 tf2 normtf data tfnorm  tfdata 
            %p=p+1;
        end; 
    end; 
    fprintf('\nDone Subject: %s\n\n',int2str(nx));
end;
save /tmp/emos1-jo.mat alltrials    % on playing   (emos 1:7)
clear alltrials
pack
save /tmp/emos2-jo.mat alltrials    % on playing   (emos 8:14)
% find mean and std for each freq for all comps
icaresults = zeros (length(freqs),length(times),clusts);
for n = 1:size(winv,2)
    tmpcomp = winv (:,n)';  %makes a totalength X 1
    tmpcomp = reshape(tmpcomp,length(freqs),length(times));  %Good
    icaresults(:,:,n) = tmpcomp;
end; 
% OR:
 load /data/common1/emotion/jo74/emXsbcptr.mat   % freqs = 4800, comps*trials = 9344
 alltrials = zeros(4800,0);
for em = 1:length(alltremo)
    tmpemo = alltremo{em};
    alltrials(:,end+1:end+size(tmpemo,2)) = tmpemo;
    fprintf('\none more emo done...%s',int2str(em));
end;
 clear alltremo
save /tmp/alltrials.mat alltrials    
floatwrite(alltrials, '/tmp/alltrials.fdt');   % playing
/data/common/matlab/ica_linux2.4 < /data/common1/emotion/jo74/ClustSnglTrEmos.sc


cat allemo1.fdt allemo2.fdt  > alltrialemos.fdt   

% run ica
/data/common/matlab/ica_linux2.4 < /data/common1/emotion/ClustSnglTrEmos.sc

% results on 9500 10192, pca to 30
% too big for matlab ica:
% activations = PCA dim X time/freq
% weights = PCA dim X comp*trials
% winv =  comp*trials X PCA dim

% find number of trials in each emotion
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
gdcomps = {subj1}; nx=1;
 load /data/common1/emotion/jo74/emXsbcptr.mat
for em = 1:length(alltremo)
    numtrials(1,em) = size(alltremo{em},2)/length(gdcomps{nx});
end; save /data/common1/emotion/jo74/numtrials.mat numtrials


% results 
wts = floatread('/data/common1/emotion/jo74/tfXcmpEmoSnglTrPC40.wts',[40 4800]);
sph = floatread('/data/common1/emotion/jo74/tfXcmpEmoSnglTrPC40.sph',[4800 4800]);

%wts = floatread('/data/common1/emotion/tfXcmpSngltrPC30.wts',[30 5500]);
%sph = floatread('/data/common1/emotion/tfXcmpSngltrPC30.sph',[5500 5500]);
%wts = floatread('/data/common1/emotion/jo74/stdmeanwts-150.wts',[25 150]);
%sph = floatread('/data/common1/emotion/jo74/stdmeansph-150.sph',[150 150]);
%wts = floatread('/data/common1/emotion/jo74/trialwts-6125-.wts',[25 6125]);
%sph = floatread('/data/common1/emotion/jo74/trialsph-6125-.sph',[6125 6125]);
alltrials = floatread('/tmp/alltrials.fdt',[4800 9344]); % rand epochs tf
load /tmp/alltrials.mat

ws = wts*sph;   clear wts sph
winv = pinv(ws);
activations = ws*alltrials;
clear alltrials

% load a dataset and run one iteration of timef for times/freqs (see beg of script)


clusts = size(ws,1); times = times(tm);
icaresults = zeros (length(freqs),length(times),clusts);
for n = 1:size(winv,2)
    tmpcomp = winv (:,n)';  %makes a totalength X 1
    tmpcomp = reshape(tmpcomp,length(freqs),length(tm));  %Good
    icaresults(:,:,n) = tmpcomp;
end; 
figure;
for n = 1:20%clusts
    subplot(4,5,n)
    imagesc(times,freqs,icaresults(:,:,n),[-5 5]);hold on;
    set(gca,'ydir','norm');
    set(gca,'ticklength',[.02 .02]);    
    if n < 9
    set(gca,'xticklabel',[]);
    end;   
    set(gca,'fontsize',14);
    title(int2str(n));
end; colorbar;

% find where each emotion starts and stops for correlated trials only run
% USE ONLY WHEN TRIALS WERE CHOSEN BY CORRELATION***
load /data/common1/emotion/Corridx.mat 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for nx = 1:length(gdcomps)
    %xx=1;
    for em = 1:length(emos)
        oneemocorrs = allcorridx{em};
       % if em == 1 & nx == 1
       % ixmap(1,xx) = size(oneemocorrs{1},2)+size(oneemocorrs{2},2)+size(oneemocorrs{3},2)+size(oneemocorrs{4},2)+size(oneemocorrs{5},2)+size(oneemocorrs{6},2)+size(oneemocorrs{7},2)+size(oneemocorrs{8},2)+size(oneemocorrs{9},2)+size(oneemocorrs{10},2)+size(oneemocorrs{11},2)+ ...
       %     size(oneemocorrs{12},2)+size(oneemocorrs{13},2)+size(oneemocorrs{14},2)+size(oneemocorrs{15},2)+size(oneemocorrs{16},2);
       % else
       % ixmap(1,xx) = size(oneemocorrs{1},2)+size(oneemocorrs{2},2)+size(oneemocorrs{3},2)+size(oneemocorrs{4},2)+size(oneemocorrs{5},2)+size(oneemocorrs{6},2)+size(oneemocorrs{7},2)+size(oneemocorrs{8},2)+size(oneemocorrs{9},2)+size(oneemocorrs{10},2)+size(oneemocorrs{11},2)+ ...
       %     size(oneemocorrs{12},2)+size(oneemocorrs{13},2)+size(oneemocorrs{14},2)+size(oneemocorrs{15},2)+size(oneemocorrs{16},2)+ixmap(1,em-1);
       % end;
        for cp = 1:length(gdcomps{nx})
        ntcomp(1,cp) = size(oneemocorrs{cp},2);
        end;
        numtrials{em} = ntcomp;
        %xx = xx+1;
    end;
end;

% Find > zero activations for clusters 1:12
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for n = 1:4%size(activations,1)
    hiact{n} = find(activations(n,:) > -std(activations(n,:)));
end;
% find sum of positive activations for all components
for clust = 1:4
    numtot = 0; clear sums
    sc = hiact{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
            if nx>1
                skipact = allsubjix{nx-1};
            end;
            for em = 1:length(emos)
                for cp = 1:length(gdcomps{nx})
                    ntrials = numtrials{em}(cp); % will have to change to cell array for more subj
                    ft = find(sc>=numtot+1 & sc<=numtot+ntrials);
                    sumcomp = sum(activations(clust,sc(ft)));
                    sums(cp) = sumcomp/ntrials*100;
                    numtot = ntrials+numtot;
                end;             
                emosumsonesubj{em} = sums;
            end;            
        end;
        allsubjsums{nx} = emosumsonesubj;
    end;
    allclustsums{clust} = allsubjsums;
end;
% find for each cluster, the sum for each emotion for each subj above a certain threshold
%cols = summer(14);
emoset = [2,4,5,7,9,11,13,1,3,6,8,10,12,14];
cols = jet(14);
figure;
for clust = 1:13
    subplot(4,4,clust); sumcell = zeros(1,0); clear clustcomps allvals
    for nx = 1:length(gdcomps)
        for em = 1:length(emos)
            sumcell(1,end+1:end+size(allclustsums{clust}{nx}{em},2)) = allclustsums{clust}{nx}{em};
        end;
    end;
    thresh = mean(sumcell);   
    for nx = 1:length(gdcomps) % only plots one subj right now
        for em = 1:length(emos)
            sums = allclustsums{clust}{nx}{em};
            if thresh < 0
                abth = find(sums < thresh);
            else
                abth = find(sums > thresh);            
            end;
            subcomps = gdcomps{nx}(abth);
            clustcomps{em} = subcomps;
            allvals{em} = sums(abth);
            ph = bar(em,sum(sums(abth))); hold on;
            set(ph,'facecolor',cols(emoset==em,:));
            ph = text(em,10,emos{em});
            set(ph,'rotation',90);
        end;
    end;
    set(gca,'xlim',[0 15]);
    set(gca,'fontsize',16);
    set(gca,'xticklabel',[]);
    title(int2str(clust));    
end;
textsc('Sums of Positive scores from all comps contributing to each emotion in each cluster; act cutoff is std(act); sum cutoff is mean of sums','title');

% make component lists for each subj with sum above certain cutoff
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
fid = fopen('/data/common1/emotion/jo74/SnglTrimTrialHiActs','a');%
for clust = 1:13
    sumcell =zeros(1,0);
        for nx = 1:length(gdcomps)
            for em = 1:length(emos)
                sumcell(1,end+1:end+size(allclustsums{clust}{nx}{em},2)) = allclustsums{clust}{nx}{em};
            end;
        end;
        thresh =  mean(sumcell)+1.5*(std(sumcell));   
        fprintf(fid,'\nICA CLUSTER: %s  POSITIVE  0-30Hz all emotions single trial input;jo74 only; activations > std(act); PCA to 40;  used non-button locked epochs selected on high correlation with other trials; threshold is mean-std of sums for cluster',int2str(clust));
        fprintf(fid,'\nCutoff SCORE: %g ',thresh);
        %NEGATIVE POSITIVE
        for nx = 1:length(gdcomps)
            for em = 1:length(emos)
                sumcell = allclustsums{clust}{nx};
                sums = sumcell{em};
                if thresh < 0
                    abth = find(sums < thresh);
                else
                    abth = find(sums > thresh);            
                end;
                subcomps = gdcomps{nx}(abth);
                clustcomps{em} = subcomps;
                allvals{em} = sums(abth);
            end;
            allsubjcomps{nx} = clustcomps;
            allsubjvals{nx} = allvals;
        end;
        for wrtemo = 1:length(emos)
            fprintf(fid,'\n%sEMOTION: %s ','% ',emos{wrtemo});
            fprintf(fid,'\njo74= [%s]; ',int2str(allsubjcomps{1}{wrtemo})); 
            fprintf(fid,'\ngdcomps%s = {jo74};',emos{wrtemo});
            % add more lines for more subj   
        end;
        allclustcomps{clust} = allsubjcomps;
        allclustvals{clust} = allsubjvals;    
    end;
fclose(fid);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For using all trials (not selected for high correlation)
% Find > zero activations for clusters 1:12
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for n = 1:4%size(activations,1)
    hiact{n} = find(activations(n,:) > -std(activations(n,:)));
end;
% find sum of positive activations for all components
for clust = 1:4
    numtot = 0; clear sums
    sc = hiact{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
            if nx>1
                skipact = allsubjix{nx-1};
            end;
            for em = 1:length(emos)
                for cp = 1:length(gdcomps{nx})
                    ntrials = numtrials(em); % will have to change to cell array for more subj
                    ft = find(sc>=numtot+1 & sc<=numtot+ntrials);
                    sumcomp = sum(activations(clust,sc(ft)));
                    sums(cp) = sumcomp/ntrials*100;
                    numtot = ntrials+numtot;
                end;             
                emosumsonesubj{em} = sums;
            end;            
        end;
        allsubjsums{nx} = emosumsonesubj;
    end;
    allclustsums{clust} = allsubjsums;
end;
% find for each cluster, the sum for each emotion for each subj above a certain threshold
%cols = summer(14);
emoset = [2,4,5,7,9,11,13,1,3,6,8,10,12,14];
cols = jet(14);
figure;
for clust = 1:4
    subplot(2,2,clust); sumcell = zeros(1,0); clear clustcomps allvals
    for nx = 1:length(gdcomps)
        for em = 1:length(emos)
            sumcell(1,end+1:end+size(allclustsums{clust}{nx}{em},2)) = allclustsums{clust}{nx}{em};
        end;
    end;
    thresh = mean(sumcell);   
    for nx = 1:length(gdcomps) % only plots one subj right now
        for em = 1:length(emos)
            sums = allclustsums{clust}{nx}{em};
            if thresh < 0
                abth = find(sums < thresh);
            else
                abth = find(sums > thresh);            
            end;
            subcomps = gdcomps{nx}(abth);
            clustcomps{em} = subcomps;
            allvals{em} = sums(abth);
            ph = bar(em,sum(sums(abth))); hold on;
            set(ph,'facecolor',cols(emoset==em,:));
            ph = text(em,10,emos{em});
            set(ph,'rotation',90);
        end;
    end;
    set(gca,'xlim',[0 15]);
    set(gca,'fontsize',16);
    set(gca,'xticklabel',[]);
    title(int2str(clust));    
end;
textsc('Sums of Positive scores from all comps contributing to each emotion in each cluster; act cutoff is std(act); sum cutoff is mean of sums','title');

% make component lists for each subj with sum above certain cutoff
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
fid = fopen('/data/common1/emotion/jo74/SnglTrimTrialHiActs','a');%
for clust = 1:13
    sumcell =zeros(1,0);
        for nx = 1:length(gdcomps)
            for em = 1:length(emos)
                sumcell(1,end+1:end+size(allclustsums{clust}{nx}{em},2)) = allclustsums{clust}{nx}{em};
            end;
        end;
        thresh =  mean(sumcell)+1.5*(std(sumcell));   
        fprintf(fid,'\nICA CLUSTER: %s  POSITIVE  0-30Hz all emotions single trial input;jo74 only; activations > std(act); PCA to 40;  used non-button locked epochs selected on high correlation with other trials; threshold is mean-std of sums for cluster',int2str(clust));
        fprintf(fid,'\nCutoff SCORE: %g ',thresh);
        %NEGATIVE POSITIVE
        for nx = 1:length(gdcomps)
            for em = 1:length(emos)
                sumcell = allclustsums{clust}{nx};
                sums = sumcell{em};
                if thresh < 0
                    abth = find(sums < thresh);
                else
                    abth = find(sums > thresh);            
                end;
                subcomps = gdcomps{nx}(abth);
                clustcomps{em} = subcomps;
                allvals{em} = sums(abth);
            end;
            allsubjcomps{nx} = clustcomps;
            allsubjvals{nx} = allvals;
        end;
        for wrtemo = 1:length(emos)
            fprintf(fid,'\n%sEMOTION: %s ','% ',emos{wrtemo});
            fprintf(fid,'\njo74= [%s]; ',int2str(allsubjcomps{1}{wrtemo})); 
            fprintf(fid,'\ngdcomps%s = {jo74};',emos{wrtemo});
            % add more lines for more subj   
        end;
        allclustcomps{clust} = allsubjcomps;
        allclustvals{clust} = allsubjvals;    
    end;
fclose(fid);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


% print component pairs to file that are highly correlated
% find best threshold and run only that one in script above to generate 'allclustcomps'
% can't do for single trial ica on only correlated trials (different trials removed for all subjects)
printnames = {'jo74'};
alpha = .01; 
fid = fopen('/data/common1/emotion/SnglTrialCorrelations','a');%
for clust = 1:10
    fprintf(fid,'\nICA CLUSTER: %s  Correlation of Activations within subject, between all components above specified threshold of activity;   ',int2str(clust));
    %NEGATIVE POSITIVE
    fprintf(fid,'\nThreshold for SUM: %g ',thresh);
    fprintf(fid,'\nCutoff is bootstrap sig of: %g ',alpha);
    clear allcorr nxcorr nxsig allsig
    for nx = 1:length(gdcomps)
        ixmap = allsubjix{nx};
        if nx>1
            skipact = allsubjix{nx-1};
        end;
        ntrials = numtrials(nx); % will have to change to cell array for more subj
        for em = 1:length(emos)
            cpofinterest = allclustcomps{clust}{nx}{em};
            nxcorr=[];   nxsig = [];
            for cp1 = 1:length(gdcomps{nx})-1
                for cp2 = cp1+1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(cp1),cpofinterest) & ismember(gdcomps{nx}(cp2),cpofinterest)
                    if nx == 1 & em == 1
                            [rsignif,accarray] = bootstat(activations(clust,(cp1-1)*ntrials+1:ntrials*cp1),activations(clust,(cp2)*ntrials+1:ntrials*(cp2+1)),'[res,indx,indy,corrs] = matcorr(arg1,arg2);','alpha',alpha,'bootside','upper');
                            nxsig(cp1,cp2) = rsignif ;                            
                            [corr,indx,indy,corrs] = matcorr(activations(clust,(cp1-1)*ntrials+1:ntrials*cp1),activations(clust,(cp2)*ntrials+1:ntrials*(cp2+1)));
                            nxcorr(cp1,cp2) = corr; 
                        elseif nx == 1 & em > 1
                            [rsignif,accarray] = bootstat(activations(clust,ixmap(em-1)+ntrials*(cp1-1)+1:ixmap(em-1)+ntrials*cp1),activations(clust,ixmap(em-1)+ntrials*(cp2-1)+1:ixmap(em-1)+ntrials*cp2),'[res,indx,indy,corrs] = matcorr(arg1,arg2);','alpha',alpha,'bootside','upper');
                            nxsig(cp1,cp2) = rsignif;
                            [corr,indx,indy,corrs] = matcorr(activations(clust,ixmap(em-1)+ntrials*(cp1-1)+1:ixmap(em-1)+ntrials*cp1),activations(clust,ixmap(em-1)+ntrials*(cp2-1)+1:ixmap(em-1)+ntrials*cp2));
                            nxcorr(cp1,cp2) = corr;
                        else
                            [rsignif,accarray] = bootstat(activations(clust,skipact+ixmap(em-1)+ntrials*(cp1-1)+1:skipact+ixmap(em-1)+ntrials*cp1),activations(clust,skipact+ixmap(em-1)+ntrials*(cp2-1)+1:skipact+ixmap(em-1)+ntrials*cp2),'[res,indx,indy,corrs] = matcorr(arg1,arg2);','alpha',alpha,'bootside','upper');
                            nxsig(cp1,cp2) = rsignif;
                            [corr,indx,indy,corrs] = matcorr(activations(clust,skipact+ixmap(em-1)+ntrials*(cp1-1)+1:skipact+ixmap(em-1)+ntrials*cp1),activations(clust,skipact+ixmap(em-1)+ntrials*(cp2-1)+1:skipact+ixmap(em-1)+ntrials*cp2));
                            nxcorr(cp1,cp2) = corr;   
                            numtot = ntrials+numtot;

                        end; % for if nx and em               
                    end; % for if ismember
                end; % for cp2
            end;  % for cp1
            emcorr{em} = nxcorr;
            emsig{em} = nxsig;
        end; % for em
        allsubjcorr{nx} =  emcorr;
        allsubjsig{nx} =  emsig;      

        for em = 1:length(emos)
            oneemocorr = emcorr{em};
            oneemosig = emsig{em};
            [x y] = find(oneemocorr > oneemosig);        
            gdcompsa{nx} = gdcomps{nx}(x');
            gdcompsb{nx} =  gdcomps{nx}(y');
            fprintf(fid,'\nEMOTION: %s   ',emos{em});
            fprintf(fid,'\n%sa = [%s]; ',printnames{nx},int2str(gdcompsa{nx}));
            fprintf(fid,'\n%sb = [%s]; ',printnames{nx},int2str(gdcompsb{nx}));
        end; 
    end;  % for nx
end;
fclose(fid);


% plot comp pairs for each emo
eeglab
EEG = pop_loadset( 'sources.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

gdcompsa = {jo74a};
gdcompsb = {jo74b};
gdcompsa = {awea,frusta,joya,angera,sada,happya,feara,lovea,jealousa,compassa,embarra,contenta,griefa, reliefa};
gdcompsb = {aweb,frustb,joyb,angerb,sadb,happyb,fearb,loveb,jealousb,compassb,embarrb,contentb,griefb,reliefb};


emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
emoset = [1,3,12,2,4,9,7,13];
figure; sb=1;
for emo = 1:8  
    subplot(3,3,sb);sb = sb+1;
    % find sources for gdcompsa set
    if ~isempty(gdcompsa{emoset(emo)})        
        allsourcesa =  EEG.dipfit.model(1);
        dipsources= EEG.dipfit.model(1);
        for w = 1:length(gdcompsa{emoset(emo)})
            dipsources(1,w)= EEG.dipfit.model(gdcompsa{emoset(emo)}(w));
        end;           
        allsourcesa(end+1:end+size(dipsources,2)) = dipsources;
        dipsources=[];    
        allsourcesa(1) = [];
    end;
    % find sources for gdcompsb set
    if ~isempty(gdcompsb{emoset(emo)})        
        allsourcesb =  EEG.dipfit.model(1);
        dipsources= EEG.dipfit.model(1);
        for w = 1:length(gdcompsb{emoset(emo)})
            dipsources(1,w)= EEG.dipfit.model(gdcompsb{emoset(emo)}(w));
        end;           
        allsourcesb(end+1:end+size(dipsources,2)) = dipsources;
        dipsources=[];
        allsourcesb(1) = [];
    end; 

    % for 3d plotting (can't use projlines currently
    allab = allsourcesa;
    allab(end+1:end+size(allsourcesb,2)) = allsourcesb;
    dipplot(allab,'dipolelength',0,'gui','off','dipolesize',25,'view',[0 0 1],'color',{'r'},'image','mri','projimg','on','spheres','on');
    
    htmp = get(gca,'children'); 
    pl=1;clear  totnum 
    for idx = 1:length(htmp)
        ctmp = get(htmp(idx),'userdata');
        if isstruct(ctmp)
           if  isempty(str2num(ctmp.rv(3:5)))
                x = str2num(ctmp.rv(3:4));            
            else
                x = str2num(ctmp.rv(3:5)); 
            end;
            totnum(pl,:) = x;pl = pl+1;        
        end;
    end;
    totnum = unique(totnum);
    totsource = size(totnum,1);  %***

    pl=1;p=1;clear newstruct
    for idx = 1:length(htmp)-1
        ctmp = get(htmp(idx),'userdata');
        if isstruct(ctmp)
             if  isempty(str2num(ctmp.rv(3:5)))
                if str2num(ctmp.rv(3:4))<=totsource/2 
                    ctmp.cnum = str2num(ctmp.rv(3:4));
                else
                    ctmp.cnum = str2num(ctmp.rv(3:4))-(totsource/2);
                end; pl = pl+1;
            else
                if str2num(ctmp.rv(3:5))<=totsource/2 
                    ctmp.cnum = str2num(ctmp.rv(3:5));
                else
                    ctmp.cnum = str2num(ctmp.rv(3:5))-(totsource/2);
                end; pl = pl+1;
            end;            
            newstruct(p) = ctmp;p=p+1;
        end;
    end;    

    ff={newstruct.cnum};
    ff=cell2mat(ff); clear allcoo 
    haf = find(ff==1);
    haf = haf(1);
    pl=1;bilines=[];
    for cc = 1:totsource/2
        clear coords
        ids = find(ff(1:haf)==cc);
        tt ={newstruct(ids).pos3d};
        tp = tt{1};
        coords(1,:) = tp(1,:);;
        if size(tt,2)>1
            tp = tt{1};
            tpbi(1,:) =  tp(1,:);
            tp = tt{2};
            tpbi(2,:) =  tp(1,:);
            bilines{pl} =  tpbi;pl=pl+1;
        end;
        ids = find(ff(haf+1:end)==cc);
        tt ={newstruct(ids+haf).pos3d};
        tp = tt{1};
        coords(2,:) = tp(1,:);;
        if size(tt,2)>1
            tp = tt{1};
            tpbi(1,:) =  tp(1,:);
            tp = tt{2};
            tpbi(2,:) =  tp(1,:);
            bilines{pl} =  tpbi;pl=pl+1;
        end;
        allcoo{cc} = coords;
    end;

    for cc = 1:size(allcoo,2)
        tp = allcoo{cc};    
        ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
        set(ph,'linewidth',2.5);
    end;
    if ~isempty(bilines)
        for cc = 1:size(bilines,2)
            tp = bilines{cc};    
            ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
            set(ph,'linestyle','--');
        end;
    end;
    camzoom(.92);
    view(37,28);
    ph=title(emos{emoset(emo)});
    set(ph,'color','y'); 
    set(ph,'fontsize',16);
end;


% plot activations of different clusters vs each other
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
emoset = [4,2,9,7,5,13,11,10,8,12,6,14,3,1]; % bad to good
emoset = [1:14]; % in temporal order
ct1 = [1,2,3,4,8];
cols = cool(14);
figure; sm = 1; clear pc
    for c1 = 1:length(ct1)-1
        clust1= ct1(c1);
        for c2 = c1+1:length(ct1)
            clust2 = ct1(c2);
            for em = 1:14
                %subplot(4,4,em)
                if emoset(em) == 1
                    %ph = plot(activations(clust1,1:ixmap(1)),activations(clust2,1:ixmap(1)),'r.');hold on;
                    pc(1,emoset(em)) = sum(activations(clust1,1:ixmap(1)));
                    pc(2,emoset(em)) = sum(activations(clust2,1:ixmap(1)));
                else
                    %ph = plot(activations(clust1,ixmap(emoset(em)-1):ixmap(emoset(em))),activations(clust2,ixmap(emoset(em)-1):ixmap(emoset(em))),'b.');hold on;
                    pc(1,emoset(em)) = sum(activations(clust1,ixmap(emoset(em)-1):ixmap(emoset(em))));
                    pc(2,emoset(em)) = sum(activations(clust2,ixmap(emoset(em)-1):ixmap(emoset(em))));
                end;
                %set(ph,'color',cols(em,:));
                %set(ph,'markersize',5);
                %plot([0 0],[-12 12],'k');
                %plot([-12 12],[0 0],'k');
                %set(gca,'xlim',[-10 10]);
                %set(gca,'ylim',[-10 10]);
                %set(gca,'fontsize',15);
                %title(emos{emoset(em)});
            end;
            subplot(3,4,sm); sm = sm+1;
            plot(pc(1,:),pc(2,:),'k.');
            for pp = 1:size(pc,2)
                text(pc(1,pp)+50,pc(2,pp),emos{pp});
            end;      
            title(['Cluster ',int2str(clust1),' vs Cluster ',int2str(clust2)]);
        end;
    end;

ph =textsc('Activations of Cluster 1 (~10 Hz;x-axis) and Cluster 2 (y-axis): Bad:Good order','title');
set(ph,'fontsize',15);

% Find cluster center and plot 
figure;silhouette(km,ones(size(km,1),1));
dst  = pdist(km, 'euclidean');
lnk = linkage(dst);

% plot only a subset of  components

% Find > zero activations for clusters 1:12
for n = 1:12%size(activations,1)
    hiact{n} = find(activations(n,:) > 0);
end;

figure; sm = 1; clear pc
clust = clust1; % to pick highly weighted comps from   
for c1 = 1:length(ct1)-1
    clust1= ct1(c1);
    for c2 = c1+1:length(ct1)
        clust2 = ct1(c2);
        sc = hiact{clust};
        for em = 1:length(emos)
            %subplot(4,4,em)
            for nx = 1:length(gdcomps)
                if ~isempty(gdcomps{nx})
                    ixmap = allsubjix{nx};
                    if nx>1
                        skipact = allsubjix{nx-1};
                    end;
                    ntrials = numtrials(nx); % will have to change to cell array for more subj
                    clear sums
                    if nx == 1 & em == 1
                        for cp = 1:length(gdcomps{nx})
                            ft = find(sc>=ntrials*(cp-1)+1 & sc<=ntrials*cp);
                            sumcomp = sum(activations(clust,sc(ft)));
                            sums(cp) = sumcomp/ntrials*100;
                            if sums(cp) > 30
                                %ph = plot(activations(clust1,ntrials*(cp-1)+1:ntrials*cp),activations(clust2,ntrials*(cp-1)+1:ntrials*cp),'b.');hold on;
                                %set(ph,'color',cols(em,:));
                                %set(ph,'markersize',5);
                                %plot([0 0],[-12 12],'k');                            
                                pc(1,emoset(em)) = sum(activations(clust1,ntrials*(cp-1)+1:ntrials*cp));
                                pc(2,emoset(em)) = sum(activations(clust2,ntrials*(cp-1)+1:ntrials*cp));
                            end;
                        end;             
                    elseif nx == 1 & em >1                 
                        for cp = 1:length(gdcomps{nx})
                            ft = find(sc>=ixmap(em-1)+ntrials*(cp-1)+1 & sc<=ixmap(em-1)+ntrials*cp);
                            sumcomp = sum(activations(clust,sc(ft)));
                            sums(cp) = sumcomp/ntrials*100;
                            if sums(cp) > 30
                                %ph = plot(activations(clust1,ixmap(em-1)+ntrials*(cp-1)+1:ixmap(em-1)+ntrials*cp),activations(clust2,ixmap(em-1)+ntrials*(cp-1)+1:ixmap(em-1)+ntrials*cp),'b.');hold on;
                                %set(ph,'color',cols(em,:));
                                %set(ph,'markersize',5);
                                %plot([0 0],[-12 12],'k');                            
                                pc(1,emoset(em)) = sum(activations(clust1,ixmap(em-1)+ntrials*(cp-1)+1:ixmap(em-1)+ntrials*cp));
                                pc(2,emoset(em)) = sum(activations(clust2,ixmap(em-1)+ntrials*(cp-1)+1:ixmap(em-1)+ntrials*cp));
                            end;
                        end;
                    else
                        for cp = 1:length(gdcomps{nx})
                            ft = find(sc>=skipact+ixmap(em-1)+ntrials*(cp-1)+1 & sc<=skipact+ixmap(em-1)+ntrials*cp);
                            sumcomp = sum(activations(clust,sc(ft)));
                            sums(cp) = sumcomp/ntrials*100;
                            if sums(cp) > 30
                                %ph = plot(activations(clust1,skipact+ixmap(em-1)+ntrials*(cp-1)+1:skipact+ixmap(em-1)+ntrials*cp),activations(clust2,skipact+ixmap(em-1)+ntrials*(cp-1)+1:skipact+ixmap(em-1)+ntrials*cp),'b.');hold on;
                                %set(ph,'color',cols(em,:));
                                %set(ph,'markersize',5);
                                %plot([0 0],[-12 12],'k');                            
                            end;
                        end;                    
                    end;   
                end;            
            end;
            %plot([0 0],[-12 12],'k');
            %plot([-12 12],[0 0],'k');
            %set(gca,'xlim',[-10 10]);
            %set(gca,'ylim',[-10 10]);
            %set(gca,'fontsize',15);
            %title(emos{emoset(em)}); 
        end;        
        subplot(3,4,sm); sm = sm+1;
        plot(pc(1,:),pc(2,:),'k.'); hold on;
        for pp = 1:size(pc,2)
            text(pc(1,pp)+2,pc(2,pp),emos{pp});
        end;      
        title(['Cluster ',int2str(clust1),' vs Cluster ',int2str(clust2)]);        
    end;
end;
ph =textsc('Activations of Cluster 1 (~10 Hz;x-axis) and Cluster 8 (~6 Hz;y-axis) Only comps highly weighted in Clust : Bad:Good order','title');
set(ph,'fontsize',15);



