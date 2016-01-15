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
        for ds = 1:length(datset)
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
                alltrials(:,end+1:end+size(comb,2)) = comb;
            end;
            clear  tf1 tf2 normtf data tfnorm  tfdata 
        end; 
    end; 
    fprintf('\nDone Subject: %s\n\n',int2str(nx));
end;
save /tmp/emos1-jo.mat alltrials    % on playing   (emos 1:7)
clear alltrials
pack
save /tmp/emos2-jo.mat alltrials    % on playing   (emos 8:14)

alltrialsm0 = alltrials;
for k=1:size(alltrials,1)
alltrialsm0(k,:) = alltrials(k,:) - mean(alltrials(k,:)')'*ones(1,size(alltrials,2));
end

% Run ICA on resulting matrix
channels = size(alltrials,1)  %  150
frames = size(alltrials,2)   % 29184
% In Linux terminal
floatwrite(alltrials, '/tmp/stdmeantrials2.fdt'); 
floatwrite(alltrialsm0, '/tmp/allstdtrials.fdt'); 

cat stdmeantrials1.fdt stdmeantrials2.fdt  > allstdtrials.fdt   
% run ica
/data/common/matlab/ica_linux2.4 < /data/common1/emotion/jo74/stdTrialICA.sc

%  pca to 25
% activations = PCA dim X time/freq
% weights = PCA dim X comp*trials
% winv =  comp*trials X PCA dim

% find number of trials in each emotion
for nx = 1:length(datset)
    EEG = pop_loadset( datset{nx},paths{1}); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    numtrials(1,nx) = size(EEG.data,3);
    ALLEEG=[];
end; save /data/common1/emotion/jo74/numtrials.mat numtrials
% results 
wts = floatread('/data/common1/emotion/jo74/stdmeanwts-150.wts',[25 150]);
sph = floatread('/data/common1/emotion/jo74/stdmeansph-150.sph',[150 150]);
alltrials = floatread('allstdtrials.fdt',[150 29184]);
ws = wts*sph;   clear wts sph
winv = pinv(ws);
activations = ws*alltrials;
clear alltrials

% load a dataset and run one iteration of timef for times/freqs (see beg of script)
subj1 =  [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46]; % jo74 for ICA on mean pwr and std
gdcomps = {subj1};

paths = {'/data/common1/emotion/jo74/'};
baseset = 'pre-post-eyes.set';
datset = {'awePress.set', 'frustrationPress.set','joyPress.set','angerPress.set','sadPress.set','happyPress.set','fearPress.set','lovePress.set' ,'jealousyPress.set','compassionPress.set','emabarrassPress.set','contentPress.set','griefPress.set','reliefPress.set'};
 eeglab
nx=1;n=1;
esph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
ewts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
EEG = pop_loadset(  'awePress.set','/data/common1/emotion/jo74/'); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG.icaweights=ewts;
EEG.icasphere=esph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);            
data = squeeze(EEG.icaact(10,:,1:10));
[tf, freqs, times, itcvals] = timefreq(data, EEG.srate,'wavelet',[3 .5], 'tlimits',[-2000 3000],'winsize',256,'freqs',[0 20],'padratio',2);
clusts = size(ws,1);

stdonly = winv(76:150,:);
stdonly = stdonly';
justmean = winv(1:75,:);
justmean = justmean';
figure;
for p = 1:25
    subplot(5,5,p)
    plot(justmean(p,:)); hold on;
    plot(justmean(p,:)-stdonly(p,:),'r');
    plot(justmean(p,:)+stdonly(p,:),'r');
    set(gca,'xlim',[0 50]);
    set(gca,'xtick',[0:10:40]);
    set(gca,'xgrid','on');
    title(int2str(p));
end;

    % find where each emotion starts and stops (no break between subjects, when applicable)
load /data/common1/emotion/jo74/numtrials.mat 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
subjname = {'jo74'};
for nx = 1:length(gdcomps)
    xx=1;
    for em = 1:length(emos)
        if em == 1 & nx == 1
            ixmap(1,xx) = numtrials(1,em)*length(gdcomps{nx});    
        else
            ixmap(1,xx) = numtrials(1,em)*length(gdcomps{nx})+ixmap(1,em-1);
        end;
        xx = xx+1;
    end;
    allsubjix{nx} = ixmap;
end;

% Find > zero activations for clusters 1:12
for n = 1:12%size(activations,1)
    %hiact{n} = find(activations(n,:));
    hiact{n} = find(activations(n,:) > 0 );
end;
% find sum of positive activations for all components
for clust = 1:12
    sc = hiact{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
            ixmap = allsubjix{nx};
            if nx>1
                skipact = allsubjix{nx-1};
            end;
            ntrials = numtrials(nx); % will have to change to cell array for more subj
            for em = 1:length(emos)
                clear sums
                if nx == 1 & em == 1
                    for cp = 1:length(gdcomps{nx})
                        ft = find(sc>=ntrials*(cp-1)+1 & sc<=ntrials*cp);
                        sumcomp = sum(activations(clust,sc(ft)));
                        sums(cp) = sumcomp/ntrials*100;
                    end;             
                elseif nx == 1 & em >1                 
                    for cp = 1:length(gdcomps{nx})
                        ft = find(sc>=ixmap(em-1)+ntrials*(cp-1)+1 & sc<=ixmap(em-1)+ntrials*cp);
                        sumcomp = sum(activations(clust,sc(ft)));
                        sums(cp) = sumcomp/ntrials*100;
                    end;
                else
                    for cp = 1:length(gdcomps{nx})
                        ft = find(sc>=skipact+ixmap(em-1)+ntrials*(cp-1)+1 & sc<=skipact+ixmap(em-1)+ntrials*cp);
                        sumcomp = sum(activations(clust,sc(ft)));
                        sums(cp) = sumcomp/ntrials*100;
                    end;                    
                end;   
                emosumsonesubj{em} = sums;
            end;            
        end;
        allsubjsums{nx} = emosumsonesubj;
    end;
    allclustsums{clust} = allsubjsums;
end;
% find for each cluster, the sum for each emotion for each subj above a certain threshold
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
emoset = [4,2,9,7,5,13,11,10,8,12,6,14,3,1]; % bad to good
%emoset = [1:14]; % in temporal order
thresh = 25;  
%cols = summer(14);
cols = cool(14);
figure;
for clust = 1:12
subplot(3,4,clust)
    for nx = 1:length(gdcomps) % only plots one subj right now
        for em = 1:length(emos)
            sums = allclustsums{clust}{nx}{emoset(em)};
                if thresh < 0
                    abth = find(sums < thresh);
                else
                    abth = find(sums > thresh);            
                end;
                subcomps = gdcomps{nx}(abth);
                clustcomps{emoset(em)} = subcomps;
                allvals{emoset(em)} = sums(abth);
                ph = bar(em,sum(sums(abth))); hold on;
                set(ph,'facecolor',cols(em,:));
                ph = text(em,10,emos{emoset(em)});
                set(ph,'rotation',90);
        end;
    end;
    set(gca,'xlim',[0 15]);
    set(gca,'xtick',[]);
    set(gca,'fontsize',16);
    title(int2str(clust));  
    allclustcomps{clust} = clustcomps;
end;
% plot activations of different clusters vs each other
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
emoset = [4,2,9,7,5,13,11,10,8,12,6,14,3,1]; % bad to good
emoset = [1:14]; % in temporal order
ct1 = [1,2,3,4,5,6,7,8,9,10,11,12];
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
        subplot(9,9,sm); sm = sm+1;
        plot(pc(1,:),pc(2,:),'k.');
        for pp = 1:size(pc,2)
            text(pc(1,pp)+50,pc(2,pp),emos{pp});
        end;      
        title(['Cluster ',int2str(clust1),' vs Cluster ',int2str(clust2)]);
    end;
end;
axcopy

%What components make up each cluster?

% make component lists for each subj with sum above certain cutoff
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};
fid = fopen('/data/common1/emotion/SnglTrialpwrOnly','a');%
for thresh =25:10:55 % sum of activations 
    for clust = 1:12
        fprintf(fid,'\nICA CLUSTER: %s POSITIVE    0-50Hz all emotions single trial input;jo74 only; activations > 0; PCA to 25; score based on SUM of > zero weightings (activations)',int2str(clust));
        fprintf(fid,'\nCutoff sum SCORE: %g ',thresh);
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
end;
fclose(fid);



figure;
for p = 1:12
    subplot(3,4,p)
hist(activations(p,:),30);
kurt(activations(p,:))
end;

figure;
bar(1,sum(abs(activations(1,:))));hold on;
bar(2,sum(abs(activations(2,:))));hold on;
bar(3,sum(abs(activations(3,:))));hold on;
bar(4,sum(abs(activations(4,:))));hold on;
bar(5,sum(abs(activations(5,:))));hold on;
bar(6,sum(abs(activations(6,:))));hold on;
bar(7,sum(abs(activations(7,:))));hold on;
bar(8,sum(abs(activations(8,:))));hold on;
bar(9,sum(abs(activations(9,:))));hold on;
bar(10,sum(abs(activations(10,:))));hold on;

