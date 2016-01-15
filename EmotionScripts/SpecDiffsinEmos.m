% goes through each dataset (emotion) and computes psd for all comps, then plots results in freq bands
eeglab
%subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];%  jo74
subj1 = [1,3:22,25:30,34,36];  % tl81  23,35 =pulse; 
subj2 = [1,2,4:9,12,13,15:19,21,22,24,26,29];% mi83  14,23 = pulse; 31,34,46,51:muscle
subj3 = [5,6,7,10,12,13,14,15,17,18,19,22];   % ms82  11,37 = pulse; 32 = muscle?
subj4 = [1:6,9:12,14:19,22:25,28,34]; % kw78
subj5 = [1:6,8,9,10,13:15,17,19,20,23,25,27,30];  % tv81  11 = pulse?? 
subj6 = [1:5,8:10,12:14,16,17,19,22,25];  % js75  7 is pulse
subj7 = [1:4,6,8:14,16:21,23,24,28];  % sr81   22,27 = pulse
subj8 = [1:7,11,13:16,18,20,21,23,25,37,38,39,40];  % an70  9,12 = pulse

gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8};

%datsets = {{'emo-1-241.set','emo-2-241.set','emo-3-241.set','emo-4-241.set','emo-5-241.set'},{'emo-1-248.set','emo-2-248.set','emo-3-248.set','emo-4-248.set','emo-5-248.set','emo-6-248.set'},{'emo-1-238.set','emo-2-238.set','emo-3-238.set','emo-4-238.set','emo-5-238.set'},{'emo-1-250.set','emo-2-250.set','emo-3-250.set','emo-4-250.set','emo-5-250.set'},{'emo-1-180.set','emo-2-180.set','emo-3-180.set','emo-4-180.set','emo-5-180.set'},{'emo-1-253.set','emo-2-253.set','emo-3-253.set','emo-4-253.set','emo-5-253.set'}};

paths = {'/data/common2/emotion/tl81/','/data/common2/emotion/mi83/','/data/common2/emotion/ms82/','/data/common2/emotion/kw78/','/data/common2/emotion/tv81/','/data/common2/emotion/js75/','/data/common2/emotion/sr81/','/data/common2/emotion/an70/'};
datpaths = {'/data/common2/emotion/tl81/ersps/','/data/common2/emotion/mi83/ersps/','/data/common2/emotion/ms82/ersps/','/data/common2/emotion/kw78/ersps/','/data/common2/emotion/tv81/ersps/','/data/common2/emotion/js75/ersps/','/data/common2/emotion/sr81/ersps/','/data/common2/emotion/an70/ersps/'};

sphs = {'/data/common2/emotion/tl81/sph241.sph','/data/common2/emotion/mi83/sph248.sph','/data/common2/emotion/ms82/sph238.sph','/data/common2/emotion/kw78/sph250pc110.sph','/data/common2/emotion/tv81/sph180-90.sph','/data/common2/emotion/js75/sph253pc100.sph','/data/common2/emotion/sr81/sph249pc100.sph','/data/common2/emotion/an70/sph250pc100.sph'};
wtss = {'/data/common2/emotion/tl81/wts241.wts','/data/common2/emotion/mi83/wts248.wts','/data/common2/emotion/ms82/wts238.wts','/data/common2/emotion/kw78/wts250pc110.wts','/data/common2/emotion/tv81/wts180-90.wts','/data/common2/emotion/js75/wts253pc100.wts','/data/common2/emotion/sr81/wts249pc100.wts','/data/common2/emotion/an70/wts250pc100.wts'};

%             tl81       mi83      ms82     kw78      tv81      js75      sr81     an70
sphsize = {[241 241],[248 248],[238 238],[250 250],[180 180],[253 253],[249 249],[250 250]};
wtssize = {[160 241],[160 248],[110 238],[110 250], [90 180],[100 253],[100 249],[100 250]};

nx=5;    
sph=floatread(sphs{nx},sphsize{nx}); 
wts=floatread(wtss{nx},wtssize{nx}); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find psd DIFFERENCES between all emos and a given test condition
nx = ;
emoset = {'prebase.set','awe.set','frustration.set','joy.set','anger.set','happy.set','sad.set','love.set','fear.set' ,'compassion.set','jealousy.set','content.set','grief.set','relief.set','disgust.set','excite.set','postbase.set'};
baseidx = 1; % 1 = prebase, etc.
cd (datpaths{nx})
load EmoSpecs.mat % loads: trialspecs freqs comment
fr = find(freqs < 40);
clear sigemo sigstdev
for cmp = 1:length(gdcomps{nx})
    basespec = trialspecs{baseidx};
    basespec = basespec(gdcomps{nx}(cmp),fr,:);
    for k = 1:length(trialspecs)
        oneemo = trialspecs{k};
        oneemo = oneemo(gdcomps{nx}(cmp),fr,:);
        for ff = 1:size(oneemo,2)
            emotest = oneemo(1,ff,:);
            basetest = basespec(1,ff,:);
             [H,P,CI,STATS] = ttest2(basetest,emotest,.05);
            if H == 0
                sigemo(cmp,:,k) = 0; 
            else                
                sigstdev(cmp,ff,k) = std(emotest)/sqrt(size(emotest,3));
                sigemo(cmp,ff,k) = mean(emotest,3)-mean(basetest,3);
            end;
        end;
        fprintf('\nOne More Emo Done: %s',emoset{k});
    end;
        fprintf('\nOne More Comp Done: %s of %s',int2str(cmp),int2str(length(gdcomps{nx})));
end;
% plot the results (by emos)
cols = jet(17);
figure; row = round(sqrt(length(gdcomps{nx}))); col = round(sqrt(length(gdcomps{nx})));
for cmp = 1:length(gdcomps{nx})
    subplot(row, col, cmp);
    for k = 1:length(emoset)
        plotdiffs = sigemo(cmp,:,k);
        ph = plot(freqs(fr),plotdiffs,'k');hold on;
        set(ph,'color',cols(k,:));
        title(int2str(gdcomps{nx}(cmp)));
    end;
end;

                
% run ANOVA on each comp for a freq or freq band
emos = {'prebase','awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief','disgust','excite','postbase'};

fr = find(freqs == 10);
fr = find(freqs >= 35 & freqs <=  45);

for y = 1:length(gdcomps{nx})
stonecomp = zeros(min(manyev),14); clear   onecomp
    for emo = 1:length(datset)
        onecomp = allsubjtrialpwr{nx}{emo};
        onecomp = onecomp(gdcomps{nx}(y),fr,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp = mean(onecomp,2); % makes a 1(comps) X 1(freqs) X trials
        onecomp = squeeze(onecomp); % makes a 1 X trials
        stonecomp(1:size(onecomp,1),emo) = onecomp;        
    end;
    statmat = stonecomp(1:min(manyev),:);
    [P,ANOVATAB,STATS] = anova1(statmat,emos','off');
figure;   [COMPARISON,MEANS,H] =multcompare(STATS,.01,'on');
end;

% plot the ratio of two freq bands
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % all comps
subj1 = [5,6,7,10,12,13,15,18,22,24,25,40,46];  % comps with 10-20 Hz peaks
subj1 = [12,13,15,16,18,24,28,46];  % comps with 10-20 Hz peaks
gdcomps = {subj1};

fr1 = find(freqs == 9);
fr1 = find(freqs >= 9 & freqs <=  11);
fr2 = find(freqs >= 18 & freqs <=  22);
cols = jet(14);
figure;
for y = 1:length(gdcomps{nx})
    subplot(3,3,y); clear stdiff  onecomp1 onecomp2
    for emo = 2:15%length(datset)
        onecomp1 = allsubjtrialpwr{nx}{emo};
        onecomp1 = onecomp1(gdcomps{nx}(y),fr1,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp1 = mean(onecomp1,2);
        onecomp2 = allsubjtrialpwr{nx}{emo};
        onecomp2 = onecomp2(gdcomps{nx}(y),fr2,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp2 = mean(onecomp2,2);
        %stdiff = sqrt((std(onecomp1)/sqrt(size(onecomp1,3))).^2 + (std(onecomp2)/sqrt(size(onecomp2,3))).^2);
        onecomp1 = squeeze(onecomp1);
        onecomp2 = squeeze(onecomp2);
        plotone = onecomp1./onecomp2;
        stdiff = std(plotone/sqrt(size(plotone,1))); hold on;
        plotone = mean(plotone,1);
        plot(mean(onecomp1),mean(onecomp2),'k.');
        ph = text(mean(onecomp1)+.01,mean(onecomp2),emos{emo});
        %ph = bar(emo,plotone);hold on;
        %if emo == 1 | emo == 16
        %    set(ph,'facecolor',[.5 .5 .5]);
        %else
        %    set(ph,'facecolor',cols(find(emoset==emo)-1,:));
        %end;     
        %plot([emo emo],[plotone-stdiff plotone+stdiff],'k');        
        %ph = text(emo,.1,emos{emo});
        %set(ph,'rotation',90);
    end;
    title(['Comp: ',int2str(gdcomps{nx}(y))]);
    %set(gca,'xlim',[0 17]);
    %set(gca,'xticklabel',[]);
end;

for y = 1:length(gdcomps{nx})
    for emo = 1:length(datset)
        onecomp1 = allsubjtrialpwr{nx}{emo};
        onecomp1 = onecomp1(gdcomps{nx}(y),fr1,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp1 = mean(onecomp1,2);
        onecomp2 = allsubjtrialpwr{nx}{emo};
        onecomp2 = onecomp2(gdcomps{nx}(y),fr2,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp2 = mean(onecomp2,2);
        onecomp1 = squeeze(onecomp1);
        onecomp2 = squeeze(onecomp2);
        statone = onecomp1./onecomp2;
        statmat(1:size(statone,1),emo) = statone;
    end;
    statmatall = statmat(1:min(manyev),:);
    [P,ANOVATAB,STATS] = anova1(statmatall,emos','off');
figure;   [COMPARISON,MEANS,H] =multcompare(STATS,.01,'on');
end;

% plot three ratios versus each other
subj1 = [12,15,18]; % list only 3
gdcomps = {subj1};
fr1 = find(freqs >= 9 & freqs <=  11);
fr2 = find(freqs >= 18 & freqs <=  22);
figure;
y = 1;
   clear stdiff  onecomp1 onecomp2
    for emo = 2:15%length(datset)
        onecomp1 = allsubjtrialpwr{nx}{emo};
        onecomp1 = onecomp1(gdcomps{nx}(y),fr1,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp1 = mean(onecomp1,2);
        onecomp2 = allsubjtrialpwr{nx}{emo};
        onecomp2 = onecomp2(gdcomps{nx}(y),fr2,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp2 = mean(onecomp2,2);
        onecomp1 = squeeze(onecomp1);
        onecomp2 = squeeze(onecomp2);
        plotone = onecomp1./onecomp2;
        plotone1 = mean(plotone,1);clear onecomp1 onecomp2
        onecomp1 = allsubjtrialpwr{nx}{emo};
        onecomp1 = onecomp1(gdcomps{nx}(y+1),fr1,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp1 = mean(onecomp1,2);
        onecomp2 = allsubjtrialpwr{nx}{emo};
        onecomp2 = onecomp2(gdcomps{nx}(y+1),fr2,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp2 = mean(onecomp2,2);
        onecomp1 = squeeze(onecomp1);
        onecomp2 = squeeze(onecomp2);
        plotone = onecomp1./onecomp2;
        plotone2 = mean(plotone,1);
        clear onecomp1 onecomp2
        onecomp1 = allsubjtrialpwr{nx}{emo};
        onecomp1 = onecomp1(gdcomps{nx}(y+2),fr1,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp1 = mean(onecomp1,2);
        onecomp2 = allsubjtrialpwr{nx}{emo};
        onecomp2 = onecomp2(gdcomps{nx}(y+2),fr2,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp2 = mean(onecomp2,2);
        onecomp1 = squeeze(onecomp1);
        onecomp2 = squeeze(onecomp2);
        plotone = onecomp1./onecomp2;
        plotone3 = mean(plotone,1);
        ph = plot3(plotone1,plotone2,plotone3,'k.');hold on;        
        ph = text(plotone1,plotone2,plotone3,emos{emo});        
    end;
    set(gca,'xgrid','on');
    set(gca,'ygrid','on');
    set(gca,'zgrid','on');
 
% plot three powers versus each other
subj1 = [7,11,18]; % list only 3
gdcomps = {subj1};
fr = find(freqs >= 18 & freqs <=  22);
figure;
y = 1;
clear stdiff  onecomp1 onecomp2
for emo = 2:15%length(datset)
    onecomp1 = allsubjtrialpwr{nx}{emo};
    onecomp1 = onecomp1(gdcomps{nx}(y),fr,:); % makes a 1(comps) X 3(freqs) X trials
    onecomp1 = mean(onecomp1,2);
    onecomp1 = squeeze(onecomp1);
    plotone1 = mean(onecomp1,1);
    clear onecomp1 onecomp2
    onecomp1 = allsubjtrialpwr{nx}{emo};
    onecomp1 = onecomp1(gdcomps{nx}(y+1),fr1,:); % makes a 1(comps) X 3(freqs) X trials
    onecomp1 = mean(onecomp1,2);
    onecomp1 = squeeze(onecomp1);
    plotone2 = mean(onecomp1,1);
    clear onecomp1 onecomp2
    onecomp1 = allsubjtrialpwr{nx}{emo};
    onecomp1 = onecomp1(gdcomps{nx}(y+2),fr1,:); % makes a 1(comps) X 3(freqs) X trials
    onecomp1 = mean(onecomp1,2);
    onecomp1 = squeeze(onecomp1);
    plotone3 = mean(onecomp1,1);
    ph = plot3(plotone1,plotone2,plotone3,'k.');hold on;        
    ph = text(plotone1,plotone2,plotone3,emos{emo});        
end;
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');


% percent change