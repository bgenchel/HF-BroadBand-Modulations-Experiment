% clusters subj-comps by power spectra (minus pre-baseline) across emos


eeglab
% input subj comp lists

%datsets = {{'emo-1-241.set','emo-2-241.set','emo-3-241.set','emo-4-241.set','emo-5-241.set'},{'emo-1-248.set','emo-2-248.set','emo-3-248.set','emo-4-248.set','emo-5-248.set','emo-6-248.set'},{'emo-1-238.set','emo-2-238.set','emo-3-238.set','emo-4-238.set','emo-5-238.set'},{'emo-1-250.set','emo-2-250.set','emo-3-250.set','emo-4-250.set','emo-5-250.set'},{'emo-1-180.set','emo-2-180.set','emo-3-180.set','emo-4-180.set','emo-5-180.set'},{'emo-1-253.set','emo-2-253.set','emo-3-253.set','emo-4-253.set','emo-5-253.set'}};

datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-243.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph243pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','sph243pc100.sph','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[243 243],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 243],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250]};

emoset = {'prebase.set','awe.set','frustration.set','joy.set','anger.set','happy.set','sad.set','love.set','fear.set' ,'compassion.set','jealousy.set','content.set','grief.set','relief.set','disgust.set','excite.set','postbase.set'};
%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nx = 1:length(gdcomps) 
nx=4;    
sph=floatread(sphs{nx},sphsize{nx}); 
wts=floatread(wtss{nx},wtssize{nx}); 

% now ready to run psd
for k = 1:length(emoset)
    EEG = pop_loadset( emoset{k}, paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    if exist('EEG.event(2)')
    ft = EEG.event(2).latency;
    else
        ft = 1;
    end;    
    for evtm = ft+128:256:size(EEG.data,2)-128  % go up by 1 sec to create non-overlapping 1 sec epochs
        EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 'fake';        
    end;
    EEG = pop_epoch( EEG,{'fake'} , [-.5 .5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-500 500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
    EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
    clear cspec;
    for cmp = 1:length(gdcomps{nx})
        for trial = 1:size(EEG.icaact,3)
            [Pxx,F] = psd(EEG.icaact(gdcomps{nx}(cmp),:,trial),256,EEG.srate,256);% much smoother
            cspec(gdcomps{nx}(cmp),:,trial) = log10(Pxx(:,1))';
        end; % to trial
    end; % to  cmp
    trialspecs{k} = cspec;
    fprintf('\none more emo done:  %s\n',emoset{k});
    ALLEEG=[];EEG=[];
end;
freqs = F;
comment = 'trials are non-overlapping 1 sec epochs from button pressing period. used psd [Pxx,F] = psd(EEG.icaact(gdcomps{nx}(cmp),:,trial),256,EEG.srate,256); and converted to log:log10(Pxx(:,1)); emo order: prebase,awe,frustration,joy,anger,happy,sad,love,fear ,compassion,jealousy,content,grief,relief,disgust,excite,postbase';
cd (datpaths{nx})
save EmoSpecs.mat trialspecs freqs comment
clear trialspecs
end;
%%%%%%%%% END COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

figure;pl = 1;
for c = 2;
    for tt = 1:size(cspec,3)
        subplot(10,10,pl)
        plot(cspec(gdcomps{nx}(c),:,tt));pl = pl+1;
        set(gca,'xlim',[1 30]);
    end;
end;
% Now plot the results
% Visualize result as bar graphs with error bars
% Each comp gets a graph displaying one freq band and all emotions
nx = 4; cd (datpaths{nx}); load  EmoSpecs.mat
subjnames = {'tl81','mi83','ms82','kw78','tv81','js75'};
emos = {'prebase','awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite','postbase'};  %
EEG = pop_loadset( 'sources.set', paths{nx});
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
%, 'frustration','disgust',
emoorder = [1,17,2,4,6,8,10,12,14,16,3,5,7,9,11,13,15];
if ceil(length(gdcomps{nx})/3) > 6
    row = 6;
else
    row = ceil(length(gdcomps{nx})/3); 
end;
col = 9;
fr = find(freqs >= 10 & freqs <= 10);
cols(1,:) = [.5 .5 .5];
cols(2,:) = [.5 .5 .5];
cols(end+1:end+length(emos),:) = jet(length(emos));

figure; pl = 1;
for y = 1:length(gdcomps{nx})
    if y == 19
xf = freqs(fr);
ph=textsc(['Subject ',subjnames{nx},'Emotion psd Power at freqs: ',num2str(xf(1)),'-',num2str(xf(end)),' Hz'],'title');
set(ph,'fontsize',14);
        figure;pl=1;
    end;    
    subplot(row,col,pl); clear sdev onecomp
    topoplot(EEG.icawinv(:,gdcomps{nx}(y)),EEG.chanlocs,  'electrodes', 'off', 'plotrad',0.5);
    title(int2str(gdcomps{nx}(y)));    
    subplot(row,col,pl+1:pl+2);
    for emo = 1:length(emos)
        onecomp = trialspecs{emo};
        onecomp = onecomp(gdcomps{nx}(y),fr,:); % makes a 1(comps) X 3(freqs) X trials
        onecomp = mean(onecomp,2);
        sdev = std(onecomp)/sqrt(size(onecomp,3));
        onecomp = mean(onecomp,3);
        onecomp = squeeze(onecomp);  % makes a 1 X 1 to plot        
        ph = bar(emo,onecomp);hold on;
        set(ph,'facecolor',cols(find(emoorder==emo),:));
        plot([emo emo],[onecomp-sdev onecomp+sdev],'k');        
        ph = text(emo,-.1,emos{emo});
        set(ph,'rotation',90);
    end;
    title(['Comp: ',int2str(gdcomps{nx}(y))]);
    set(gca,'xlim',[0 length(emos)+1]);
    set(gca,'xticklabel',[]);
    pl = pl+3;
end;
xf = freqs(fr)
ph=textsc(['Subject ',subjnames{nx},'Emotion psd Power at freqs: ',num2str(xf(1)),'-',num2str(xf(end)),' Hz'],'title');
set(ph,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use spectral diffs from baseline to cluster across subjects
% collect a matrix of freqs(emo-prebase) X subj*components*emotions
 clear sph wts cspec
nx = 9; cd (datpaths{nx});load EmoSpecs.mat
alldat = zeros(length(freqs),0);
%for nx = 1:length(gdcomps)
%    cd (datpaths{nx})
%    load EmoSpecs.mat
    for cp= 1:length(gdcomps{nx})
        basespec = trialspecs{1}; % use prebase as baseline
        basespec = mean(basespec,3);
        for e = 1:length(trialspecs)   
            tmpmat = trialspecs{e};
                lengthemos(e) = size(tmpmat,3);
            for tr = 1:size(tmpmat,3)
                tmpinput = tmpmat(gdcomps{nx}(cp),:,tr); % makes 1 X 129 (length(freqs))
                tmpinput = tmpinput-basespec(gdcomps{nx}(cp),:);%  prebase will be subtracted from all
                alldat(:,end+1) = tmpinput';
            end;
        end;
    fprintf('\nOne more Component done: %s of %s\n',int2str(cp),int2str(length(gdcomps{nx})));
    end;clear tmpmat tmpinput basespec
    fprintf('\nOne more subj done: %s of %s\n',int2str(nx),int2str(length(gdcomps)));
%end;
floatwrite(alldat,'/data/common2/emotion/Subj4Specs.fdt');
save trialnums.mat lengthemos freqs cp
size(alldat)
lengthemos
% to find all lengthemos
for nx = 1:length(paths)
    cd (paths{nx})
    load trialnums.mat
    numtrials{nx} = lengthemos;
end;
save /data/common2/emotion/numtrials.mat numtrials

%numbers of trials in each emotion (in order) from each subject
nx=1  [87   142   189   130   140   146   139   119   136   133    95   103   157 102   147   119    85]   
nx=2  [101   212   201   186   162   213   238   207   212   189   178   159   178 159   209   236    99]   
nx=3   [110   147   164   156   139   193   186   212   261   208   187   198   203 186   183   196  116] 
nx=4   [107   249   211   166   195   177   153   200   220   200   199   213   203 196   215   205   110]   
nx=5   [95   211   161   160   191   164   146   165   130   183   190   134   212 128   152   156   101]   
nx=6   [103    92   120    78   167   128    99   125    80   132   119   141   200 109   141    92   105]
nx=7   [94   129   129   123   238   133   127   138   143   184   134   142   228  183   161   119   101]
nx=8   [102   156   133    82    81   111    97   189    85    73    87    76   115 68    66    91   107]
nx=9   [108   181   233   215   213   250   159   233   162   145   143   181   161 126   110   157   112]

numtrials = {[87   142   189   130   140   146   139   119   136   133    95   103   157 102   147   119    85],[101   212   201   186   162   213   238   207   212   189   178   159   178 159   209   236    99],[110   147   164   156   139   193   186   212   261   208   187   198   203 186   183  196  116],[107   249   211   166   195   177   153   200   220   200   199   213   203 196   215   205   110],[95   211   161   160   191   164   146   165   130   183   190   134   212 128   152   156   101],[103    92   120    78   167   128    99   125    80   132   119   141   200 109   141    92   105],[94   129   129   123   238   133   127   138   143   184   134   142   228  183   161   119   101],[102   156   133    82    81   111    97   189    85    73    87    76   115 68    66    91   107],[108   181   233   215   213   250   159   233   162   145   143   181   161 126   110   157   112]};

%%%%%%%%%%%%%  RUN ICA ON RESULTS  %%%%%%%%%%%5
cd /data/common2/emotion/
cat Subj1Specs.fdt  Subj2Specs.fdt  Subj3Specs.fdt  Subj4Specs.fdt  Subj5Specs.fdt  Subj6Specs.fdt  Subj7Specs.fdt Subj8Specs.fdt > allspecs.fdt  
%size of 1:5 = 62901  + 62780  + 36540  + 51504  + 58938  +38589  + 52626 + 36099 + 40446 = 
% 129 X 368057
% run ica in linux
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/ClustPwrICA.sc

wts = floatread('/data/common2/emotion/wts129pc50.wts',[50 129]);
sph = floatread('/data/common2/emotion/sph129pc50.sph',[129 129]);
alldat = floatread('/data/common2/emotion/allspecs.fdt',[129 364888]);
ws = wts*sph;
winv = pinv(ws);
activations = ws*alldat;
load /data/common2/emotion/numtrials.mat  % cell of numbers of trials (numtrials), freqs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% INDIVIDUAL ICA CLUSTERING FOR EACH SUBJECT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj1 = [1,3:22,25:30,34,36];  % tl81  23,35 =pulse; 
subj2 = [1,2,4:9,12,13,15:19,21,22,24,26,29];% mi83  14,23 = pulse; 31,34,46,51:muscle
subj3 = [5,6,7,10,12,13,14,15,17,18,19,22];   % ms82  11,37 = pulse; 32 = muscle?
subj4 = [1:5,8:10,12:14,16,17,19,22,25];  % js75  7 is pulse
% all before are button presses
subj5 = [1:6,9:12,14:19,22:25,28,34]; % kw78
subj6 = [1:6,8,9,10,13:15,17,19,20,23,25,27,30];  % tv81  11 = pulse?? 
subj7 = [1:4,6,8:14,16:21,23,24,28];  % sr81   22,27 = pulse
subj8 = [1:7,11,13:16,18,20,21,23,25,37,38,39,40];  % an70  9,12 = pulse
subj9 = [2,3,4,5,7,8,9,10,11,15,16,17,18,25];  % sg75  14 = pulse

gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9};
load /data/common2/emotion/numtrials.mat  % cell of numbers of trials (numtrials), freqs

subjdat = {'Subj1Specs.fdt' , 'Subj2Specs.fdt' , 'Subj3Specs.fdt',  'Subj4Specs.fdt' , 'Subj5Specs.fdt', 'Subj6Specs.fdt' , 'Subj7Specs.fdt', 'Subj8Specs.fdt', 'Subj9Specs.fdt'};  
subjwts = {'Subj1pc20.wts' , 'Subj2pc20.wts' , 'Subj3pc20.wts',  'Subj4pc20.wts' , 'Subj5pc20.wts', 'Subj6pc20.wts' , 'Subj7pc20.wts', 'Subj8pc20.wts', 'Subj9pc20.wts'};  
subjsph = {'Subj1pc20.sph' , 'Subj2pc20.sph' , 'Subj3pc20.sph',  'Subj4pc20.sph' , 'Subj5pc20.sph', 'Subj6pc20.sph' , 'Subj7pc20.sph', 'Subj8pc20.sph', 'Subj9pc20.sph'};  
datsize = [ 62901    62780   36540   51504   58938  38589   52626  36099  40446];% = 440423

toi = {[1:5,7,8,9,10,11,13],[1:3,5:8],[3,4,5,8,9,10,11,15],[1,3,4,6,7,8,9],[1:8],[1:4,6,7,9,10],[1:5,7,8,9,11],[1:3,5,6,7],[1,2,3,6,7,8,10,12]};
for nx = 1:length(gdcomps)
wts = floatread(['/data/common2/emotion/',subjwts{nx}],[20 129]);
sph = floatread(['/data/common2/emotion/',subjsph{nx}],[129 129]);
alldat = floatread(['/data/common2/emotion/',subjdat{nx}],[129 datsize(nx)]);
ws = wts*sph;
winv = pinv(ws);
activations = ws*alldat;
figure;
xmat = 100/size(activations(1,420:550),2);
ph = plot([xmat:xmat:100],sort(activations(3,420:550)));hold on;
xmat = 100/size(activations(1,1:87),2);
ph = plot([xmat:xmat:100],sort(activations(3,1:87)),'r');hold on;
plot([0 100],[0 0],'k');

figure;
for n = 1:length(toi{nx})
    subplot(round(sqrt((length(toi{nx}))))+1,round(sqrt((length(toi{nx})))),n)
    ph = plot(freqs,winv(:,toi{nx}(n)));hold on;
    set(ph,'linewidth',2);
    set(gca,'xlim',[0 50]);
    set(gca,'xgrid','on');
    set(gca,'xtick',[0:10:50]);
    set(gca,'xticklabel',[0:10:50]);
    plot([get(gca,'xlim')],[0 0],'k-');
    set(gca,'ticklength',[.02 .02]); 
    set(gca,'fontsize',14);
%    if n < size(winv,2)-3
%        set(gca,'xticklabel',[]);
%    end;   
    title(int2str(n));
end;axcopy 
ph= textsc(['Templates of Spectral Differences from Pre-baseline in Subject: ',int2str(nx)],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
end;


toi = {[1:5,7,8,9,10,11,13],[1:3,5:8],[3,4,5,8,9,10,11,15],[1,3,4,6,7,8,9],[1:8],[1:4,6,7,9,10],[1:5,7,8,9,11],[1:3,5,6,7],[1,2,3,6,7,8,10,12]};
emos = {'prebase','awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
% Find > zero activations for clusters 1:12
nclust = length(toi{nx});  fac = 1.5;
cut = 1; % 1 for positive wts, 0 for neg
for n = 1:nclust  %size(activations,1)
    if cut == 1
        tpact = sort(activations(n,:));
        cutnum = tpact(round(length(tpact)*.8));
        hiact{n} = find(activations(n,:) >= cutnum); % cutnum
       %hiact{n} = find(activations(toi{nx}(n),:) > fac*std(activations(toi{nx}(n),:)));
    else
        tpact = sort(activations(n,:));
        cutnum = tpact(round(length(tpact)*.2));
        hiact{n} = find(activations(n,:) <= cutnum); % cutnum
        hiact{n} = find(activations(toi{nx}(n),:) < -fac*std(activations(toi{nx}(n),:)));    
    end;
end;
% find sum of positive activations for all components
for clust = 1:nclust
    numtot = 0; clear sums
    sc = hiact{clust};
    ntrials = numtrials{nx};
    for em = 1:length(emos)
        for cp = 1:length(gdcomps{nx})
            emotrials = ntrials(em); 
            ft = find(sc>=numtot+1 & sc<=numtot+emotrials);
            sumcomp = sum(activations(clust,sc(ft)));
            sums(cp) = sumcomp/emotrials*100;
            numtot = emotrials+numtot;
        end;             
        emosumsonesubj{em} = sums;  clear sums
    end;            
    allclustsums{toi{nx}(clust)} = emosumsonesubj;
    fprintf('\n One More Cluster done: %s of %s',int2str(clust),int2str(nclust));
end;

% find for each cluster, the sum for each emotion for each subj above a certain threshold
%cols = summer(14);
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
cols(1,:) = [.5 .5 .5];
cols(2:16,:) = jet(15);
cols(17,:) = [.5 .5 .5];

figure; col = 4; row = ceil(nclust/2); pl = 1;
for clust = 1:nclust
    subplot(row,col,pl)
    if cut==1
        plot(freqs,winv(:,toi{nx}(clust)));hold on;pl = pl+1;
    else
        plot(freqs,-winv(:,toi{nx}(clust)));hold on;pl = pl+1;    
    end;
    set(gca,'xtick',[0:10:50]);
    set(gca,'xticklabel',[0:10:50]);
    set(gca,'xlim',[0 50]);
    set(gca,'xgrid','on');
    plot([get(gca,'xlim')],[0 0],'k-');
    set(gca,'ticklength',[.02 .02]); 
    set(gca,'fontsize',14);
    %title(int2str(toi{nx}(clust)));
    subplot(row,col,pl); sumcell = zeros(1,0); clear clustcomps allvals
    for em = 1:length(emos)
        onesum = allclustsums{toi{nx}(clust)}{em};
        sumcell(1,end+1:end+size(onesum,2)) = onesum;
    end;
    if cut == 1
        thresh = mean(sumcell(find(sumcell)))+std(sumcell);
    else
        thresh = mean(sumcell(find(sumcell)))-std(sumcell);    
    end;
    for em = 1:length(emos)
        plotsum = zeros(0);sumcell = zeros(0);
        onesum = allclustsums{toi{nx}(clust)}{emoset(em)};
        sumcell(1,end+1:end+size(onesum,2)) = onesum;
        if thresh < 0
            tmpsum = sumcell(sumcell <  thresh); 
            plotsum(end+1:end+length(tmpsum)) = tmpsum;
        ph = bar(em,-sum(plotsum)); hold on;
        else
            tmpsum = sumcell(sumcell >  thresh); 
            plotsum(end+1:end+length(tmpsum)) = tmpsum;
        ph = bar(em,sum(plotsum)); hold on;
        end;
        set(ph,'facecolor',cols(em,:));
            ph = text(em,2,emos{emoset(em)});
        set(ph,'rotation',90);
    end;
    set(gca,'xlim',[0 18]);
    set(gca,'fontsize',14);
    set(gca,'xticklabel',[]);
    title(int2str(toi{nx}(clust)));   pl = pl+1; 
end;
if cut == 1
    cutname = 'Positive'; perc = 'top';
else
    cutname = 'Negative';perc = 'bottom';
end;
%ph = textsc(['Sums of ',cutname,' scores; All comps; Subj: ',int2str(nx),', for each emotion; act cutoff is ',num2str(fac),'*std(act); sum cutoff: mean(sum)+std(sums)'],'title');
ph = textsc(['Sums of ',cutname,' scores; All comps; Subj ',int2str(nx),', for each emotion; act cutoff is ',perc,' 20%; sum cutoff: mean(sum)+std(sums)'],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Which components contribute to each cluster
printname = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9'};
realname = {'tl81','mi83','ms82','kw78','tv81','js75','sr81','an70','sg75'};
clear clustcomps subjcomps emocomps
fid = fopen('/data/common2/emotion/PwrDiffClusts','a');%
for nx = 1:length(gdcomps)
    cut = 0; % 1 for positive wts, 0 for neg
if cut == 1
    cutname = 'Positive'; perc = 'top';else
    cutname = 'Negative';perc = 'bottom';end;

    load /data/common2/emotion/numtrials.mat  % cell of numbers of trials (numtrials), freqs
    wts = floatread(['/data/common2/emotion/',subjwts{nx}],[20 129]);
    sph = floatread(['/data/common2/emotion/',subjsph{nx}],[129 129]);
    alldat = floatread(['/data/common2/emotion/',subjdat{nx}],[129 datsize(nx)]);
    ws = wts*sph;winv = pinv(ws);activations = ws*alldat;
    % Find > zero activations for clusters 1:12
    nclust = length(toi{nx}); 
    for n = 1:nclust  %size(activations,1)
        if cut == 1
            tpact = sort(activations(n,:));
            cutnum = tpact(round(length(tpact)*.8));
            hiact{n} = find(activations(n,:) >= cutnum); % cutnum
        else
            tpact = sort(activations(n,:));
            cutnum = tpact(round(length(tpact)*.2));
            hiact{n} = find(activations(n,:) <= cutnum); % cutnum
        end;
    end;
    % find sum of positive activations for all components
    for clust = 1:nclust
        numtot = 0; clear sums
        sc = hiact{clust};
        ntrials = numtrials{nx};
        for em = 1:length(emos)
            for cp = 1:length(gdcomps{nx})
                emotrials = ntrials(em); 
                ft = find(sc>=numtot+1 & sc<=numtot+emotrials);
                sumcomp = sum(activations(clust,sc(ft)));
                sums(cp) = sumcomp/emotrials*100;
                numtot = emotrials+numtot;
            end;             
            emosumsonesubj{em} = sums;  clear sums
        end;            
        allclustsums{toi{nx}(clust)} = emosumsonesubj;
        fprintf('\n One More Cluster done: %s of %s',int2str(clust),int2str(nclust));
    end;
    for clust = 1:nclust
        fprintf(fid,'\n\nICA CLUSTER: %s, %s    0-50Hz; Power Diffs from prebase spectrum clustered by ICA; PCA to 20; Activation cut off is %s * std(activations(clust,:));sum cutoff is mean of sums +- std(sums)',int2str(toi{nx}(clust)),cutname,num2str(fac));
        fprintf(fid,'\n%s  SUBJECT: %s  ','%    ',printname{nx});
        sumcell = zeros(1,0); 
        for em = 1:length(emos)
            onesum = allclustsums{toi{nx}(clust)}{em};
            sumcell(1,end+1:end+size(onesum,2)) = onesum;
        end;
        if cut == 1
            thresh = mean(sumcell(find(sumcell)))+std(sumcell);
        else
            thresh = mean(sumcell(find(sumcell)))-std(sumcell);    
        end;
        fprintf(fid,'\n%s  Threshold (mean+-std): %s  ','%    ',int2str(thresh));
        clear subjcomps emocomps
        for em = 1:length(emos)
            fprintf(fid,'\n%s  EMOTION: %s  ','%    ',emos{em});
            onesum = allclustsums{toi{nx}(clust)}{em};
            if thresh < 0
                abth = find(onesum < thresh);
            else
                abth = find(onesum > thresh);            
            end;
            subjcomps = gdcomps{nx}(abth);
            fprintf(fid,'\n%s= [%s]; %s   %s',emos{em},int2str(subjcomps),'   % ',realname{nx});
        end;
        fprintf(fid,'\n%sgdcomps = {prebase,awe,frustration,joy,anger,happy,sad,love,fear,compassion,jealousy,content,grief,relief,disgust,excite,postbase};',printname{nx});       
    end;
end;
fclose(fid);




%%%%%  Find correlated activation pairs
% print component pairs to file that are highly correlated
% eeglab % need eeglab for bootstat!!
toi = {[1:5,7,8,9,10,11,13],[1:3,5:8],[3,4,5,8,9,10,11,15],[1,3,4,6,7,8,9],[1:8],[1:4,6,7,9,10],[1:5,7,8,9,11],[1:3,5,6,7],[1,2,3,6,7,8,10,12]};

    load /data/common2/emotion/numtrials.mat  % cell of numbers of trials (numtrials), freqs
printname = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9'};
realname = {'tl81','mi83','ms82','kw78','tv81','js75','sr81','an70','sg75'};
alpha = .00001; 
pos = 1;  % 1 indicates taking positive correlation, 0 negative corr
if pos == 1
    cutname = 'Positive'; else
    cutname = 'Negative'; end;
for nx = 1:length(gdcomps)
fid = fopen('/data/common2/emotion/SnglTrialPwrCorrel','a');%
    wts = floatread(['/data/common2/emotion/',subjwts{nx}],[20 129]);
    sph = floatread(['/data/common2/emotion/',subjsph{nx}],[129 129]);
    alldat = floatread(['/data/common2/emotion/',subjdat{nx}],[129 datsize(nx)]);
    ws = wts*sph;winv = pinv(ws);activations = ws*alldat;
    % Find > zero activations for clusters 1:12
    nclust = length(toi{nx}); 
    fprintf(fid,'\n%s  SUBJECT: %s,     (%s)  ','%    ',printname{nx},realname{nx});
    for clust = 1:nclust
        fprintf(fid,'%s \nICA CLUSTER: %s; %s Weights; Correlation of Activations within subject, significant using bootstat/distfit to get p val indicated below; ALL comps tested ',' %  ',int2str(clust),cutname);
    fprintf(fid,'\n%s      Significance level (alpha):   %0.5f,   ','%    ',alpha);
        nxcorr=[]; possig = [];negsig = []; printcorr = [];
        ntrials = numtrials{nx}; 
        for cp1 = 1:length(gdcomps{nx})-1
        numtot = 0;   
            for em = 1:length(emos)
                emotrials = ntrials(em);
                for cp2 = cp1+1:length(gdcomps{nx})
                    dat1 = activations(clust,numtot+(cp1-1)*sum(ntrials)+1:numtot+(cp1-1)*sum(ntrials)+emotrials);
                    dat2 = activations(clust,numtot+(cp2-1)*sum(ntrials)+1:numtot+(cp2-1)*sum(ntrials)+emotrials);
                    [rsignif,accarray] = bootstat(dat1,dat2,'[res,indx,indy,corrs] = matcorr(arg1,arg2);','alpha',alpha,'distfit','on');
                    possig(cp1,cp2,em) = rsignif(1,2) ;
                    negsig(cp1,cp2,em) = rsignif(1,1) ;                        
                    [corr,indx,indy,corrs] = matcorr(dat1,dat2);
                    nxcorr(cp1,cp2,em) = corr; 
                end; % for cp2
                numtot = numtot + emotrials;
            end;% end for em          
        end;  % for cp1 
        possig(find(possig<0))=.7;
        negsig(find(negsig>0))=-.7;
        if pos == 1
            for em = 1:length(emos)
                printcorr = [];
                %fprintf(fid,'\n%s  Emotion: %s     ','%    ',emos{em});               
                [x y] = find(nxcorr(:,:,em) > possig(:,:,em));  
                gdcompsa{nx} = gdcomps{nx}(x');
                gdcompsb{nx} =  gdcomps{nx}(y');
                for r = 1: length(x)
                    printcorr(1,r) = nxcorr(x(r),y(r),em);     
                end;                    
                fprintf(fid,'\n%sa= [%s]; ',emos{em},int2str(gdcompsa{nx}));
                fprintf(fid,'\n%sb= [%s]; ',emos{em},int2str(gdcompsb{nx}));
                if ~isempty(printcorr)
                    fprintf(fid,'\n%s correl:  ',' % ');
                    for r = 1:length(printcorr)
                        fprintf(fid,'\t %1.2f, ',printcorr(1,r));
                    end;          
                end;
            end;
        else 
            [x y] = find(nxcorr < negsig);        
            gdcompsa{nx} = gdcomps{nx}(x');
            gdcompsb{nx} =  gdcomps{nx}(y');            
            for r = 1: length(gdcompsa{nx})
                printcorr(1,r) = nxcorr(gdcompsa{nx}(r),gdcompsb{nx}(r));     
            end;   
        end;       
        clustcorr{clust} = nxcorr;
        clustpossig{clust} = possig;
        clustnegsig{clust} = negsig;
        fprintf(fid,'\n%sgdcompsa = {prebasea,awea,frustrationa,joya,angera,happya,sada,lovea,feara,compassiona,jealousya,contenta,griefa,reliefa,disgusta,excitea,postbasea};',printname{nx});       
        fprintf(fid,'\n%sgdcompsb = {prebaseb,aweb,frustrationb,joyb,angerb,happyb,sadb,loveb,fearb,compassionb,jealousyb,contentb,griefb,reliefb,disgustb,exciteb,postbaseb};',printname{nx});       
        fprintf('\none more subject done: %s of %s',int2str(nx),int2str(length(gdcomps)));
    end;   % for clust  
    allcorr{nx} = clustcorr;
    allpossig{nx} = clustpossig;
    allnegsig{nx} = clustnegsig;
fclose(fid);
end;  % for nx


%  Find mean and variance for each component in each template
emos = {'prebase','awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
nclust = 10;clear allsubjmean allsubjvar
for nx = 1:length(gdcomps) 
    clear emomean emovar allclustmean allclustvars 
    wts = floatread(['/data/common2/emotion/',subjwts{nx}],[20 129]);
    sph = floatread(['/data/common2/emotion/',subjsph{nx}],[129 129]);
    alldat = floatread(['/data/common2/emotion/',subjdat{nx}],[129 datsize(nx)]);
    ws = wts*sph;winv = pinv(ws);activations = ws*alldat;  % alldat includes pre and post base
    ntrials = numtrials{nx};
    for clust = 1:nclust
        numtot = 0;
        for cp = 1:length(gdcomps{nx})
            for em = 1:length(emos)
                emotrials = ntrials(em);                
                meancomp = mean(activations(clust,numtot+1:numtot+emotrials));
                %if cp == 5 & em == 11
                %    figure;
                %    hist(activations(clust,numtot+1:numtot+emotrials),100);
                %end;            
                onevar = var(activations(clust,numtot+1:numtot+emotrials));
                allmeans(em) = meancomp;
                allvar(em) = onevar;
                numtot = emotrials+numtot;
            end;             
            emomean{gdcomps{nx}(cp)} = allmeans;  clear allmeans
            emovar{gdcomps{nx}(cp)} = allvar;  clear allvar
        end;            
        allclustmean{clust} = emomean;
        allclustvars{clust} = emovar;
    end;
    allsubjmean{nx} = allclustmean;
    allsubjvar{nx} = allclustvars;
end;

% Plot results
toi = {[1,2,3,4],[1,2,7],[1,2,3,4,8],[1,2,3,6],[1,2,5],[1,3,6,7],[1,3,7,8],[1,2,5,6,7],[1,3,6,7]};
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
cols(1,:) = [.5 .5 .5];
cols(2:16,:) = jet(15);
cols(17,:) = [.5 .5 .5];
%cps = [1,4,15,16,17,19,21,24,26,29];
for nx = 1:length(gdcomps)
cps = gdcomps{nx};
    for clust = 1:length(toi{nx})
        figure;
        for cc = 1:length(cps)
            subplot(round(sqrt(length(cps)))+1,round(sqrt(length(cps))),cc)
            cp = cps(cc);
            for em = 2:length(emos)
                plotmean = allsubjmean{nx}{toi{nx}(clust)}{cp}(emoset(em));
                %plotmean = allsubjvar{nx}{toi{nx}(clust)}{cp}(emoset(em));
                ph = bar(em-1,plotmean); hold on;
                set(ph,'facecolor',cols(em,:));
                ph = text(em-1,0,emos{emoset(em)});
                set(ph,'rotation',90);
            end;
            set(gca,'xlim',[0 17]);
            set(gca,'fontsize',14);
            set(gca,'xticklabel',[]);
            title(int2str(cp));     
        end;
        ph = textsc(['Mean of activations for each component; Subject ',int2str(nx),'; Cluster ',int2str(clust)],'title');
        set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);        
    end;
end;

    
% plot 3d mean activations vs each other
nx= 2;  clust1 = 1; clust2 = 7; clust3 = 2; 
cols(1,:) = [.5 .5 .5];
cols(2:16,:) = jet(15);
cols(17,:) = [.5 .5 .5];

for cc = 1:20
    cp = gdcomps{nx}(cc);
cp = 26;
    figure;
for em = 2:length(emos)
plotmean1 = allsubjmean{nx}{clust1}{cp}(emoset(em));
plotmean2 = allsubjmean{nx}{clust2}{cp}(emoset(em));
plotmean3 = allsubjmean{nx}{clust3}{cp}(emoset(em));
ph = plot3(plotmean1,plotmean2,plotmean3,'k.');hold on;
set(ph,'markersize',25);
set(ph,'color',cols(em,:));
ph = text(plotmean1+.01,plotmean2,plotmean3,emos{emoset(em)});
set(ph,'fontsize',14);
end;
set(gca,'fontsize',14);
xlabel(['Template ',int2str(clust1),' Mean Activation']);
ylabel(['Template ',int2str(clust2),' Mean Activation']);
zlabel(['Template ',int2str(clust3),' Mean Activation']);
ph = textsc(['Mean Activations of Template ',int2str(clust1),',  ',int2str(clust2),' and ',int2str(clust3),' of Component ',int2str(cp)],'title');
set(ph,'fontsize',14);
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');

end;


% plot 2d mean activations vs each other
nx= 2;  clust1 = 1; clust2 = 7; 
cols(1,:) = [.5 .5 .5];
cols(2:16,:) = jet(15);
cols(17,:) = [.5 .5 .5];
for cc = 1:20
    cp = gdcomps{nx}(cc);
cp = 7;
    figure;
for em = 2:length(emos)
plotmean1 = allsubjmean{nx}{clust1}{cp}(emoset(em));
plotmean2 = allsubjmean{nx}{clust2}{cp}(emoset(em));
ph = plot(plotmean1,plotmean2,'k.');hold on;
set(ph,'markersize',25);
set(ph,'color',cols(em,:));
ph = text(plotmean1+.01,plotmean2,emos{emoset(em)});
set(ph,'fontsize',14);
end;
set(gca,'fontsize',14);
xlabel(['Template ',int2str(clust1),' Mean Activation']);
ylabel(['Template ',int2str(clust2),' Mean Activation']);
ph = textsc(['Mean Activations of Template ',int2str(clust1),',  ',int2str(clust2),' of Component ',int2str(cp)],'title');
set(ph,'fontsize',14);
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');

end;

% Plot activations of three components from different templates vs each other
cp1 = 7;     clust1 = 2;
cp2 = 8;    clust2 = 2;
cp3 = 15;    clust3 = 2;
    figure;
for em = 2:length(emos)
plotmean1 = allsubjmean{nx}{clust1}{cp1}(emoset(em));
plotmean2 = allsubjmean{nx}{clust2}{cp2}(emoset(em));
plotmean3 = allsubjmean{nx}{clust3}{cp3}(emoset(em));
ph = plot3(plotmean1,plotmean2,plotmean3,'k.');hold on;
set(ph,'markersize',25);
set(ph,'color',cols(em,:));
ph = text(plotmean1+.01,plotmean2,plotmean3,emos{emoset(em)});
set(ph,'fontsize',14);
end;
set(gca,'fontsize',14);
xlabel(['Template ',int2str(clust1),', Component ',int2str(cp1)]);
ylabel(['Template ',int2str(clust2),', Component ',int2str(cp2)]);
zlabel(['Template ',int2str(clust3),', Component ',int2str(cp3)]);
ph = textsc(['Mean Activations of Templates ',int2str(clust1),',  ',int2str(clust2),' and ',int2str(clust3)],'title');
set(ph,'fontsize',14);
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');


% run ica on activation means
% put subj means into one matrix
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
nx = 4; ica2mat = zeros(length(emos),0);pl = 1;
for clust = 1:20
    for cp = 1:length(gdcomps{nx})
        for em = 1:length(emos)
            onemean = allsubjmean{nx}{clust}{gdcomps{nx}(cp)}(emoset(em));
            ica2mat(em,pl) = onemean;             
        end;
        %ica2mat(:,pl) = ica2mat(:,pl) - mean(ica2mat(:,pl) );
        pl = pl+1;
    end;
end;

 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(ica2mat,'stop',1e-8);
 ws = weights*sphere;
 winv = pinv(ws);
 
 figure;
 for newclust = 1:size(winv,2)
     subplot(5,4,newclust)
    ph= plot(winv(:,newclust));hold on;
    set(ph,'linewidth',2);
    set(gca,'xlim',[1 16]);
    set(gca,'xgrid','on');
    %plot([get(gca,'xlim')],[0 0],'k-');
    set(gca,'ticklength',[.02 .02]); 
    set(gca,'fontsize',14);
    title(int2str(newclust));
end;axcopy 

ph= textsc(['Templates of Spectral Differences from Pre-baseline in Subject: ',int2str(nx)],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
 
figure;
plot(activations(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot activations for each template for each subject

% Plot results
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
cols(1,:) = [.5 .5 .5];
cols(2:16,:) = jet(15);
cols(17,:) = [.5 .5 .5];
%cps = [1,4,15,16,17,19,21,24,26,29];
for nx = 1:length(gdcomps)
cps = gdcomps{nx};
    clear emomean emovar allclustmean allclustvars 
    wts = floatread(['/data/common2/emotion/',subjwts{nx}],[20 129]);
    sph = floatread(['/data/common2/emotion/',subjsph{nx}],[129 129]);
    alldat = floatread(['/data/common2/emotion/',subjdat{nx}],[129 datsize(nx)]);
    ws = wts*sph;winv = pinv(ws);activations = ws*alldat;  % alldat includes pre and post base
    ntrials = numtrials{nx};
figure;
    for clust = 1:8
        numtot = 0;
        for cc = 1:length(cps)
            subplot(8,1,clust)
            cp = cps(cc);
            for em = 2:length(emos)
                emotrials = ntrials(em);                
                onecomp = activations(clust,numtot+1:numtot+emotrials);
                ph = plot(onecomp,'k-'); hold on;
                set(ph,'color',cols(em,:));                
                numtot = emotrials+numtot;
            end;
            title(int2str(clust));     
        end;
    end;
end;


