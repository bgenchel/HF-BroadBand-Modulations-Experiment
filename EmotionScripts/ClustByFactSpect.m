%Clusters intra-subject template spectra across subj by kmeans then looks for emo distributions

% early button press sessions (not instructed to only press on surge)
subj1 = [1,3:22,25:30,34,36];  % tl81  23,35 =pulse; 
subj2 = [1,2,4:9,12,13,15:19,21,22,24,26,29];% mi83  14,23 = pulse; 31,34,46,51:muscle
subj3 = [5,6,7,10,12,13,14,15,17,18,19,22];   % ms82  11,37 = pulse; 32 = muscle?
subj4 = [1:5,8:10,12:14,16,17,19,22,25];  % js75  7 is pulse
subj5 = [1:6,9:12,14:19,22:25,28,34]; % kw78
%  instructed to press only at peak of emotion
subj6 = [1:3,5,6,8:11,14:17,20,23];% jo82 
subj7 = [1,2,5:10,15,22];% kl80 
subj8 = [1,2,6,8,9,11,12];% ar81 
subj9 = [1,2,3,5:11,13,14,16:20,22:24,29,30,32,36,37,42]; % eb79  But no button presses actually recorded
subj10 = [2,3,5:8,10:12,16:18,21,23,24,27,28]; % dg75  
subj11 = [1,2,4:18,21:23,29,32,34];  % an82  29 maybe
subj12 = [2,5,7,8,12,14,15,17,20,21,25,26,29];  % jw84
% NO button press sessions
subj13 = [1:6,8,9,10,13:15,17,19,20,23,25,27,30];  % tv81  11 = pulse?? 
subj14 = [1:4,6,8:14,16:21,23,24,28];  % sr81   22,27 = pulse  19 maybe
subj15 = [1:7,11,13:16,18,20,21,23,25,37,38,40];  % an70  9,12 = pulse
subj16 = [2:5,7:11,15:18,25];  % sg75  14 = pulse
subj17 = [1:3,6:8,10:14,16,19];% mr72
subj18 = [2,3,5,7:12,14,16:20];% dk74
subj19 = [4:8,10:15,17:24,26:28,30,33];% dn86  
subj20 = [1,3,6,7,10:16,18,19,21];% mr71
subj21 = [1:6,8:12,14:24,32];  % md85
%  instructed to press only at peak of emotion
subj22 = [1:3,6,7,9,10,11,13:15,17:25,28,29];  % mr72, second session
subj23 = [1:3,7:9,11:15,18:20,22,25]; %cj82   19 is 24%, 20 is 19%,22 outside brain
subj24 = [2:4,6:13,15,17,21,24];  % kc66
subj25 = [9,11,13,14,16,18,19,21,23,24];  % ts79
subj26 = [1,4,6,7,8,11,12,14:20];  % es76
gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21, subj22, subj23, subj24, subj25, subj26};
paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/','/mr72-2/','/cj82/','/kc66/','/ts79/','/es76/'};

datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-243.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set','emo-1-246.set','emo-1-248.set','emo-1-240.set','emo-1-254.set','emo-1-246.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/','/mr72-2/','/cj82/','/kc66/','/ts79/','/es76/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph243pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph','sph246pc100.sph','sph248pc100.sph','sph240pc100.sph','sph254pc100.sph','sph246pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','wts243pc100.wts','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts','wts246pc100.wts','wts248pc100.wts','wts240pc100.wts','wts254pc100.wts','wts246pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[243 243],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250],[246 246],[248 248],[240 240],[254 254],[246 246]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 243],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250],[100 246],[100 248],[100 240],[100 254],[100 246]};
 
subjspecs = {'Subj1trials.fdt','Subj2trials.fdt','Subj3trials.fdt','Subj4trials.fdt','Subj5trials.fdt','Subj6trials.fdt','Subj7trials.fdt','Subj8trials.fdt','Subj9trials.fdt','Subj10trials.fdt','Subj11trials.fdt','Subj12trials.fdt','Subj13trials.fdt','Subj14trials.fdt','Subj15trials.fdt','Subj16trials.fdt','Subj17trials.fdt','Subj18trials.fdt','Subj19trials.fdt','Subj20trials.fdt','Subj21trials.fdt','Subj22trials.fdt','Subj23trials.fdt','Subj24trials.fdt','Subj25trials.fdt','Subj26trials.fdt'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};

%sphfile = {'Subj1pc10.sph','Subj2pc10.sph','Subj3pc10.sph','Subj4pc10.sph','Subj5pc10.sph','Subj6pc10.sph','Subj7pc10.sph','Subj8pc10.sph','Subj9pc10.sph','Subj10pc10.sph','Subj11pc10.sph','Subj12pc10.sph','Subj13pc10.sph','Subj14pc10.sph','Subj15pc10.sph','Subj16pc10.sph','Subj17pc10.sph','Subj18pc10.sph','Subj19pc10.sph','Subj20pc10.sph','Subj21pc10.sph','Subj22pc10.sph','Subj23pc10.sph','Subj24pc10.sph','Subj25pc10.sph','Subj26pc10.sph'};
%wtsfile = {'Subj1pc10.wts','Subj2pc10.wts','Subj3pc10.wts','Subj4pc10.wts','Subj5pc10.wts','Subj6pc10.wts','Subj7pc10.wts','Subj8pc10.wts','Subj9pc10.wts','Subj10pc10.wts','Subj11pc10.wts','Subj12pc10.wts','Subj13pc10.wts','Subj14pc10.wts','Subj15pc10.wts','Subj16pc10.wts','Subj17pc10.wts','Subj18pc10.wts','Subj19pc10.wts','Subj20pc10.wts','Subj21pc10.wts','Subj22pc10.wts','Subj23pc10.wts','Subj24pc10.wts','Subj25pc10.wts','Subj26pc10.wts'};
sphfile = {'Subj1pc15.sph','Subj2pc15.sph','Subj3pc15.sph','Subj4pc15.sph','Subj5pc15.sph','Subj6pc15.sph','Subj7pc15.sph','Subj8pc15.sph','Subj9pc15.sph','Subj10pc15.sph','Subj11pc15.sph','Subj12pc15.sph','Subj13pc15.sph','Subj14pc15.sph','Subj15pc15.sph','Subj16pc15.sph','Subj17pc15.sph','Subj18pc15.sph','Subj19pc15.sph','Subj20pc15.sph','Subj21pc15.sph','Subj22pc15.sph','Subj23pc15.sph','Subj24pc15.sph','Subj25pc15.sph','Subj26pc15.sph'};
wtsfile = {'Subj1pc15.wts','Subj2pc15.wts','Subj3pc15.wts','Subj4pc15.wts','Subj5pc15.wts','Subj6pc15.wts','Subj7pc15.wts','Subj8pc15.wts','Subj9pc15.wts','Subj10pc15.wts','Subj11pc15.wts','Subj12pc15.wts','Subj13pc15.wts','Subj14pc15.wts','Subj15pc15.wts','Subj16pc15.wts','Subj17pc15.wts','Subj18pc15.wts','Subj19pc15.wts','Subj20pc15.wts','Subj21pc15.wts','Subj22pc15.wts','Subj23pc15.wts','Subj24pc15.wts','Subj25pc15.wts','Subj26pc15.wts'};
load /data/common2/emotion/clusters/subjdims.mat % subjdims  numtrials comment
freqs = [1:.5:50]; 
% the following is subject subjective emotion scores with zeros where people couldn't get into it
tl81=[8,8,8,8,8,8,8,8,8,8,0,8,0,8,8];
mi83=[7,8,7,7,8,8,8,8,0,7,8,6,7,0,5];
ms82=[5,6,4,4,9,5,6,0,5,5,7,8,5,6,7];
js75=[6,7,6,7,7,9,8,8,9,7,0,9,6,7,7]; % zero out jealousy
%js75=[6,7,6,7,7,9,8,8,9,7,5,9,6,7,7];
kw78=[8.5,9,7,6,6,9,9,7,8,5,0,4,0,4,9];
jo82=[8,8,7,7,6,8,8,9,5,9,9,5,9,8,9];
kl80=[8,8,8,9,9,9,8,9,9,9,9,7,8,9,9];
ar81=[7,7,7,8,8,8,8,8,6,6,8,8,8,8,8];
eb79=[9,9,9,9,9,9,9,9,9,9,7,9,9,9,9];
dg75=[6,9,4,7,4,9,8,5,3,4,3,3,4,4,4];
an82=[7,7,9,7,8,7,7,8,7,9,7,9,8,5,7];
jw84=[7,6,8,6,8,6,7,8,8,7,5,8,6,8,7];
tv81=[5,7,8,7,9,8,0,6,0,8,8,5,7,8,6];
sr81=[7,8,8,9,6,7,8,7,7,8,5,8,6,7,7];
an70=[8,7,7,7,7,7,7,8,8,7,7,6,7,7,7];
sg75=[6,5,7,7,8,5,8,6,0,7,5,6,6,8,7];
mr72=[8,7,5,5,6,7,8,8,2,8,0,7,7,6,6];
dk74=[7,7,9,9,5,7,8,7,6,6,0,5,7,7,8];
dn86=[7,6,8,8,6,7,8,7,9,8,9,5,7,8,9];
mr71=[7,0,6,8,6,5,6,8,6,6,0,7,5,8,7];
md85=[6,7,8,8,8,7,6,8,9,6,6,9,7,8,7];
mr722=[6,4,4,0,5,6,5,7,4,6,0,6,4,4,6];
cj82=[7,9,7,6,9,8,5,9,9,5,6,8,6,7,8];
kc66=[8,6,8,8,0,7,7,8,6,6,0,0,7,8,7];
ts79 =[5,5,7,8,9,6,6,6,6,0,6,5,5,6,0];
es76 =[6,7,8,7,8,7,7,10,8,7,0,7,7,8,6]; % substituted some values for blank ratings

allsubj = {tl81,mi83,ms82,js75,kw78,jo82,kl80,ar81,eb79,dg75,an82,jw84,tv81,sr81,an70,sg75,mr72,dk74,dn86,mr71,md85,mr722,cj82,kc66,ts79,es76};

% reorder emo scores
cats = {'awe','joy','frustration','anger','sad','happy','content','love','fear','compassion','jealousy','grief','relief','excite','disgust'}; % this is the order of the scores
% pick ONE of the following :   **********************8
%realorder = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'}; % this is the order of data (spectra)
realorder = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % this is the reordered order in terms of valence/activity
%   *********************************
for w = 1:length(allsubj)
    for q = 1:length(allsubj{w})
        newsubj{w}(q) = allsubj{w}(find(ismember(cats,realorder{q})));
    end;
end;% I checked and this works

for w = 1:length(allsubj)
    newsubj{w} = newsubj{w}/max(newsubj{w}); % divides all scores by max score to make internally consistent
end;
%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);
clear newact
    for tp = 1:size(winv,2)
        for cp = 1:length(gdcomps{nx})
            newact(tp,cp,:) = activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp);
        end;      
    end;
    for tp = 1:size(winv,2)
        for ff = 1:length(freqs)
            tmpact = newact(tp,:,ff);
            tvar(tp,ff) = var(tmpact);
        end;
    end;
    % to Normalize by lowest variance freqencies for each template
    for tp = 1:size(winv,2)
        mvar = min(tvar(tp,:));
        lowvar{tp} = find(tvar(tp,:) >= mvar & tvar(tp,:) <= mvar+mvar*.4); % find freqs with low var
        tmpcps = newact(tp,:,lowvar{tp}); tmpcps = squeeze(tmpcps); % makes a cp X freqs
        newact(tp,:,:) = newact(tp,:,:)/mean(std(tmpcps'));   % mean of std over all comps        
    end;
    
    % find most variant components by zscore cutoff
    clear tplists
    for tp = 1:size(winv,2)
        cplist = zeros(0);forstd = newact(tp,:,:); forstd = squeeze(forstd);% makes a cp X freqs
        cut = mean(std(forstd'))+std(std(forstd'));  % includes comps from mean to .3 *std
        for cp = 1:length(gdcomps{nx})
            if std(forstd(cp,:)) > cut
                cplist(end+1) = gdcomps{nx}(cp);
            end;
        end;
        tplists{tp} = cplist;
    end;
    nxlists{nx} = tplists;
    
    fprintf('\n One More SUBJECT Done: %i',nx);
end;
eeglab
%%%%%%% DO ABOVE ONLY ONCE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE ONE OF THE FOLLOWING (SELECTS SUBJECTS OF INTEREST)
button = [2:12,21:26]; % all button presses, early and 'only when you feel it' subjects
button = [13:20]; % no button press (apart from the first one)
button = [2,4:6,8:12,14,17:21,23,25,26]; % all 'good' subjects (ones that said they got into it)
button = [2:21,23:26];  % not mr72-2 or tl81 (weird spectra)
button = [1:26];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; % bad to good with passive inbetween
pl = 1;  alltpspecs = zeros(0,length(freqs));  
tt = 1;clear clustspecs factmap tpemoscores dipolelist
for nx = 1:length(button)
    subjmap(nx) = pl;
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{button(nx)}],[15 subjdims{button(nx)}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);pp=1; clear clustspecs
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{button(nx)}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    emomap = ones(1,1); %%%%% for finding emoscores
    for e = 2:length(numtrials{button(nx)})+1
        emomap(1,e) = emomap(e-1) + numtrials{button(nx)}(e-1); % marks where each emotion STARTS
    end;
    for tp = 1:size(winv,2)
        coi = nxlists{button(nx)}{tp};
        factmap(tp) = length(coi)*2; % times 2 because of both orientations
        tpwts =  winv(:,tp)'; 
        for cp = 1:length(coi)
            rcp = find(coi(cp) == gdcomps{button(nx)});
            % Specs are collected twice (once is *-1)
            clustspecs(pp,:) = activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp); 
       
            for e = 2:length(numtrials{button(nx)})-1  % start with 2 for straight nums (not diffs)
                e=emoorder(e);
                tempmat = tpwts(emomap(e):emomap(e+1)-1); 
                if ~isempty(tempmat(find(tempmat>0)))
               tpemoscores(e-1,pl) = sum(tempmat(find(tempmat>0)))/length(tempmat);%normalize by # trials
                else
                    tpemoscores(e-1,pl) = 0;
                end;                
            end; ;pp=pp+1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            trackscores(pl,:) = [nx,tp];
            tpemoscores(find(newsubj{nx}==0),pl) = 0; % zero out 'bad' emotions
            if max(tpemoscores(:,pl)) > 0
            tpemoscores(:,pl) = tpemoscores(:,pl)/max(tpemoscores(:,pl)); % normalize to max score
            end;
            dipolelist(pl) =  EEG.dipfit.model(nxlists{button(nx)}{tp}(cp));pl = pl+1; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clustspecs(pp,:) = activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp)*-1;
            for e = 2:length(numtrials{button(nx)})-1  % start with 2 for straight nums (not diffs)
                e=emoorder(e);
                tempmat = tpwts(emomap(e):emomap(e+1)-1); 
                if ~isempty(tempmat(find(tempmat<0)))
                    tpemoscores(e-1,pl) =  abs(sum(tempmat(find(tempmat<0))))/length(tempmat);
                else
                    tpemoscores(e-1,pl) = 0;
                end;                
            end; pp=pp+1;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            trackscores(pl,:) = [nx,tp];
            tpemoscores(find(newsubj{nx}==0),pl) = 0; % zero out 'bad' emotions
            if max(abs(tpemoscores(:,pl))) > 0
            tpemoscores(:,pl) = tpemoscores(:,pl)/max(abs(tpemoscores(:,pl))); % normalize to max score
            end;
            dipolelist(pl) =  EEG.dipfit.model(nxlists{button(nx)}{tp}(cp));pl = pl+1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end;
        %tpspecs(tp,:) = mean(clustspecs,1);
    end;
    stpmap{nx} = factmap;
    alltpspecs(end+1:end+size(clustspecs,1),:) = clustspecs;
    %alltpspecs(end+1:end+size(tpspecs,1),:) = tpspecs;
    fprintf('\n One More SUBJECT Done: %i',nx);
    ALLEEG=[];EEG=[];
end;
subjmap(nx+1) = pl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = find(freqs<51);alltpspecs = alltpspecs(:,fr);
%reps = 2;[nclust,centr,clst,Cg] =Kmeangrp(alltpspecs,50,reps,1); % find optimal cluster number, don't plot
nclust = 37;% All subjects bu #1 and 22
nclust = 35;% for "only good subjects"
[kout, C,sumd, allsums] = kmeans(alltpspecs,nclust,'replicates',5);
% zero out outliers (> 3std)
newpole = dipolelist;
for q = 1:size(allsums,2)
    fout = allsums(:,q)/std(allsums(:,q));
    thisclust = find(kout==q);
    for w = 1:length(thisclust)
        if ismember(thisclust(w),find(fout>3))
            kout(thisclust(w))=0;
            newpole(thisclust(w))=[];
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Optional plotting of all clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for cl = 1:nclust
    oneclust = find(kout == cl);
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),cl)
    x=mean(alltpspecs(oneclust,:),1); plot(freqs,x);hold on;        
    for w = 1:length(oneclust)
        plot(freqs(fr),alltpspecs(oneclust(w),:)');hold on;
    end;
    ph = plot(freqs(fr),x,'r'); set(ph,'linewidth',2);
    set(gca,'xlim',[1 50]);
    title(['Cluster ',int2str(cl)]);set(gca,'xgrid','on');
end;axcopy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next find emo contributions, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear  collectscores collectscores tmpscores tmpscore  % get mnemodiff from PCAforEmos.m nx is rel to 'button' list
for cl = 1:nclust
    oneclust = find(kout == cl);
    collectscores = tpemoscores(:,oneclust)';
    subjtp = trackscores(oneclust,:);
    for w = 1:length(subjtp(:,1))
        onesubj = subjtp(w,1);
        if length(find(onesubj(:,1) == subjtp(:,1))) > 1
            likesubj = find(onesubj == subjtp(:,1));
            for q = 2:length(likesubj)
                if subjtp(likesubj(q),2) == subjtp(likesubj(1:q-1),2)
                    collectscores(likesubj(q),:) = 0;
                end;
            end;
        end;
    end; % all above ensures that each template only counted once for each cluster   
    for e = 1:size(collectscores,2)
        tmpscore = collectscores(:,e);
        tmpscore = tmpscore(find(tmpscore~=0));
        tmpscores(e) = mean(tmpscore);
    end;    
    clustscores(cl,:) = tmpscores;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot mean clusters and emo bars together
emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'};
cols = jet(15);
figure;pp = 1; 
%row = round(sqrt(nclust*2)); col = ceil(sqrt(nclust*2));
row = round(sqrt(nclust*2))+1; col = floor(sqrt(nclust*2));
for cl = 1:nclust
    oneclust = find(kout == cl);
    subplot(row,col,pp)
    x=mean(alltpspecs(oneclust,:),1); plot(freqs(fr),x);hold on;   pp = pp+1;      
    set(gca,'xlim',[1 50]);set(gca,'xtick',[5:5:50]);set(gca,'xgrid','on'); title(['Cluster ',int2str(cl)]);
    set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40 [] 50});
    subplot(row,col,pp)
    for e = 1:size(clustscores,2)
        ph =bar(e,clustscores(cl,e));hold on;
        set(ph,'facecolor',cols(e,:));
        ph=text(e,0,emo2{e}); set(ph,'fontsize',7);
        set(ph,'rotation',90);title(['Cluster ',int2str(cl)]);set(gca,'box','off');
    end;  set(gca,'xlim',[0 16]); set(gca,'ylim',[0 max(clustscores(cl,:))+max(clustscores(cl,:))*.01]);
    set(gca,'xticklabel',[]);     pp = pp+1;    
end;axcopy
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    print  -dpsc2 -Pcoloring 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run pca on emoscores  TO 3<<<<
 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(clustscores','pca',3);
 ws = weights*sphere; winv = pinv(ws);
cols = jet(15);
 figure; subplot(2,1,1)
 for e = 1:size(winv,1)
 ph =plot3(winv(e,1),winv(e,2),winv(e,3),'.');hold on;set(gca,'fontsize',14);
 set(ph,'markersize',20);
 set(ph,'color',cols(e,:));
 ph = text(winv(e,1),winv(e,2),winv(e,3),emo2{e});
 set(ph,'color',cols(e,:)); 
 pl =plot3([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)],[-.04  winv(e,3)]);
 set(pl,'color',cols(e,:))             
 end;
 set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             xlabel('Winv 1'); ylabel('Winv 2'); zlabel('Winv 3'); subplot(2,1,2)
 ph = plot(activations(1,:),'b');hold on;set(ph,'linewidth',2);set(gca,'fontsize',14);
 ph = plot(activations(2,:),'g');set(ph,'linewidth',2);
 ph = plot(activations(3,:),'r');set(ph,'linewidth',2);
xlabel('Spectral Cluster Number'); ylabel('ICA Weight');set(gca,'xtick',[1:1:10]);
set(gca,'xticklabel',[1:1:10]); set(gca,'xgrid','on');set(gca,'xlim',[1 10]);
legend('Winv 1','Winv 2','Winv 3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run pca on emoscores  TO ANYTHING ELSE<<<<

[weights,sphere,compvars,bias,signs,lrates,activations] = runica(clustscores','pca',5);
 ws = weights*sphere; winv = pinv(ws);cols = jet(15);
figure;  row = ceil(sqrt(size(winv,2)*2)); col = round(sqrt(size(winv,2)*2)); pp = 1;
for cb = 1:size(winv,2)-2
    for bc = cb+1:size(winv,2)-1
        for cbc = bc+1:size(winv,2)
            subplot(row,col,pp)
            for e = 1:size(winv,1)
                ph =plot3(winv(e,cb),winv(e,bc),winv(e,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e,cb),winv(e,bc),winv(e,cbc),emo2{e});
                set(ph,'color',cols(e,:)); 
            end;
            mx =  max(winv(:,cb));mnx = min(winv(:,cb));
            my =  max(winv(:,bc));mny = min(winv(:,bc));
            mz =  max(winv(:,cbc));mnz = min(winv(:,cbc));
            set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
            zl = get(gca,'zlim');
            for e = 1:size(winv,1)
                pl =plot3([winv(e,cb) winv(e,cb)],[winv(e,bc) winv(e,bc)],[zl(1)  winv(e,cbc)]);
                set(pl,'color',cols(e,:))             
            end;
            set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);
            
            xlabel(['Winv ',int2str(cb)]); ylabel(['Winv ',int2str(bc)]); zlabel(['Winv ',int2str(cbc)]); 
            pp = pp+1;
        end;
    end;
end;axcopy
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster pairs that are opposites
cpairs = {[8,9],[18,1],[18,9],[8,1],[7,13],[14,13],[2,11],[10,12],[3,16],[4,6],[5,17],[15,19]};
row = 3; col = 4; cols(1,:) = [1 0 0]; cols(2,:) = [0 0 1];
barcols = jet(15);
for cl = 1:length(cpairs)
    if cl == 1 | cl == 7 |  cl == 13
figure;pl = 1;
end;
   for pr = 1:length(cpairs{cl})       
    subplot(row,col,pl); 
    oneclust = find(kout == cpairs{cl}(pr));
    x=mean(alltpspecs(oneclust,:),1);  ph = plot(freqs,x);hold on; set(ph,'color',cols(pr,:));
    set(ph,'linewidth',2);    set(gca,'xlim',[1 50]);set(gca,'xtick',[5:5:50]);
    set(gca,'xgrid','on');title(['Clusters ',int2str(cpairs{cl}(1)),' and ',int2str(cpairs{cl}(2))]);
    set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40 [] 50});
    legend(['Cluster ',int2str(cpairs{cl}(1))],['Cluster ',int2str(cpairs{cl}(2))],4);
    for e = 1:size(clustscores,2)
        ecompare(pr,e) = clustscores(cpairs{cl}(pr),e);
    end; 
   end;   pl = pl+1;
    subplot(row,col,pl)
    for e = 1:size(clustscores,2)
        ph =bar(e,ecompare(1,e) - ecompare(2,e));hold on;
        set(ph,'facecolor',barcols(e,:));
    end;    
    title(['Emotion Differences:  Clusters ',int2str(cpairs{cl}(1)),' minus ',int2str(cpairs{cl}(2))]);
    set(gca,'box','off');  set(gca,'xlim',[0 16]); set(gca,'xtick',[]);set(gca,'xticklabel',[]);
    slo = min(ecompare(1,:) - ecompare(2,:));shi = max(ecompare(1,:) - ecompare(2,:));
    if slo < 0 & shi < 0
        set(gca,'ylim',[slo- abs(slo)*.05 0]);
    for e = 1:size(clustscores,2)
        ph=text(e,slo+abs(slo)*.05,emo2{e}); set(ph,'fontsize',12);  set(ph,'rotation',90);  
    end;
    elseif slo> 0 & shi > 0 
        set(gca,'ylim',[0 shi+ abs(shi)*.05]);
    for e = 1:size(clustscores,2)
        ph=text(e,0,emo2{e}); set(ph,'fontsize',12);  set(ph,'rotation',90);  
    end;
    else
        set(gca,'ylim',[slo- abs(slo)*.05 shi+ abs(shi)*.05]);
    for e = 1:size(clustscores,2)
        ph=text(e,slo+abs(slo)*.05,emo2{e}); set(ph,'fontsize',12);  set(ph,'rotation',90);  
    end;
    end;    
    pl = pl+1;  
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
end;
print  -dpsc2 -Pcoloring 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clist = [2,11,7,13,14,10,12,3,16,4,6,5,17,15,19,18,1,8,9];
clist = [1:nclust];
figure;dipcols = jet(nclust);
for cl = 1:nclust
    pc = clist((nclust+1)-cl);
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),(nclust+1)-cl)
    plotdips = newpole(find(kout == pc));pl = 1;clear cutdips % newpole replaced dipolelist
    for w= 1:length(plotdips)-1
        for q = 1:w:length(plotdips)
            if plotdips(w).posxyz(1,:) == plotdips(q).posxyz(1,:)
                cutdips(pl) = w;pl = pl+1;
            end;
        end;
    end;cutdips = unique(cutdips);plotdips(cutdips) = [];
    if ~isempty(plotdips)
    dipplot(plotdips,'image','mri','gui','off','normlen','on','dipolelength',0,'dipolesize',35,'spheres','on','color',{dipcols(cl,:)},'projlines','on','projimg','on'); view(60,20);  camzoom(1.2)
    end;
end;
ph=title(int2str(pc));set(ph,'fontsize',14);set(ph,'color','r');
set(gcf,'color','w');
ph = textsc('Dipole locations for spectral template clusters across subjects; All "good" subjects included','title');set(ph,'fontsize',14);set(ph,'color','r');
set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
print  -dpsc2 -Pcoloring 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do stats on emo wts differences
% Cluster pairs that are opposites
cpairs = {[8,9],[18,1],[18,9],[8,1],[7,13],[14,13],[2,11],[10,12],[3,16],[4,6],[5,17],[15,19]};
for cl = 1:length(cpairs)
    oneclust = find(kout == cpairs{cl}(1));
    collectscores1 = tpemoscores(:,oneclust)';
    oneclust = find(kout == cpairs{cl}(2));
    collectscores2 = tpemoscores(:,oneclust)';
    for e = 1:size(collectscores1,2)
        [H,stats(1,e),CI,st] = ttest2(collectscores1(:,e),collectscores2(:,e),.01,0);
    end;
    allstats(cl,:) = stats;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find sig diffs between avg emotion scores for each cluster
alpha = .05; clear alldiffs allmasks sigtp
for cl = 1:nclust
alls = zeros(0,1);grpvar = zeros(0,1);    oneclust = find(kout == cl);
    collectscores1 = tpemoscores(:,oneclust)';    
    subjtp = trackscores(oneclust,:);
    for w = 1:length(subjtp(:,1))
        onesubj = subjtp(w,1);
        if length(find(onesubj(:,1) == subjtp(:,1))) > 1
            likesubj = find(onesubj == subjtp(:,1));
            for q = 2:length(likesubj)
                if subjtp(likesubj(q),2) == subjtp(likesubj(1:q-1),2)
                    collectscores1(likesubj(q),:) = 0;
                end;
            end;
        end;
    end; % all above ensures that each template only counted once for each cluster   
    for e = 1:size(collectscores1,2)
        cs = collectscores1(:,e);
        % for matlab anova
        cs(find(cs==0)) = NaN;
        collectscores1(:,e) = cs;
        % for Strauss anova:
        cs(find(cs==0))=[];
        alls(end+1:end+length(cs),1) = cs;
        grpvar(end+1:end+length(cs),1) = e;
    end;
     [F,pr,df,ss,ms,varcomp,varprop] = anova(alls,grpvar);
        [P,ANOVATAB,STATS] = anova1(collectscores1,emo2,'off');
[compare,means,H] = multcompare(STATS,alpha,'off');
    for q = 1:size(compare,1)
        diffmat(compare(q,1),compare(q,2)) = compare(q,4);
        if  compare(q,5)>0 & compare(q,3)<0
            mask(compare(q,1),compare(q,2)) = 0;
        elseif compare(q,5)<0 & compare(q,3)>0
            mask(compare(q,1),compare(q,2)) = 0;    
        else        
            mask(compare(q,1),compare(q,2)) = 1;    
        end;    
    end;
    alldiffs{cl} = diffmat;
    allmasks{cl} = mask;
    sigtp(cl) = pr;
end;
figure;
for cl = 1:length(alldiffs)
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),cl)
    plotdiffs = alldiffs{cl};
    plotdiffs(find(allmasks{cl}==0))=0;
    imagesc(plotdiffs,[-.25 .25]);
    for w = 1:size(alldiffs{cl},1)
        ph =text(1,w,emo2{w});set(ph,'fontsize',7);
    end;  
    set(gca,'xtick',[2:2:14]); set(gca,'xgrid','on');title(['Cluster ',int2str(cl)]);
end;axcopy
textsc(['Mulitple Comparisons between emotions; Only SIGNIFICANT voxels (p < ',num2str(alpha),')'],'title');

