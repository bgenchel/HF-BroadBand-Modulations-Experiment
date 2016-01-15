% looks at the magnitude of power differences in different emotions (using weights from ICA)
eeglab
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
subj24 = [];  % kc66
subj25 = [9,11,13,14,16,18,19,21,23,24];  % ts79
subj26 = [1,4,6,7,8,11,12,14,16:19];  % es76

gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21, subj22, subj23, subj24, subj25, subj26};
%************   Frontal/temporal Componenents ONLY  *******************************
subj1 = [3,11,13,15,18,19,21,22,26,27,28,29,30,34];  % tl81  23,35 =pulse; 
subj2 = [2,6,9,12,13,15,19,26];% mi83  14,23 = pulse; 31,34,46,51:muscle
subj3 = [5,12,13,14,22];   % ms82  11,37 = pulse; 32 = muscle?
subj4 = [9,13,16,17,19,25];  % js75  7 is pulse
subj5 = [11,14,16,18,19,24,28,34]; % kw78

%  instructed to press only at peak of emotion
subj6 = [1,6,8,11,14,15,16,17,20];% jo82 
subj7 = [1,6,7,8];% kl80 
subj8 = [1,9,11,12];% ar81 
subj9 = [8,9,11,13,16,18,20,29,36]; % eb79  But no button presses actually recorded
subj10 = [2,6,16,17,18,21,23,24,27]; % dg75  
subj11 = [12,14,21,22,23,29,34];  % an82  29 maybe
subj12 = [2,8,12,14,15,20,25];  % jw84

% NO button press sessions
subj13 = [3,4,6,13,15,17,19,20,23,27,30];  % tv81  11 = pulse?? 
subj14 = [1,2,4,9,11,12,14,17,19,20,28];  % sr81   22,27 = pulse  19 maybe
subj15 = [1,11,13,20,38];  % an70  9,12 = pulse
subj16 = [11,15,16,18,25];  % sg75  14 = pulse
subj17 = [2,11,12,3,14,16,19];% mr72
subj18 = [2,10,11,16,17,18,20];% dk74
subj19 = [8,11,18,20,23,24,26,27,28,33];% dn86  
subj20 = [1,3,10,13,14,15,16,18,19,21];% mr71
subj21 = [2,6,15,20,32];  % md85

%  instructed to press only at peak of emotion
subj22 = [3,10,15,17,18,19,20,21,28];  % mr72, second session
subj23 = [1,9,12,18]; %cj82   19 is 24%, 20 is 19%,22 outside brain
subj24 = [];  % kc66
subj25 = [9,14,18,19,21,23,24];  % ts79
subj26 = [1,7,11,14,15,18,19];  % es76
gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21, subj22, subj23, subj24, subj25, subj26};

 
datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-243.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set','emo-1-246.set','emo-1-248.set','emo-1-240.set','emo-1-254.set','emo-1-246.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/','/mr72-2/','/cj82/','/kc66/','/ts79/','/es76/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph243pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph','sph246pc100.sph','sph248pc100.sph','sph240pc100.sph','sph254pc100.sph','sph246pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','wts243pc100.wts','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts','wts246pc100.wts','wts248pc100.wts','wts240pc100.wts','wts254pc100.wts','wts246pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[243 243],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250],[246 246],[248 248],[240 240],[254 254],[246 246]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 243],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250],[100 246],[100 248],[100 240],[100 254],[100 246]};
subjspecs = {'Subj1trials.fdt','Subj2trials.fdt','Subj3trials.fdt','Subj4trials.fdt','Subj5trials.fdt','Subj6trials.fdt','Subj7trials.fdt','Subj8trials.fdt','Subj9trials.fdt','Subj10trials.fdt','Subj11trials.fdt','Subj12trials.fdt','Subj13trials.fdt','Subj14trials.fdt','Subj15trials.fdt','Subj16trials.fdt','Subj17trials.fdt','Subj18trials.fdt','Subj19trials.fdt','Subj20trials.fdt','Subj21trials.fdt','Subj22trials.fdt','Subj23trials.fdt'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};

% analyze output of ICA
sphfile = {'Subj1pc10.sph','Subj2pc10.sph','Subj3pc10.sph','Subj4pc10.sph','Subj5pc10.sph','Subj6pc10.sph','Subj7pc10.sph','Subj8pc10.sph','Subj9pc10.sph','Subj10pc10.sph','Subj11pc10.sph','Subj12pc10.sph','Subj13pc10.sph','Subj14pc10.sph','Subj15pc10.sph','Subj16pc10.sph','Subj17pc10.sph','Subj18pc10.sph','Subj19pc10.sph','Subj20pc10.sph','Subj21pc10.sph','Subj22pc10.sph','Subj23pc10.sph','Subj24pc10.sph'};
wtsfile = {'Subj1pc10.wts','Subj2pc10.wts','Subj3pc10.wts','Subj4pc10.wts','Subj5pc10.wts','Subj6pc10.wts','Subj7pc10.wts','Subj8pc10.wts','Subj9pc10.wts','Subj10pc10.wts','Subj11pc10.wts','Subj12pc10.wts','Subj13pc10.wts','Subj14pc10.wts','Subj15pc10.wts','Subj16pc10.wts','Subj17pc10.wts','Subj18pc10.wts','Subj19pc10.wts','Subj20pc10.wts','Subj21pc10.wts','Subj22pc10.wts','Subj23pc10.wts','Subj24pc10.wts'};

emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
load /data/common2/emotion/subjdims.mat % subjdims  numtrials comment
freqs = [1:.5:50]; 

%  for tl81 ONLY based on self rating
happyemos = [4,16];% only "active" happy emotions: joy,  excitement
happypass = [6,2,12,14]; % only "passive" happy emotions: awe, content, ,happy (bad at relief)

sademos = [7,13]; % sad and grief

mademos = [3,5]; % ang,frust only  (bad at jealousy)

baseemos = [1];
happyemos = [4];% only "active" happy emotions: joy,  excitement
happypass = [16]; % only "passive" happy emotions: awe, content, ,happy (bad at relief)

sademos = [6]; % sad and grief

mademos = [2]; % ang,frust only  (bad at jealousy)

baseemos = [12];
% unincorporated: compassion, disgust, love, fear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first part uses raw trial spectra for averaging
for nx=1:length(paths);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_loadset( datsets{nx},['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);c=1;w=c+round(sqrt(length(gdcomps{nx})))+1;
    sbsph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    sbwts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
    EEG.icaweights=sbwts;  EEG.icasphere=sbsph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for cc = 1:length(gdcomps{nx}) 
        if c == round(sqrt(length(gdcomps{nx})))+2 | c==(round(sqrt(length(gdcomps{nx})))+2)*2| c==(round(sqrt(length(gdcomps{nx})))+2)*3| c==(round(sqrt(length(gdcomps{nx})))+2)*4| c==(round(sqrt(length(gdcomps{nx})))+2)*5| c==(round(sqrt(length(gdcomps{nx})))+2)*6
            c = c+round(sqrt(length(gdcomps{nx})))+1; w = c+round(sqrt(length(gdcomps{nx})))+1;
        end;        
        sbplot((ceil(sqrt(length(gdcomps{nx})))+2)*2,round(sqrt(length(gdcomps{nx})))+1,c);
        topoplot(EEG.icawinv(:,subjorder(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5);
        sbplot((ceil(sqrt(length(gdcomps{nx})))+2)*2,round(sqrt(length(gdcomps{nx})))+1,w);
        psd(EEG.icaact(gdcomps{nx}(cc),:),512,EEG.srate,512);hold on;w = w+2;
    end;
    ALLEEG=[];EEG=[];
    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;emomap(end) = emomap(end);

    % icamatall is emotions*trials X freqs*comps
    sph=floatread(['/data/common2/emotion/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;
    activations = ws*icamatall;
    winv = pinv(ws);    
    %%%%%%%
    for e = 1:length(numtrials{nx})
        realspec = icamatall(emomap(e):emomap(e+1)-1,:);
        emospecs{e} = realspec;
    end;   
    figure;pl = 1;
%%%%%%%%%%%%%%%%%%% >>  OR << %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same thing but using ersps and a moving average
ls = 1;% linewidth for spectral diffs
for nx = 1:length(gdcomps)
    if iseven(round(sqrt(length(gdcomps{nx}))))
        col = round(sqrt(length(gdcomps{nx})))+2;
    else
        col = round(sqrt(length(gdcomps{nx})))+3;
    end;    
    row = ceil(length(gdcomps{nx})/(col/2))*2;
    if iseven(row)
        row=row;
    else
        row = row+1;
    end;    
    figure;pl = 2;
    EEG = pop_loadset( datsets{nx},['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    sbsph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    sbwts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
    EEG.icaweights=sbwts;  EEG.icasphere=sbsph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
c=1;w=c+col;
    for cc = 1:length(gdcomps{nx}) 
        if c == col+1| c== (col*2)+1| c== (col*3)+1| c== (col*4)+1| c== (col*5)+1| c== (col*6)+1| c== (col*7)+1| c== (col*8)+1| c== (col*9)+1| c== (col*10)+1| c== (col*11)+1| c== (col*12)+1| c== (col*13)+1
            c = c+col; w = w+col;
        end;        
        sbplot(row,col,c);
        topoplot(EEG.icawinv(:,gdcomps{nx}(cc)),EEG.chanlocs,'electrodes','off','plotrad',.5);c=c+2;
        sbplot(row,col,w);
        psd(EEG.icaact(gdcomps{nx}(cc),:),512,EEG.srate,512);hold on;w = w+2;
        set(gca,'xlim',[1 35]); xlabel([]);ylabel([]); set(gca,'fontsize',7);
        title(['Comp ',int2str(gdcomps{nx}(cc))]);
    end;
    ALLEEG=[];EEG=[];

    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    icamatall = zeros(0,length(freqs)*length(gdcomps{nx}));
    for k = 1:length(Alllongersps)
        oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
        for cmp = 1:length(gdcomps{nx})
            [outdata,outx] = movav(oneemo(:,:,gdcomps{nx}(cmp)),[1:size(oneemo,2)],40);
            trialmat = zeros(0,size(outdata,1));
            for trl = 1:5:size(outdata,2)-11
                onetrl = mean(outdata(:,trl:trl+10),2);
                trialmat(end+1,:) = onetrl';
            end;
            if cmp == 1
                icamat = zeros(size(trialmat,1),0);
            end;        
            icamat(:,end+1:end+size(trialmat,2)) = trialmat;
        end;
        icamatall(end+1:end+size(icamat,1),:) = icamat;
        fprintf('\n One More Emotion Done: %i',k);
        subjtrials(k) = size(trialmat,1);    
    end;
    %floatwrite(icamatall, ['/data/common2/emotion/',subjspecs{nx}]); 
    subjdims{nx} = size(icamatall);
    numtrials{nx} = subjtrials;
    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    % icamatall is emotions*trials X freqs*comps
    for e = 1:length(numtrials{nx})
        realspec = icamatall(emomap(e):emomap(e+1)-1,:);
        emospecs{e} = realspec;
    end;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cmp = 1:length(gdcomps{nx})
        happyspecs = zeros(0,length(freqs));sadspecs = zeros(0,length(freqs));madspecs = zeros(0,length(freqs));
        basespecs = zeros(0,length(freqs));happypassspecs = zeros(0,length(freqs));
        for em = 1:length(happyemos)
            happyspecs(end+1:end+size(emospecs{happyemos(em)},1),:) = emospecs{happyemos(em)}(:,length(freqs)*(cmp-1)+1:length(freqs)*(cmp));
        end;
        for em = 1:length(happypass)
            happypassspecs(end+1:end+size(emospecs{happypass(em)},1),:) = emospecs{happypass(em)}(:,length(freqs)*(cmp-1)+1:length(freqs)*(cmp));
        end;
        for em = 1:length(sademos)
            sadspecs(end+1:end+size(emospecs{sademos(em)},1),:) = emospecs{sademos(em)}(:,length(freqs)*(cmp-1)+1:length(freqs)*(cmp));
        end;
        for em = 1:length(mademos)
            madspecs(end+1:end+size(emospecs{mademos(em)},1),:) = emospecs{mademos(em)}(:,length(freqs)*(cmp-1)+1:length(freqs)*(cmp));
        end;
        for em = 1:length(baseemos)
            basespecs(end+1:end+size(emospecs{baseemos(em)},1),:) = emospecs{baseemos(em)}(:,length(freqs)*(cmp-1)+1:length(freqs)*(cmp));
        end;
        for fr = 1:length(freqs)
            happystd(fr) = std(happyspecs(:,fr))/sqrt(size(happyspecs,1));
            happypassstd(fr) = std(happypassspecs(:,fr))/sqrt(size(happypassspecs,1));
            sadstd(fr) = std(sadspecs(:,fr))/sqrt(size(sadspecs,1));
            madstd(fr) = std(madspecs(:,fr))/sqrt(size(madspecs,1));
            basestd(fr) = std(basespecs(:,fr))/sqrt(size(basespecs,1));
        end;    
        happyspecs = mean(happyspecs,1);
        happypassspecs = mean(happypassspecs,1);
        sadspecs = mean(sadspecs,1);
        madspecs = mean(madspecs,1);
        basespecs = mean(basespecs,1);
        sbplot(ceil(row/2),col,pl);        
        ph = plot(freqs,happyspecs,'r-');hold on;pl = pl+2;
        set(ph,'linewidth',ls);
        ph = plot(freqs,happypassspecs,'m-');hold on;
        set(ph,'linewidth',ls);
        ph = plot(freqs,sadspecs,'b-');hold on;
        set(ph,'linewidth',ls);
        ph = plot(freqs,madspecs,'k-');hold on;
        set(ph,'linewidth',ls);
        ph = plot(freqs,basespecs,'k-');hold on;
        set(ph,'color',[.5 .5 .5]);     set(ph,'linewidth',ls-.5);
        set(gca,'xlim',[1 50]);    set(gca,'fontsize',7);    
        for ff = 1:length(freqs)/2
            xf = ff*2; % mult fact must be the same as /by in line above
            ph=plot([freqs(xf)  freqs(xf)],[happyspecs(xf)-happystd(xf) happyspecs(xf)+happystd(xf)]);
            set(ph,'color','r');            set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[happyspecs(xf)+happystd(xf) happyspecs(xf)+happystd(xf)]);
            set(ph,'color','r');             set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[happyspecs(xf)-happystd(xf) happyspecs(xf)-happystd(xf)]);
            set(ph,'color','r');             set(ph,'linewidth',ls);              
            ph=plot([freqs(xf)  freqs(xf)],[happypassspecs(xf)-happypassstd(xf) happypassspecs(xf)+happypassstd(xf)]);
            set(ph,'color','m');            set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[happypassspecs(xf)+happypassstd(xf) happypassspecs(xf)+happypassstd(xf)]);
            set(ph,'color','m');             set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[happypassspecs(xf)-happypassstd(xf) happypassspecs(xf)-happypassstd(xf)]);
            set(ph,'color','m');             set(ph,'linewidth',ls);              
            ph=plot([freqs(xf)  freqs(xf)],[sadspecs(xf)-sadstd(xf) sadspecs(xf)+sadstd(xf)]);
            set(ph,'color','b');            set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[sadspecs(xf)+sadstd(xf) sadspecs(xf)+sadstd(xf)]);
            set(ph,'color','b');             set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[sadspecs(xf)-sadstd(xf) sadspecs(xf)-sadstd(xf)]);
            set(ph,'color','b');           set(ph,'linewidth',ls);    
            ph=plot([freqs(xf)  freqs(xf)],[madspecs(xf)-madstd(xf) madspecs(xf)+madstd(xf)]);
            set(ph,'color','k');           set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[madspecs(xf)+madstd(xf) madspecs(xf)+madstd(xf)]);
            set(ph,'color','k');             set(ph,'linewidth',ls);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[madspecs(xf)-madstd(xf) madspecs(xf)-madstd(xf)]);
            set(ph,'color','k');             set(ph,'linewidth',ls);    
            ph=plot([freqs(xf)  freqs(xf)],[basespecs(xf)-basestd(xf) basespecs(xf)+basestd(xf)]);
            set(ph,'color',[.5 .5 .5]);            set(ph,'linewidth',ls-.5);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[basespecs(xf)+basestd(xf) basespecs(xf)+basestd(xf)]);
            set(ph,'color',[.5 .5 .5]);             set(ph,'linewidth',ls-.5);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[basespecs(xf)-basestd(xf) basespecs(xf)-basestd(xf)]);
            set(ph,'color',[.5 .5 .5]);             set(ph,'linewidth',ls-.5);    
        end;
        fprintf('\n One More Componnent Done: %i',cmp);
        fmax = [happyspecs happypassspecs sadspecs madspecs basespecs];
        axmax = max(fmax); axmin = min(fmax); set(gca,'ylim',[axmin-(abs(axmin*.05)) axmax+(abs(axmax*.05))]);
        set(gca,'xtick',[10:10:50]);set(gca,'xgrid','on');    
    end;% to cmp
    axcopy
    ph =textsc(['Subject ',int2str(nx),': Averaged Spectra among ACTIVE-Happy(red), PASSIVE-Happy(mag), Sad(blue) and Mad(black) Emotions; prebase(grey)'],'title');set(ph,'fontsize',14);
    set(ph,'fontsize',14); set(gcf,'PaperOrientation','landscape'); set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    clear Alllongersps  sbwts sbsph icamatall icamat oneemo outdata emospecs  realspec
    fprintf('\n One More SUBJECT Done: %i (of %i)\n',nx,length(gdcomps));

end; % to nx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  look at imagesc of emos with prebase baseline. Why the offsets? Esp with prebase??
figure; pl = 1;
for nx = 9:9%length(gdcomps)
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    %for k = 1:length(Alllongersps)
        k = 2;
figure; pl = 1;
    oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
        for cmp = 1:8
        subplot(4,2,pl)
        [outdata,outx] = movav(oneemo(:,:,gdcomps{nx}(cmp)),[1:size(oneemo,2)],100);
        %mnsp = mean(outdata,2);
        %outdata = rmbase(outdata);
        imagesc(outdata,[-10 10]);pl = pl+1;
        %imagesc(oneemo(:,:,gdcomps{nx}(cmp)),[-10 10]);pl = pl+1;
        end;
    %end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search for emotions with greatest variance from baseline, then group        
for nx = 25:length(gdcomps)  
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    base=  Alllongersps{1};
    clear allcmpdiffs plotdiffs allcmpH allemodiffs allcmpP
    for cmp = 1:length(gdcomps{nx})
        clear allemodiffs
        %figure; 
        for k = 2:length(Alllongersps)
            oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
            basecmp = base(:,:,gdcomps{nx}(cmp));
            avgbase = mean(basecmp,2);
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            avgcmp = mean(onecmp,2);
            for ff = 1:size(onecmp,1)
                [H(ff),P(ff),CI,STATS] = ttest2(onecmp(ff,:),basecmp(ff,:),.01);
                stddiffs(ff) = sqrt(std(basecmp(ff,:)).^2 - std(onecmp(ff,:)).^2);
            end;
            allH(k-1,:) = H;  allP(k-1,:) = P;
            diffs = avgbase - avgcmp; diffs = diffs';  % 1 X freqs
            allemodiffs(k-1,:) = diffs;
            %ph = plot(freqs,diffs,'k-'); hold on;
            %set(ph,'color',cols(k,:));     set(ph,'linewidth',1);
            %set(gca,'xlim',[1 50]);    set(gca,'fontsize',7);              
    fprintf('\n emotion %i Done',k);
        end;
        allcmpdiffs(:,:,cmp) = allemodiffs; % makes a Emo X Freqs X Cmp
        allcmpH(:,:,cmp) = allH; allcmpP(:,:,cmp) = allP; 
    end;

    plotdiffs = allcmpdiffs;
    plotdiffs(allcmpH==0) = 0;

    % look at distance between all emos at each freq
    keepdiffs = plotdiffs;     
    plotdiffs(16,:,:) = []; minf = 1; mf = 29; fqs = 2; % >fqs of max or min diff freqs needed to be included
    clear cmppairs savepairs allplotidx plotleg
    for cmp = 1:length(gdcomps{nx});
        mxidxs = zeros(1,length(freqs));mnidxs= zeros(1,length(freqs));
        for ff = minf:mf% length(freqs) % 29 = 15 Hz; 39 = 20 Hz; 49 = 25 Hz; 59 = 30 hz; 69=35Hz;79=40Hz;89=45Hz;
            [mxdiff mxidx] = max(plotdiffs(:,ff,cmp));
            %if mxdiff ~= 0
            mxidxs(ff) = mxidx;
            %end;
            [mndiff mnidx] = min(plotdiffs(:,ff,cmp));
            %if mndiff ~= 0
            mnidxs(ff) = mnidx;
            %end;
            if mxdiff == 0 & mndiff == 0
                mnidxs(ff)=0; mxidxs(ff) = 0;
            end;        
            savepairs{ff} = [int2str(freqs(ff)),'; ',emos{mxidxs(ff)+1},'; ' emos{mnidxs(ff)+1}];
        end;
        pl = 1; plotidxs =[]; savepairs=[];
        for w = 1:size(plotdiffs,1)
            if find(mxidxs == w) > fqs
                plotidxs(pl) = w;pl = pl+1;  % ensures more than two freq diff
            end;
        end;
        if ~isempty(plotidxs)
        plotleg{cmp} = emos(plotidxs+1); allplotidx{cmp} = plotidxs; 
        cmppairs {cmp} = savepairs;
        end;        
    end;
    fprintf('\n Comp Search for Subject: %i',nx);
    %PLOT
    emo2 ={  'anger'   'frustration'    'jealousy'    'fear'   'disgust'    'grief'  'sad'    'compassion'    'love'    'relief'     'content'   'awe'   'happy'    'joy'       'excite' }; % same as emoorder

    figure;cols = jet(15);pl = 2; c=1; % for one subject
    if round(sqrt(length(gdcomps{nx})))*2 == 9  | round(sqrt(length(gdcomps{nx})))*2 == 12| round(sqrt(length(gdcomps{nx})))*2 == 15
        col = round(sqrt(length(gdcomps{nx})))*2;
    elseif round(sqrt(length(gdcomps{nx})))*2 == 8  | round(sqrt(length(gdcomps{nx})))*2 == 11| round(sqrt(length(gdcomps{nx})))*2 == 14
        col = round(sqrt(length(gdcomps{nx})))*2+1;
    elseif round(sqrt(length(gdcomps{nx})))*2 == 7 | round(sqrt(length(gdcomps{nx})))*2 == 10| round(sqrt(length(gdcomps{nx})))*2 == 13
        col = round(sqrt(length(gdcomps{nx})))*2+2;        
    end;    
    row = ceil(length(gdcomps{nx})/(col/3));
    EEG = pop_loadset( 'sources.set',['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    for cc = 1:length(gdcomps{nx}) 
        sbplot(row,col,c);
        topoplot(EEG.icawinv(:,gdcomps{nx}(cc)),EEG.chanlocs,'electrodes','off','plotrad',.5);c=c+3;
        title(['Comp ',int2str(gdcomps{nx}(cc))]);
    end;
    ALLEEG=[];EEG=[];

    emoorder = [4,2,10,8,14,12,6,9,7,13,11,1,5,3,15,16]; 
    for cmp = 1:length(gdcomps{nx})
        subplot(row,col,pl:pl+1)
        for r = 1:length(allplotidx{cmp})
            if ~isempty(allplotidx{cmp}(r))
            ph = plot(freqs(1:mf),plotdiffs(allplotidx{cmp}(r),1:mf,cmp)); hold on;
            set(ph,'linewidth',2);
            set(ph,'color',cols(emoorder==allplotidx{cmp}(r),:));set(gca,'fontsize',14);
            set(gca,'box','off'); set(gca,'xgrid','on');set(gca,'xlim',[1 freqs(mf)]);
            %title(['Comp ',int2str(gdcomps{nx}(cmp))]);
            end;
        end;pl = pl+3;
        %legend(plotleg{cmp},3);
    end;
    axcopy
    ph=textsc(['Subject ',int2str(nx),'; Emotions chosen by variance between ',int2str(freqs(minf)),' Hz and ',int2str(freqs(mf)),' Hz; > ',int2str(fqs),' timepnts'],'title');
    set(ph,'fontsize',12);
    fprintf('\n One More Subject Done: %i',nx);
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
print  -dpsc2 -Pcoloring 
close
end;

% for legend
figure;
for r = 1:15
    ph = plot(freqs,plotdiffs(r,:,cmp)); hold on;
    set(ph,'color',cols(r,:));set(gca,'fontsize',14);
    set(gca,'box','off');
end;
ph=legend(emo2,3);set(ph,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'emoorder' based on following order of matrix:
             'awe'    'frustration'    'joy'    'anger'    'happy'    'sad'

  Columns 7 through 12

    'love'    'fear'    'compassion'    'jealousy'    'content'    'grief'

  Columns 13 through 16

    'relief'    'disgust'    'excite'    'postbase'

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plot only specific clusters of components
gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21, subj22, subj23, subj24, subj25, subj26};
clust = 8;
minf = 1; mf = 29; fqs = 8; % >fqs of max or min diff freqs needed to be included
% 29 = 15 Hz; 39 = 20 Hz; 49 = 25 Hz; 59 = 30 hz; 69=35Hz;79=40Hz;89=45Hz;
ls = 2;% linewidth for spectral diffs
totalnum = length(gdcomps{1})+length(gdcomps{2})+length(gdcomps{3})+length(gdcomps{4})+length(gdcomps{5})+length(gdcomps{6})+length(gdcomps{7})+length(gdcomps{8})+length(gdcomps{9})+length(gdcomps{10})+length(gdcomps{11})+length(gdcomps{12})+length(gdcomps{13})+length(gdcomps{14})+length(gdcomps{15})+length(gdcomps{16})+length(gdcomps{17})+length(gdcomps{18})+length(gdcomps{19})+length(gdcomps{20})+length(gdcomps{21})+length(gdcomps{22})+length(gdcomps{23})+length(gdcomps{24})+length(gdcomps{25})+length(gdcomps{26});
figure;cols = jet(15);
    if round(sqrt(totalnum))*2 == 9  | round(sqrt(totalnum))*2 == 12| round(sqrt(totalnum))*2 == 15
        col = round(sqrt(totalnum))*2;
    elseif round(sqrt(totalnum))*2 == 8  | round(sqrt(totalnum))*2 == 11| round(sqrt(totalnum))*2 == 14
        col = round(sqrt(length(gdcomps{nx})))*2+1;
    elseif round(sqrt(totalnum))*2 == 7 | round(sqrt(totalnum))*2 == 10| round(sqrt(totalnum))*2 == 13
        col = round(sqrt(totalnum))*2+2;        
    end;    
    row = ceil(totalnum/(col/3));
    
c=1;pl = 2;
for nx = 1:length(gdcomps)
if ~isempty(gdcomps{nx})
    EEG = pop_loadset( 'sources.set',['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cc = 1:length(gdcomps{nx}) 
        sbplot(row,col,c);
        topoplot(EEG.icawinv(:,gdcomps{nx}(cc)),EEG.chanlocs,'electrodes','off','plotrad',.5);c=c+3;
        title(['Comp ',int2str(gdcomps{nx}(cc))]);
    end;
    ALLEEG=[];EEG=[];

    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    base=  Alllongersps{1};
    clear allcmpdiffs plotdiffs allcmpH allemodiffs allcmpP
    for cmp = 1:length(gdcomps{nx})
        clear allemodiffs
        for k = 2:length(Alllongersps)
            oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
            basecmp = base(:,:,gdcomps{nx}(cmp));
            avgbase = mean(basecmp,2);
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            avgcmp = mean(onecmp,2);
            for ff = 1:size(onecmp,1)
                [H(ff),P(ff),CI,STATS] = ttest2(onecmp(ff,:),basecmp(ff,:),.01);
                stddiffs(ff) = sqrt(std(basecmp(ff,:)).^2 - std(onecmp(ff,:)).^2);
            end;
            allH(k-1,:) = H;  allP(k-1,:) = P;
            diffs = avgbase - avgcmp; diffs = diffs';  % 1 X freqs
            allemodiffs(k-1,:) = diffs;
            fprintf('\n emotion %i Done',k);
        end;
        allcmpdiffs(:,:,cmp) = allemodiffs; % makes a Emo X Freqs X Cmp
        allcmpH(:,:,cmp) = allH; allcmpP(:,:,cmp) = allP; 
    end;

    plotdiffs = allcmpdiffs;
    plotdiffs(allcmpH==0) = 0;

    % look at distance between all emos at each freq
    keepdiffs = plotdiffs;     
    plotdiffs(16,:,:) = []; 
    clear cmppairs savepairs allplotidx plotleg
    for cmp = 1:length(gdcomps{nx});
        mxidxs = zeros(1,length(freqs));mnidxs= zeros(1,length(freqs));
        for ff = minf:mf% length(freqs) 
            [mxdiff mxidx] = max(plotdiffs(:,ff,cmp));
            mxidxs(ff) = mxidx;
            [mndiff mnidx] = min(plotdiffs(:,ff,cmp));
            mnidxs(ff) = mnidx;
            if mxdiff == 0 & mndiff == 0
                mnidxs(ff)=0; mxidxs(ff) = 0;
            end;        
            savepairs{ff} = [int2str(freqs(ff)),'; ',emos{mxidxs(ff)+1},'; ' emos{mnidxs(ff)+1}];
        end;
        p = 1; plotidxs =[]; savepairs=[];
        for w = 1:size(plotdiffs,1)
            if find(mxidxs == w) > fqs
                plotidxs(p) = w;p = p+1;  % ensures more than two freq diff
            end;
        end;
        if ~isempty(plotidxs)
        plotleg{cmp} = emos(plotidxs+1); allplotidx{cmp} = plotidxs; 
        cmppairs {cmp} = savepairs;
        end;        
    end;
    fprintf('\n Comp Search for Subject: %i',nx);
    %PLOT
    emo2 ={  'anger'   'frustration'    'jealousy'    'fear'   'disgust'    'grief'  'sad'    'compassion'    'love'    'relief'     'content'   'awe'   'happy'    'joy'       'excite' }; % same as emoorder

    emoorder = [4,2,10,8,14,12,6,9,7,13,11,1,5,3,15,16]; 
    for cmp = 1:length(gdcomps{nx})
        subplot(row,col,pl:pl+1)
        for r = 1:length(allplotidx{cmp})
            if ~isempty(allplotidx{cmp}(r))
            ph = plot(freqs(1:mf),plotdiffs(allplotidx{cmp}(r),1:mf,cmp)); hold on;
            set(ph,'linewidth',2);
            set(ph,'color',cols(emoorder==allplotidx{cmp}(r),:));set(gca,'fontsize',14);
            set(gca,'box','off'); set(gca,'xgrid','on');set(gca,'xlim',[1 freqs(mf)]);
            title(['S-',int2str(nx),'; C-',int2str(gdcomps{nx}(cmp))]);
            end;
        end;pl = pl+3;
    end;
    fprintf('\n One More Subject Done: %i',nx);
end;
end;
ph=textsc(['Cluster # ',int2str(clust),'; Emotions chosen by variance between ',int2str(freqs(minf)),' Hz and ',int2str(freqs(mf)),' Hz; > ',int2str(fqs),' timepnts'],'title');
set(ph,'fontsize',12);
axcopy
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    
print  -dpsc2 -Pcoloring 
