% takes single trial output of newcrossf and clusters between subjects

eeglab
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
subj10 = [2:4,6:8,11:13,16,18,22,23,24,27,35]; % dg75  
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

gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21};

 
datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-243.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph243pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','sph243pc100.sph','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[243 243],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 243],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250]};

subjspecs = {'Subj1trials.fdt','Subj2trials.fdt','Subj3trials.fdt','Subj4trials.fdt','Subj5trials.fdt','Subj6trials.fdt','Subj7trials.fdt','Subj8trials.fdt','Subj9trials.fdt','Subj10trials.fdt','Subj11trials.fdt','Subj12trials.fdt','Subj13trials.fdt','Subj14trials.fdt','Subj15trials.fdt','Subj16trials.fdt','Subj17trials.fdt','Subj18trials.fdt','Subj19trials.fdt','Subj20trials.fdt','Subj21trials.fdt'};
emoset = {'prebase.set','awe.set','frustration.set','joy.set','anger.set','happy.set','sad.set','love.set','fear.set' ,'compassion.set','jealousy.set','content.set','grief.set','relief.set','disgust.set','excite.set','postbase.set'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nx = 1:length(gdcomps) 
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 

    % now ready to run psd
    for k = 1:length(emoset)
        EEG = pop_loadset( emoset{k}, ['/data/common2/emotion/',paths{nx}]);
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
        EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),size(EEG.icaact,2)*size(EEG.icaact,3));
        for cmp1 = 1:length(gdcomps{nx})-1
            for cmp2 = cmp1+1:length(gdcomps{nx})
figure; [coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher,alltfX,alltfY] =newcrossf(EEG.icaact(gdcomps{nx}(cmp1),ft:end),EEG.icaact(gdcomps{nx}(cmp2),ft:end),size(EEG.icaact,2)-(ft-1),[0 ((size(EEG.icaact,2)-(ft-1))/256)*1000],EEG.srate,[3 .5],'type','phasecoher','timesout',((size(EEG.icaact,2)-ft)/256)*4,'freqs',[1:.5:50]);
            longcrs(:,:,gdcomps{nx}(cmp2)) =coh ;  % freqsXtimesXcomponent   
            clf
            end;
            allcrs{gdcomps{nx}(cmp1)} = longcrs;
        end;
        Alllongcrs{k} = allcrs;
        ALLEEG=[];EEG=[];
    end;
    save ContDataCRS.mat Alllongcrs freqsout timesout
end;

