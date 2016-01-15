% finds mean power before and after button press begins for each emotion

eeglab
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % jo74 emot
%%%****  get pressevents first from EmotionPresses256.m  (copied below)
ALLEEG=[];
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bemabarrass','bcontent','bgrief','brelief'};
 EEG = pop_loadset( 'ButtonOnly.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

for k = 1:length(emos)
    if k == 14
        endt = 250;
    else
        endt=300;
    end;    
    EEG = pop_epoch( EEG, {  emos{k}  }, [-.5  endt], 'newname', emos{k}, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname', emobuts{k});
    EEG.data = rmbase(EEG.data,EEG.pnts,1:512); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
end;
 ALLEEG(1)=[];
%figure;plot(EEG.data); %to make sure you got the right thing
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    
    for ev = 1:length(EEG.event)
        if EEG.event(ev).Event_Type(1) == 'S'
            seltime = EEG.event(ev).latency;
        end;
    end;
    EEG = pop_select( EEG, 'point',[1 seltime] );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, k);
end;
%figure;plot(EEG.data); %to make sure you got the right thing
clear pressevents
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    maxdata = max(EEG.data);
    thresh = maxdata*.1;  % sets threshold at 10% of max value
    r=1; 
    for g = 1: size(EEG.data,2)
        if EEG.data(1,g)>thresh & EEG.data(1,g-1)<thresh
            EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
            EEG.event(end).type = 'press';
            EEG.event(end).Event_Type = 'Response';
        end;    
    end;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
pressevents{k} = EEG.event(end-(r-2):end);
end;
ALLEEG=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
gdcomps = {subj1};
paths = {'/data/common1/emotion/jo74/'};
datset = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};
lotime = -500;
hitime = 500;
maxfreq = 30;
wsize = 256;
tlim = [-1500 1500];
prat = 2;
wave = [3 .5];
pl = 1;
nx=1;n=1;
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
for em = 1:length(datset)
    EEG = pop_loadset( datset{em},paths{nx}); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.event = pressevents{em};  % get from EmotionPresses256.m        
    ft = EEG.event(1).latency;    % this part makes fake epochs to have equal num of epochs for each emotion
    for evtm = 384:256:ft-512  % go up by 1.5 sec to create 3 sec epochs til 2 sec before button
        EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 'fake';        
    end;
    EEG = pop_epoch( EEG,{'fake'} , [-1.5 1.5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1500 1500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,0,subj1 ,5,5,0,1);        
    EEG = pop_jointprob(EEG,0,subj1 ,5,5,0,1);
    for cmp = 1:length(subj1)
        for tr = 1:size(EEG.icaact,3)
            [Pxx,F] = psd(EEG.icaact(subj1(cmp),:,tr),256,EEG.srate,256);% much smoother
            cspec(subj1(cmp),:,tr) = log10(Pxx(:,1))';
        end;        
    end;
    allemopwrb4{1,em} = cspec;
    ALLEEG=[];EEG=[];
    EEG = pop_loadset( datset{em},paths{nx}); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.event = pressevents{em};  % get from EmotionPresses256.m        
    ft = EEG.event(1).latency;    % this part makes fake epochs to have equal num of epochs for each emotion
    for evtm = ft+384:256:ft+38400  % go up by 1.5 sec to create 3 sec epochs
        EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 'fake';        
    end;
    EEG = pop_epoch( EEG,{'fake'} , [-1.5 1.5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1500 1500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,0,subj1 ,5,5,0,1);        
    EEG = pop_jointprob(EEG,0,subj1 ,5,5,0,1);
    for cmp = 1:length(subj1)
        for tr = 1:size(EEG.icaact,3)
            [Pxx,F] = psd(EEG.icaact(subj1(cmp),:,tr),256,EEG.srate,256);% much smoother
            cspec(subj1(cmp),:,tr) = log10(Pxx(:,1))';
        end;        
    end;
    allemopwr{1,em} = cspec;
    ALLEEG=[];EEG=[];
end;
freqs = F;

% plot results
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
% page per component
pl=1; %for error bar spacing
for cp = 1:length(subj1)
    figure;
    for em = 1:length(datset)
        b4oneemo = allemopwrb4{em};
        oneemo = allemopwr{em};
        for f = 1:size(oneemo,2)
            sdemo(subj1(cp),f) = std(oneemo(subj1(cp),f,:))/sqrt(size(oneemo,3));
        end;
        for f = 1:size(oneemo,2)
            sdemob4(subj1(cp),f) = std(b4oneemo(subj1(cp),f,:))/sqrt(size(b4oneemo,3));
        end;
        subplot(4,4,em)
        ph=plot(freqs,mean(oneemo(subj1(cp),:,:),3),'r');hold on;
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf)]);
            set(ph,'color','r');
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf)]);
            set(ph,'color','r'); 
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf)]);
            set(ph,'color','r'); 
            set(ph,'linewidth',1);    
        end;
        ph=plot(freqs,mean(b4oneemo(subj1(cp),:,:),3),'b');hold on;
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf)]);
            set(ph,'color','b');
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf)]);
            set(ph,'color','b'); 
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf)]);
            set(ph,'color','b'); 
            set(ph,'linewidth',1);    
        end;
        set(gca,'xlim',[2 45]);
        set(gca,'fontsize',16);
        set(gca,'box','off');
        title(emos{em}); 
        fprintf('One More Emo done: %s',emos{em});
    end;
textsc(['Component: ',int2str(subj1(cp))],'title');
end;


%One page per emotion
figure;
pl=1; %for error bar spacing
for em = 1:length(datset)
    figure;
    for cp = 1:length(subj1)
        b4oneemo = allemopwrb4{em};
        oneemo = allemopwr{em};
        subplot(4,4,cp)
         ph=plot(freqs,mean(oneemo(subj1(cp),:,:),3),'r');hold on;
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf)]);
            set(ph,'color','r');
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)+sdemo(subj1(cp),xf)]);
            set(ph,'color','r'); 
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf) mean(oneemo(subj1(cp),xf,:),3)-sdemo(subj1(cp),xf)]);
            set(ph,'color','r'); 
            set(ph,'linewidth',1);    
        end;
        ph=plot(freqs,mean(b4oneemo(subj1(cp),:,:),3),'b');hold on;
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf)]);
            set(ph,'color','b');
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)+sdemob4(subj1(cp),xf)]);
            set(ph,'color','b'); 
            set(ph,'linewidth',1);    
            ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf) mean(b4oneemo(subj1(cp),xf,:),3)-sdemob4(subj1(cp),xf)]);
            set(ph,'color','b'); 
            set(ph,'linewidth',1);    
        end;
    set(gca,'xlim',[2 45]);
    set(gca,'fontsize',16);
    set(gca,'box','off');
          title(int2str(subj1(cp))); 
    end;
textsc(['Emotion: ',emos{em}],'title');
end;

