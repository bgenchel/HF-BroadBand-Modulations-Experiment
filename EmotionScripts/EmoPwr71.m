pop_topoplot(EEG,0, [37:71] , 'rr80 Emotion Components',[6 6] ,0, 'electrodes', 'off', 'plotrad',0.5, 'masksurf', 'on');

eeglab

datset = {'w-scenes.set','emot1.set'};
subj1 = [7:12,14,17:20,22,24,26];  %ap80
subj2 = [5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,31,42,47];  % rr80
gdcomps = {subj1 subj2};

paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/rr80/'}
emos = {'awe', 'frust','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief'};


for nx = 1:length(paths)
    cd (paths{nx})
    load wts
    load sph
    EEG = pop_loadset( datset{nx},paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    for k = 1:length(emos)
        if k == 13
            endt = (EEG.event(end).latency-EEG.event(end-3).latency)/256;
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
            if EEG.event(ev).type(1) == '2'
                seltime = EEG.event(ev).latency;
                break
            end;
        end;
        EEG = pop_select( EEG, 'point',[1 seltime] );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, k);
    end;
    %figure;plot(EEG.data); %to make sure you got the right thing
    clear pressevents
    for k = 1:length(emos)
        EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
        maxdata = max(EEG.data(size(EEG.data,1),:));
        thresh = maxdata*.1;  % sets threshold at 10% of max value
        r=1; 
        for g = 2: size(EEG.data,2)
            if EEG.data(size(EEG.data,1),g)>thresh & EEG.data(size(EEG.data,1),g-1)<thresh
                EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
                EEG.event(end).type = 'press';
                EEG.event(end).Event_Type = 'Response';
            end;    
        end;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);  
        pressevents{k} = EEG.event(end-(r-2):end);
    end;
    
    for em = 1:length(emos)
        EEG = eeg_retrieve(ALLEEG, em); CURRENTSET = em;
        EEG = pop_select( EEG, 'nochannel',72);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
        % now recalculate ica.act
        EEG.icaweights = wts;
        EEG.icasphere = sph;    EEG.icaact = [];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

        ft = pressevents{em}(1).latency;
        for evtm = ft+500:1500:size(EEG.data,2)  % go up by 1.5 sec to create 3 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';        
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-1.5 1.5]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
        EEG = pop_rmbase( EEG,[-1500 1500]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
        EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
        
        for cmp = 1:length(gdcomps{nx})
            for tr = 1:size(EEG.icaact,3)
                [Pxx,F] = psd(EEG.icaact(gdcomps{nx}(cmp),:,tr),256,EEG.srate,256);% much smoother
                cspec(gdcomps{nx}(cmp),:,tr) = log10(Pxx(:,1))';
            end;
            
        end;
        allemopwr{1,em} = cspec;
    end;
    ALLEEG=[];EEG=[];
    comment = 'psd power spectra in good comps during each emotion (the different cell arrays) in this order: awe, frustration,joy,anger,sad,happy,fear,love ,jealousy,compassion,content,grief,relief; starts at first button press to end of emotin(space bar); 3 sec non-overlapping windows to get "single trials" ';
    freqs = F;
    cd (paths{nx})
    save SnglTrSpecPwr.mat allemopwr freqs comment
end;


% visualize results with error bars
nx=1;
cd (paths{nx})
load SnglTrSpecPwr.mat
           
% plot error manually 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
plotemos = [1:13];
posset = [1,3,6,8,10,11,13];
negset = [2,4,5,7,9,12];
colneg = winter(6);
colpos = autumn(7);
figure;n=1; p=1;  pl = 1; %(for error bar spacing)
for c = 1:length(gdcomps{nx})
    n=1; p=1;
    subplot(5,4,c)
    for em=1:length(plotemos)
        em = plotemos(em);
        oneemo = allemopwr{em};
        for f = 1:size(oneemo,2)
            sdemo(gdcomps{nx}(c),f) = std(oneemo(gdcomps{nx}(c),f,:))/sqrt(size(oneemo,3));
        end;
        ph=plot(freqs,mean(oneemo(gdcomps{nx}(c),:,:),3));hold on;
        if ismember(em,negset)
            set(ph,'color',colneg(n,:));
            currcol = colneg(n,:);n=n+1;
        elseif ismember(em,posset)
            set(ph,'color',colpos(p,:));
            currcol = colpos(p,:);p=p+1;
        else
            set(ph,'color','k');
            currcol = 'k';%p=p+1;            
        end;
        
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(oneemo(gdcomps{nx}(c),xf,:),3)-sdemo(gdcomps{nx}(c),xf) mean(oneemo(gdcomps{nx}(c),xf,:),3)+sdemo(gdcomps{nx}(c),xf)]);
            set(ph,'color',currcol);
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(gdcomps{nx}(c),xf,:),3)+sdemo(gdcomps{nx}(c),xf) mean(oneemo(gdcomps{nx}(c),xf,:),3)+sdemo(gdcomps{nx}(c),xf)]);
            set(ph,'color',currcol); 
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(gdcomps{nx}(c),xf,:),3)-sdemo(gdcomps{nx}(c),xf) mean(oneemo(gdcomps{nx}(c),xf,:),3)-sdemo(gdcomps{nx}(c),xf)]);
            set(ph,'color',currcol); 
            set(ph,'linewidth',1);    
        end;
    end;
    set(gca,'xlim',[2 45]);
    set(gca,'fontsize',16);
    set(gca,'box','off');
    title(int2str(gdcomps{nx}(c)));
end;
legend({'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'});
axcopy

figure; % colorbar
for c = 1:length(posset)
ph=barh(c,1);hold on;
set(ph, 'FaceColor',colpos(c,1:3));
end;
set(gca,'xticklabel',[]);
set(gca,'ytick',[1:7]);
set(gca,'yticklabel', emos(posset));
figure;
for c = 1:length(negset)
ph=barh(c,1);hold on;
set(ph, 'FaceColor',colneg(c,1:3));
end;
set(gca,'xticklabel',[]);
set(gca,'ytick',[1:6]);
set(gca,'yticklabel', emos(negset));

% 
fr = find(freqs >= 8 & freqs <=  9);
emoset = [2,4,5,7,9,12,1,3,6,8,10,11,13];
cols = jet(13);
figure;
for cp = 1:length(gdcomps{nx})
    subplot(5,5,cp); clear sdev plotone
    for em = 1:13
        plotone = allemopwr{em};
        plotone = plotone(gdcomps{nx}(cp),fr,:);
        plotone = mean(plotone,2);
        plotone = squeeze(plotone);
        sdev = std(plotone)/sqrt(length(plotone));
        plotone = mean(plotone);
        ph = bar(em,plotone);hold on;
        set(ph,'facecolor',cols(find(emoset==em),:));
        plot([em em],[plotone-sdev plotone+sdev],'k');        
        ph = text(em,.1,emos{em});
        set(ph,'rotation',90);
    end;
    title(['Comp: ',int2str(gdcomps{nx}(cp))]);
    set(gca,'xlim',[0 15]);
    set(gca,'xticklabel',[]);
end;

