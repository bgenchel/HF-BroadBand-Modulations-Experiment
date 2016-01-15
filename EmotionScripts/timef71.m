% run timef on old 71 channel emotion data to make it congruent with 256 data 

eeglab
subj1 = [7:12,14,17:20,22,24,26];  %ap80
subj2 = [5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,31,42,47];  % rr80
subj2 = [5,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,31,42,47];  % rr80
gdcomps = {subj1 subj2};
paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/rr80/'};
datpaths = {'/data/common1/emotion/ap80/ersps/','/data/common1/emotion/rr80/ersps/'};
baseset = {'w-scenes.set','emot1.set'};
% get baseline spectra
comment = ' Emotion pilot ERSPs.  epoched randomly across all emos (-1.5 1.5s) and cleaned by autorejection (4,4). taken for baseline spectrum for emotion ersps. padratio:, winsize, baseline: whole epoch (-1.5 1.5s)';
figure;
%for nx = 1:length(paths)
nx=2;
cd (paths{nx})
    load wts
    load sph
    EEG = pop_loadset( baseset{nx}, paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    % now recalculate ica.act
    EEG.icaweights = wts;
    EEG.icasphere = sph;    EEG.icaact = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
    for evtm = 500:1500:size(EEG.data,2)  % go up by 1.5 sec to create 3 sec epochs
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

    comp_ersp=zeros(65,200,gdcomps{nx}(end));
    ersp_boot = zeros(2,65,gdcomps{nx}(end));
    baseline = zeros(gdcomps{nx}(end),65);
    for n=1:length(gdcomps{nx})
        [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(gdcomps{nx}(n),:), EEG.pnts,[EEG.xmin*1000 EEG.xmax*1000], EEG.srate,[3 .5], 'winsize',256,'padratio',4,'baseline',1500,'plotitc','off','alpha',.01);
        comp_ersp(:,:,gdcomps{nx}(n))=ersp;
        ersp_boot(:,:,gdcomps{nx}(n)) =  erspboot;
        baseline(gdcomps{nx}(n),:) = powbase;
        clear ersp itc powbase  erspboot itcboot
        clf
    end;
    cd (datpaths{nx})
    save baseline.mat baseline comp_ersp ersp_boot  times freqs comment ;
    ALLEEG=[]; EEG=[];
end;


    
    % Now run timefs using this baseline on randomly epoched emotion periods beginning with button press
nx=2;
paths = {'/data/common/emotion/ap80/','/data/common/emotion/rr80/'};
datpaths = {'/data/common1/emotion/ap80/ersps/','/data/common1/emotion/rr80/ersps/'};
emos = {'awe', 'frust','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
emobuts = {'bawe', 'bfrustration','bjoy','banger','bsad','bhappy','bfear','blove','bjealousy','bcompassion','bcontent','bgrief','brelief'};
datset = {'w-scenes.set','emot1.set'};
EEG = pop_loadset( datset{nx},paths{nx} );
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

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
    maxdata = max(EEG.data(72,:));
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

comment = ' Emotion pilot ERSPs.  epoched on button press within button-pressing period. taken for baseline spectrum for emotion ersps. padratio:, winsize, coher, baseline: whole session (-1.5 1.5s)';
cd (datpaths{nx})
load baseline.mat
cd (paths{nx})
load wts
load sph
figure;
for k = 1:length(emos)
    EEG = eeg_retrieve(ALLEEG, k); CURRENTSET = k;
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'overwrite', 'on');
    EEG.icaweights = wts;
    EEG.icasphere = sph;    EEG.icaact = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
    ft = pressevents{k}; ft = ft(1).latency;
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
 
    comp_ersp=zeros(65,200,gdcomps{nx}(end));
    ersp_boot = zeros(2,65,gdcomps{nx}(end));
    %itc_boot= zeros(65,gdcomps{nx}(end));
    %comp_itc=zeros(65,200,gdcomps{nx}(end));
    for n=1:length(gdcomps{nx})
        [ersp,itc,powbase,times,freqs,erspboot,itcboot]=timef( EEG.icaact(gdcomps{nx}(n),:), EEG.pnts,[EEG.xmin EEG.xmax], 250,[3 .5], 'winsize',256,'padratio',4,'powbase',baseline(gdcomps{nx}(n),:),'alpha',.01,'plotitc','off');
        comp_ersp(:,:,gdcomps{nx}(n))=ersp;
        %comp_itc( :,:,gdcomps{nx}(n))=itc;
        ersp_boot(:,:,gdcomps{nx}(n)) =  erspboot;
        %itc_boot(:,gdcomps{nx}(n)) = itcboot';
        clear ersp itc powbase  erspboot itcboot
        clf
    end;
    erspcell{k} = comp_ersp;
    %itccell{k} = comp_itc;
    ebootcell{k} = ersp_boot;
    %ibootcell{k} = itc_boot;
end;

cd (datpaths{nx})
save RandEpochEmoERSPs.mat baseline erspcell ebootcell  times freqs comment ;


% plot the results

%  Plot each component, all emotions
% load in chanlocs and adjust   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab  % create sources for plotting maps and sources,done for rr80
nx=2;
paths = {'/data/common/emotion/ap80/','/data/common/emotion/rr80/'};
datpaths = {'/data/common1/emotion/ap80/ersps/','/data/common1/emotion/rr80/ersps/'};
emos = {'awe', 'frust','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
datset = {'w-scenes.set','emot1.set'};
cd (paths{nx})
load wts
load sph
EEG = pop_loadset( datset{nx},paths{nx} );
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
EEG = pop_select( EEG, 'nochannel',72);
%EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/ap80/ap80.elp', 'filetype',''}, 'forcelocs',{0, 'X', 'Cz'}, 'convert',{ 'chancenter',[],1});
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'load',{ '/data/common1/emotion/rr80/rr80.elp', 'filetype', 'autodetect'}, 'forcelocs',{0, 'X', 'Cz'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[]);');
    EEG.icaweights = wts;
    EEG.icasphere = sph;    EEG.icaact = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_select( EEG, 'point',[1 512] );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset( EEG,'sources.set' ,paths{nx});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
subj1 = [7:12,14,17:20,22,24,26];  %ap80
subj2 = [5,9,10,11,12,13,14,16,17,18,19,20,21,22,42];  % rr80
gdcomps = {subj1 subj2};
paths = {'/data/common/emotion/ap80/','/data/common/emotion/rr80/'};
datpaths = {'/data/common1/emotion/ap80/ersps/','/data/common1/emotion/rr80/ersps/'};
emos = {'awe', 'frust','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','content','grief','relief'};
datset = {'w-scenes.set','sources.set'};

nx=2;
cd (paths{nx})
load wts
load sph
EEG = pop_loadset( datset{nx},paths{nx} );
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
cd (datpaths{nx})
load RandEpochEmoERSPs.mat
fr = find(freqs<30);

row = length(gdcomps{nx});
col = length(erspcell)+2;
figure;pl=3;tp = 1;  lim = 7;
for cp = 1:length(gdcomps{nx})
    if cp == 18
        figure;pl=3;tp=1;
    end;    
    subplot(row,col,tp)
    topoplot(EEG.icawinv(:,gdcomps{nx}(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5);
    title(int2str(gdcomps{nx}(cp)));tp=tp+1;
    subplot(row,col,tp)
    plot(freqs(fr),baseline(gdcomps{nx}(cp),fr));
    hold on
    set(gca,'xlim',[2 30]);
    set(gca,'yticklabel',[]);
    set(gca,'xgrid','on');
    tp=tp+col-1;   

    for p = 1:length(erspcell)
        tmpemo =  erspcell{p};  % 1 is ersp, 2 is itc
        subplot(row,col,pl)        
        imagesc(times,freqs(fr),tmpemo(fr,:,gdcomps{nx}(cp)),[-lim lim]);
        pl=pl+1;
        set(gca,'ydir','norm');        
        set(gca,'fontsize',5);
        set(gca,'ticklength',[.04 .04]);
        if cp == 1
            ph =  title(emos{p});
            set(ph,'rotation',45);
            set(ph,'fontsize',14);
        end;
        set(gca,'xticklabel',[]);
        
    end;
    pl=pl+2;
end; 
colorbar;
textsc('rr80 Emotion Time/Freq vs Whole Trial Baseline','title');

