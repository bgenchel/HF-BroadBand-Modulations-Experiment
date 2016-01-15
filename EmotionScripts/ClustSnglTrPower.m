% runs ica clustering on indiv subject power spectra vs fixation spectra. 

% input subject lists from /data/common2/emotion/Complist
 
eeglab

%emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones

%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datsets = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
for nx = 29:length(gdcomps) 
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
    for ds = 1:length(datsets)%numsets(nx)
        % now ready to run timef
        %EEG = pop_loadset( ['emo-',int2str(ds),'-',int2str(gdchan(nx)),'.set'], ['/data/common2/emotion/',paths{nx}]);
        EEG = pop_loadset( datsets{ds}, ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG.icaweights=wts; EEG.icasphere=sph;EEG.icawinv=[];
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_saveset(EEG, datsets{ds},['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        ALLEEG=[];EEG=[];
    end;
end;

% use newtimef to find power over long time periods, scaling wavelet decomp 
%done subj :1:2
comment = 'each emotion continuous dataset epoched into non-overlapping trials and noise removed by autorejection. then data reshaped back to continuous data and submitted to newtimef for continuous data; baseline from all emotions, not any baseline. epoched in overlapping 4s epochs with noise removed.emo order:anger,frustration,jealousy,fear ,disgust,grief,sad,compassion,love,relief,content,awe,happy,joy,excite';
figure;
for nx = 35:35%length(gdcomps) 
    cd (['/data/common2/emotion/',paths{nx},'/ersps/']);load EmoBaseline.mat
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 

    % now ready to run timef
    for k = 1:length(emos)
        EEG = pop_loadset( [emos{k},'.set'], ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %if exist('EEG.event(2)')  % for baselines
            ft = EEG.event(2).latency;
        %else
        %    ft = 1;
        %end; 
        for evtm = ft+128:256:size(EEG.data,2)-128  % go up by 1 sec to create non-overlapping 1 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';        
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-.5 .5]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
        EEG = pop_rmbase( EEG,[-500 500]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG.icaweights=wts; EEG.icasphere=sph;EEG.icawinv=[];
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_eegthresh(EEG,0,gdcomps{nx} ,-400,400,EEG.xmin,EEG.xmax,0,1);
        EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
        EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
        runemo = reshape(EEG.icaact,size(EEG.icaact,1),size(EEG.icaact,2)*size(EEG.icaact,3));
        for cmp = 1:length(gdcomps{nx})
           [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(runemo(gdcomps{nx}(cmp),:),size(runemo,2),[0 (size(runemo,2)/256)*1000],EEG.srate,[3 .5],'timesout',(size(runemo,2)/256)*4,'plotitc','off','freqs',[1:.5:50],'powbase',baseline(gdcomps{nx}(cmp),:), 'verbose','off','freqs',[1:.5:50]);
            longersp(:,:,gdcomps{nx}(cmp)) = ersp;  % freqsXtimesXcomponent   
            clf
            fprintf('\nOne More Component Done: %i of %i\n',cmp,length(gdcomps{nx}));  clf  
        end;
        Alllongersps{k} = longersp;clear longersp
        ALLEEG=[];EEG=[];
            fprintf('\nOne More EMOTION Done: %s\n',emos{k});  clf  
    end;
    save ContDataERSPs.mat Alllongersps freqs times
    clear Alllongersps longersp ersp itc powbase times freqs erspboot itcboot
    fprintf('\nOne More Subj Done: %i\n',nx);   
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now pull in data for clustering subj by subj
%load /data/common2/emotion/FrontalComps.mat gdcomps
% need to do frntal: nx = 33 30 29 25 24 
% need to do reg: nx = 24
ICA_LINUX = '/data/common/matlab/ica_linux2.4';
for nx = 1:34%length(gdcomps)
    ICA_SCRIPT = ['/data/common2/emotion/clusters/Subj',int2str(nx),'ClustPwr.sc'];
    cd (['/data/common2/emotion/',paths{nx},'/ersps/']);load ContDataERSPs.mat
    icamatall = zeros(0,length(freqs)*length(gdcomps{nx}));clear subjtrials
    for k = 1:length(Alllongersps)
        oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
        for cmp = 1:length(gdcomps{nx})
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            trialmat = zeros(0,size(oneemo,1));
            if nx == 10 | nx==23 | nx == 33
                increm = 2;
            else increm = 1;
            end; % start not at beginning for nx = 25,31,33 because too long, but not 10,23 long           
            for trl = round(size(oneemo,2)*.05):increm:size(onecmp,2)-round(size(oneemo,2)*.05) % 1,3,5,... avg over 4 (1sec)
                onetrl = mean(onecmp(:,trl:trl),2);
                onetrl = onetrl/sqrt(mean(onetrl.^2));
                                   % ersp = ersp/sqrt(mean(ersp.^2)); % normalize by div by RMS
                %savemean(k,cmp,trl) = mean(onetrl);
                %onetrl = onetrl-mean(onetrl); % subtract the mean
                trialmat(end+1,:) = onetrl';
            end;
            if cmp == 1
                icamat = trialmat;  icamat(:,1:end) =[];
            end;
            icamat(:,end+1:end+size(trialmat,2)) = trialmat;
        end;
        icamatall(end+1:end+size(icamat,1),:) = icamat;
        subjtrials(k) = size(trialmat,1); clear icamat   
        fprintf('\n One More EMOTION Done: %i of %i',k,length(Alllongersps));
    end;
    %SubjSpecMeans{nx} = savemean;
    % write out data in proper format for ICA computations
    if  length(gdcomps{1}) < 20
        floatwrite(icamatall, ['/data/common2/emotion/clusters/',Frontsubjspecs{nx}]);
    else               
        floatwrite(icamatall, ['/data/common2/emotion/clusters/',subjspecs{nx}]);    
    end;
    % generate an ICA script 
    pcadims{nx} = round(sqrt(size(icamatall,2)/8));% div by # = desired ratio frames/dims^2
    fid = fopen(ICA_SCRIPT, 'w');
    if  length(gdcomps{1}) < 20
        fprintf(fid, 'DataFile %s\n', ['/data/common2/emotion/clusters/',Frontsubjspecs{nx}]);
    else
        fprintf(fid, 'DataFile %s\n', ['/data/common2/emotion/clusters/',subjspecs{nx}]);% Front or not%
    end;    
    fprintf(fid, 'chans %d\n', size(icamatall,1));
    fprintf(fid, 'frames %d\n',size(icamatall,2));
    if  length(gdcomps{1}) < 20
        fprintf(fid, 'WeightsOutFile %s\n', ['/data/common2/emotion/clusters/FRSubj',int2str(nx),'pc',int2str(pcadims{nx}),'.wts']);
        fprintf(fid, 'SphereFile %s\n', ['/data/common2/emotion/clusters/FRSubj',int2str(nx),'pc',int2str(pcadims{nx}),'.sph']);  
    else
        fprintf(fid, 'WeightsOutFile %s\n', ['/data/common2/emotion/clusters/Subj',int2str(nx),'pc',int2str(pcadims{nx}),'.wts']); 
        fprintf(fid, 'SphereFile %s\n', ['/data/common2/emotion/clusters/Subj',int2str(nx),'pc',int2str(pcadims{nx}),'.sph']);  
    end;        
    fprintf(fid, 'sphering on\n');
    fprintf(fid, 'bias on\n');
    fprintf(fid, 'extended 1\n');
    fprintf(fid, 'pca %d\n', pcadims{nx}); 
    fprintf(fid, 'lrate 2.0e-4\n');
    fprintf(fid, 'blocksize 0\n');
    fprintf(fid, 'stop 1e-07\n');
    fprintf(fid, 'maxsteps 2000\n');
    fprintf(fid, 'posact on\n');
    fprintf(fid, 'annealstep 0.98\n');
    fprintf(fid, 'annealdeg 60\n');
    fprintf(fid, 'momentum 0\n');
    fprintf(fid, 'verbose on\n');
    fclose(fid);
    % run the ICA program
    run_ica_str = [ICA_LINUX ' < ' ICA_SCRIPT];
    [status, result] = system(run_ica_str);
    %%%%%
    pcadims{nx} = round(sqrt(size(icamatall,2)/8));% div by # = desired ratio frames/dims^2
    subjdims{nx} = size(icamatall);
    numtrials{nx} = subjtrials;   
    if  length(gdcomps{1}) < 20
        wtsfile{nx} = ['FRSubj',int2str(nx),'pc',int2str(pcadims{nx}),'.wts'];
        sphfile{nx} = ['FRSubj',int2str(nx),'pc',int2str(pcadims{nx}),'.sph'];    
    else
        wtsfile{nx} = ['Subj',int2str(nx),'pc',int2str(pcadims{nx}),'.wts'];
        sphfile{nx} = ['Subj',int2str(nx),'pc',int2str(pcadims{nx}),'.sph'];
    end;    
    fprintf('\n One More SUBJECT Done: %i\n',nx);
end;
comment = 'all  components; cells are subjects; dimension is the icamatall matrix with trials X freqs*comps; to be used for linux ica; numtrials tells how many trials, or timepoints, are included for each subj/emotion; NO smoothing across trials. Trial 1:1:end-40 (to eliminate end) number of pca dims adjusted to ~ 8 pnts/wt';
load /data/common2/emotion/clusters/subjdims.mat subjdims numtrials pcadims wtsfile sphfile comment
load /data/common2/emotion/clusters/Frontsubjdims.mat subjdims numtrials pcadims wtsfile sphfile comment gdcomps

%save /data/common2/emotion/clusters/subjdims.mat subjdims numtrials pcadims wtsfile sphfile comment
save /data/common2/emotion/clusters/Frontsubjdims.mat subjdims numtrials pcadims wtsfile sphfile comment gdcomps


% run ica from linux:
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/clusters/ClustPwrICA.sc
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/clusters/ClustPwrICA1.sc
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/clusters/ClustPwrICA2.sc
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/clusters/ClustPwrICA3.sc
ex
%%%%%%%%%%%%%%%%  PLOT ICA RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots a page of spectral templates for each subject
% enter grandgdcomps from DipoleClusters
% CHOOSE ONE OF THE FOLLOWING:
load /data/common2/emotion/clusters/subjdims.mat subjdims numtrials pcadims wtsfile sphfile comment
load /data/common2/emotion/clusters/Frontsubjdims.mat subjdims numtrials pcadims wtsfile sphfile gdcomps comment
%%%%%%%%% 
load /data/common2/emotion/clusters/KmeanClustDips.mat clustcps kout C sumd allsums keeptrack clustorder 
load /data/common2/emotion/clusters/KmeanClustFrontDips.mat clustcps kout C sumd allsums keeptrack clustorder
%load /data/common2/emotion/clusters/KmeanClustDipsSpec.mat clustcps kout C sumd allsums keeptrack clustorder  
load /data/common2/emotion/clusters/allPreBasespecs.mat freqs allspec comment
spfreqs = freqs;fr = find(spfreqs > 2 & spfreqs < 35);
freqs = [1:.5:50]; 
addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
for nx = 1:10%length(gdcomps)
    clear sph wts 
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);
    EEG = pop_loadset('sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    
    figure;row = size(activations,1)+2; col = length(gdcomps{nx})+1;
    subjorder = zeros(0);clear  clustnum
    for ct = 1:length(clustorder)
        rct = clustorder(ct);
        if ~isempty(clustcps{rct}{nx})
             match=ismember(clustcps{rct}{nx},gdcomps{nx});
             if ~isempty(find(match))
                subjorder(end+1:end+length(match)) = clustcps{rct}{nx}(match ) ;
             end;             
        end;
    end;   pl = 2;  
    for cp = 1:length(subjorder)
        subplot(row,col,pl)
        topoplot(EEG.icawinv(:,subjorder(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
        set(gca,'fontsize',7);
        title(int2str(subjorder(cp)));
    end;
    for cp = 1:length(subjorder)
        subplot(row,col,pl+1)
        rcp = find(allgdcomps{nx} == subjorder(cp));
        ph = plot(spfreqs(fr),allspec{nx}(rcp,fr),'k-');hold on;
        set(gca,'fontsize',7);set(gca,'xgrid','off');pl = pl+1;
        set(gca,'xlim',[2 40]);   set(gca,'xtick',[5:5:40]);
        set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40});  
        set(gca,'yticklabel',[]);set(gca,'xgrid','on');
    end;pl = pl+1;
    for tp = 1:row-2
        subplot(row,col,pl)
        hist(winv(:,tp),75);pl = pl+1;hold on;
        set(gca,'fontsize',7);set(gca,'xlim',[-1.5 1.5]);
        plot([0 0],[get(gca,'ylim')],'r-');
        set(gca,'yticklabel',[]);
        for cmp = 1:length(subjorder)
            rcp = find(gdcomps{nx}== subjorder(cmp) );
            subplot(row,col,pl)
            plot(freqs,activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp)); pl = pl+1;hold on;
            plot([0 0],[get(gca,'ylim')],'g');
            set(gca,'fontsize',7);set(gca,'box','off');
            set(gca,'xlim',[0 40]); set(gca,'ylim',[-3.5 10.5]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] [] [] 30 [] []}); 
            set(gca,'xgrid','on');
            %set(gca,'yticklabel',[]);  set(gca,'xticklabel',[]);            
        end;
    end;
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    ph=textsc(['Subject ',int2str(nx),'; Spectral Templates for all comps (pos orientation)'],'title');
    set(ph,'fontsize',14);ALLEEG=[];EEG=[];    
end;axcopy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  plot several factors on same plot (of same component)
load /data/common2/emotion/clusters/subjdims.mat % subjdims  numtrials comment
freqs = [1:.5:50];
%%%%%%%%%%%%%%%%  PLOT ICA RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots a page of spectral templates for each subject
load /data/common2/emotion/clusters/KmeanClustDips.mat clustcps kout C sumd allsums keeptrack clustorder 
%load /data/common2/emotion/clusters/KmeanClustFrontDips.mat clustcps kout C sumd allsums keeptrack clustorder
%load /data/common2/emotion/clusters/KmeanClustDipsSpec.mat clustcps kout C sumd allsums keeptrack clustorder  
load /data/common2/emotion/clusters/allPreBasespecs.mat freqs allspec comment
spfreqs = freqs;fr = find(spfreqs > 2 & spfreqs < 35);
freqs = [1:.5:50]; 
facs = [1:10]; cols = jet(length(facs));
for nx = 1:1%length(gdcomps)
    if round(sqrt(length(gdcomps{nx})))*2 == 9  | round(sqrt(length(gdcomps{nx})))*2 == 12| round(sqrt(length(gdcomps{nx})))*2 == 15
        col = round(sqrt(length(gdcomps{nx})))*2;
    elseif round(sqrt(length(gdcomps{nx})))*2 == 8  | round(sqrt(length(gdcomps{nx})))*2 == 11| round(sqrt(length(gdcomps{nx})))*2 == 14
        col = round(sqrt(length(gdcomps{nx})))*2+1;
    elseif round(sqrt(length(gdcomps{nx})))*2 == 7 | round(sqrt(length(gdcomps{nx})))*2 == 10| round(sqrt(length(gdcomps{nx})))*2 == 13
        col = round(sqrt(length(gdcomps{nx})))*2+2;        
    end;    
    row = ceil(length(gdcomps{nx})/(col/3));
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;
    activations = ws*icamatall;
    winv = pinv(ws);
    sbsph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    sbwts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
    EEG = pop_loadset('sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    EEG.icaweights=sbwts;  EEG.icasphere=sbsph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    figure;pl = 2;cc=1;
    p = 1;clear subjorder clustnum
    for ct = 1:length(clustorder)
        rct = clustorder(ct);
        if ~isempty(clustcps{rct}{nx})
            for cp = 1:length(clustcps{rct}{nx})
                subjorder(p) = clustcps{rct}{nx}(cp);
                clustnum(p) = rct; p = p+1;
            end;
        end;
    end;   
    for cp = 1:length(subjorder)
        subplot(row,col,cc)
        topoplot(EEG.icawinv(:,subjorder(cp)),EEG.chanlocs,'electrodes','off'); cc = cc+3;
        set(gca,'fontsize',7);
        title(int2str(subjorder(cp)));
        subplot(row,col,pl:pl+1)
        rcp = find(gdcomps{nx}== subjorder(cp) );
        for tp = 1:length(facs)
           ph=plot(freqs,activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp));hold on;
            set(ph,'color',cols(tp,:)); set(ph,'linewidth',1.5); set(gca,'fontsize',7);
            set(gca,'xlim',[0 40]);     %set(gca,'ylim',[-3.5 6.5]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on');
        end; pl = pl+3;       
    end;
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    ph=textsc(['Subject ',int2str(nx),'; Spectral Templates for all comps (pos orientation)','; Factors ',int2str(facs)],'title');
    set(ph,'fontsize',14);ALLEEG=[];EEG=[];    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot by Factor, plot comp spectra on top of each other
freqs = [1:.5:50]; % MUST PUT IN GRANDGDCMPS FIRST (DIPOLE CLUSTERS)
row = length(paths)-3; col = 10;  cols = cool(10);
figure; pl = 1;
for nx = 1:length(paths)-3
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;
    activations = ws*icamatall;
    for tp = 1:10
        subplot(row,col,pl)
        for cp = 1:length(gdcomps{nx})
            ph=plot(freqs,activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp));hold on;
            set(ph,'color',cols(tp,:)); set(gca,'fontsize',7);
            set(gca,'xlim',[0 40]);     set(gca,'ylim',[-3.5 6.5]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on');            
        end;pl = pl+1;
        if nx == 1
            ph =title(['Factor ',int2str(tp)]);set(ph,'fontsize',14);
        end;        
    end;
    fprintf('\n One More SUBJECT Done: %i\n',nx);
end;
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
ph=textsc(['All Subjects,'; Spectral Templates for all comps (pos orientation)'],'title');
set(ph,'fontsize',14);    

                
%%%%%%%%%%%%%%%%****************************%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot percentile graph to show variations between emmotions
emo2 ={  'anger'   'frustration'    'jealousy'    'fear'   'disgust'    'grief'  'sad'    'compassion'    'love'    'relief'     'content'   'awe'   'happy'    'joy'       'excite' }; % same as emoorder
emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17]; clear cols
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
row = 5; col = 4; 
for nx = 1:1
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    perc = {'10%','20%','30%','40%','50%','60%','70%','80%','90%'};    
    figure;pl = 1;  cpcols = hsv(length(gdcomps{nx})); p = 1;
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)';
        if mean(tpwts < 0)
            tpwts = tpwts*-1;
            activations(tp,:) = activations(tp,:)*-1;
        end;
        
        subplot(row,col,p)
        mxcp = max(activations(tp,:)); mncp = min(activations(tp,:));
        for cp = 1:length(gdcomps{nx})
            ph=plot(freqs,activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp));hold on;
            set(ph,'color',cpcols(cp,:)); set(ph,'linewidth',1.5); 
        end;p = p+1;
            set(gca,'xlim',[0 40]);     set(gca,'ylim',[mncp+mncp*.01 mxcp+mxcp*.01]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on'); set(gca,'fontsize',10); title(['Factor ',int2str(tp)]);          
        if nx == 1
            ph =title(['Factor ',int2str(tp)]);set(ph,'fontsize',14);
        end;        
        
        subplot(row,col,p); 
        clear basemat basemat1
        for e = 1:length(numtrials{nx})
            e=emoorder(e);
            tempmat = tpwts(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);clear newmat
            if e == 1
                pl=1; for pk = .1:.1:.9
                    basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat1,basemat1-basemat1,'.');hold on;
                set(ph,'color',[.5 .5 .5]);        set(ph,'markersize',18);
                pl=1; 
                for pk = .03:.01:.98
                    basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat,basemat-basemat);hold on;
                handvec(1,e) = ph;
                set(ph,'color',[.5 .5 .5]);        set(ph,'linewidth',1.75); 
                mx = max(basemat-basemat); mn = min(basemat-basemat);
            else
                pl=1;
                for pk = .1:.1:.9
                    newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat1,newmat1-basemat1,'.');hold on;
                set(ph,'color',cols(e,:));        set(ph,'markersize',18);pl=1;
                for pk = .03:.01:.98
                    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat,newmat-basemat);hold on;
                handvec(1,e) = ph;
                set(ph,'color',cols(e,:));         set(ph,'linewidth',1.75); 
                %if e == length(emos)
                %    yl = get(gca,'ylim');
                %    for tx = 1:length(basemat1)
                %        text(basemat1(tx),yl(1)-yl(1)*.1,perc{tx});
                %    end;
                %end;         
                mx1 = max(newmat-basemat); mn1 = min(newmat-basemat);
                if mx1 > mx
                    mx = mx1;
                end;
                if mn1 < mn
                    mn = mn1;
                end;                
            end; 
            set(gca,'xlim',[basemat(1)-basemat(1)*.02 basemat(end)+basemat(end)*.02]);
            set(gca,'box','off'); set(gca,'fontsize',10); title(['Factor ',int2str(tp)]);
        end;
        set(gca,'ylim',[mn+mn*.01 mx+mx*.01]); 
        p = p+1;
    end;
    ph=xlabel('Prebaseline Weights');        set(ph,'fontsize',14);
    ph=ylabel('Emotion - Prebase Weights');        set(ph,'fontsize',14); 
    ph=textsc(['Subject ',int2str(nx),'; All Comp Templates and Percentile Comparisons with Prebaseline Weights for all Factors'],'title');set(ph,'fontsize',12);
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    %print  -dpsc2 -Pcoloring 
    %close
end;

%g = legend(handvec,{'prebase'  'anger'   'frustration'    'jealousy'    'fear'   'disgust'    'grief'  'sad'    'compassion'    'love'    'relief'     'content'   'awe'   'happy'    'joy'       'excite' 'postbase'},-1);set(g,'fontsize',14);
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  FIND WHICH EMOTIONS WEIGH HIGHLY IN EACH TEMPLATE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% activations =  pca dims X freqs*comps
% winv = trials*emos X pca dims
% need to scan winv for num of trials above or below a threshold for each emotion.
clustorder = [2,8,1,13,4,11,14,3,12,10,9,6,15,5,7]; % order of dipole clusters to image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 4;ALLEEG=[];EEG=[];
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 

p = 1;clear subjorder clustnum
    for ct = 1:length(clustorder)
        rct = clustorder(ct);
        if ~isempty(grandgdcomps{rct}{nx})
            for cp = 1:length(grandgdcomps{rct}{nx})
                subjorder(p) = grandgdcomps{rct}{nx}(cp); 
                clustnum(p) = rct; p = p+1;
            end;
        end;
    end;  
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;
    activations = ws*icamatall;
    winv = pinv(ws);
emomap = ones(1,1);
for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
end;emomap(end) = emomap(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear templatewts
    
pos = 1; % 1 for positive weights, 0 otherwise
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
emoorder = [1,5,3,11,9,15,13,7,10,8,14,12,2,6,4,16,17]; clear cols
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [.5 .5 .5];

for tp = 1:size(winv,2)
    tpwts =  winv(:,tp);
    if pos == 0
        for e = 1:length(numtrials{nx})
            wtsdist = tpwts(emomap(e):emomap(e+1)-1);
            wtsdist = sort(wtsdist);
            thresh = wtsdist(round(length(wtsdist)*.25));
            emowts(e) = mean(tpwts(tpwts(emomap(e):(emomap(e+1)-1))< thresh));
            %emowts(e) = sum(tpwts(tpwts(emomap(e):(emomap(e+1)-1))> thresh)/numtrials{nx}(e));
        end;        
    else
        for e = 1:length(numtrials{nx})
            wtsdist = tpwts(emomap(e):emomap(e+1)-1);
            wtsdist = sort(wtsdist);
            thresh = wtsdist(round(length(wtsdist)*.75));
            emowts(e) = mean(tpwts(tpwts(emomap(e):(emomap(e+1)-1))> thresh));
            %emowts(e) = sum(tpwts(tpwts(emomap(e):(emomap(e+1)-1))> thresh)/numtrials{nx}(e));
        end;
    end;
    templatewts{tp} = emowts;
end;

% figure with bar graph, scalp maps and templates
% for 9-13 comps: 4 X 5
% for 14-18 comps 5 x 5
% for 19-22 comps 5 X 6
% for 23-28 comps 6 X 6
if length(gdcomps{nx}) < 11
    row = 4; col = 5;
elseif    length(gdcomps{nx}) > 10 & length(gdcomps{nx}) < 16
    row = 6; col = 5;
elseif    length(gdcomps{nx}) > 15 & length(gdcomps{nx}) < 19
    row = 6; col = 6;
elseif    length(gdcomps{nx}) > 18 & length(gdcomps{nx}) < 22
    row = 6; col = 7;
elseif    length(gdcomps{nx}) > 21 & length(gdcomps{nx}) < 29
    row = 8; col = 7;
end;

%  REQUIRES subjorder(see above), emos and emoorder and dataset with EEG.icawinv up
tp = 1;  
figure;pl = 3;
sbplot(row-3,col,1:2)

figure;for tp = 1:15
    subplot(4,4,tp)
for e = 1:length(emos)
    ph = bar(e,templatewts{tp}(emoorder(e))); hold on;
    set(ph,'facecolor',cols(e,:));
    ph = text(e,0,emos{emoorder(e)});
    set(ph,'rotation',90);
end;set(gca,'xlim',[0 18]);
title(['Template ',int2str(tp)]);
end;
figure;for e =2:length(emos)-1
    ph = bar(e-1,templatewts{tp}(emoorder(e))); hold on;
    set(ph,'facecolor',cols(e,:));
    ph = text(e-1,0,emos{emoorder(e)});
    set(ph,'rotation',90);set(ph,'fontsize',16);set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
end;
set(gca,'xlim',[0 16]);

for cmp = 1:col-2
    sbplot(row,col,pl);
    topoplot(EEG.icawinv(:,subjorder(cmp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
    set(gca,'fontsize',7);
    title(int2str(subjorder(cmp)));
end;pl = pl+2;
for cmp = 1:col-2
    sbplot(row,col,pl);
    if pos == 1
    plot(freqs,activations(tp,length(freqs)*(cmp-1)+1:length(freqs)*(cmp))); pl = pl+1;
    set(gca,'ylim',[-3.5 6.5]);   
    else
    plot(freqs,activations(tp,length(freqs)*(cmp-1)+1:length(freqs)*(cmp))*-1); pl = pl+1;
    set(gca,'ylim',[-6.5 3.5]);   
    end;        
    set(gca,'fontsize',7);
    set(gca,'xlim',[0 40]);
    set(gca,'xtick',[5:5:40]);
    set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
    set(gca,'xgrid','on');
end;
for cp = cmp+1:cmp+col
    sbplot(row,col,pl);
    topoplot(EEG.icawinv(:,subjorder(cp)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
    set(gca,'fontsize',7);
    title(int2str(subjorder(cp)));
end;
for cp = cmp+1:cmp+col
    sbplot(row,col,pl);
    if pos == 1
    plot(freqs,activations(tp,length(freqs)*(cp-1)+1:length(freqs)*(cp))); pl = pl+1;
    set(gca,'ylim',[-3.5 6.5]);   
    else
    plot(freqs,activations(tp,length(freqs)*(cp-1)+1:length(freqs)*(cp))*-1); pl = pl+1;
    set(gca,'ylim',[-6.5 3.5]);   
    end;
    set(gca,'fontsize',7);
    set(gca,'xlim',[0 40]);
    set(gca,'xtick',[5:5:40]);
        set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
    set(gca,'xgrid','on');
end;
if row > 4
    for cc = cp+1:cp+col
        sbplot(row,col,pl);
        topoplot(EEG.icawinv(:,subjorder(cc)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
        set(gca,'fontsize',7);
        title(int2str(subjorder(cc)));
    end;
end;
if row > 4
    pl = (cp+2)*2+col+1;
    for cc = cp+1:cp+col
        sbplot(row,col,pl);
    if pos == 1
    plot(freqs,activations(tp,length(freqs)*(cc-1)+1:length(freqs)*(cc))); pl = pl+1;
    set(gca,'ylim',[-3.5 6.5]);   
    else
    plot(freqs,activations(tp,length(freqs)*(cc-1)+1:length(freqs)*(cc))*-1); pl = pl+1;
    set(gca,'ylim',[-6.5 3.5]);   
    end;
        set(gca,'fontsize',7);
        set(gca,'xlim',[0 40]);
        set(gca,'xtick',[5:5:40]);
        set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
        set(gca,'xgrid','on');
    end;
end;
if row > 6
    for ccc = cc+1:cc+col
        sbplot(row,col,pl);
        topoplot(EEG.icawinv(:,subjorder(ccc)),EEG.chanlocs,'electrodes','off','plotrad',.5); pl = pl+1;
        set(gca,'fontsize',7);
        title(int2str(subjorder(ccc)));
    end;
end;    
if row > 6
    pl = (cc+2)*2+col+1;
    for ccc = cc+1:cc+col
        sbplot(row,col,pl);
        if pos == 1
            plot(freqs,activations(tp,length(freqs)*(ccc-1)+1:length(freqs)*(ccc))); pl = pl+1;
            set(gca,'ylim',[-3.5 6.5]);   
        else
            plot(freqs,activations(tp,length(freqs)*(ccc-1)+1:length(freqs)*(ccc))*-1); pl = pl+1;
            set(gca,'ylim',[-6.5 3.5]);   
        end;
        set(gca,'fontsize',7);
        set(gca,'xlim',[0 40]);
        set(gca,'xtick',[5:5:40]);
        set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
        set(gca,'xgrid','on');
    end;
end;    

%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,14,6,4,16,17];
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [.5 .5 .5];
figure;row = 3;  col = 4; pl = 1;
for tp = 1:length(templatewts)
    subplot(row, col, tp)
    for e = 1:length(emos)
        ph = bar(e,templatewts{tp}(emoorder(e))); hold on;
        set(ph,'facecolor',cols(e,:));
        ph = text(e,0,emos{emoorder(e)});
        set(ph,'rotation',90);
    end;set(gca,'xlim',[0 18]);
    title(['Template ',int2str(tp)]);
end;

% plot all emotion distributions on single graph for each template
emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,14,6,4,16,17];
clear cols; cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [.5 .5 .5];
figure;row = 3;  col = 4; pl = 1;
for tp = 1:size(activations,1)
    subplot(row, col, tp)
    plot([0 100],[0 0],'k-');hold on;
    for e = 1:length(emos)
        xax = [100/numtrials{nx}(emoorder(e)):(100/numtrials{nx}(emoorder(e))):100];
        ph = plot(xax,sort(winv(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1),tp)));hold on;
        set(ph,'color',cols(emoorder(e),:));set(ph,'linewidth',2);
    end;
    title(['Subj: ',int2str(nx),'; Template ',int2str(tp)]);
end;
legend('prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase');


% for one template, show wt distribution for each emotion
emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,14,6,4,16,17];
for tp = 1:size(winv,2)
    figure; row = 4; col = 5;
    for e = 1:length(emos)
        subplot(row,col,e)
        mnwt = mean(winv(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1),tp));
        hist(winv(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1),tp),50);hold on;
        set(gca,'xlim',[min(winv(:,tp)) max(winv(:,tp))]);
        set(gca,'ylim',[0 30]);
        %set(gca,'edgecolor',cols(emoorder(e),:));
        ph=plot([mnwt mnwt],[get(gca,'ylim')],'g-');
        set(ph,'linewidth',1.5);
        ph=plot([0 0],[get(gca,'ylim')],'r-');
        set(ph,'linewidth',1.5);ph =title(emos(emoorder(e)));
        set(ph,'fontsize',14);
    end;
    ph =textsc(['Template ',int2str(tp),' for Subject ',int2str(nx),': Distribution of Weights across Emotions'],'title');
    set(ph,'fontsize',14);
end;

% possible metric:
% every emotion can have one score representing the mean of the wts for a given template 
% Then, across templates you could make a visual representation of the constellation of "weights"
figure; row = 3;col = 3;pl = 1;
for nx = 1:9%length(gdcomps)
    if pl == 10
        ph =textsc(['Subject ',int2str(nx),' Mean Wts Across Emotions'],'title');
        set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        figure; pl = 1;
    end;
    p = 1;clear subjorder clustnum
    for ct = 1:length(clustorder)
        rct = clustorder(ct);
        if ~isempty(grandgdcomps{rct}{nx})
            for cp = 1:length(grandgdcomps{rct}{nx})
                subjorder(p) = grandgdcomps{rct}{nx}(cp); 
                clustnum(p) = rct; p = p+1;
            end;
        end;
    end;  
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);
    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;emomap(end) = emomap(end);
    clear emoscore
    for tp = 1:size(winv,2)
        for e = 1:length(emos)
            emoscore(e,tp) = mean(winv(emomap(e):(emomap(e+1)-1),tp));
        end;
    end;
    subplot(row,col,pl);
    imagesc(emoscore,[-2 2]);hold on;pl = pl+1;
    %imagesc(emostd,[-2 2]);hold on;pl = pl+1;
    for e = 1:length(emos)
        ph = text(.75,e,emo2{e});
        set(ph,'fontsize',12);set(ph,'color','w');
    end;
    set(gca,'yticklabel',[]);set(gca,'fontsize',14);
    xlabel('Templates');ylabel('Emotions');colorbar;
    title(['Subject ',int2str(nx)]);
end;
ph =textsc(['Mean Weights Across Emotions and Templates'],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

% plot in 3 dimensions to see if emotions separate
clear cols; cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [.5 .5 .5];
tp1 = 1;  tp2 = 2;  tp3 = 3;  
nx = 1; 
    sph=floatread(['/data/common2/emotion/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;
    activations = ws*icamatall;
    winv = pinv(ws);
    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;emomap(end) = emomap(end);
    clear emoscore emostd
    for tp = 1:size(winv,2)
        for e = 1:length(emos)
            emoscore(e,tp) = mean(winv(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1),tp));
            emostd(e,tp) = std(winv(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1),tp));
        end;
    end;
figure;
for k = 1:size(emoscore,1)
ph =plot3(emoscore(k,tp1),emoscore(k,tp2),emoscore(k,tp3),'k.');hold on;
set(ph,'color',cols(k,:)); set(ph,'markersize',35);
set(gca,'xgrid','on');set(gca,'ygrid','on');set(gca,'zgrid','on');
text(emoscore(k,tp1)+.01,emoscore(k,tp2)+.01,emoscore(k,tp3)+.01,emos{emoorder(k)});
end;
xlabel(['Template ',int2str(tp1)]);ylabel(['Template ',int2str(tp2)]);zlabel(['Template ',int2str(tp3)]);
ph =textsc(['Subject ',int2str(nx)],'title');
set(ph,'fontsize',15);

% k means clustering of weights
emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,14,6,4,16,17];
idx = kmeans(winv,10);
figure;
for e=1:length(emos)
    subplot(5,4,e)
    plotclust{e} = hist(idx(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1)),3)/numtrials{nx}(e);
hist(idx(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1)));
set(gca,'ylim',[0 200]);set(gca,'xlim',[0 11]);
title(emos{emoorder(e)});
end;
clear cols; cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [.5 .5 .5];
figure;
for e=1:length(plotclust)
ph=plot3(plotclust{e}(1),plotclust{e}(2),plotclust{e}(3),'k.');hold on;
set(ph,'color',cols(e,:)); set(ph,'markersize',15);
text(plotclust{e}(1)+.01,plotclust{e}(2)+.01,plotclust{e}(3)+.01,emos{emoorder(e)});
set(gca,'xgrid','on');set(gca,'ygrid','on');set(gca,'zgrid','on');
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instead of subj by subj, put in freqs X trials from all subjects to find max indep across emos
freqs = [1:.5:50];
icamatall = zeros(length(freqs),0);
for nx = 11:length(gdcomps)
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    for k = 1:length(Alllongersps)
        oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
        icamat = zeros(length(freqs),0);
        for cmp = 1:length(gdcomps{nx})
            [outdata,outx] = movav(oneemo(:,:,gdcomps{nx}(cmp)),[1:size(oneemo,2)],20);
            icamat(:,end+1:end+size(outdata,2)) = outdata;
        end;
        icamatall(:,end+1:end+size(icamat,2)) = icamat;
        fprintf('\n One More Emotion Done: %i',k);
        subjtrials(k) = size(outdata,2);
    end;
    fprintf('\n One More SUBJECT Done: %i\n',nx);
    numtrials{nx} = subjtrials;  clear Alllongersps  totrials icamat  onecmp oneemo
end;
floatwrite(icamatall, '/data/common2/emotion/freqsXtrials.fdt'); 
size(icamatall)
save /data/common2/emotion/numtrials-frqsXtrials.mat numtrials freqs
/data/common/matlab/ica_linux2.4 < /data/common2/emotion/ClustPwrICA.sc

% run ica on Freqs X Comps*Trials for one subject
[wts,sph,compvars,bias,signs,lrates,activations]  = runica(icamatall,'pca',25,'extended',1,'stop',1e-8);
ws = wts*sph;
activations = ws*icamatall;
winv = pinv(ws);
figure;
for tp = 1:25
    subplot(5,5,tp)
    ph =plot(freqs,winv(:,tp),'k-');
    set(ph,'linewidth',1.5);
    set(gca,'xlim',[0 50]);
    set(gca,'ylim',[-5 6]);
    set(gca,'xgrid','on');
    title(['Template ',int2str(tp)]);
end;
% k means clustering of weights
idx = kmeans(winv,4);
figure;
for e=1:length(emos)
    subplot(5,4,e)
hist(idx(emomap(emoorder(e)):(emomap(emoorder(e)+1)-1)));
title(emos{emoorder(e)});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze output of ICA
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
load /data/common2/emotion/numtrials-frqsXtrials.mat numtrials freqs
ixmap = zeros(1,1); % ixmap marks where each subject BEGINS (+1)
for nx = 2:length(numtrials)
    ixmap(1,nx) = ixmap(1,nx-1) + sum(numtrials{nx-1})*length(gdcomps{nx-1});
end;

sph=floatread('/data/common2/emotion/freqsXtrialsPC25.sph',[99 99]); 
wts=floatread('/data/common2/emotion/freqsXtrialsPC25.wts',[25 99]); 
icamatall = floatread('/data/common2/emotion/freqsXtrials.fdt',[99 765581]);
ws = wts*sph;
activations = ws*icamatall;
winv = pinv(ws);


figure;
for tp = 1:25
    subplot(5,5,tp)
    ph =plot(freqs,winv(:,tp),'k-');
    set(ph,'linewidth',1.5);
    set(gca,'xlim',[0 50]);
    set(gca,'ylim',[-5 6]);
    set(gca,'xgrid','on');
    title(['Template ',int2str(tp)]);
end;

% find weightings for each emotion/subject
% activations are PCA dims X trials (subj,comps,emotion,trials)
% emoscore gives center of gravity weight for each emotion
% plots for each subject and image of scores: emotions X template #

emoorder = [1,3,5,11,9,15,13,7,10,8,2,12,14,6,4,16,17];
figure; row = 4;  col = 4;  pl =1;
for nx = 18:21%length(numtrials)
    clear plotscore varscore kurtscore skewscore
    for tp = 1:12%size(activations,1)
        for k = 1:length(numtrials{nx})
            emo = emoorder(k);
            if emo == 1
                skewscore(k,tp) = skewness(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx}))); % tells how far and in which direction distribution skew is
                kurtscore(k,tp) = kurtosis(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})));
                varscore(k,tp) = var(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})));
                plotscore(k,tp) = mean(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})));
            else
                skewscore(k,tp) = skewness(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx})));
                kurtscore(k,tp) = kurtosis(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx})));
                varscore(k,tp) = var(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx})));
                plotscore(k,tp) = mean(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx})));
            end;
        end;
    end;
    subplot(row,col,pl)
    imagesc(plotscore,[min(min(plotscore)) abs(min(min(plotscore)))]); pl = pl+1;
    for e = 1:length(emos)
        ph = text(1,e,emos{emoorder(e)});
        set(ph,'fontsize',8);set(ph,'color','w');
    end;colorbar;
    set(gca,'yticklabel',[]);set(gca,'fontsize',14);
    xlabel('Templates');ylabel('Emotions');
    title(['Subject ',int2str(nx),' Center of Gravity']);
    subplot(row,col,pl)
    imagesc(skewscore,[-max(max(skewscore)) max(max(skewscore))]); pl = pl+1;
    for e = 1:length(emos)
        ph = text(1,e,emos{emoorder(e)});
        set(ph,'fontsize',8);set(ph,'color','w');
    end;colorbar;
    set(gca,'yticklabel',[]);set(gca,'fontsize',14);
    title(['Skewness']);
    subplot(row,col,pl)
    imagesc(kurtscore,[0 max(max(kurtscore))]); pl = pl+1;
    for e = 1:length(emos)
        ph = text(1,e,emos{emoorder(e)});
        set(ph,'fontsize',8);set(ph,'color','w');
    end;colorbar;
    set(gca,'yticklabel',[]);set(gca,'fontsize',14);
    title(['Kurtosis']);
    subplot(row,col,pl)
    imagesc(varscore,[0 max(max(varscore))]); pl = pl+1;
    for e = 1:length(emos)
        ph = text(1,e,emos{emoorder(e)});
        set(ph,'fontsize',8);set(ph,'color','w');
    end;colorbar;
    set(gca,'yticklabel',[]);set(gca,'fontsize',14);
    title(['Variance']);
end;
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

% look at actual distributions
nx=4; tp =4;
figure; row = 4;  col = 4;  pl =1;
for k = 1:16    
    subplot(row,col,k);
    emo = emoorder(k);
    if emo == 1
        hist(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})),25);hold on;
    else
        hist(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx})),25);hold on;       
    end;
    plot([0 0],[get(gca,'ylim')],'r-');
    title(emos{emo});
end;

% need to see what components within distribution are on fringes
load /data/common2/emotion/subjdims.mat % subjdims  numtrials comment
fac = .1;  % percent of top or bottom of distribution to sample
for nx = 1:length(gdcomps)
    for tp = 1:12
        for k = 1:length(emos)
            emo = emoorder(k);
            if emo == 1
                dist = sort(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})));
                valcut = dist(end-round(length(dist)*fac));
                hiscores = find(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx}))>valcut);
            else
                dist = sort(activations(tp,ixmap(nx)+1:ixmap(nx)+sum(numtrials{nx}(emo))*length(gdcomps{nx})));
                valcut = dist(end-round(length(dist)*fac));
                hiscores = find(activations(tp,ixmap(nx)+sum(numtrials{nx}(1:emo-1))*length(gdcomps{nx})+1:ixmap(nx)+sum(numtrials{nx}(1:emo))*length(gdcomps{nx}))>valcut);
            end;
           for cp = 1:length(gdcomps{nx})
               if ~isempty(find(hiscores>(cp-1)*numtrials{nx}(k)+1 & hiscores < cp*numtrials{nx}(k)))
                   cpmat(k,cp) = length(find(hiscores>(cp-1)*numtrials{nx}(k)+1 & hiscores < cp*numtrials{nx}(k))); % tells number of trials included in top fac%
               end;
           end; % for cp
        end; % for k
    end;
end;

           