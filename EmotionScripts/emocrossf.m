% runs crossf on button-pressing emotion subjects

eeglab
subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];%  jo74
subj2 = [1,3:23,25:30,34,35,36];  % tl81
subj3 = [1,2,4:9,12:24,26,30,31,33,34,46,51];   % mi83
subj4 = [10,12,13,15,17,18,20,21,22,23,28,33,36,45,58];   % ms82
gdcomps = {subj1, subj2, subj3, subj4};

paths = {'/data/common2/emotion/jo74/','/data/common2/emotion/tl81/','/data/common2/emotion/mi83/','/data/common2/emotion/ms82/'};
datpaths = {'/data/common2/emotion/jo74/crosscomps/','/data/common2/emotion/tl81/crosscomps/','/data/common2/emotion/mi83/crosscomps/','/data/common2/emotion/ms82/crosscomps/'};

sphs = {'/data/common2/emotion/jo74/sph252-160.sph','/data/common2/emotion/tl81/sph241.sph','/data/common2/emotion/mi83/sph248.sph','/data/common2/emotion/ms82/sph241.sph'};
wtss = {'/data/common2/emotion/jo74/wts252-160.wts','/data/common2/emotion/tl81/wts241.wts','/data/common2/emotion/mi83/wts248.wts','/data/common2/emotion/ms82/wts241.wts'};


sphsize = {[252 252],[241 241],[248 248],[241 241],};
wtssize = {[160 252],[160 241],[160 248],[160 241]};

nx = 1 ; % input desired subj number
sph=floatread(sphs{nx},sphsize{nx}); 
wts=floatread(wtss{nx},wtssize{nx}); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Run timef on button presses of each emotion
emoset = {'prebase.set','awe.set','frustration.set','joy.set','anger.set','happy.set','sad.set','love.set','fear.set' ,'compassion.set','jealousy.set','content.set','grief.set','relief.set','disgust.set','excite.set','postbase.set'};

comment = 'Crossf on epochs created by non-overlapping time periods during button press periods(while experiencing emotion) ;crossf(EEG.icaact(gdcomps{nx}(index1),:),EEG.icaact(gdcomps{nx}(index2),:), EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000], EEG.srate, [3 .5], alpha, 0.01,padratio, 4,winsize,256);saved in order: prebase,awe.set,frustration.set,joy.set,anger.set,happy.set,sad.set,love.set,fear.set ,compassion.set,jealousy.set,content.set,grief.set,relief.set,disgust.set,excite.set,postbase.set'; 

figure;  
for e = 17:length(emoset)
    EEG = pop_loadset( emoset{e},paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);   
  % make fake events every 3 sec throughout  
      if exist('EEG.event(2)')
    ft = EEG.event(2).latency;
    else
        ft = 1;
    end;    
    for evtm = ft+384:768:size(EEG.data,2)-384  
        EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 'fake';        
    end;
    EEG = pop_epoch( EEG,{'fake'} , [-1.5 1.5]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1500 1500]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% load wts
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% reject noise
    EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
    EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);

    coher=zeros(63,200,gdcomps{nx}(end));
    crsboot=zeros(63,2,gdcomps{nx}(end));
    crsangle=zeros(63,200,gdcomps{nx}(end));
    cohercell = cell (1,gdcomps{nx}(end));
    crsbootcell = cell (1,gdcomps{nx}(end));
    crsanglecell = cell (1,gdcomps{nx}(end));
    %-- always list lowest to highest comp numbers-------------
    %done =
    for index1=1:length(gdcomps{nx})-1
        for index2=index1+1:length(gdcomps{nx})
            [coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher, alltfX, alltfY] = crossf(EEG.icaact(gdcomps{nx}(index1),:),EEG.icaact(gdcomps{nx}(index2),:), EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000], EEG.srate, [3 .5], 'alpha', 0.01,'padratio', 4,'winsize',256);
            coher(:,:,gdcomps{nx}(index2)) = coh;
            cohercell{1,gdcomps{nx}(index1)} = coher;
            crsboot(:,:,gdcomps{nx}(index2)) = cohboot;
            crsbootcell{1,gdcomps{nx}(index1)} = crsboot;
            crsangle(:,:,gdcomps{nx}(index2)) = cohangles;
            crsanglecell{1,gdcomps{nx}(index1)}= crsangle;
            clear coh mcoh cohboot cohangles allcoher alltfX alltfY
            clf
        fprintf('\nShowing Comp: %i vs Comp %i\n',index1,index2);
        fprintf('\n..... in dataset: %s\n',emoset{e});
        end
    end;
    emocoh{e} = cohercell;
    emoboot{e} = crsbootcell;
    emoangles{e} = crsanglecell;
    clear coher cohercell crsboot crsbootcell crsangle crsanglecell
    ALLEEG=[]; EEG=[];
end;
cd (datpaths{nx})    
save RndEpochsXcoh.mat emocoh freqsout timesout comment emoboot
save RndEpochsXang.mat emoangles freqsout timesout comment emoboot
clear emocoh emoboot emoangles    
ALLEEG=[]; EEG=[];
pack

load RndEpochsXcoh.mat
 oneboot = emoboot{1};                                      
 onecoh = emocoh{1};  
 save prebase.mat  onecoh freqsout timesout comment oneboot 
 oneboot = emoboot{2};                                      
 onecoh = emocoh{2};  
 save awe.mat  onecoh freqsout timesout comment oneboot    
 oneboot = emoboot{3};                                 
 onecoh = emocoh{3};                                   
 save frustration.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{4};                                         
 onecoh = emocoh{4};                                           
 save joy.mat  onecoh freqsout timesout comment oneboot        
 oneboot = emoboot{5};                                 
 onecoh = emocoh{5};                                   
 save anger.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{6};                                   
 onecoh = emocoh{6};                                     
 save happy.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{7};                                   
 onecoh = emocoh{7};                                     
 save sad.mat  onecoh freqsout timesout comment oneboot  
 oneboot = emoboot{8};                                 
 onecoh = emocoh{8};                                   
 save love.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{9};                                  
 onecoh = emocoh{9};                                    
 save fear.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{10};                                  
 onecoh = emocoh{10};                                   
 save compassion.mat  onecoh freqsout timesout comment oneboot
 oneboot = emoboot{11};                                       
 onecoh = emocoh{11};                                         
 save jealousy.mat  onecoh freqsout timesout comment oneboot  
 oneboot = emoboot{12};                                     
 onecoh = emocoh{12};                                       
 save content.mat  onecoh freqsout timesout comment oneboot 
 oneboot = emoboot{13};                                    
 onecoh = emocoh{13};                                    
 save grief.mat  onecoh freqsout timesout comment oneboot  
 oneboot = emoboot{14};                                    
 onecoh = emocoh{14};                                  
 save relief.mat  onecoh freqsout timesout comment oneboot 
 oneboot = emoboot{15};                                   
 onecoh = emocoh{15};                                  
 save disgust.mat  onecoh freqsout timesout comment oneboot 
 oneboot = emoboot{16};                                   
 onecoh = emocoh{16};                                  
 save excite.mat  onecoh freqsout timesout comment oneboot 
 oneboot = emoboot{17};                                   
 onecoh = emocoh{17};                                  
 save postbase.mat  onecoh freqsout timesout comment oneboot 

% Visualize Crosses
emos = {'PREBASELINE Cross coherence',' AWE -- Randomly epoched during button presses', ' FRUSTRATION -- Randomly epoched during button presses',' JOY -- Randomly epoched during button presses',' ANGER -- Randomly epoched during button presses',' HAPPY -- Randomly epoched during button presses',[' ' ...
                    'SAD -- Randomly epoched during button presses'],' LOVE -- Randomly epoched during button presses' ,' FEAR -- Randomly epoched during button presses',' COMPASSION -- Randomly epoched during button presses',[' JEALOUSY -- ' ...
                    'Randomly epoched during button presses'],' CONTENT -- Randomly epoched during button presses',' GRIEF -- Randomly epoched during button presses','RELIEF -- Randomly epoched during button presses','DISGUST -- Randomly epoched during button presses','EXCITEMENT -- Randomly epoched during button presses','POSTBASELINE Cross coherence'};
nx= 4;
EEG = pop_loadset( 'sources.set', paths{nx});
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
compsrun = gdcomps{nx}; ; 
lim = .7;
f = find (freqsout<30);
t = find (timesout>-1500 & timesout<1500);
for e = 1:length(emos)
    figure; clear forsum
    % first plot scalp maps in top row
    for mp = 2:length(compsrun)   
        subplot(length(compsrun),length(compsrun),mp)
    topoplot(EEG.icawinv(:,compsrun(mp)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
        title(int2str(compsrun(mp)));
    end;
    subplot(length(compsrun),length(compsrun),length(compsrun)+1)
    topoplot(EEG.icawinv(:,compsrun(mp)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
    title(int2str(compsrun(1)));
    for kp = 2:length(compsrun)-1
        subplot(length(compsrun),length(compsrun),length(compsrun)*kp+kp)
        topoplot(EEG.icawinv(:,compsrun(kp)),EEG.chanlocs,'electrodes','off');
        title(int2str(compsrun(kp)));
    end;    
    p = ((length(compsrun)-1)/2)*length(compsrun);
    for w=1:length(compsrun)-1
        a = compsrun(w);
        b = compsrun(w+1:end);
        cohercell = emocoh{e};
        Diff = cohercell{1,a};
        crsbootcell = emoboot{e};
        bootmat1 = crsbootcell {1,a};   %chooses a 3D array from cell array
        minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
        minmask = repmat (minboot,[1,size(timesout,2),1]); %makes (63,200,71) matrix of minboot
        Diff(find(Diff <= minmask)) = 0;%finds and zeros cohmats less than minmask
        dd = Diff(f,t,:);
        forsum(:,:,w) = sum(dd,3);
        for m=w+1:length(b)+w
            % subplot(length(compsrun)-1,length(compsrun)-1,((length(compsrun)-1)*(w-1))+(m-1))
            subplot(length(compsrun),length(compsrun),((length(compsrun))*(w-1))+(m-1)+length(compsrun)+1)
            imagesc(timesout(t),freqsout(f),Diff(f,t,b(m-w)),[-lim lim]);
            hold on
            set(gca,'ticklength',[0.05 0.05]); 
            set(gca,'tickdir','out');      
            set(gca,'xtick',[-500:250:500]); 
            set(gca,'xticklabel',[]);    
            set(gca,'fontsize',5);    
            set(gca,'ytick',[5:5:30]);
            set(gca,'yticklabel',[]); 
            %set(gca,'xlim',[-lenset(1,e)/2 lenset(1,e)]);
            %plot([0 0],[0 30],'k-');
            title(int2str(b(m-w)))
        end;
    end;
   %set(gca,'xticklabel',[{-2 [] -1 [] 0 [] 1 [] 2}]);    
    set(gca,'yticklabel',[5:5:30]);    
    colorbar
    textsc(emos{e},'title');
    %axcopy
    newlim = max(max(sum(forsum)))-.4;
    sbplot(5,10,32:33)
    imagesc(timesout(t),freqsout(f),sum(forsum,3),[-newlim newlim]); colorbar
            set(gca,'ticklength',[0.03 0.03]); 
     title('Sum of Cross Cohs');
end;
