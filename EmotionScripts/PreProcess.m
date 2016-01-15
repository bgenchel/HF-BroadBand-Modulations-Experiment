% new preprocessing script for emotion expt.

ssh juggling
qlogin
matlab


addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
newpaths{22} = [newpaths{22}(1:22),'mr74/'];
%newpaths{22} = [newpaths{22}(1:22),'mr72-2/'];


badchan{38} = {'A12'    'B1'    'B2'    'B3'    'C1'    'C6'    'C19'    'C20'    'D3' 'D6'    'D23'    'D27'    'D32'    'E5'    'F4'    'F32'    'G9'    'G21'  'H17'};


% ImportEmoRaw(bdf,elp,savepath,savename,saveto,blocks,getpress);
EEG = pop_loadset('filename','ButtonEvents.set','filepath',[newpaths{nx}]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
load /data/common4/emotion/ICstuff.mat b l m 
  emos = [emos 'prebase', 'postbase'];
savename = 'NoFilt';
easysubjs = [4:8,10:34];% 1:3,9,35 by hand
% subjs 1,2,3 need exit for excite
% subj 3: can't use E3 ref and:
%          EEG.event(47).type = 'fear';
%          EEG.event(54).type = 'jealousy';
%          EEG.event(68).type = 'grief';
% subj 2:
%          EEG.event(129).type = 'fear';
%          EEG.event(149).type = 'jealousy';
%          EEG.event(177).type = 'grief';
% subj 9,35 has no button (all datsets done)
% subj 1 has two '15' codes, but no '20'.change 2nd 'love' to 'grief'...  (EEG.event(139) = 'grief';?)
% gotta check Subj 1:3 events carefully...subj 1,2,3 done.
% fix subjs 11(problem?)    15(suspect events)   23?(2 bdf files)
      x=openbdf([newpaths{nx},'.bdf']);
redosubjs = [];
newdir = '/data/projects/julieo/emotion';
for nxx = 19:length(easysubjs)
  nx=easysubjs(nxx);
  ALLEEG=[];EEG=[];
  EEG = pop_loadset('sources.set' ,newpaths{nx});
  xold = {EEG.chanlocs.labels};
  wts = EEG.icaweights; sph = EEG.icasphere; 
  [EEG]=ImportEmoRaw([newpaths{nx}(end-4:end-1),'.bdf'],[newpaths{nx}(end-4:end-1),'-256.elp'],fullpaths{nx},savename,[newdir,newpaths{nx}(end-5:end)],[],1); %  get presses
  [ALLEEG EEG index] = eeg_store(ALLEEG, EEG, 1);
  %EEG = pop_select( EEG, 'nochannel',badchan{nx} ); 
  EEG.gdcomps = gdcomps{nx};
  EEG.blink = b{nx};
  EEG.lateyes = l{nx};
  EEG.muscle = m{nx};
  delmore = setdiff({EEG.chanlocs.labels},xold);
  EEG = pop_select( EEG, 'nochannel',delmore ); 
  EEG.icaweights=wts; EEG.icasphere=sph; 
  EEG = eeg_checkset(EEG);
  [ALLEEG EEG index] = eeg_store(ALLEEG, EEG, 1);
  
  stev= zeros(1,length(emos));
  sptm= zeros(1,length(emos));
  spev= zeros(1,length(emos));
  for em = 1:length(emos)
    try
      if em < 16
        %ev = find(ismember({EEG.event.type},emos{em}));
        %tms = (EEG.event(ev+1).latency-EEG.event(ev).latency)/256;
        %tm = tms/3;
        %nd = (EEG.event(ev+1).latency)/256;
        for ev = 1:length(EEG.event)
          if strcmp(emos{em},EEG.event(ev).type)
            emevent(1,em) = ev; 
            for evv = ev+1:length(EEG.event)
              if strcmp('press1',EEG.event(evv).type)
                stev(1,em) = evv;
                for evvv = evv+1:length(EEG.event)
                  if strcmp('exit',EEG.event(evvv).type)
                    sptm(1,em) = (EEG.event(evvv).latency-EEG.event(evv).latency)/EEG.srate;
                    spev(1,em) = evvv; 
                    break;
                  end;
                end;  
                break;
              end;
            end;
            break;
          end;
        end;
        EEG = pop_epoch( EEG,{} , [-.5 sptm(1,em)+.5],'eventindices',stev(1,em), 'newname', [fullpaths{nx}(end-4:end-1),'-',emos{em}]);            
        EEG = pop_rmbase( EEG, [], [1:size(EEG.data,2)] );
        EEG.comments = pop_comments('', '', strvcat(EEG.comments,'Removed baseline after epoching on emotion'));
      else
        %EEG = pop_epoch( EEG,{emos{em}} , [round(tm) round(tm)+floor(tm)*2],'newname', [fullpaths{nx}(end-4:end-1),'-',emos{em}]);    % subj 9 and 35  
        EEG = pop_epoch( EEG,{emos{em}} , [-.5 120],'newname', [fullpaths{nx}(end-4:end-1),'-',emos{em}]);            
        EEG = pop_rmbase( EEG, [], [1:size(EEG.data,2)] );
        EEG.comments = pop_comments('', '', strvcat(EEG.comments,'Removed baseline after epoching on emotion'));
      end;
      EEG = pop_saveset( EEG, [emos{em},'_NF.set'], [newdir,newpaths{nx}(end-5:end)],'savemode','twofiles');  
    catch
      %redosubjs = [redosubjs,nx];
    end;
    EEG = ALLEEG(1);
  end;
end;     
datset = {'anger_NF.set','frustration_NF.set','jealousy_NF.set','fear_NF.set' ,'disgust_NF.set','grief_NF.set','sad_NF.set','compassion_NF.set','love_NF.set','relief_NF.set','content_NF.set','awe_NF.set','happy_NF.set','joy_NF.set','excite_NF.set','prebase_NF.set','postbase_NF.set'}; % for all new ones

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% # datasets/subj: made by (manually input subjs 3(6),22(5),23(9)
for nx = 1:21 % 22 by hand
    x=openbdf([fullpaths{nx},fullpaths{nx}(end-4:end-1),'.bdf']);
    numframes = x.Head.NRec;
    nsets(nx) = ceil(numframes/1000);
    fid = fopen([fullpaths{1}(1:end-5),'Nsets'],'a');
    fprintf(fid, '\nnsets(%s) = %s;', int2str(nx),int2str(nsets(nx)));
    fclose(fid);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% filtered above .5 Hz, no lowpass, 137/212 reference
savename = 'Emo-HP';
for nx=1:35
    clear blocks
    x=openbdf([fullpaths{nx},fullpaths{nx}(end-4:end-1),'.bdf']);
    numframes = x.Head.NRec;
    strt = 0;
    for bl = 1:ceil(numframes/1000)
        if bl ~= ceil(numframes/1000)
            blocks{bl} = [strt strt+999];
        else
            blocks{bl} = [strt numframes];
        end;
        strt = strt + 1000;
    end;
    %ImportEmoRaw([fullpaths{nx}(end-4:end-1),'.bdf'],[fullpaths{nx}(end-4:end-1),'-256.elp'],fullpaths{nx},savename,blocks,1); % get presses
    ImportEmoRaw([fullpaths{nx}(end-4:end-1),'.bdf'],[fullpaths{nx}(end-4:end-1),'-256.elp'],fullpaths{nx},savename,blocks,0); % don't get presses
    %ImportEmoRawNew([fullpaths{nx}(end-4:end-1),'.bdf'],[fullpaths{nx}(end-4:end-1),'.elp'],fullpaths{nx},savename,blocks); %For all data post hf45 (and including)
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bl = 1:nsets(nx)
    EEG = pop_loadset([savename,'-',int2str(bl),'.set'] ,fullpaths{nx});
    EEG=pop_chanedit(EEG,  'load',{['/data/common4/emotion/jc66/jc66.elp'],'filetype','autodetect'});
    EEG=pop_chanedit(EEG,  'eval','', 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[]);');
    EEG = pop_saveset( EEG,[savename,'-',int2str(bl),'.set'],fullpaths{nx} ); 
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of bad channels for each subject in DataInfo.m
%nx = 38;
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed

numchans = 254; pcs = 0;
zfac = [];  ALLEEG=[];EEG=[];
savename = 'Emo-HP';
wtstem = 'EMO';rejthresh = [1000 5];
for nx = 1:35
   clear datsets
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'.set'];
    end;

%[EEG,delchans,nchan] = Cut2ChanX(datsets,wtstem,fullpaths{nx},numchans,zfac,['E3' 'G23',badchan{nx}],[],[5 5 5]);
[EEG,delchans,nchan] = Cut2ChanX(datsets,wtstem,fullpaths{nx},numchans,zfac,{'E3' 'G23'},[],[5 5 5]);
    %[status, result] = system(['\rm ',fullpaths{nx},'Emo-HP-?.set']);  % delete orig datasets
    
    fid = fopen([fullpaths{1}(1:end-5),'Nchan'],'a'); % save # of chans for each subj
    fprintf(fid, '\nnchans(%s) = %s;', int2str(nx),int2str(nchan));
    fclose(fid);
    fid = fopen([fullpaths{1}(1:end-5),'DelChan'],'a'); % save deleted chans as 'bad'
    fprintf(fid, '\nbadchan{%s} = {', int2str(nx));
    for d = 1:length(delchans)
        fprintf(fid, ' ''%s'' ', delchans{d});
    end;
    fprintf(fid, '};');    fclose(fid);
clear datsets
    for bl = 1:nsets(nx)
        %datsets{bl} = ['Emo-',int2str(bl),'.set'];
        datsets{bl} = [wtstem,'-',int2str(nchan),'-',int2str(bl),'.set'];
    end;
    [EEG,numframes,nchan] = MakeFloats(tmpdatsets,fullpaths{nx},wtstem,[],[],rejthresh,0);    
    [run_ica_str] = LinuxICA(wtstem,'ICAinstruct',fullpaths{nx},sum(numframes),nchan,0,1);
end;

dat = floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.fdt'],[nchans(nx) inf],[],0);
dat = floatread([fullpaths{nx},wtstem,int2str(nchan),'.fdt'],[nchan inf],[],0);

/data/common/matlab/ica_linux2.4 < /data/common4/emotion/hf45/ICAinstruct.sc
/data/common/matlab/ica_linux2.4 < /data/common4/emotion/jc66/ICAinstruct.sc

%%%%% Use Jason's Formica:------------------

formica0([fullpaths{nx},wtstem,int2str(nchans(nx)),'.fdt'],['/home/julie/EmotICA/',paths{nx}],nchans(nx),sum(numframes),'num_models',1);
% call in weights

dat = floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.fdt'],[nchans(nx) inf],[],0);
nframes = size(dat,2); clear dat
[A,varord,LLt,LL,c,W,S,gm,alpha,mu,sbeta,rho,mn,nd,svar] = loaddat5(['/home/julie/EmotICA/',paths{nx}],nchans(nx),nchans(nx),1,3,10000,nframes);

EEG = pop_loadset( datsets{ds},fullpaths{nx},'all' );
EEG.icawinv = A(:,varord);
EEG.icaweights = W(varord,:);
EEG.icasphere = S;
pop_topoplot(EEG,EEG.chanlocs,[1:36],'','electrodes','off');



%%%%-------------------------------------------------------
% after cut dataset is created:  
for ds = 1:nsets(nx)
    EEG = pop_loadset( datsets{ds},fullpaths{nx},'all' );
    EEG = pop_select(EEG,'nochannel',[24,106]);
    EEG = pop_saveset( EEG,datsets{ds},fullpaths{nx} ); 
end;

% check weights in progress:
chlocs = EEG.chanlocs;
% then wait til ica starts:
watchcomps([fullpaths{nx},'binica.wts'],[fullpaths{nx},wtstem,int2str(nchans(nx)),'.sph'],chlocs,length(chlocs),10,10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% use a different elpfile
badchan{30} = { 'A12'  'A13'  'A14'  'A20'  'A21'  'A22'  'A27'  'A32'  'A6'  'A7'  'B13'  'B14'  'B17'  'B20'  'B26'  'C16'  'C19'  'C26'  'C7'  'D1'  'D11'  'D23'  'D25'  'D4'  'D5'  'E1'  'E2'  'E3'  'E32'  'EXG1'  'EXG4'  'EXG5'  'F1'  'F11'  'F14'  'F18'  'F31'  'G11'  'G23'  'G5'  'G9'  'H10'  'H24'  'H7'  'H9' 'C5'};
badchan{33} = { 'A1'  'A18'  'A31'  'A6'  'A7'  'C18'  'C8'  'D1'  'D13'  'D17'  'D18'  'D19'  'D24'  'D25'  'E19'  'E20'  'E23'  'E24'  'E26'  'E29'  'E3'  'E30'  'E5'  'E7'  'EXG1'  'EXG4'  'EXG5'  'F23'  'F24'  'F29'  'F30'  'F32'  'G10'  'G11'  'G13'  'G23'  'G27'  'H11'  'H16'  'H18'  'H19' 'H20'  'H21'  'H23'  'H24' 'F4' 'F22' 'F28' 'C16' 'F15' 'H9' 'H7' 'B19' 'E4' 'F19' 'D20'};
nx = 33; % or 30,33
clear datset
for bl = 1:nsets(nx)
    datset = [savename,'-',int2str(bl),'.set'];
    EEG = pop_loadset(datset ,fullpaths{nx});
    dels = [];
    for ch = 1:length(EEG.chanlocs)
        if ~isempty(find(strcmp(EEG.chanlocs(ch).labels,badchan{nx})))
            dels = [dels ch];
        end;
    end
    EEG = pop_select(EEG,'nochannel',dels);    
    sph=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.sph'],[nchans(nx) nchans(nx)],[],0); 
    wts=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.wts'],[nchans(nx) nchans(nx)],[],0); 
    EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[]; EEG.icaact = [];
    EEG = eeg_checkset(EEG);
        
    datset = [savename,'-',int2str(bl),'-',int2str(size(EEG.data,1)),'.set'];
    EEG = pop_saveset( EEG,datset,fullpaths{nx} ); 
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check out comps and choose eye blink comps
 load /data/common4/emotion/complists.mat complist15 complist20
savename = 'Emo-HP';wtstem = 'EMO';ds=2;
savename = 'NoFilt';fullpaths=newpaths;
for nx = 1:length(fullpaths)
    EEG = pop_loadset('sources.set' ,newpaths{nx});
    df = EEG.dipfit;
  for em = 1:length(emos)
    %clear datsets
    %for bl = 1:nsets(nx)
    %    datsets{bl} = [wtstem,'-',int2str(nchans(nx)),'-',int2str(bl),'.set'];
    %end;
    ALLEEG=[];EEG=[];
    EEG = pop_loadset([emos{em},'_NF.set'] ,[newdir,fullpaths{nx}(end-5:end)]); %EEG.dipfit = df;
    %EEG = pop_loadset(datsets{ds} ,fullpaths{nx}); %EEG.dipfit = df;
    sph=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.sph'],[nchans(nx) nchans(nx)],[],0); 
    wts=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.wts'],[nchans(nx) nchans(nx)],[],0); 
    EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[]; EEG.icaact = [];
    EEG = eeg_checkset(EEG);
    frlim = [0 70]; 
    figure; pl=1;  row = 10; col = 9;
    for c = 1:30%length(complist15{nx})
        %cc = complist15{nx}(c);
        cc=c;
        if c == 41 | c == 81
            textsc([fullpaths{nx}(end-4:end),' Emotion ICs/Spectra'],'title');axcopy
            figure; pl=1;
        end;        
        sbplot(row,col,pl);
        topoplot(EEG.icawinv(:,cc),EEG.chanlocs(1:nchans(nx)),  'electrodes', 'off');
        title(int2str(cc));    
        
        sbplot(row,col,[pl+1 pl+2]);     
        [Pxx freqs] = pwelch(EEG.icaact(cc,:),EEG.srate*2,EEG.srate,EEG.srate*3,EEG.srate);hold on; 
        Pxx = 10*log10(Pxx);
        ph = plot(freqs(find(freqs>frlim(1)&freqs<frlim(2))),Pxx(find(freqs>frlim(1)&freqs<frlim(2))),'b-','linewidth',2);
        %title(['rv ',num2str(round(EEG.dipfit.model(cc).rv*100)/100)]);    
        title(['Cp ',int2str(cc)]);    
        set(gca,'xtick',[0:10:frlim(2)]);
        set(gca,'xlim',frlim);    set(gca,'xgrid','on');    pl=pl+3;
        plot([10 10],[get(gca,'ylim')],'r-');
    end;    
    textsc([fullpaths{nx}(end-4:end),' Emotion ICs/Spectra'],'title');axcopy
    gdcomps{nx} = input('List Good Comps: ');
    %b{nx} = input('Eye blink component number: ');
    %l{nx} = input('Lateral eye mvment number: ');
    %m{nx} = input('Clear muscle components: ');
    h{nx} = input('Clear pulse component: ');
    close; close;
end;
pop_topoplot(EEG,0, [4,13,17,20] ,[fullpaths{nx}(end-4:end-1),' emotion'],[] ,0, 'electrodes', 'off', 'masksurf', 'on');
save /data/common4/emotion/gdcomps.mat gdcomps 
save /data/common4/emotion/ICstuff.mat b l m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save weights into datasets if they are good:----
for ds = 1:nsets(nx)
    EEG = pop_loadset([wtstem,'-',int2str(nchans(nx)),'-',int2str(ds),'.set'] ,fullpaths{nx});
    EEG.icaweights=wts;EEG.icasphere=sph;
    EEG.icawinv=[]; EEG.icaact = [];
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset( EEG,[wtstem,'-',int2str(nchans(nx)),'-',int2str(ds),'.set'],fullpaths{nx} );    
end
%------------------------------
% find dipoles and complist:---------------------------------
transform = [0 0 0 0 0 0 85 85 85];% spherical
ALLEEG=[];EEG=[];[complist] = FitDipoles(datsets{ds},fullpaths{nx},'sources.set','sphere',transform,0);
EEG = pop_loadset('sources.set' ,fullpaths{nx}); 
[gdcomps{nx}] = AutoICselect(EEG,15,-25,-45,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call in data and weights-- and save with blink comp(s)
load /data/common4/emotion/ICstuff.mat b l m
load /data/common4/emotion/gdcomps.mat gdcomps 
% only need to do [33]; when they are done!
for nx = 1:35
    clear datsets
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    ALLEEG=[];EEG=[];
    %sph=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.sph'],[nchans(nx) nchans(nx)],[],0); 
    %wts=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.wts'],[nchans(nx) nchans(nx)],[],0); 
    for ds = 1:length(datsets)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset(datsets{ds} ,fullpaths{nx});
        %EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[]; EEG.icaact = [];
        %EEG = eeg_checkset(EEG);
        EEG.gdcomps = gc;
        %EEG.gdcomps = gdcomps{nx};
        %EEG.blink = b{nx};
        %EEG.lateyes = l{nx};
        %EEG.muscle = m{nx};
        EEG = pop_saveset( EEG,datsets{ds},fullpaths{nx} ); 
    end;
    for em = 1:length(emos)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset([emos{em},'.set'] ,newpaths{nx});
        EEG.gdcomps = gc;
        EEG = pop_saveset( EEG,[emos{em},'.set'] ,newpaths{nx} ); 
    end;
end;
frlim = [3 128];
for nx = 1:35
    ALLEEG=[];EEG=[];
    EEG = pop_loadset('sources.set' ,newpaths{nx});
    muscles = EEG.muscle;
    if length(muscles) > 5
    dipmod = EEG.dipfit.model;
    EEG = pop_loadset('anger.set' ,newpaths{nx});
    figure; row = round(sqrt(length(muscles)));col = ceil(sqrt(length(muscles))); pl =1;
    for cc = 1:length(muscles)
        [Pxx freqs] = pwelch(EEG.icaact(muscles(cc),:),EEG.srate*2,EEG.srate,EEG.srate*3,EEG.srate);
        Pxx = 10*log10(Pxx);
        sbplot(row,col,pl); pl = pl+1;
        topoplot(EEG.icawinv(:,muscles(cc)),EEG.chanlocs,  'electrodes', 'off');
        title(int2str(muscles(cc)));    
        
        sbplot(row,col,pl); pl = pl+1;
        ph = plot(freqs(find(freqs>frlim(1)&freqs<frlim(2))),Pxx(find(freqs>frlim(1)&freqs<frlim(2))),'b-','linewidth',2);hold on;
         set(gca,'xtick',[0:10:frlim(2)]);set(gca,'xticklabel',[]);
        set(gca,'xlim',frlim); set(gca,'xscale','log');  
        plot([10 10],[get(gca,'ylim')],'r-');
       title(['rv ',num2str(round(dipmod(muscles(cc)).rv*100)/100)]);    
    end;
    ALLEEG=[];EEG=[];    EEG = pop_loadset('sources.set' ,newpaths{nx});
    EEG.muscle
    EEG.muscle = input('new muscles: ');
    EEG = pop_saveset( EEG,'sources.set',newpaths{nx} ); 
    close
    end;
end;
SpecCoModPlot('sources.set',newpaths{nx},[],[1:14],savedat,[3 125],'n',0,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
EEG = pop_loadset('sources.set' ,newpaths{nx});
pop_dipplot( EEG,[21,26], 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.2/standard_BESA/avg152t1.mat', 'projimg', 'on', 'projlines', 'on', 'normlen', 'on');

    pop_dipplot( EEG,35, 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.2/standard_BEM/standard_mri.mat', 'projimg', 'on', 'projlines', 'on', 'normlen', 'on');

    SpecCoModPlot('sources.set',newpaths{nx},[],[1:15],savedat,[3 125],'n',0,[]);

    
% check to see that gdcomps dipoles are in the head
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nx = 1:35
    ALLEEG=[];EEG=[];
    EEG = pop_loadset('sources.set' ,newpaths{nx});
    %pop_dipplot( EEG,EEG.blink , 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.1/standard_BESA/avg152t1.mat', 'projlines', 'on', 'normlen', 'on');view(0,90)
figure;pop_topoplot(EEG,0, EEG.blink , fullpaths{nx}(end-4:end-1),[] ,0, 'electrodes', 'off', 'masksurf', 'on');
eeglab redraw
    %view(50,30)
    %gc = EEG.gdcomps;
    %EEG.ventfrontal = input('inf frontal comp: ');
    EEG = pop_saveset( EEG,'sources.set',newpaths{nx} ); 
    close;close;
end;
    clear datsets
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    for ds = 1:length(datsets)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset(datsets{ds} ,fullpaths{nx});
        EEG.gdcomps = gc;
        EEG = pop_saveset( EEG,datsets{ds},fullpaths{nx} ); 
    end;
    for em = 1:length(emos)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset([emos{em},'.set'] ,newpaths{nx});
        EEG.gdcomps = gc;
        EEG = pop_saveset( EEG,[emos{em},'.set'] ,newpaths{nx} ); 
    end;
end;
eeglab redraw

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert press events:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
butpeeps = [1:2,4:8,10:33];% nx=9 has no button
for nxx = 29:30%length(butpeeps) 
    nx = butpeeps(nxx);
    searchnext = 0;    
    clear datsets
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    for ds = 1:length(datsets)
        ALLEEG=[];EEG=[];
        EEG = pop_loadset(datsets{ds} ,fullpaths{nx}); clrevs = [];
        for ev = 1:length(EEG.event)
            if strcmp(EEG.event(ev).type,'press') | strcmp(EEG.event(ev).type,'press1')
                clrevs = [clrevs ev];
            end;
        end;
        EEG = pop_editeventvals(EEG, 'delete',clrevs);
        EEG = eeg_checkset(EEG,'eventconsistency');
        [EEG,searchnext] = GetPresses(EEG,searchnext);
        EEG = pop_saveset( EEG,datsets{ds},fullpaths{nx} ); 
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try a decomp on back-proj of lower comps (no improvement on later ICs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=4; subcp = 20;
savename = 'Emo-HP';
for bl = 1:nsets(nx)
    datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
end;
for ds = 1:length(datsets)
    ALLEEG=[];EEG=[];
    EEG = pop_loadset(datsets{ds} ,fullpaths{nx});
    sph=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.sph'],[nchans(nx) nchans(nx)],[],0); 
    wts=floatread([fullpaths{nx},wtstem,int2str(nchans(nx)),'.wts'],[nchans(nx) nchans(nx)],[],0); 
    EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[]; EEG.icaact = [];
    EEG = eeg_checkset(EEG);
    EEG = pop_subcomp( EEG, [1:subcp], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on'); 
    EEG.icaweights=[];EEG.icasphere=[];EEG.icawinv=[]; EEG.icaact = [];
    EEG = pop_saveset( EEG,[datsets{ds}(1:end-8),'-subcomp.set'],fullpaths{nx} ); 
end;
% run ICA on result:
ALLEEG=[];EEG=[];
wtstem = 'SUB';rejthresh = [1000 5]; clear datsets
for bl = 1:nsets(nx)
    datsets{bl} = [savename,'-',int2str(bl),'-subcomp.set'];
end;
[EEG,numframes,nchan] = MakeFloats(datsets,fullpaths{nx},wtstem,0,[],rejthresh); 
pcs = nchans(nx) - subcp;
LinuxICA(wtstem,fullpaths{nx},sum(numframes),nchan,pcs,1);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find dipoles
   for ic = 1:length(EEG.dipfit.model)
       xx(1,ic) =EEG.dipfit.model(ic).rv;
   end;
   complist15{nx} = find(xx < .15);
   complist20{nx} = find(xx < .2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds=2; % nx1:2,4:20 done with co-reg (18:19 needs optimization)
load /data/common4/emotion/complists.mat complist15 complist20
dnsubj = [1:19,27:35];% DONE with manual coreg
[20:26]; % done

for nxx = 1:length(dnsubj)
    nx = dnsubj(nxx);
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end; ALLEEG=[];EEG=[];
    EEG = pop_loadset(datsets{ds}, fullpaths{nx});
    eeglab redraw; 
    % co-register electrodes with model by hand*********************  
    EEG = pop_saveset( EEG,'sources.set' ,fullpaths{nx}); 
    % now run autofit (this cuts down channels after using for spectra)
    ALLEEG=[];EEG=[];
    [complist15{nx} complist20{nx}] = FitDipoles('sources.set',fullpaths{nx});
end;
 save /data/common4/emotion/complists.mat complist15 complist20
 
% you can use this complist to find 'gdcomps', they are all comps with < 15% rv

EEG = pop_loadset('sources.set', '/data/common2/emotion/js75/');
pop_topoplot(EEG,0, [3,5,7,9,13,15,18,19,22,24] ,'',[] ,0, 'electrodes', 'off', 'masksurf', 'on');
figure; eegplot(EEG.icaact([3,5,7,9,13,15,18,19,22,24],:),'spacing',12,'winlength',20, 'dispchans',50)
figure; eegplot(EEG.icaact([1:20],:),'spacing',12,'winlength',20, 'dispchans',50)

nx=33;
ALLEEG = pop_delset( ALLEEG, [1] );
EEG = pop_loadset( 'filename', 'Emo-HP-7-189.set', 'filepath', fullpaths{nx});
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
ALLEEG = pop_delset( ALLEEG, [1] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create emotion datasets from press1 to exit
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[emoorders] = FindEmoOrder(fullpaths,emos)

    emos = [emos 'prebase', 'postbase'];
savename = 'Emo-HP';wtstem = 'EMO';
% subjs 9,34,35 without button presses had 'press1' inserted 10% through time
{EEG.event.type}
ev1 = 5;   ev2 = 6;
x=.1*(EEG.event(ev2).latency - EEG.event(ev1).latency);
EEG.event(end+1).type = 'press1';
EEG.event(end).latency = EEG.event(ev1).latency+x;
EEG = eeg_checkset(EEG,'eventconsistency');
{EEG.event.type}
EEG = pop_saveset( EEG,tmpdatsets{ds},fullpaths{nx} ); 
    clear tmpdatsets
    for bl = 1:nsets(nx)
        tmpdatsets{bl} = ['tmp',savename,'-',int2str(bl),'-223.set'];
    end;
    
for nx = 1:16%length(fullpaths)
    clear datsets
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    ALLEEG=[];EEG=[];
    for ds = 2:length(datsets)
        EEG = pop_loadset( tmpdatsets{ds},fullpaths{nx});
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );

        concat = zeros(1,length(emos));
        stev= zeros(1,length(emos));
        sptm= zeros(1,length(emos));
        spev= zeros(1,length(emos));
        for em = 1:length(emos)
            oops = 1;
            for ev = 1:length(EEG.event)
                if strcmp(emos{em},EEG.event(ev).type)
                    emevent(1,em) = ev; 
                    if em < 16
                        if ev == length(EEG.event)                           
                            oops = 0;  % if emotion is last event
                        else                            
                            for evv = ev+1:length(EEG.event)
                                if strcmp('press1',EEG.event(evv).type)
                                    stev(1,em) = evv; oops = 1;
                                    if evv == length(EEG.event)
                                        oops = 0;
                                    else 
                                        for evvv = evv+1:length(EEG.event)
                                            if strcmp('exit',EEG.event(evvv).type)
                                                sptm(1,em) = (EEG.event(evvv).latency-EEG.event(evv).latency)/EEG.srate;
                                                spev(1,em) = evvv; oops = 1; 
                                                break;
                                            end;
                                        end;  
                                        break;
                                    end;
                                end;
                            end; 
                            break;
                        end;
                    else
                        if ev == length(EEG.event) 
                            oops = 0;
                        else 
                            stev(1,em) = ev; sptm(1,em) = 119;% pre/post-baseline, 2 min
                        end;
                    end;
                end;
            end;     
            if oops == 0
                concat(1,em) = 1;
            end;   
            if concat(1,em) == 0 & stev(1,em) ~= 0
                EEG = pop_epoch( EEG,{} , [-1 sptm(1,em)+.1],'eventindices',stev(1,em), 'newname', [fullpaths{nx}(end-4:end-1),'-',emos{em}]);        
                EEG = eeg_checkset( EEG );
                EEG = pop_rmbase( EEG, [], [1:EEG.pnts] );
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
                fprintf('\nSaving dataset for %s\n',emos{em});
                EEG = pop_saveset( EEG,[emos{em},'.set'],fullpaths{nx} ); 
                ALLEEG = pop_delset( ALLEEG, [2] );
                [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
            end;    
        end;
        if ~isempty(find(concat))
            EEG = pop_select(EEG,'point',[EEG.event(emevent(find(concat))).latency-512 size(EEG.data,2)]);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
            EEG = pop_loadset( tmpdatsets{ds+1},fullpaths{nx});
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
            [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
            EEG = pop_mergeset( ALLEEG, [1  2], 1);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on'); 
            % remove boundary events and rename numbered events
            for ev = length(EEG.event):-1:1   % done
                if EEG.event(ev).type(1) == 'b'
                    EEG.event(ev) = [];
                end;
            end; 
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
            em = find(concat);
            concat = zeros(1,length(emos));
            stev= zeros(1,length(emos));
            sptm= zeros(1,length(emos));
            spev= zeros(1,length(emos));
            oops = 1;
            for ev = 1:length(EEG.event)
                if strcmp(emos{em},EEG.event(ev).type)
                    emevent(1,em) = ev; 
                    if em < 16
                        if ev == length(EEG.event)                           
                            oops = 0;  % if emotion is last event
                        else                            
                            for evv = ev+1:length(EEG.event)
                                if strcmp('press1',EEG.event(evv).type)
                                    stev(1,em) = evv; oops = 1;
                                    if evv == length(EEG.event)
                                        oops = 0;
                                    else 
                                        for evvv = evv+1:length(EEG.event)
                                            if strcmp('exit',EEG.event(evvv).type)
                                                sptm(1,em) = (EEG.event(evvv).latency-EEG.event(evv).latency)/EEG.srate;
                                                spev(1,em) = evvv; oops = 1; 
                                                break;
                                            end;
                                        end;  
                                        break;
                                    end;
                                end;
                            end; 
                            break;
                        end;
                    else
                        stev(1,em) = ev; sptm(1,em) = 119;% pre/post-baseline, 2 min
                    end;
                end;
            end;  oops      
            if oops == 0
                concat(1,em) = 1;
            end;
            if concat(1,em) == 0 & stev(1,em) ~= 0
                EEG = pop_epoch( EEG,{} , [-1 sptm(1,em)+1],'eventindices',stev(1,em), 'newname', [fullpaths{nx}(end-4:end-1),'-',emos{em}]);            
                EEG = eeg_checkset( EEG );
                EEG = pop_rmbase( EEG, [], [1:EEG.pnts] );
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
                fprintf('\nSaving dataset for %s\n',emos{em});
                EEG = pop_saveset( EEG,[emos{em},'.set'],fullpaths{nx} ); 
                ALLEEG = pop_delset( ALLEEG, [1 2] );
            else
                fprintf('\nERROR: could not complete dataset for emotion %s\n',emos{em});
            end;    
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Run spectral modulator decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssh juggling
qlogin
matlab
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab file loads all subject info needed
newdir = '/data/projects/julieo/emotion/';
fullpaths = cell(1,length(newpaths));
for nx =1:length(newpaths)
  fullpaths{nx} = [newdir,newpaths{nx}(end-4:end)];
end;
fullpaths{22} = [fullpaths{22}(1:30),'mr74/'];


datset = {'anger_NF.set','frustration_NF.set','jealousy_NF.set','fear_NF.set' ,'disgust_NF.set','grief_NF.set','sad_NF.set','compassion_NF.set','love_NF.set','relief_NF.set','content_NF.set','awe_NF.set','happy_NF.set','joy_NF.set','excite_NF.set'}; % for all new ones


datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set','prebase.set','postbase.set'}; % for all new ones

savedat = 'SpecCoModNoFiltWave'; frqlim = [2 128];overlap = 1;pcfac = 2; nfreqs = 370;freqscale = 'quad'; wsize = 1; pwrdecomp=[frqlim(1) frqlim(2) 256*wsize+1];% wsize=epsize in sec
savedat = 'SpecCoModNoFiltTest2'; frqlim = [.5 80];overlap = 1;pcfac = 2; nfreqs = 500;freqscale = 'quad'; %mi83,3 sec
savedat = 'SpecCoModNoFiltTest'; frqlim = [1.4 128];overlap = 1;pcfac = 2; nfreqs = 370;freqscale = 'quad'; wsize = 2;
savedat = 'SpecCoModNoFilt'; frqlim = [1.4 128];overlap = 1;pcfac = 2; nfreqs = 400;freqscale = 'quad'; wsize = 2;
savedat = 'SpecCoModMuscle'; frqlim = [3 125];overlap = 2;pcfac = 2; nfreqs = 370;freqscale = 'quad'; 
savedat = 'SpecCoModMoreFreqs'; frqlim = [3 125];overlap = 2;pcfac = 2;nfreqs = 370;freqscale = 'quad';
savedat = 'SpecCoModFFT'; frqlim = [3 125];overlap = 2;pcfac = 2;nfreqs = 370;freqscale = 'quad';
savedat = 'SpecCoModNoOvrlap'; frqlim = [3 125];overlap = 1;pcfac = 2;nfreqs = 370;freqscale = 'quad';
savedat = 'SpecCoModMoreFreqsmovav'; frqlim = [3 125];overlap = 8;pcfac = 2;nfreqs = 370;freqscale = 'quad'; smth = 6;

savedat = 'SpecCoModWaveTest';frqlim = [3 125];overlap = 2;pcfac = 2;nfreqs = 370; freqscale = 'quad'; % wavelet = [3 6] ?
savedat = 'SpecCoModWaveTest2'; frqlim = [6 125]; % all wavelet = 6
savedat = 'SpecCoModWaveTest3'; frqlim = [6 125]; % wavelet = [6 12]
savedat = 'SpecCoModWaveTest4'; frqlim = [3 125]; % wavelet = [3 3]
savedat = 'SpecCoModWaveTest5'; frqlim = [8 125]; % wavelet = [8 8]
savedat = 'SpecCoModWaveTest6'; frqlim = [30 125]; % wavelet = [10 10]
savedat = 'SpecCoModWaveTest7'; frqlim = [30 125]; % wavelet = [20 20]
savedat = 'SpecCoModWaveTest8'; frqlim = [3 35]; % wavelet = [6 6] (win 513)
savedat = 'SpecCoModWaveTest9'; frqlim = [30 125]; % wavelet = [50 50]

savedat = 'SpecCoModWaveTest10'; frqlim = [30 125]; % wavelet = [100 100]
savedat = 'SpecCoModWaveTest11'; frqlim = [3 35]; % wavelet = [9 12] (win 641: 2.5 sec)
savedat = 'SpecCoModWaveTest12'; frqlim = [3 125]; % wavelet = [9 12] (win 641: 2.5 sec)

savedat = 'SpecCoModLongWave'; frqlim = [3 125]; pcfac = 100; freqscale = 'linear'; nfreqs = 200;overlap = 2;% constant number of pcs% wavelet
savedat = 'SpecCoModWave';frqlim = [3 125];pcfac = 2;freqscale = 'quad';nfreqs = 370;overlap = 2; % wavelet [3 125]
savedat = 'SpecCoModWavePrePost';frqlim = [3 125];pcfac = 2;freqscale = 'quad';nfreqs = 370;overlap = 2; % wavelet [3 125]
savedat = 'SpecCoModTight';frqlim = [3 125];pcfac = 75;freqscale = 'quad';nfreqs = 300;overlap = 70; % wavelet [3 125]
 savedat = ['mi83decimationtest',int2str(dec)]; frqlim = [3 125];overlap = 2;nfreqs = 370;freqscale = 'quad';
 smth = 4;savedat = ['mi83decimationtestmovav',int2str(smth)]; frqlim = [3 125];overlap = 2;nfreqs = 370;freqscale = 'quad';% 15 datsets

nx=2; comps = gdcomps{nx};
dec = [];
%for dec = 2:10:42
%for dec = 1.1:.4:6
    savedat = ['mi83decimationtest',num2str(dec)];pcfac = 2; % results in 54 PCs
     savedat = ['mi83decimationtest',int2str(dec),'PC25'];pcfac = 25;
     savedat = ['mi83decimationtest',num2str(dec),'PC100'];pcfac = 100;
     savedat = ['mi83decimationDblOvrlp',int2str(dec)]; overlap = 4; pcfac = 100;datset = {'anger.set','jealousy.set','sad.set','love.set','content.set','awe.set','excite.set'};% half datasets
     savedat = ['mi83decimationTrplOvrlp',int2str(dec)]; overlap = 8;pcfac = 100; datset = {'jealousy.set','sad.set','awe.set','excite.set'};% quarter datasets
     savedat = ['mi83decimation4set',int2str(dec)]; overlap = 2;pcfac = 100; datset = {'jealousy.set','sad.set','awe.set','excite.set'};% quarter datasets
     savedat = ['mi83decimation2set',int2str(dec)]; overlap = 2;pcfac = 2; datset = {'jealousy.set','awe.set'};% quarter datasets
     savedat = ['mi832setOvrlp8',int2str(dec)]; overlap = 8;pcfac = 2; datset = {'jealousy.set','awe.set'};% quarter datasets
     savedat = ['mi83decimation7set',int2str(dec)]; overlap = 2;pcfac = 100; datset = {'anger.set','jealousy.set','sad.set','love.set','content.set','awe.set','excite.set'};% half datasets
    SpecCoModAnal(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,[frqlim 257],dec); % wavelet
end;
savedat = 'mi83ChannelTest';frqlim = [3 125];overlap = 2;pcfac = 2;nfreqs = 370;freqscale = 'quad';
freqs = linspace(.01, 5, 50); 
dec=[];nogo = []; % also 5,20 is messy, low # trials
for nx =1:35 %7 no disgust;  11?    15    23 didn't run  
  fprintf('\nLoading subject %s.\n',int2str(nx));
  %EEG = pop_loadset('sources.set' ,newpaths{nx});  
  %comps = [EEG.gdcomps EEG.ventfrontal EEG.muscle];
  comps = gdcomps{nx};
    %ALLEEG=[];EEG=[];
    try
      %clear allspec
      %for ds = 1:length(datset)
      %[allspec(:,:,ds), freqs,comment] = CalcSpectra(datset{ds},['/data/projects/julieo/emotion/',newpaths{nx}(end-4:end)],gdcomps{nx},[],freqs,[],[],'db','raw');
      %end;
      %subjspec{nx} = allspec;
    SpecCoModAnal(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,pwrdecomp); % wavelet
    %SpecCoModAnal(datset,['/data/projects/julieo/emotion/',fullpaths{nx}(end-4:end)],comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,pwrdecomp); % wavelet
    %SpecCoModAnal(datset,['/data/projects/julieo/emotion/',newpaths{nx}(end-4:end)],comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,wsize); % FFT
    %SpecCoModAnal(datset,fullpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,wsize); % FFT
    catch
      nogo = [nogo,nx];
      end;
end;
nogo
CoModDecomp(fullpaths{nx},savedat,savedat ,[],pcfac,[],'raw',{256,freqscale,1,frqlim(2),nfreqs},[frqlim(1) frqlim(2) 257]); %  

save /data/projects/julieo/emotion/LowFreqSpecs.mat subjspec freqs
lfemos = [];subjinfo = [];
for nx=1:35
  if isempty(find(subjspec{nx}==0))
  for ic = 1:size(subjspec{nx},1)
    lfemos = [lfemos;squeeze(mean(subjspec{nx}(ic,:,:),2))'];
    subjinfo = [subjinfo;[nx ic]];
  end
  end
end;
dd = pdist(lfemos,'correlation');
links = linkage(dd,'complete');
clss = cluster(links,'cutoff',.5);
[H, T] = dendrogram(links,6);
for r=1:length(T)
  cls{r} = subjinfo(find(T==r),:);
  clspattern{r} = lfemos(find(T==r),:);
end;

    %SpecCoModAnal(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,1); % FFT
    %SpecCoModAnal(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,[frqlim 257],dec); % wavelet
    %SpecCoModAnalSmooth(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,[],nfreqs,[frqlim 257],[],smth); % wavelet

    %CoModDecomp(newpaths{nx},'SpecCoMod',savedat ,[],pcfac,[],'raw',{256,freqscale,frqlim(1),frqlim(2),nfreqs},0); % FFT
    %CoModDecomp(newpaths{nx},'SpecCoMod',savedat ,[],pcfac,[],'raw',{256,freqscale,frqlim(1),frqlim(2),nfreqs},[frqlim(1) frqlim(2)-1 129]); % high freq decomp 
   % CoModDecomp(newpaths{nx},'SpecCoMod',savedat ,[],pcfac,[],'raw',{256,freqscale,frqlim(1),frqlim(2),nfreqs},[1.3 frqlim(2) 129]); % all freq wavelet (timefreqs)
    %CoModDecomp(newpaths{nx},'SpecCoModWaveTest11',savedat ,[],pcfac,[],'raw',{256,freqscale,frqlim(1),frqlim(2),nfreqs},[3 6 129]); % Rey's wavelet 
end;
%%%% Plot the results:--------
s = load([fullpaths{nx},savedat,'.mat']);
wts = floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);
sph = floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0);
ws = wts*sph;winv = pinv(ws);
speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) inf],[],0);
specwts = speceig*winv; winv = specwts;     
data = floatread([fullpaths{nx},savedat,'.fdt'],[s.numrows s.numframes],[],0); acts = ws*data;
data = floatread([fullpaths{nx},savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);
alldat = floatread([fullpaths{nx},savedat,'RAW.fdt'],[length(s.complist) 256*wsize length(s.rowmeans)],[],0);

for nx = subjlist
SpecCoModPlot('awe_NF.set',fullpaths{nx},gdcomps{nx},[],savedat,frqlim,'ims',0,[],'off',[]);
end;
figure; [newpl] = PlotIMweights(newpaths{nx},savedat,emos,[],[1:24],'n',3,4,1)

SpecCoModPlot('sources.set',newpaths{nx},gdcomps{nx},[],savedat,frqlim,'n',0,[],'off',[]);

SpecCoModPCA('awe_NF.set',fullpaths{nx},savedat,[],[]); % just the pca templates

for nx = subjlist
    str = ['load ',newpaths{nx},savedat,'.mat'];eval(str);  
    rowmeans = zeros(numtrials,1); numrows = pcs;
    str = ['save ',newpaths{nx},savedat,'.mat meanpwr numrows numframes freqs freqscale keeptrack rmepochs dstrials  pcs rowmeans complist overlap datset eigfile'];eval(str);    
%SpecCoModPlot('sources.set',newpaths{nx},[],[],savedat,[35 128],'n',0,[]);
end;
tmpSpecCoModPlot('sources.set',newpaths{nx},[],[],savedat,[3 125],'n',0,[]);
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/ChanDecomp/mi83pg1.jpg -djpeg
SpecCoModPlot('sources.set',newpaths{nx},gdcomps{nx},[1:30],savedat,frqlim,'ims',0,[],'off',[]);
SpecCoModPlot('sources.set',newpaths{nx},[],[],savedat,[3 125],'ims',0,[]);
SpecCoModPlot('sources.set',newpaths{nx},[],[31:45],savedat,[3 125],'n',0,[]);
SpecCoModPlot('sources.set',newpaths{nx},ic,im,savedat,[3 125],'n',0,[]);

str = ['print /data/common2/emotion/Figs/GammaIMacts2.eps -depsc -painters -adobe'];eval(str)
str = ['print /data/common2/emotion/Figs/GammaIMacts2.jpg -djpeg'];eval(str)
templs = 22;
PlotSpecFacEnv('sources.set',savedat,newpaths{nx},facs,cps,[],[],frqlim,1,.99,0);
textsc(['IM ',int2str(templs)],'title');
str = ['print /data/common2/emotion/Figs/GammaAlpha3.eps -depsc -painters -adobe'];eval(str)
set(gca,'ylim',[-52 -24]);
set(gca,'ylim',[-55 -24]);

%for dec = 2:10:42
for dec = 1.2:.2:6
     %savedat = ['mi83decimationtest',int2str(dec),'PC25'];
     savedat = ['mi83decimationtest',int2str(dec)];
     %savedat = ['mi83decimationtest',int2str(dec),'PC100'];
     %savedat = ['mi83decimationDblOvrlp',int2str(dec)]; 
     SpecCoModPlot('sources.set',newpaths{nx},[],[39:45],savedat,[3 125],'n',0,[]);
SpecCoModPlot('sources.set',newpaths{nx},gdcomps{nx},[],savedat,[3 125],'ims',0,[],'off',[]);
SpecCoModPlot('sources.set',newpaths{nx},ic,im,savedat,[3 125],'ims',0,[],'off',[]);
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83dec12IMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83nodecIMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83nodecset4IMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83nodecset4DblOvrlpIMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83nodecset4TrplOvrlpIMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83-2setOvrlp50IMs.jpg -djpeg'];eval(str)
     str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83-2setOvrlp87IMs.jpg -djpeg'];eval(str)
     str = ['print /data/common1/emotion/mi83/',savedat,'.eps -depsc2'];eval(str)
     close all
end;

s = load([newpaths{nx},savedat,'.mat']);
wts = floatread([newpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);
sph = floatread([newpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0);
ws = wts*sph;winv = pinv(ws);
speceig = floatread([newpaths{nx},s.eigfile],[length(s.rowmeans) inf],[],0);
specwts = speceig*winv; winv = specwts;     
data = floatread([newpaths{nx},savedat,'DAT.fdt'],[s.numrows s.numframes],[],0);
alldat = floatread([newpaths{nx},savedat,'RAW.fdt'],[length(s.complist) srate*2 length(s.rowmeans)],[],0);
alldat = floatread([newpaths{nx},savedat,'RAW.fdt'],[length(s.rowmeans) length(s.complist)*length(s.freqs)],[],0);

[pc,eigvec,sv] = runpca(data,100);
figure; plot(max(sv))

ntrials = 125;pcfac = 1;auxpath = [];
savedat = 'SpecHPModEmos'; 
CoModDecomp(newpaths{nx},savedat, [savedat,'SUB'],ntrials,pcfac,auxpath,'spec',{});
SpecCoModPlot('sources.set',newpaths{nx},[],[],[savedat],[3 128],'ims',0,[]);

ims = [1,4]; comp = [2,5];frqlim=[3 125];
SpecCoModPlot('sources.set',newpaths{nx},comp,ims,savedat,frqlim,'ims',0,[],[],[]);ph=textsc(['Subject mi83'],'title');% templates

SpecCoModPlot('sources.set',newpaths{nx},comp,ims,savedat,frqlim,'off',1,[],[],[]);ph=textsc(['Subject mi83'],'title');% back proj    

print /home/julie/Talks/mi83IMs4Movies.eps -depsc2 -adobe -painters
[pcared pv] = SpecCoModPvaf(savedat,newpaths{nx},0,[],1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare decomp tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 2;pl=1; clear corr
figure; row = 6; col = 6;
for dec = 12:10:32
for dec = 1.1:.4:6
%for nx = 1:35
    %savedat1 = 'SpecCoModWave';
    savedat1 = 'SpecCoModMoreFreqs';
    %savedat1 = ['mi83decimation4set',int2str(dec)]; 
    %savedat2 = ['mi83decimation7set',int2str(dec)]; 
    savedat2 = ['mi83decimationtest',num2str(dec)];% pc 54
    %savedat2 = ['mi83decimationtest',int2str(dec),'PC25'];
    %savedat2 = ['mi83decimationtest',num2str(dec),'PC100'];
    %savedat2 = ['mi83decimationDblOvrlp',int2str(dec)];
    %savedat2 = ['mi83decimationTrplOvrlp',int2str(dec)];
    %savedat2 = 'SpecCoModMoreFreqsmovav';
    %savedat2 = 'SpecCoModWave';
    [corr(pl,:)] = CompareIMdecomps(savedat2,savedat1,newpaths{nx});pl = pl+1;
    [corr] = CompareIMdecomps(savedat1,savedat2,newpaths{nx}); mncorrs(1,nx) = mean(corr);
    %sbplot(row,col,pl); pl = pl+1;
    plot(sort(abs(corr'))); hold on;
end;
figure; cols = jet(16);set(gca,'fontsize',15);
for c = 1:size(corr,1)
    ph = plot(sort(abs(corr(c,:)')),'linewidth',2,'color',cols(c,:));   
    hold on;
end;
legend({'4494 windows','3295 windows','2602 windows','2149 windows','1831 windows','1595 windows','1412 windows','1267 windows','1150 windows','1052 windows','969 windows','899 windows','838 windows','419 windows','228 windows','157 windows'},'location','southeast');
%legend({'2511 windows','419 windows','228 windows','157 windows'},'location','northwest');
set(gca,'xlim',[1 size(corr,2)]);xlabel('IMs');ylabel('Absolute correlation between best-matching IMs');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Manuscripts/Gamma/Frontiers/mi83decCorrs.jpg -djpeg'];eval(str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many data points per subject?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nx=1:length(gdcomps)
    if ~isempty(gdcomps{nx})
        s = load([newpaths{nx},savedat,'.mat']);  
        npnts(1,nx) = length(s.rowmeans);
    end;
end;
mean(npnts(find(npnts)))
mean = 3,567

  Columns 1 through 6

        3365        5022        3056        4192         865        2878

  Columns 7 through 12

        1873        2151        9500        4613        2793        2076

  Columns 13 through 18

        1878        3369        1461        3177        2540        4435

  Columns 19 through 24

        1281        1761        3311        2207        6125        3976

  Columns 25 through 30

        3081        3949        3324        2277        4040        5069

  Columns 31 through 35

        7133        2845        4826        3785        6596
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newsave = 'TESTrunpca';
newsave = 'TESTdecomp';

allsaves = {'TESTpcfac.01','TESTpcfac.05','TESTpcfac.1','TESTpcfac.3','TESTpcfac.5','TESTpcfac1','TESTpcfac2','TESTpcfac4','TESTpcfac6'};

t = 7;
CoModDecomp(newpaths{nx},savedat, allsaves{t},[],.05,[],'raw',{256,'log',3,125,150});
SpecCoModPlot('sources.set',newpaths{nx},[],[1:14],savedat,[3 125],'n',0,[]);
SpecCoModPlot('sources.set',newpaths{nx},[],[15:29],savedat,[3 125],'n',0,[]);
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Manuscripts/Gamma/pcfactests/',newsave,'.jpg -djpeg'];  eval(str);

clear topcorr corr
for t = 1:length(allsaves)-1
    s = load([newpaths{nx},allsaves{t},'.mat']);
wts = floatread([newpaths{nx},allsaves{t},'.wts'],[s.pcs s.numtrials],[],0);
sph = floatread([newpaths{nx},allsaves{t},'.sph'],[s.numtrials s.numtrials],[],0);
ws = wts*sph;winv = pinv(ws);
icamatall = floatread([newpaths{nx},allsaves{t},'.fdt'],[s.numtrials s.numframes],[],0);
    %sph=floatread([newpaths{nx},allsaves{t},'.sph'],[s.pcs s.pcs],[],0); 
    %wts=floatread([newpaths{nx},allsaves{t},'.wts'],[s.pcs s.pcs],[],0); 
    %icamatall = floatread([newpaths{nx},allsaves{t},'.fdt'],[s.pcs s.numframes],[],0);    
    ws = wts*sph;    act1 = ws*icamatall; 
    for tt = t+1:length(allsaves)
        s = load([newpaths{nx},allsaves{tt},'.mat']);
        sph=floatread([newpaths{nx},allsaves{tt},'.sph'],[s.pcs s.pcs],[],0); 
        wts=floatread([newpaths{nx},allsaves{tt},'.wts'],[s.pcs s.pcs],[],0); 
        icamatall = floatread([newpaths{nx},allsaves{tt},'.fdt'],[s.pcs s.numframes],[],0);    
        ws = wts*sph;    act2 = ws*icamatall; 
        [corr,indx,indy,corrs] = matcorr(act2,act1);
        topcorr(t,tt,:) = corr(1:17);
    end    
end;

figure; row = round(sqrt(size(topcorr,3))); col = ceil(sqrt(size(topcorr,3)));
for t = 1:size(topcorr,3)
    sbplot(row,col,t)
    imagesc(topcorr(:,:,t),[-1 1]);
end
cbar;
    
    
%s = load([newpaths{nx},newsave,'.mat']);
%pc = floatread([newpaths{nx},newsave,'DAT.fdt'],[s.pcs s.numframes],[],0);% final pc b4 ICA
s = load([newpaths{nx},savedat,'.mat']);pcs = 50;
pwr = floatread([newpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);
    [pc1,eigvec,sv] = runpca(pwr,pcs);
    [eigvec,sv,pc2] = pcsquash(pwr,pcs); % eigvec is the wts matrix in this version
 [corr,indx,indy,corrs] = matcorr(pc1,pc2);

activations = pc;
EEG = pop_loadset('sources.set' ,newpaths{nx});  
figure; plotcomps = gdcomps{nx};
row = 12; pl = 1;
col = length(plotcomps);
for cp = 1:length(plotcomps)
    sbplot(row,col,pl)
    topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
    set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
end;
for tp = 1:11%size(activations,1)
    for cp = 1:length(plotcomps)
        if pl>row*col
            figure; pl = 1;
            for cp = 1:length(plotcomps)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
                set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
            end;
        end;
        rcp = find(plotcomps(cp) == s.complist);    
        sbplot(row,col,pl)
        if strcmp(s.freqscale,'quad') % quadratic spacing
            tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
            quadplot(s.freqs(fr),tpact(:,fr),1.75,'b');
            pl = pl+1;hold on;
        elseif strcmp(s.freqscale,'log')  % log spacing
            logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),1.75,'b'); pl = pl+1;hold on;
        else % otherwise linear
            plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',1.75); pl = pl+1;hold on;
            set(gca,'xtick',[10:10:frqlim(2)]);
            set(gca,'xticklabel',{10 [] 30 [] []});
            set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
        end;  
        set(gca,'ylim',[-.18 .18]);
        set(gca,'ylim',[-250 250]);
        set(gca,'xgrid','on');
        set(gca,'fontsize',7);set(gca,'box','off');
        set(gca,'ticklength',[.03 .03]);
        plot([get(gca,'xlim')],[0 0],'r-');
        if pl <= (row-1)*col+1
            set(gca,'xticklabel',[]);
        end; 
        if cp == 1
            title(['IM ',int2str(tp)]);
        end;
    end;
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find kurtosis of all IMs

s = load([fullpaths{nx},savedat,'.mat']);
wts = floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0);
sph = floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0);  
ws = wts*sph;winv = pinv(ws);
dat = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);
acts = ws*dat;  clear dat
for im = 1:size(acts,1)
    kurts(1,im) = kurt(acts(im,:));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run multi-IC TW spectral decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssh juggling
qlogin
matlab
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
   
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones

savedat = 'SpecCoModTW'; 
savedat = 'SpecCoModTWmusc'; 
%pcs=[];  frqlim = [3 128];freqscale = 'log';
pcs=[];  frqlim = [4 128];freqscale = 'log';
overlap = 2;pcfac = 2;
auxpath = [];
for nx = 1:35
   ALLEEG=[];EEG=[];
   EEG = pop_loadset('sources.set' ,newpaths{nx});  
   comps = [EEG.gdcomps EEG.ventfrontal EEG.muscle];
   %comps = EEG.gdcomps;
   ALLEEG=[];EEG=[];
   
   SpecCoModTW(datset,newpaths{nx},comps,savedat,frqlim,freqscale,pcfac,overlap,auxpath,65);
end;



SpecCoModTWPlot('sources.set',newpaths{nx},[],[],savedat,[0 128],'n',0,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run single-IC TW spectral decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones

savestem = 'SpecModTWic';
pcs=[];  frqlim = [3.5 125];freqscale = 'log';
overlap = 2;pcfac = 2;
auxpath = [];
for nx = 1:35
   ALLEEG=[];EEG=[];
   for ic = 1:length(gdcomps{nx})
       savedat = [savestem,int2str(gdcomps{nx}(ic))];; 
   
       SpecCoModTW(datset,newpaths{nx},gdcomps{nx}(ic),savedat,frqlim,freqscale,pcfac,overlap,auxpath,94);
   end;
end;
SpecCoModTWPlot('sources.set',newpaths{nx},[],[],savedat,[0 128],'n',0,[]);

nx=1;ic = 1;
name = [savestem,int2str(gdcomps{nx}(ic))];  
PlotTWspecPCs(name,fullpaths{nx});
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a single-channel dataset to know emo order (subjs 36 and on)
savename = 'Emo-HP';wtstem = 'EMO'; 
%%% for subjs 36 and on:
for nx = 36:length(fullpaths)
    clear datsets    
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    for ds = 1:length(datsets)
        EEG = pop_loadset( datsets{ds},fullpaths{nx},'all');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, ds );
        EEG = pop_select(EEG,'channel',1);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, ds );
    end;
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, ds, 'retrieve',1, 'study',0);
    EEG = pop_mergeset( ALLEEG, [1:length(datsets)], 0);
    % remove boundary events and rename numbered events
    for ev = length(EEG.event):-1:1   % done
        if EEG.event(ev).type(1) == 'b'
            EEG.event(ev) = [];
        end;
    end; % will exceed matrix dimensions if it takes any out
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [fullpaths{nx}(end-4:end-1),' Emotion Order']);
    %EEG = pop_saveset( EEG,'ButtonOnly.set',fullpaths{nx} ); 
    EEG = pop_saveset( EEG,'ButtonEvents.set',fullpaths{nx} ); 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a single-channel dataset to know emo order (subjs 1:35)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
butpeeps = [1:2,4:8,10:34];% nx=9,35 has no button
for nxx = 1:length(butpeeps) % nx=9 has no button
    nx = butpeeps(nxx);
    clear datsets    
    for bl = 1:nsets(nx)
        datsets{bl} = [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'];
    end;
    ALLEEG=[];EEG=[];
    for ds = 1:length(datsets)
        EEG = pop_loadset( datsets{ds},fullpaths{nx},'all');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, ds );
        EEG.icaweights=[];EEG.icasphere=[];EEG.icawinv=[]; EEG.icaact = [];
        EEG.icachansind = [];EEG.chaninfo = [];
        EEG = pop_select(EEG,'nochannel',[2:size(EEG.data,1)]);
        EEG.data(2,:) = EEG.button; EEG.nbchan = 2;
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG,'nochannel',[1]);
        EEG = rmfield(EEG,'button');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, ds );
    end;
    EEG = pop_mergeset( ALLEEG, [1:length(datsets)], 0);
    delevs = [];
    for ev = 1:length(EEG.event)   % remove boundaries
        if strcmp(EEG.event(ev).type,'boundary')
            delevs = [delevs ev];
        end;
    end; 
    EEG = pop_editeventvals(EEG, 'delete',delevs);
    EEG = eeg_checkset(EEG);    
    EEG = pop_saveset( EEG,'ButtonEvents.set',newpaths{nx} ); 
end;
%EEG = pop_loadset( 'ButtonOnly.set',fullpaths{nx});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import press events to 
butpeeps = [1:2,4:8,10:34];% nx=9 has no button, nor does 35
for nxx = 1:length(butpeeps) 
    nx = butpeeps(nxx); 
    ALLEEG=[];EEG=[];
    clear evs
    EEG = pop_loadset( 'ButtonEvents.set',newpaths{nx});
    evs = EEG.event;    datlen = size(EEG.data,2);
    EEG = pop_loadset( 'Heart2.set',newpaths{nx}); 
    if datlen == size(EEG.data,2)
        EEG.event = evs;
    end;    
    EEG = pop_saveset( EEG,'Heart2.set',newpaths{nx} );   
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

