%%%%  script to run auto-rejection ica on all emotion subjects.

zfactor = -.055;
alllocs = {EEG.chanlocs.Z}; alllocs = cell2mat(alllocs);
lowlocs = find(alllocs<zfactor);
EEG = pop_select( EEG, 'nochannel',lowlocs );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    
numsets = [5,6,5,6,5, 5,4,4,9,7,  7,4,5,5,4,  5,5,5,4,5,  5,5,9,6,6,  6,6,6,6,7,  8,5,9,5,8];
gdchan = [241,248,238,253,250, 244,248,231,250,243, 251,251,180,249,250,  246,237,250,250,251, 250,246,248,240,254, 246,203,252,249,249, 253,236,205,247,254];

for nx = 29:length(gdcomps) 
    sph=floatread(['/data/common4/emotion/',paths{nx},sphs{nx}],sphsize{nx},[],0); 
    wts=floatread(['/data/common4/emotion/',paths{nx},wtss{nx}],wtssize{nx},[],0); 
    for ds = 1:numsets(nx)
        EEG = pop_loadset( ['emo-',int2str(ds),'-',int2str(gdchan(nx)),'.set'], ['/data/common4/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG.icaweights=wts; EEG.icasphere=sph;EEG.icawinv=[];
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_saveset(EEG, ['emo-',int2str(ds),'-',int2str(size(EEG.data,1)),'.set'],['/data/common4/emotion/',paths{nx}]);
        %EEG = pop_saveset(EEG, datsets{ds},['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        ALLEEG=[];EEG=[];
    end;
end;

eeglab
numchans = 110;
for nx = 1:length(paths)
    for ds = 1:numsets(nx)
        datset{ds} = ['emo-',int2str(ds),'-',int2str(gdchan(nx)),'.set'];
        newset{ds} = ['cut-',int2str(numchans),'-',int2str(ds),'.set'];
    end;
    
    datpath = ['/data/common4/emotion/',paths{nx}];
    zfactor = -.055;
    rmref = {'E11','G20'};
    preschan = [];rmchans = [];
    chkchans = 1;  corrsel = 1;   
    [EEG,delchans] = Cut2ChanX2(datset,datpath,numchans,zfactor,rmref,rmchans,preschan,chkchans,corrsel);

    wtstem = 'cutchan';pcs = numchans;chancomp = [];
    [EEG,numframes] = MakeFloats2(newset,datpath,wtstem,numchans,[],chancomp);
    LinuxICA2(wtstem,datpath,sum(numframes),numchans,pcs,1);
  
end;
