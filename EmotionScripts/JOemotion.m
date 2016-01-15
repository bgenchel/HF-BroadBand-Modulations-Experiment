% Julie Onton emotion expt
fullpaths{1} = '/data/common4/emotion/jo74/';
datset = {'JOnew1.set','JOnew2.set','JOnew3.set','JOnew4.set'};
newset = {'Continuous1.set','Continuous2.set','Continuous3.set','Continuous4.set'};
nchan = 245;ALLEEG=[];EEG=[];
for ds = 1:length(datset)
    EEG = pop_loadset( 'filename',datset{ds} , 'filepath',fullpaths{1});
    EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
    EEG = pop_saveset( EEG,newset{ds}, fullpaths{1});
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wtstem = 'ICA';rejthresh = [1000 5];ALLEEG=[];EEG=[];
tmpstem = [wtstem,int2str(day),'_'];
[EEG,numframes,nchan] = MakeFloats(newset,fullpaths{1},wtstem,[],[],rejthresh,0);    
[run_ica_str] = LinuxICA(wtstem,['ICAinstruct'],fullpaths{1},sum(numframes),nchan,0,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_loadset( 'filename',newset{1} , 'filepath',fullpaths{1}); 
sph=floatread([fullpaths{1},wtstem,int2str(nchan),'.sph'],[nchan nchan],[],0); 
wts=floatread([fullpaths{1},wtstem,int2str(nchan),'.wts'],[nchan nchan],[],0);       
EEG.icaweights = wts; EEG.icasphere = sph;
EEG.icaacts = []; EEG.icawinv = [];
EEG = eeg_checkset(EEG);
pop_topoplot(EEG,0, [1:24] , 'jo74 - 1',[] ,0, 'electrodes', 'off');
