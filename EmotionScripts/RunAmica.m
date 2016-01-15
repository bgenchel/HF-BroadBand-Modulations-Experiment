% run Amica on emotion datasets
eeglab
addpath('/home/julie/MatlabScripts/emotion/')
DataInfo

nx = 2; 

%-----------------------------------------------------------
%% run AMICA: running nx=15, 19, 20, 6 (filtered)
%-----------------------------------------------------------
  
ALLEEG=[];EEG=[];
EEG = pop_loadset('filename','Emo-HP-Merged-232.set','filepath',fullpaths{nx});
nchan = size(EEG.data,1);
wtsphname = ['AMICA',int2str(nchan)];% 
floatwrite(EEG.data,[fullpaths{nx},wtsphname,int2str(nchan),'.fdt']);numframes = size(EEG.data,2); 
runamica([fullpaths{nx},wtsphname,int2str(nchan),'.fdt'],[fullpaths{nx},wtsphname,'/'],nchan,numframes,'qsub','on','do_reject',1,'max_iter',5000,'do_newton',1,'numrej',3,'numprocs',3,'num_models',5);%

% Call in results:
EEG = pop_loadset('filename',['Emo-HP-Merged-232.set'],'filepath',fullpaths{nx});nchan = size(EEG.data,1);
wtsphname = ['AMICA',int2str(nchan)];% 
mod = loadmodout([fullpaths{nx},wtsphname,'/']);
% transfer weights to EEG structure:
for m = 1:size(mod.W,3)
    EEG.icaweights = mod.W(:,:,m); EEG.icasphere = mod.S; EEG.icawinv = [];EEG.icaact = [];EEG.icaalgorithm = 'amica';      
    EEG = eeg_checkset(EEG); 
    winvs(:,:,m) = EEG.icawinv;

    pop_topoplot(EEG,0, [1:25] ,['Model Number ',int2str(m)] ,[5 5] ,0, 'electrodes', 'off', 'masksurf', 'on');
end;
for m = 1:size(mod.W,3)-1
    for mm = m+1:size(mod.W,3)
    [corr] = PlotWinvCorr(winvs(:,:,m),winvs(:,:,mm),EEG.chanlocs,EEG.chanlocs,96);
    textsc(['Model ',int2str(m),' vs Model ',int2str(mm)],'title');
    end;
end;

for n = 0:10000:size(mod.LLt,2)
figure; plot(mod.LLt(:,n+1:n+10000)');legend({'1','2','3','4','5'});
end;
% plot 2 decomps vs each other:
datset1 = 'AllDat.set';EEG = pop_loadset( 'filename', datset1, 'filepath', fullpaths{nx});nchan = size(EEG.data,1);
wts1 = ['amicaclean',int2str(nchan)]; 
mod = loadmodout([fullpaths{nx},wts1,'/']);
EEG.icaweights = mod.W; EEG.icasphere = mod.S; EEG.icawinv = [];EEG.icaact = [];
EEG = eeg_checkset(EEG); wv1 = EEG.icawinv; clocs1 = EEG.chanlocs;

datset2 = 'AllDat.set'; EEG = pop_loadset( 'filename', datset2, 'filepath', fullpaths{nx});nchan = size(EEG.data,1);
wts2 = ['AMICAall',int2str(nchan)];
mod = loadmodout([fullpaths{nx},wts2,'/']);
EEG.icaweights = mod.W; EEG.icasphere = mod.S; EEG.icawinv = [];EEG.icaact = [];
EEG = eeg_checkset(EEG);wv2 = EEG.icawinv; clocs2 = EEG.chanlocs;

[corr] = PlotWinvCorr(wv1,wv2,clocs1,clocs2,96);

%-------------------------------------------
% look at mod.LLt activations throughout expt
%-------------------------------------------
EEG = pop_loadset('filename',['Emo-HP-Merged-232.set'],'filepath',fullpaths{nx});nchan = size(EEG.data,1);
bigevents = [];
for ev = 1:length(EEG.event)
    if ismember(EEG.event(ev).type,emos)
        bigevents = [bigevents,ev];
    end;
end;
marktimes = cell2mat({EEG.event(bigevents).latency});
evvec = zeros(1,0);
for p = 1:length(marktimes)
    evvec(end+1:marktimes(p)) = p;
end;
evvec(end+1:size(mod.LLt,2)) = p+1;

deltimes = find(mod.LLt==0);
plotmod = mod.LLt; plotmod(deltimes) = [];
evvec(deltimes) = [];newmarks =[];
for p = 1:max(evvec)
    [v i] = find(evvec==p);
    newmarks = [newmarks,i(end)];
end;

cols = {'r','b','g','m','c'};
for m = 1:size(mod.LLt,1)
    figure;
    plot(mod.LLt(m,:),cols{m});hold on;
    for n = 1:length(marktimes)
        plot([marktimes(n) marktimes(n)],[get(gca,'ylim')],'k-');
    end;
end;


    

