% clusters components across subjects by 3D dippole locations

%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Get dipole locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nx = 10:10%length(gdcomps)
    EEG = pop_loadset(datsets{nx} ,['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'point',[1:100] );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    sph=floatread(['/data/common2/emotion',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion',paths{nx},wtss{nx}],wtssize{nx}); 
    EEG.icaweights=wts;
    EEG.icasphere=sph;
    EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_multifit(EEG, gdcomps{nx}, 'settings',{}, 'threshold',40, 'plotopt',{ 'normlen', 'on', 'image', 'fullmri'});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    EEG = pop_saveset( EEG,'sources.set' , ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    ALLEEG=[];EEG=[];
end;
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /data/common2/emotion/clusters/allPreBasespecs.mat freqs allspec comment
fr = find(freqs > 2 & freqs < 15);clustspecs = zeros(0,length(fr));
allloc = zeros(0,6);keeptrack = zeros(0,2);
for nx = 1:35%length(gdcomps)
    if ~isempty(clustcomps{nx})
addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
        EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
        for cp = 1:length(clustcomps{nx})
            subjwinv = EEG.dipfit.model(clustcomps{nx}(cp)).posxyz;
            if size(subjwinv,1) > 1
                subjwinv(1,end+1:end+size(subjwinv,2)) = subjwinv(2,:); subjwinv(2,:) = [];
            else
                subjwinv(1,end+1:end+3) = zeros(1,3);
            end;            
            allloc(end+1,:) = zscore(subjwinv)';
            clustspecs(end+1,:) = zscore(allspec{nx}(find(gdcomps{nx} == clustcomps{nx}(cp)),fr))';
            keeptrack(end+1,:) = [nx,clustcomps{nx}(cp)];
        end;
        ALLEEG=[];EEG=[];
    end;
end;
clear acts wts
    [pc,eigvec,sv] = runpca(clustspecs,16);
acts{1} = eigvec; %wts{1} = pc;
    [pc,eigvec,sv] = runpca(allloc);
acts{2} = eigvec; %wts{2} = pc;
kmeansmat  = [acts{1} acts{2}];

kmeansmat  = [allloc clustspecs];

% try kmean
reps = 2;
        [optk,centr,clst,Cg] =Kmeangrp(kmeansmat,30,reps,1); % find optimal cluster number, don't plot
cnum =16;  % 10 is too few(larger clusters, no fml); 15 is pretty good, 16 no better, 17 gets a fml, but splits sup par into left/right, 18 really good: gies right and left occ (not just right)

[kout, C,sumd, allsums] = kmeans(allloc,cnum,'replicates',5);

clear clustcps
for clust = 1:size(C,1)
    cpidx = find(kout == clust); 
    %cpidx(find(abs(zscore(allsums(cpidx,clust)))>1)) = [];
             % deleting > ? stds cleans up signif 
    relcomps = keeptrack(cpidx,:); 
    cpoi = cell(1,length(gdcomps));   
    for w = 1:length(cpidx)
    cpoi{relcomps(w,1)}(end+1) = relcomps(w,2);
    end;
    clustcps{clust} = cpoi;
end;
 
clustorder = [1,13,6,14,15,10,11,2,16,12,8,4,7,5,9,3]; % order of dipole clusters to image
save /data/common2/emotion/clusters/KmeanClustDips.mat clustcps kout C sumd allsums keeptrack clustorder 
clustorder = [8,14,15,16,10,6,4,3,2,1,13,9,7,12,11,5]; % order of dipole clusters to image
save /data/common2/emotion/clusters/KmeanClustDipsSpec.mat clustcps kout C sumd allsums keeptrack clustorder 
clustorder = [3,2,5,4,10,9,1,7,6,8]; % order of dipole clusters to image
save /data/common2/emotion/clusters/KmeanClustFrontDips.mat clustcps kout C sumd allsums keeptrack clustorder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize mean specs from all clusters
figure; row = 4; col = 4;
for clust = 1:size(C,1)
    cpidx = find(kout == clust); 
    subplot(row,col,clust)
ph = plot(freqs(fr),clustspecs(cpidx,:),'c');hold on;
ph = plot(freqs(fr),mean(clustspecs(cpidx,:)),'r');
set(gca,'xlim',[2 15]); set(gca,'xgrid','on');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
for nx = 1:length(paths)  
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/scott/matlab');
rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
figure; row =4; col=4; cols = jet(17);
for clust = 1:cnum   
    numsubj = 0;numdips = 0;
    for x = 1:length(gdcomps)
        if ~isempty(clustcps{clust}{x})
            numsubj = numsubj+1;
            numdips = numdips+length(clustcps{clust}{x});
        end;
    end;    
            EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    allbesa =  EEG.dipfit.model(1); %colset = zeros(1,length(clustcps{clust}));
    for ss = 1:length(clustcps{clust})
        if ~isempty(clustcps{clust}{ss})        
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            dipsources = EEG.dipfit.model(clustcps{clust}{ss}(1));
            for w = 1:length(clustcps{clust}{ss})
                dipsources(1,w) = EEG.dipfit.model(clustcps{clust}{ss}(w));
            end;     
            if isfield(dipsources,'diffmap')
                dipsources = rmfield(dipsources,'diffmap');
            end;
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            %colset(ss) = length(dipsources);
            dipsources = [];
        end; 
    end;
    allbesa(1) = [];
    %% Assigns each subject a different color; and makes subjidx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dcols = jet(length(clustcps{clust})); 
    %plotcol = cell(0); pl = 0;clear subjidxmat
    %for ss = 1:length(clustcps{clust})
    %    if colset(ss) ~= 0
    %        for nn = 1:colset(ss)
    %            plotcol{pl+nn} = dcols(ss,:);
    %            subjidxmat(pl+nn) = ss;
    %        end;
    %        pl = pl+nn;
    %    end;
    %end;
    %optdipplot = {allbesa,'dipolelength',0,'gui','off','dipolesize',18,'image','mri','spheres','on','color',plotcol};
    %figure; mydipoleentropy( optdipplot, subjidxmat, 'distance',10);  cbar;
    %figure; dipoledensity( optdipplot, 'subjind',subjidxmat, 'method','distance','methodparam',10, 'mri_view' ,'top', 'gridsize',30); 
    %ph = textsc(['Subj dipole density; Dipole Cluster: ',int2str(clust),'; ',int2str(numdips),' from ',int2str(numsubj),' subj'],'title');   
    %set(ph,'fontsize',16); set(ph,'color','r');
    %savename = ['print -djpeg /home/julie/Scott/Emotion/DipoleClustDens-',int2str(clust),'.jpg'];
    %eval(savename)
    %close
    subplot(row,col,clust);
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(clust,:)},'projlines','on','projimg','on');
    %if orient == 1
    ph=text(60,-100,145,['Clust: ',int2str(clust)]);
    %else
    %    ph=text(60,-100,145,['Clust: -',int2str(clust)]);    
    %end;
    set(ph,'color','y'); set(ph,'fontsize',14);
    view(60,20);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cluster dipoles after ICA clustering of spectra/emotion weights
allloc = zeros(0,6);dp=1;
for nx = 1:length(subjfactors)
    if ~isempty(subjfactors{nx})
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    for tp = 1:length(subjfactors{nx})
        for cp = 1:length(nxlists{nx}{subjfactors{nx}(tp)})
            rcp = find(nxlists{nx}{subjfactors{nx}(tp)}(cp) == gdcomps{nx});
            subjwinv = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp)).posxyz;
            if size(subjwinv,1) > 1
                subjwinv(1,end+1:end+size(subjwinv,2)) = subjwinv(2,:); subjwinv(2,:) = [];
            else
                subjwinv(1,end+1:end+3) = zeros(1,3);
            end;            
            allloc(end+1,:) = subjwinv;
            dipsource{dp} = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));dp= dp+1;
        end;        
        oneact = activations;
        zerout= [1:15];  zerout(subjfactors{nx}(fac)) = [];
        oneact(zerout,:) = 0;
        backprojdat = winv * oneact;  %figure;pl = 1;
        for k = 1:length(bde)
            rk = bde(k);
            badspecs(k,:) = mean(backprojdat(emomap(rk):emomap(rk+1)-1,:),1);
        end;
        for k = 1:length(gde)
            rk = gde(k);
            goodspecs(k,:) = mean(backprojdat(emomap(rk):emomap(rk+1)-1,:),1);
        end;
        badspecs = mean(badspecs,1);goodspecs = mean(goodspecs,1);
        
    end;
    end;
    ALLEEG=[];EEG=[];
end;
(length(freqs)*(rcp-1)+1:length(freqs)*rcp)
[optk,centr,clst,Cg] =Kmeangrp(allloc,4,reps,1); % find optimal cluster number, don't plot

% try kmeans
cnum = 3;  
ids = kmeans(allloc,optk,'replicates',25);


figure;
for clust = 1:max(ids)
idxs = find(ids == clust);
subplot(2,2,clust)
for w=  1:length(idxs)
      dipplot(dipsource{idxs(w)},'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','projlines','on','projimg','on'); camzoom(.7)
      end;view(60,20); 
end;
