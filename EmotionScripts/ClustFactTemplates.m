% cluster spectral templates from ICA clustering of single trial spectra


button = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
button = [13:20]; % no button press (apart from the first one)
button = [1,2,4:6,8:12,14,17:21,23,25:30,31,33]; % all 'good' subjects (ones that said they got into it)
button = [1:21,23:27];  % not mr72-2
button = [1,3:9,12,14,16,17,19,21,22,23,24,26,27]; % females
button = [1,4:6,8,9,12,14,17,19,21,23,26,27]; % 'good' females
button = [2,10,11,13,15,18,20,25]; % males
button = [2,10,11,18,20,25]; % 'good' males
eeglab
pl = 1;  alltpspecs = zeros(0,length(freqs)+3); keeptrack = []; clear clustspecs factmap subj_comps
for nx = 1:length(button)
    subjmap(nx) = pl;
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{button(nx)}],[pcadims{button(nx)} subjdims{button(nx)}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);pp=1; clear clustspecs factmap
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{button(nx)}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    clustspecs = zeros(0,length(freqs)+3);
    for tp = 1:size(winv,2)
        coi = nxlists{button(nx)}{tp};        
        for cp = 1:length(coi)
            rcp = find(coi(cp) == gdcomps{button(nx)});
            clustspecs(pp,1:length(freqs)) = activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp); 
            clustspecs(pp,1:length(freqs)) = clustspecs(pp,1:length(freqs))/max(abs(clustspecs(pp,1:length(freqs))));
            clustspecs(pp+1,1:length(freqs)) = clustspecs(pp,1:length(freqs))/max(abs(clustspecs(pp,1:length(freqs))))*-1;
            keeptrack(pl,:) = [button(nx),tp,gdcomps{button(nx)}(rcp)]; 
            keeptrack(pl+1,:) = [button(nx),tp,gdcomps{button(nx)}(rcp)]; 
            xyz = EEG.dipfit.model(gdcomps{button(nx)}(rcp)).posxyz(1,:);
            clustspecs(pp,length(freqs)+1:end) = xyz/max(abs(xyz));pl = pl+2;pp=pp+2;
            if  ~isempty(find(EEG.dipfit.model(gdcomps{button(nx)}(rcp)).posxyz(2,:)))                
                clustspecs(pp,1:length(freqs)) = activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp); 
                clustspecs(pp,1:length(freqs)) = clustspecs(pp,1:length(freqs))/max(abs(clustspecs(pp,1:length(freqs))));
                clustspecs(pp+1,1:length(freqs)) = clustspecs(pp,1:length(freqs))/max(abs(clustspecs(pp,1:length(freqs))))*-1;
                keeptrack(pl,:) = [button(nx),tp,gdcomps{button(nx)}(rcp)]; 
                keeptrack(pl+1,:) = [button(nx),tp,gdcomps{button(nx)}(rcp)]; 
                xyz = EEG.dipfit.model(gdcomps{button(nx)}(rcp)).posxyz(2,:);
                clustspecs(pp,length(freqs)+1:end) = xyz/max(abs(xyz));pl = pl+2;pp=pp+2;
            end;            
        end;
        %for cp = 1:length(gdcomps{button(nx)})
        %    clustspecs(pp,:) = activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp); 
        %    keeptrack(pl,:) = [button(nx),tp,gdcomps{button(nx)}(cp)]; pl = pl+1;pp=pp+1;
        %end; % for all comps, all factors,not just high variance ones
        %tpspecs(tp,:) = mean(clustspecs,1); % for mean specs over whole factor
    end;
    alltpspecs(end+1:end+size(clustspecs,1),:) = clustspecs;
    %alltpspecs(end+1:end+size(tpspecs,1),:) = tpspecs;
    fprintf('\n One More SUBJECT Done: %i\n',nx);ALLEEG=[];EEG=[];
end;
%%%%%%%% run kmeans clustering  %%%%%%%%%%%%%%%%%%%%%%%%%
% cut off freqs 0-3Hz
pcdims = 102;
[pc,eigvec,sv] = runpca(alltpspecs,pcdims);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcdims = round(sqrt(size(alltpspecs,1)/10)); % div by # is pnts/wts desired
 [weights,sphere,compvars,bias,signs,lrates,activations] = runica (alltpspecs','PCA',pcdims,'extended',1,'stop',1e-7,'maxsteps',2000);
winv = pinv(weights*sphere);
figure;
for clust = 1:size(winv,2)
    subplot(ceil(sqrt(size(winv,2))),round(sqrt(size(winv,2))),clust)
    plot(freqs,winv(1:length(freqs),clust),'b');
    set(gca,'xlim',[1 50]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nclusts = 40;
[kout, C,sumd, allsums] = kmeans(eigvec,nclusts,'replicates',5 );%eigvec
% test clustering by collecting distances to centroids
% ,'start',startk % better without seeded locations
% find highest weighted components in each cluster:
clear clustcomps
for clust = 1:size(C,1)
    cp = find(kout == clust); 
    cp(find(abs(zscore(allsums(cp,clust)))>1)) = [];
             % deleting > ? stds cleans up signif 
    cpidx{clust} = cp;
    relcomps = keeptrack(cp,:); 
    cpoi = cell(1,length(gdcomps));   
    for w = 1:length(cp)
    cpoi{relcomps(w,1)}(end+1) = relcomps(w,3);
    end;
    clustcomps{clust} = cpoi;
end;
     
% determine templates for each cluster:
figure; col = ceil(sqrt(nclusts)); row = round(sqrt(nclusts)); % number of clusters
for comp = 1:size(C,1)
    subplot(row,col,comp); clear pic
    for template = 1:size(pc,1)
        pic(template,:) = pc(template,1:length(freqs))*C(comp,template);
    end;        
    pic = mean(pic,1);
        ph =plot(freqs,pic,'b');  set(ph,'linewidth',2);          
        hold on; plot([0 0],[0 50],'k-');   title(['Kmeans Cluster ',int2str(comp)]);
        set(gca,'xlim',[1 50]);  set(gca,'xgrid','on');  set(gca,'xticklabel',[10:10:50]);
end;
axcopy; clear pic
ph =textsc(['Emotion Spectral Clusters (Trials X Freqs); 24 Subj; PCA Dims Retained: ',int2str(pcdims)],'title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

figure;cols = jet(size(C,1));
for comp = 1:size(C,1)
    for template = 1:size(pc,1)
        pic(template,:) = pc(template,length(freqs)+1:end)*C(comp,template);
    end;        
    pic = mean(pic,1);
ph = plot3(pic(1,1),pic(1,2),pic(1,3),'k.');hold on;
set(ph,'color',cols(comp,:));set(ph,'markersize',25);set(gca,'xgrid','on');set(gca,'ygrid','on');set(gca,'zgrid','on');

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Plot Dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
for nx = 1:length(paths)  
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
% look at dipole locations;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/scott/matlab');
%figure; row = round(sqrt(length(clustcomps))); col = ceil(sqrt(length(clustcomps))); cols = jet(length(clustcomps));
    tl = ['Emotion Spectral Clusters; PCA to ',int2str(pcdims),'; RMS weights > mean+1.5*std of each template '];
for clust =1:length(clustcomps) 
    allbesa =  EEG.dipfit.model(1); colset = zeros(1,length(clustcomps{clust}));
    for ss = 1:length(clustcomps{clust})
        if ~isempty(clustcomps{clust}{ss})        
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            dipsources = EEG.dipfit.model(clustcomps{clust}{ss}(1));
            for w = 1:length(clustcomps{clust}{ss})
                dipsources(1,w) = EEG.dipfit.model(clustcomps{clust}{ss}(w));
            end;           
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            colset(ss) = length(dipsources);
            dipsources = [];
        end; 
    end;
    allbesa(1) = [];
    %% Assigns each subject a different color; and makes subjidx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcols = jet(length(clustcomps{clust})); subjidxmat = zeros(0);
    plotcol = cell(0); pl = 0;
    for ss = 1:length(clustcomps{clust})
        if colset(ss) ~= 0
            for nn = 1:colset(ss)
                plotcol{pl+nn} = dcols(ss,:);
                subjidxmat(pl+nn) = ss;
            end;
            pl = pl+nn;
        end;
    end;
    optdipplot = {allbesa,'dipolelength',0,'gui','off','dipolesize',18,'image','mri','spheres','on','color',plotcol};
    figure; dipoleentropy( optdipplot, subjidxmat, 'distance',10);  cbar;
    set(gcf,'color','w');    ph = textsc(['Subject dipole density (10); Cluster: ',int2str(clust)],'title');   set(ph,'fontsize',14);
    savename = ['print -djpeg /home/julie/Scott/Spectral_Clusters/AllSubjClustFactDensity-',int2str(clust),'.jpg'];
    eval(savename)
    close
end;


dens3d = dipoleentropy( optdipplot, subjidxmat, 'distance',10,ones(1,1e5),10);  
figure;dipoleentropy( optdipplot, subjidxmat, 'distance',10,ones(1,length(allbesa)),70);  

    subplot(row,col,clust);
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(clust,:)},'projlines','on','projimg','on');pl = pl+1;
    ph=text(60,-100,145,['Clust: ',int2str(clust)]);
    set(ph,'color','y'); set(ph,'fontsize',14);
    view(60,20);
end;
set(gcf,'color','w');
ph = textsc(tl,'title');
set(ph,'color','r');set(ph,'fontsize',14);














[weights,sphere,compvars,bias,signs,lrates,activations] = runica(forclust','pca',30);pc = activations;
figure;plot(pc,'k-');
 for q = 1:size(pc,2)
     pc(:,q) = pc(:,q)-mean(pc(:,q));
 end;
 for q = 1:size(pc,1)
     pc(q,:) = pc(q,:)-mean(pc(q,:));
 end;
nclust =10 ; 
  [idx,C,sumd,D,outliers] = robust_kmeans(pc,nclust,3);
%The input data has to be in the form: spectra x comp, like 10x160 (10 PCA comp to represent the spectra on 160 ICA components).

kout = idx;
figure;
for cl = 1:nclust
    oneclust = find(kout == cl);
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),cl)
    for w = 1:length(oneclust)
        plot(freqs(fr),forclust(oneclust(w),:)');hold on;        
    end;
    set(gca,'xlim',[1 50]);
    title(['Cluster ',int2str(cl)]);
end;axcopy
%%%%%%%% >> OR <<%%%%%%%%%%%%% 
reps = 5;
%[optk,centr,clst,Cg] =Kmeangrp(clustspecs,20,reps,1); % find optimal cluster number, don't plot
%[kout, C,sumd, allsums] = kmeans(clustspecs,20,'replicates',reps);

[optk,centr,clst,Cg] =Kmeangrp(alltpspecs,11,reps,1); % find optimal cluster number, don't plot
nclust = 20;
[kout, C,sumd, allsums] = kmeans(forclust,nclust,'replicates',2);
%%%%%%%%%%%%%%%%%%%%% 
figure;
for cl = 1:nclust
    oneclust = find(kout == cl);
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),cl)
    for w = 1:length(oneclust)
        plot(freqs,alltpspecs(oneclust(w),:)','k-');hold on;
    end;
    x=mean(alltpspecs(oneclust,:),1); ph=plot(freqs,x,'r-');hold on;        
    set(gca,'xlim',[1 50]);    title(['Cluster ',int2str(cl)]);
    set(gca,'xtick',[10:10:50]); set(gca,'xgrid','on');  
end;axcopy
%%%%%%%%%%%%%%%%%%%%%
%%%   USE ICA TO CLUSTER SPECTRAL TEMPLATES  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 
nclust = 25;
fr = find(freqs>3 & freqs < 35);
forclust = alltpspecs(:,fr);
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(forclust','pca',nclust,'stop',1e-7,'maxsteps',2000);
ws = weights*sphere; winv = pinv(ws);
%%%%%%%%%%%%%%%%%%%%% 
 figure;% for WT
 for clust = 1:size(activations,1)
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),clust) 
     plot(freqs,activations(clust,:),'k-');
    set(gca,'xlim',[1 max(freqs)]);    title(['Cluster ',int2str(cl)]);
    set(gca,'xtick',[10:10:max(freqs(fr))]); set(gca,'xgrid','on');      
 end;
 figure; % for TW
 for clust = 1:size(winv,2)
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),clust) 
     plot(freqs(fr),winv(:,clust)','k-');
    set(gca,'xlim',[1 max(freqs(fr))]);    title(['Cluster ',int2str(clust)]);
    set(gca,'xtick',[10:10:max(freqs(fr))]); set(gca,'xgrid','on');      
 end;
%   find absolute wts for each component and choose cluster members   %%%%%%%%%%
cut = 6;
for clust = 1:size(activations,1)
hiwts = find(abs(activations(clust,:)) > cut);
cc{clust} = subj_comps(hiwts,:);
end; clustcomps = cell(1,size(activations,1));
for clust = 1:length(cc)
    for nx = 1:length(button)
        comps = zeros(1,0);
        for x= 1:length(cc{clust})
            if cc{clust}(x,1) == button(nx)
                comps(end+1) = cc{clust}(x,2);
            end;
        end;comps = unique(comps);
        clustcomps{clust}{nx} = comps ;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
for nx = 1:length(paths)  -1
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
% look at dipole locations;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcdims = 25;  addpath('/home/scott/matlab');
figure; row = round(sqrt(length(clustcomps))); col = ceil(sqrt(length(clustcomps))); cols = jet(length(clustcomps));
    tl = ['Emotion Spectral Clusters; PCA to ',int2str(pcdims),'; abs weights > 3 '];
for clust =1:6%length(clustcomps) 
    allbesa =  EEG.dipfit.model(1); colset = zeros(1,length(clustcomps{clust}));
    for ss = 1:length(clustcomps{clust})
        if ~isempty(clustcomps{clust}{ss})        
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            dipsources = EEG.dipfit.model(clustcomps{clust}{ss});
            for w = 1:length(clustcomps{clust}{ss})
                dipsources(1,w) = EEG.dipfit.model(clustcomps{clust}{ss}(w));
            end;           
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            colset(ss) = length(dipsources);
            dipsources = [];
        end; 
    end;
    allbesa(1) = [];
    %% Assigns each subject a different color; and makes subjidx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcols = jet(length(clustcomps{clust})); subjidxmat = zeros(0);
    plotcol = cell(0); pl = 0;
    for ss = 1:length(clustcomps{clust})
        if colset(ss) ~= 0
            for nn = 1:colset(ss)
                plotcol{pl+nn} = dcols(ss,:);
                subjidxmat(pl+nn) = ss;
            end;
            pl = pl+nn;
        end;
    end;
    optdipplot = {allbesa,'dipolelength',0,'gui','off','dipolesize',18,'image','mri','spheres','on','color',plotcol};
    figure; dipoleentropy( optdipplot, subjidxmat, 'distance',10);  cbar;
    set(gcf,'color','w');    ph = textsc(['Subject dipole density (10); Cluster: ',int2str(clust)],'title');   set(ph,'fontsize',14);

    subplot(row,col,clust);
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(clust,:)},'projlines','on','projimg','on');pl = pl+1;
    ph=text(60,-100,145,['Clust: ',int2str(clust)]);
    set(ph,'color','y'); set(ph,'fontsize',14);
    view(60,20);
end;
set(gcf,'color','w');
ph = textsc(tl,'title');
set(ph,'color','r');set(ph,'fontsize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear emoscores collectscores clustscores factmap  % get mnemodiff from PCAforEmos.m nx is rel to 'button' list
for cl = 1:nclust
    oneclust = find(kout == cl);
    for cat = 1:length(oneclust)
        clear emoscores
        for nx = 1:length(button)
            factmap(1,1) = 0;
            factmap(2:11) = stpmap{nx}(1:10);
            if oneclust(cat) >= subjmap(nx) & oneclust(cat) < subjmap(nx+1)
                for tp = 1:length(factmap)-1
                    if oneclust(cat)>=subjmap(nx)+sum(factmap(1:tp))&oneclust(cat)<subjmap(nx)+sum(factmap(1:tp+1))
                        for e = 1:size(mnemodiff,2)
                            emoscores(1,e) = mnemodiff(tp,e,nx);
                        end;
                    end;
                end;
            end;
        end;        
        collectscores(cat,:) = emoscores;
    end;
    clustscores(cl,:) = mean(collectscores,1);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%emoorder = [4,2,10,8,14,12,6,9,7,13,11,1,5,3,15]; clear cols  % no pre or post baseline
emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'};
cols = jet(length(emo2));
figure;
for clust = 1:size(clustscores,1)
    subplot(ceil(sqrt(nclust)),round(sqrt(nclust)),clust)
    for e = 1:size(clustscores,2)
        ph =bar(e,clustscores(clust,e));hold on;
        set(ph,'facecolor',cols(e,:));
        ph=text(e,0,emo2{e}); set(ph,'fontsize',7);
        set(ph,'rotation',90);
    end;
    set(gca,'xlim',[0 16]);
end;
% plot mean clusters and emo bars together
figure;pp = 1; row = ceil(sqrt(nclust*2)); col = round(sqrt(nclust*2));
for cl = 1:nclust
    oneclust = find(kout == cl);
    subplot(row,col,pp)
    x=mean(alltpspecs(oneclust,:),1); plot(freqs,x);hold on;   pp = pp+1;      
    set(gca,'xlim',[1 50]);set(gca,'xtick',[5:5:50]);set(gca,'xgrid','on'); title(['Cluster ',int2str(cl)]);
    set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40 [] 50});
    subplot(row,col,pp)
    for e = 1:size(clustscores,2)
        ph =bar(e,clustscores(cl,e));hold on;
        set(ph,'facecolor',cols(e,:));
        ph=text(e,0,emo2{e}); set(ph,'fontsize',7);
        set(ph,'rotation',90);title(['Cluster ',int2str(cl)]);set(gca,'box','off');
    end;  set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]); pp = pp+1;    
end;axcopy
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    print  -dpsc2 -Pcoloring 

% run pca on emoscores
 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(clustscores','pca',3);
 ws = weights*sphere; winv = pinv(ws);
cols = jet(15);
 figure; 
 for e = 1:size(winv,1)
 ph =plot3(winv(e,1),winv(e,2),winv(e,3),'.');hold on;
 set(ph,'markersize',20);
 set(ph,'color',cols(e,:));
 ph = text(winv(e,1),winv(e,2),winv(e,3),emo2{e});
 set(ph,'color',cols(e,:)); 
 pl =plot3([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)],[-.04  winv(e,3)]);
 set(pl,'color',cols(e,:))             
 end;
 set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             xlabel('Winv 1'); ylabel('Winv 2'); zlabel('Winv 3'); 


% find dipoles of clustered spectra
figure;cols = jet(20);
for cl = 1:20
    dipcol = cols(cl,:); clear plotpoles
    oneclust = find(kout == cl);
    for cat = 1:length(oneclust)
        for nx = 1:length(gdcomps)           
            factmap(1,1) = 0;
            factmap(2:11) = stpmap{nx}(1:10);
            if oneclust(cat) >= subjmap(nx) & oneclust(cat) < subjmap(nx+1)
                for tp = 1:length(factmap)-1
                    if oneclust(cat)>=subjmap(nx)+sum(factmap(1:tp))&oneclust(cat)<subjmap(nx)+sum(factmap(1:tp+1))
                        EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
                        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
                        if cat ==1 
                            plotpoles =  EEG.dipfit.model(1);
                            plotpoles(end+1) = EEG.dipfit.model(dipolelist(cat));
                            plotpoles(1) = [];
                        else
                            plotpoles(end+1) = EEG.dipfit.model(dipolelist(cat));
                        end;  
                        ALLEEG=[];EEG=[];
                    end;
                end;
            end;
        end;        
    end;
    subplot(5,4,cl)                        
    dipplot(plotpoles,'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',{dipcol},'projlines','on','projimg','on');hold on;     
    view(60,20); title(int2str(cl));    
end;
set(gcf,'color','w');
