% clusters spectra across subjects to find typical spectral shifts. 

cd (['/data/common2/emotion',paths{3},'/ersps/']);load ContDataERSPs.mat
icamatall = zeros(length(freqs),0);clear subjtrials
index = 1;
for nx = 2:length(gdcomps)
%index = 1; % for within subj
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    %icamatall = zeros(length(freqs),0);clear subjtrials  % for within subj decomp
   pl = 1;clear savemean
    for cmp = 1:length(gdcomps{nx})
        for k = 1:length(Alllongersps)
            oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            trialmat = zeros(size(oneemo,1),0);
            for trl = 1:5:size(onecmp,2)-6  % moved up by 3 each time except for nx=23 (4)
                onetrl = mean(onecmp(:,trl:trl),2);%  avg over 4 (1 sec)
                trialmat(:,end+1) = onetrl;
                keeptrack(index,:) = [nx,gdcomps{nx}(cmp),k];index = index+1;
            end;
            icamatall(:,end+1:end+size(trialmat,2)) = trialmat;
            subjtrials(k) = size(trialmat,2); clear icamat   
        end;
        numtrials{nx} = subjtrials;%clear subjtrials
        fprintf('\n One More COMPONENT Done: %i of %i',cmp,length(gdcomps{nx}));
    end;
    fprintf('\n One More SUBJECT Done: %i\n',nx);
%floatwrite(icamatall, ['/data/common2/emotion',paths{nx},'ersps/','FrqsXTrls.fdt']);    % for within subj decomp
end;
dims = size(icamatall);
pcs = 40% round(sqrt(size(icamatall,2)/8));% div by # = desired ratio frames/dims^2
%save /data/common2/emotion/clusters/FreqsXTrialsClusts.mat numtrials gdcomps dims freqs pcs keeptrack 
save /data/common2/emotion/clusters/FrqsXTrlsNoSmoothClusts.mat numtrials gdcomps dims freqs pcs keeptrack 

floatwrite(icamatall, '/data/common2/emotion/clusters/FreqsXTrialsNoSmooth.fdt');    % for across subject decomp
ICA_LINUX = '/data/common/matlab/ica_linux2.4';
ICA_SCRIPT = ['/data/common2/emotion/clusters/ClustPwr.sc'];
%ICA_SCRIPT = ['/data/common2/emotion',paths{nx},'ClustPwr.sc'];% for within subj decomp
% generate an ICA script 
fid = fopen(ICA_SCRIPT, 'w');
fprintf(fid, 'DataFile %s\n', '/data/common2/emotion/clusters/FreqsXTrials.fdt');
%fprintf(fid, 'DataFile %s\n', ['/data/common2/emotion',paths{nx},'ersps/','FrqsXTrls.fdt']);% for within subj decomp
fprintf(fid, 'chans %d\n', size(icamatall,1));
fprintf(fid, 'frames %d\n',size(icamatall,2));
fprintf(fid, 'WeightsOutFile %s\n', '/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC40.wts');
fprintf(fid, 'SphereFile %s\n', '/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC40.sph');  
%fprintf(fid, 'WeightsOutFile %s\n',['/data/common2/emotion',paths{nx},'ersps/','FrqsXTrlsPC',int2str(pcs),'.wts']);% for within subj decomp
%fprintf(fid, 'SphereFile %s\n',['/data/common2/emotion',paths{nx},'ersps/','FrqsXTrlsPC',int2str(pcs),'.sph']);% for within subj decomp  
fprintf(fid, 'sphering on\n');
fprintf(fid, 'bias on\n');
fprintf(fid, 'extended 1\n');
fprintf(fid, 'pca %d\n', pcs); %pcs ********
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
%save FreqsXTrialsClusts.mat subjtrials gdcomps dims freqs pcs keeptrack   % for within subj decomp
end;

/data/common/matlab/ica_linux2.4 < /data/common2/emotion/clusters/ClustPwr.sc

% Call in wts and sph
load /data/common2/emotion/clusters/FreqsXTrialsClusts.mat  numtrials gdcomps dims freqs pcs keeptrack 
sph=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC50.sph',[dims(1) dims(1)]); 
wts=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC50.wts',[pcs dims(1)]); 
sph=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC40.sph',[dims(1) dims(1)]); 
wts=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC40.wts',[pcs dims(1)]); 
%sph=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC30.sph',[dims(1) dims(1)]); 
%wts=floatread('/data/common2/emotion/clusters/FreqsXTrialsPC30.wts',[30 dims(1)]); 

load /data/common2/emotion/clusters/FrqsXTrlsNoSmoothClusts.mat 
%sph=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC50.sph',[dims(1) dims(1)]); 
%wts=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC50.wts',[pcs dims(1)]); 
%sph=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC30.sph',[dims(1) dims(1)]); 
%wts=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC30.wts',[pcs dims(1)]); 
sph=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC40.sph',[dims(1) dims(1)]); 
wts=floatread('/data/common2/emotion/clusters/FrqsXTrlsNoSmoothPC40.wts',[pcs dims(1)]); 
icamatall = floatread('/data/common2/emotion/clusters/FreqsXTrials.fdt',[dims(1) dims(2)]);
ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize spectral templates:
figure;row  = 5; col = 6;
for sp = 1:30
subplot(row,col,sp)
%plot(activations(sp,:));
ph = plot(freqs,winv(:,sp));
set(gca,'xtick',[0:10:50]);set(ph,'linewidth',2);set(gca,'xlim',[0 50]);
set(gca,'xgrid','on');set(gca,'ylim',[-4 4]);
title(int2str(sp));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  call in wts for each subj and plot 10 spectra:
figure;row  = 10; col = 10;
ntemps = 10;ktrack = zeros(0);kmrms = zeros(0,ntemps);  emon = 0; % 1 collects rms for each emotion, 0 just by comp
ss = 1; pl = 1;
for nx = 1:length(paths)
    if nx == 11 | nx == 21 | nx == 31
        figure; pl = 1;
    end;    
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);clear indivcprms
    load FreqsXTrialsClusts.mat subjtrials gdcomps dims freqs pcs keeptrack
    sph=floatread(['/data/common2/emotion',paths{nx},'/ersps/','FrqsXTrlsPC',int2str(pcs),'.sph'],[dims(1) dims(1)]); 
    wts=floatread(['/data/common2/emotion',paths{nx},'/ersps/','FrqsXTrlsPC',int2str(pcs),'.wts'],[pcs dims(1)]); 
    icamatall = floatread('/data/common2/emotion/clusters/FreqsXTrials.fdt',[dims(1) dims(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);
    for sp = 1:ntemps
        collectwinv(ss,:) = winv(:,sp)';subjassign(ss,:) = [nx,sp];ss = ss+1;
        subplot(row,col,pl);ph = plot(freqs,winv(:,sp));pl = pl+1;set(gca,'xtick',[0:10:50]);
        set(ph,'linewidth',2);set(gca,'xlim',[0 50]);set(gca,'xgrid','on');set(gca,'ylim',[-4 3]);
    end;
    kmap = ones(1,1); clear emorms indivcprms
    for p = 2:length(subjtrials)+1
        kmap(1,p) = subjtrials(p-1)+kmap(p-1);
    end;
    for clust = 1:ntemps
    startcp = 0;
        for cp = 1:length(gdcomps{nx})
            if emon == 1
                for k = 1:length(subjtrials)
                    %emorms(k,cp,clust) = sqrt(mean(activations(clust,startcp+kmap(k):startcp+kmap(k+1)-1).^2));
                    emorms(k,cp,clust) = mean(activations(clust,startcp+kmap(k):startcp+kmap(k+1)-1));
                end;
            else
                indivcprms(cp,clust) = sqrt(mean(activations(clust,startcp+kmap(1):startcp+kmap(end)-1).^2));
            end;
            startcp = startcp + kmap(end)-1;
        end;
        %ktrack(end+1,:) = [nx,gdcomps{nx}(cp)];
    end;    
    if emon == 1
        subjemorms{nx} = emorms;        
    else
        kmrms(end+1:end+size(indivcprms,1),:) = indivcprms; 
        indivrms{nx} = indivcprms;
    end;
    fprintf('\n One More SUBJECT Done: %i',nx);
end;
save /data/common2/emotion/clusters/ClustTemplates.mat collectwinv subjassign freqs
load /data/common2/emotion/clusters/ClustTemplates.mat collectwinv subjassign freqs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcdims = 5;
[weights,sphere,compvars,bias,signs,lrates,activations] = runica (collectwinv,'PCA',pcdims,'extended',1,'stop',1e-7,'maxsteps',2000); winv = pinv(weights*sphere);
%%%% 'WT' ICA:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; for x=1:size(activations,1)            
    subplot(ceil(sqrt(pcdims)),round(sqrt(pcdims)),x)
    plot(freqs,activations(x,:));
    set(gca,'xlim',[1 50]);    title(['Cluster ',int2str(x)]);set(gca,'xgrid','on');
end;axcopy
clear clusttemps
for clust = 1:pcdims
    oneclust = winv(:,clust);
    cutval = zscore(winv(:,clust));
    clusttemps{clust} = subjassign(find(abs(cutval) > 2),:); % determines zscore cut for weights
end;
% find components for each cluster by going back to individual ICAs and RMS
clear clustcomps
for clust = 1:pcdims
    clustcps = cell(1,length(gdcomps));
    nxtemp = clusttemps{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(find(nxtemp(:,1) == nx))
            toi = nxtemp(find(nxtemp(:,1) == nx),2)';
            indivcprms = indivrms{nx}; % gives cp,clust
            for tp = 1:length(toi)
                distrib = indivcprms(:,toi(tp));
                tpcomps = gdcomps{nx}(abs(zscore(distrib))>1.5);
                clustcps{nx}(end+1:end+length(tpcomps)) = tpcomps;
            end;
        end;
    end;
    clustcomps{clust} = clustcps;
end;
save /data/common2/emotion/clusters/misc2.mat clustcomps
load /data/common2/emotion/clusters/misc2.mat clustcomps
%%%% 'TW' ICA:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure;for x=1:size(winv,2)
subplot(ceil(sqrt(pcdims)),round(sqrt(pcdims)),x)
plot(freqs,winv(:,x))
    set(gca,'xlim',[1 50]);    title(['Cluster ',int2str(x)]);set(gca,'xgrid','on');
end;axcopy
clear clusttemps
for clust = 1:pcdims
    cutval = zscore(activations(clust,:));
    nxtemp = keeptrack(find(abs(cutval) > 1),:);  % determines zscore cut for weights
    % keeptrack is [nx,gdcomps{nx}(cp),emo]
    clustcps = cell(1,length(gdcomps));
    for nx = 1:length(gdcomps)
        if ~isempty(find(nxtemp(:,1) == nx))
            coi = nxtemp(find(nxtemp(:,1) == nx),2)';% gives cp
            for cp = 1:length(toi)
                distrib = activations(clust,toi(cp));
                tpcomps = gdcomps{nx}(find(distrib > mean(distrib) + 1.5 * std(distrib)));
                clustcps{nx}(end+1:end+length(tpcomps)) = tpcomps;
            end;
        end;
    end;
    clustcomps{clust} = clustcps;
end;
% find components for each cluster by going back to individual ICAs and RMS
clear clustcomps
for clust = 1:pcdims
    clustcps = cell(1,length(gdcomps));
    nxtemp = clusttemps{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(find(nxtemp(:,1) == nx))
            toi = nxtemp(find(nxtemp(:,1) == nx),2)';
            indivcprms = indivrms{nx}; % gives cp,clust
            for tp = 1:length(toi)
                distrib = indivcprms(:,toi(tp));
                tpcomps = gdcomps{nx}(find(distrib > mean(distrib) + 1.5 * std(distrib)));
                clustcps{nx}(end+1:end+length(tpcomps)) = tpcomps;
            end;
        end;
    end;
    clustcomps{clust} = clustcps;
end;
save /data/common2/emotion/clusters/misc2.mat clustcomps
load /data/common2/emotion/clusters/misc2.mat clustcomps

%  Cluster component spectral templates across subjects using kmeans  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
nclusts = 20;
[ids, C, SUMD, D] = kmeans(collectwinv,nclusts,'replicates',10);
save /data/common2/emotion/clusters/misc.mat ids C SUMD D
load /data/common2/emotion/clusters/misc.mat ids C SUMD D
figure;
for clust = 1:nclusts
    oneclust = find(ids == clust);
    cutval = zscore(D(oneclust,clust));
    oneclust(find(abs(cutval) > .75)) = [];
    subplot(ceil(sqrt(nclusts)),round(sqrt(nclusts)),clust)
    for w = 1:length(oneclust)
        plot(freqs,collectwinv(oneclust(w),:),'b');hold on;        
    end;
    plot(freqs,mean(collectwinv(oneclust,:),1),'r');
    set(gca,'xlim',[1 50]);    title(['Cluster ',int2str(clust)]);
    clusttemps{clust} = subjassign(oneclust,:);
end;axcopy
save /data/common2/emotion/clusters/misc1.mat clusttemps
load /data/common2/emotion/clusters/misc1.mat clusttemps

% find components for each cluster by going back to individual ICAs and RMS
for kclust = 1:nclusts
    clustcomps = cell(1,length(gdcomps));
    nxtemp = clusttemps{kclust};
    for nx = 1:length(gdcomps)
        if ~isempty(find(nxtemp(:,1) == nx))
            toi = nxtemp(find(nxtemp(:,1) == nx),2)';
            indivcprms = indivrms{nx}; % gives cp,clust
            for tp = 1:length(toi)
                distrib = indivcprms(:,toi(tp));
                tpcomps = gdcomps{nx}(find(distrib > mean(distrib) + 1.5 * std(distrib)));
                clustcomps{nx}(end+1:end+length(tpcomps)) = tpcomps;
            end;
        end;
    end;
    kclustcomps{kclust} = clustcomps;
end;
save /data/common2/emotion/clusters/misc2.mat kclustcomps
load /data/common2/emotion/clusters/misc2.mat kclustcomps

%%%%% decompose emotion means with ica within subj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcdims = 3;
for nx = 1:6%length(gdcomps)
    emorms = subjemorms{nx};
    [weights,sphere,compvars,bias,signs,lrates,activations] = runica (reshape(emorms,size(emorms,1),size(emorms,2)*size(emorms,3)),'PCA',pcdims,'extended',1,'stop',1e-7,'maxsteps',2000);
    winv = pinv(weights*sphere);
    cols = jet(15);
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
    figure;  % for one point per emotion (mean)
    for wv = 1:size(winv,2)
        subplot(round(sqrt(size(winv,2))),ceil(sqrt(size(winv,2))),wv)
        for e = 1:size(winv,1)
            ph=plot(e,winv(e,wv),'.');hold on;
            set(ph,'markersize',20);set(ph,'color',cols(e,:));
            ph = text(e,winv(e,wv),emo2{e});
            set(ph,'color',cols(e,:)); 
        end;
    end;
    textsc(['Subject ',int2str(nx),'; PCA to ',int2str(pcdims),'; ICA on Emos X Comps*Templates'],'title');
    cpacts{nx} = reshape(activations,size(winv,2),size(emorms,2),size(emorms,3));
    nclusts = 3;[ids, C, SUMD, D] = kmeans(winv,nclusts,'replicates',10);
    for r = 1:size(winv,2)
        emgrps{r} = {emos{find(ids==r)}};
    end;
    allegrps{nx} = emgrps;
end;


% cluster emoscores by kmeans:
nclusts = 25;
[ids, C, SUMD, D] = kmeans(kmrms,nclusts,'replicates',10);
% determine templates for each cluster:
figure;pl = 1; col = ceil(sqrt(nclusts)); row = round(sqrt(nclusts)); % number of clusters
for comp = 1:size(C,1)
    subplot(row,col,comp); clear pic
    for template = 1:ntemps % number of dimensions clustered over
        pic(:,template) = winv(:,template)*C(comp,template);
    end;        
    pic = mean(pic,2);    plot(freqs,pic,'k-');  pl = pl+1;        
    hold on;   set(gca,'ylim',[-4 2]); title(['Kmeans Cluster ',int2str(comp)]);
    set(gca,'xtick',[0:10:50]);set(gca,'xticklabel',[]);set(gca,'xgrid','on');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%**********************************%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%**********************************%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  All subject spectral decomposition analysis:
% Use 'keeptrack' to find components and emotion stats
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
for clust = 1:24
    for nx = 2:length(gdcomps) % subject 1 not included
        %cpemorms = zeros(length(gdcomps{nx}),length(emos));
        %cpemovar = zeros(length(gdcomps{nx}),length(emos));
        clear cpemorms cpemovar
        for emo = 1:length(emos)
            for cp = 1:length(gdcomps{nx})
                trialidx = find(keeptrack(:,1) == nx & keeptrack(:,2) == gdcomps{nx}(cp) & keeptrack(:,3) == emo);
                cpemoacts = activations(clust,trialidx);
                cpemorms(cp,emo) = sqrt(mean(cpemoacts.^2));
                cpemovar(cp,emo) = var(cpemoacts);
            end;
        end;
        actstats{nx}{1} = cpemorms;
        actstats{nx}{2} = cpemovar;
        fprintf('\n one more SUBJECT done: %s',int2str(nx));
    end;
    cluststats{clust} = actstats;
end;
save /data/common2/emotion/clusters/SpecClusteringStats.mat cluststats
%%%%%%%%%%%
% Find interesting trends
% 1) find if any emotions have higher rms or variance in a cluster

for clust = 1:24
    for nx = 2:length(cluststats{clust})
        for emo = 1:length(emos)
            emorms(nx,emo) = max(cluststats{clust}{nx}{1}(:,emo)); % rms
            emovar(nx,emo) = max(cluststats{clust}{nx}{2}(:,emo)); % var        
        end;
        emorms(nx,:) = zscore(emorms(nx,:))';
        emovar(nx,:) = zscore(emovar(nx,:))';        
        for cp = 1:length(gdcomps{nx})  % make this across subject
            colmeans(1,cp) = mean(cluststats{clust}{nx}{1}(cp,:));
        end;
        clustcomps{clust}{nx} = gdcomps{nx}(find(zscore(colmeans) > 1.25));
    end;
    emostats{clust}{1} = emorms;
    emostats{clust}{2} = emovar; clear colmeans
end;
save /data/common2/emotion/clusters/SpecClusteringStats.mat cluststats clustcomps emostats
                
%%%%  2) cluster emos according to scores in each cluster
rv = 1;clear colemo 
for clust = 1:24
    for nx = 2:length(emostats)
        for emo = 1:size(emostats{nx}{rv},2)
            colemo(emo,clust) = emostats{clust}{rv}(nx,emo);
        end;
    end;
end;
[idx, C, sumd, d] = kmeans(colemo,3,'replicates',5);
% plot the results
emos = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excite'}; % for all new ones
figure; cols = jet(15);
for emo = 1:length(emos)
    ph = plot3(d(emo,1),d(emo,2),d(emo,3),'k.');hold on;
    set(ph,'color',cols(emo,:));set(ph,'markersize',25);
    ph = text(d(emo,1),d(emo,2),d(emo,3),emos{emo});
    set(ph,'color',cols(emo,:)); zl = get(gca,'zlim');
end;set(gca,'xgrid','on'); set(gca,'ygrid','on'); set(gca,'zgrid','on');
for emo = 1:length(emos)
    pl =plot3([d(emo,1) d(emo,1)],[d(emo,2) d(emo,2)],[zl(1) d(emo,3)]);
    set(pl,'color',cols(emo,:))             
end;

%%%%%  3) Find comps with similar patterns over emotions
%%!!!  need to select only components with some weighting in the given cluster
rv = 2; 
figure; pp = 1;
for clust = 1:4
    pl =1; clear concat
    for nx = 2:length(cluststats{clust})
        for cp = 1:length(gdcomps{nx})  % make this across subject
            colmeans(1,cp) = mean(cluststats{clust}{nx}{rv}(cp,:));
        end;
        for cp = 1:size(cluststats{clust}{nx}{rv},1)
            if ~isempty(zscore(colmeans>1.5)) & ismember(cp,find(zscore(colmeans)>1.25))
                concat(pl,:) = cluststats{clust}{nx}{rv}(cp,:);
                minitrack(pl,:) = [nx,gdcomps{nx}(cp)];pl  = pl+1;
            end;
        end;
    end;
    [idx, C, sumd, d] = kmeans(concat,3,'replicates',5);
    % what comps with whom?
    tpcp = cell(1,length(gdcomps));clear clustcomps
    for kclust = 1:max(idx)
        precutidx = find(idx == kclust);
        forcut = d(precutidx,kclust);    
        oneclust = minitrack(precutidx(find(forcut < mean(forcut) + std(forcut))),:);
        for addon = 1:size(oneclust,1)
            tpcp{oneclust(addon,1)}(end+1) = oneclust(addon,2);        
        end;
        clustcomps{kclust} = tpcp;
    end;
    %save /data/common2/emotion/clusters/tempkclusts.mat clustcomps
    % Plot the patterns across emotions
    %figure;
    for kclust = 1:max(idx)
        subplot(4,3,pp)
        %subplot(2,2,kclust)
        bar(mean(concat(find(idx == kclust),:),1));pp = pp+1;
        set(gca,'xlim',[0 16]);title(['Kclust: ',int2str(kclust)]);
    end;
end;
%%%%%%%%%  4) Cluster comps based on max RMS over all templates
pl = 1; = 1; % RMS
clustcomps = cell(1,24);
for mk = 1:length(clustcomps)
clustcomps{mk} = cell(1,length(gdcomps));
end;
for nx = 2:length(cluststats{clust})
    clear concat colmeans
    for cp = 1:size(cluststats{clust}{nx}{rv},1)
        for cc = 1:length(gdcomps{nx})  % make this across subject
            colmeans(1,cc) = mean(cluststats{clust}{nx}{rv}(cc,:));
        end; 
        for clust = 1:24
            concat(1,clust) = mean(cluststats{clust}{nx}{rv}(cp,:));
        end;        
        [val clust] = max(concat);
        if ~isempty(zscore(colmeans)>1.5) 
            clustcomps{clust}{nx}(end+1) = gdcomps{nx}(cp);
        end;
    end;
end;
save /data/common2/emotion/clusters/tempkclusts.mat clustcomps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find how many clusters each comp is in. 
whatclusts = cell(1,length(gdcomps));
for nx = 1:length(gdcomps)
    for cp = 1:length(gdcomps{nx})
    count = 0; whichclust = zeros(1,0);
        for clust = 1:length(clustcomps)
            if ~isempty(find(clustcomps{clust}{nx}==gdcomps{nx}(cp)))
                whichclust(end+1) = clust;
                count = count+1;
            end;
        end;
        compcount{nx}(cp) = count;
        whatclusts{nx}{cp} = whichclust;
    end;
end;
% plot results
collect = zeros(1,0);
for nx = 1:length(gdcomps)
    tp = compcount{nx}(find(compcount{nx}~=0));
    collect(end+1:end+length(tp)) = tp;
end;
figure;  hist(collect);
        
% which clusters go together?
q=1;qq=1;qqq=1;a=1;clear common2 common3 allcommon1  allcommon2 allcommon3
plot2 = zeros(9,9);plot2too = zeros(9,9);
for nx = 1:length(gdcomps)
p = 1;  pp = 1; ppp = 1; 
    for cp = 1:length(whatclusts{nx})
        if length(whatclusts{nx}{cp}) == 1
            allcommon1(a,:) = whatclusts{nx}{cp}; a =  a+1;
        end;        
        if length(whatclusts{nx}{cp}) == 2
            allcommon2(q,:) = whatclusts{nx}{cp}; q=  q+1;
            plot2(whatclusts{nx}{cp}(1),whatclusts{nx}{cp}(2))=plot2(whatclusts{nx}{cp}(1),whatclusts{nx}{cp}(2))+1;
            common2{nx}(p,:) = whatclusts{nx}{cp};p = p+1;            
        end;
        if length(whatclusts{nx}{cp}) == 3
            allcommon3(qq,:) = whatclusts{nx}{cp};qq=  qq+1;
            for x = 1:2
                for y=x+1:3
                plot2too(whatclusts{nx}{cp}(1,x),whatclusts{nx}{cp}(1,y)) = plot2too(whatclusts{nx}{cp}(1,x),whatclusts{nx}{cp}(1,y))+1;
                end;
            end;
            
            common3{nx}(pp,:) = whatclusts{nx}{cp}; pp = pp+1;
        end;
        if length(whatclusts{nx}{cp}) == 4
            allcommon4(qqq,:) = whatclusts{nx}{cp};qqq=  qqq+1;
            common4{nx}(ppp,:) = whatclusts{nx}{cp}; ppp = ppp+1;
        end;
        
    end;
end;
figure; imagesc(plot2);colorbar;
figure; imagesc(plot2too);colorbar;
[std,idx] = sort(allcommon2(:,1)); % puts pairs in order to see pattern
allcommon2 = allcommon2(idx,:);
[std,idx] = sort(allcommon3(:,1));
allcommon3 = allcommon3(idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Plot Dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find mean weight for each component for each cluster for density plot
clear compval
for clust=1:size(activations,1)
    for nx = 2:length(gdcomps) % no subject 1
        for cp = 1:length(gdcomps{nx})
            onecp = activations(clust,find(keeptrack(:,1) == nx&keeptrack(:,2) == gdcomps{nx}(cp)));  
            compval{clust}{nx}(cp) = mean(onecp);
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
for nx = 1:length(paths)  
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
% look at dipole locations;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/scott/matlab');
pcdims = 40;
figure; pp = 1;
row = 3; col = 4;
%row = round(sqrt(length(clustcomps))); col = ceil(sqrt(length(clustcomps))); 
cols = jet(12);%length(clustcomps));
    tl = ['Emotion Spectral Clusters; PCA to ',int2str(pcdims),'; RMS weights > mean+1.5*std of each template '];
for clust =1:12%length(clustcomps) 
    allbesa =  EEG.dipfit.model(1); colset = zeros(1,length(clustcomps{clust}));
    for ss = 2:length(clustcomps{clust})
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
    plotcol = cell(0); pl = 0;clear plotwts
    for ss = 2:length(clustcomps{clust})
        if colset(ss) ~= 0
            for nn = 1:colset(ss)
                plotcol{pl+nn} = dcols(ss,:);
                subjidxmat(pl+nn) = ss;
                %plotwts(pl+nn) = compval{clust}{ss}(nn);
            end;
            pl = pl+nn;
        end;
    end;
    %optdipplot = {allbesa,'dipolelength',0,'gui','off','dipolesize',18,'image','mri','spheres','on','color',plotcol};
    %figure; mydipoleentropy( optdipplot, subjidxmat, 'distance',10);  cbar;
    %figure; dipoledensity( optdipplot, 'subjind',subjidxmat, 'method','distance','methodparam',10, 'mri_view' ,'top', 'gridsize',16,'weights',plotwts); 
    %set(gcf,'color','w');    ph = textsc(['Subject dipole density (10); Cluster: ',int2str(clust)],'title');   set(ph,'fontsize',14);
    %savename = ['print -djpeg /home/julie/Scott/Spectral_Clusters/AllSubjClusterDensity-',int2str(clust),'.jpg'];
    %eval(savename)
    %close
%end;
    subplot(row,col,pp);
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(clust,:)},'projlines','on','projimg','on');pp = pp+1; camzoom(.9);
    ph=text(60,-100,145,['Clust: ',int2str(clust)]);
    set(ph,'color','y'); set(ph,'fontsize',14);
    view(60,20);
end;
set(gcf,'color','w');
ph = textsc(tl,'title');
set(ph,'color','r');set(ph,'fontsize',14);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjidxmat =

  Columns 1 through 13 

     4     4     5     6     7     9     9    11    11    11    11    11    11

  Columns 14 through 26 

    14    14    14    14    14    14    14    14    14    15    15    18    18

  Columns 27 through 39 

    18    18    18    18    19    19    20    20    20    20    20    20    20

  Columns 40 through 52 

    20    20    20    21    22    23    23    24    24    24    25    26    26

  Columns 53 through 61 

    29    30    30    30    30    30    31    31    31

    optdipplot = 

  Columns 1 through 6

    [1x61 struct]    'dipolelength'    [0]    'gui'    'off'    'dipolesize'

  Columns 7 through 13

    [18]    'image'    'mri'    'spheres'    'on'    'color'    {1x61 cell}
