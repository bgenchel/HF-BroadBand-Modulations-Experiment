% new approaches to clustering Co-modulation patterns
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
%str = ['load /data/common4/emotion/GoodComps.mat ']; eval(str);
savedat = 'SpecCoModMuscle';  % includes muscle and inf frontal
savedat = 'SpecCoModWave';
savedat = 'SpecCoModMoreFreqs';
savedat = 'SpecCoModWavePrePost'
savedat = 'SpecCoModNoFilt'; 
for nx =1:12
SpecCoModPlot('sources.set',newpaths{nx},[],[1:15],savedat,[3 128],'n',0,[]);
end;
nx4,ic2
nx21,ic1,3
nx=29,6

newdir = '/data/projects/julieo/emotion/';
for nx = 1:length(newpaths)
  npaths{nx} = [newdir,newpaths{nx}(end-4:end)];
end;
for nx =24:35
SpecCoModPlot('awe_NF.set',npaths{nx},[],[1:15],savedat,[0 128],'n',0,[]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cluster spectral co-mod templates across subjects
% use maxval > 3.25 and/or allrms > 1.25
clear data activations winv alldip kptk allbigs bigwts clustfacs mnspecs orivec
clustfacs = [];  mnspecs = []; kptk = [];
subjlist = [2:21,23:31,33:35];subjlist = [1:14,16:18];
freqscale = 'quad';
percmax =.5; % percent of max template to take as a 'comod'
 [clustfacs,templcell,mnspecs,kptk,allbigs,bigwts,onebig,orivec,modcorr,freqs,outcomods,rawrms,rawcorr] = CollectCoModTempls(savedat,npaths,subjlist,percmax);
% [clustfacs,templcell,mnspecs,kptk,allbigs,bigwts,onebig,orivec,modcorr,freqs,outcomods,rawrms,rawcorr] = CollectCoModTempls(savedat,fullpaths,subjlist,percmax);
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    [pcared{nx} pv{nx}] = SpecCoModPvaf(savedat,newpaths{nx},[],[],0,0);
end;
save /home/julie/tmppvaf.mat pcared pv

figure; set(gca,'fontsize',16);
badims = [];allrms = [];
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    allrms = [allrms,reshape(rawrms{nx},[1 size(rawrms{nx},1)*size(rawrms{nx},2)])];
    pvs = repmat(pv{nx}',[1 size(rawrms{nx},2)]);
    plot(rawrms{nx},pvs,'k.');hold on;
    %x=rawrms{nx} - rawcorr{nx} ;
    %x = max(x,[],2)';    
    %badims{nx} = find(x>1);  
end;
figure; hist(allrms,100);set(gca,'fontsize',16); set(gca,'ylim',[0 1750]);hold on;
xlabel('Relative RMS'); ylabel('Number of Templates');
plot([.5 .5],[get(gca,'ylim')],'r-','linewidth',3);
print /home/julie/Manuscripts/Gamma/Frontiers/RMS-Hist.jpg -djpeg

xlabel('Percent variance accounted for'); ylabel('Relative RMS');
%xlabel('Correlation with largest RMS template'); ylabel('Relative RMS');
plot([get(gca,'xlim')],[.5 .5],'r-','linewidth',3);
print /home/julie/Manuscripts/Gamma/Frontiers/RMS-CorrScatter.jpg -djpeg
     
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
SpecCoModPlot('sources.set',newpaths{nx},[],badims{nx},savedat,[3 125],'n',0,[]);
end;
p=1; clear outcomods
for percmax = .1:.1:.9
[clustfacs,templcell,mnspecs,kptk,allbigs,bigwts,onebig,orivec,modcorr,freqs,outcomods{p},rawrms,rawcorr] = CollectCoModTempls(savedat,newpaths,subjlist,percmax);p=p+1;
end;
origfreqs = freqs;
figure; cols = jet(length(outcomods)); 
set(gca,'fontsize',16);
for p = 1:length(outcomods)
    ph = plot(sort(outcomods{p}),'linewidth',2); hold on;
    set(ph,'color',cols(p,:));
end;
set(gca,'xlim',[0 length(outcomods{1})+1]);
legend({'10%','20%','30%','40%','50%','60%','70%','80%','90%'},'location','northwest');
ylabel('Percent of  "Comodulated" ICs per IM'); xlabel('IMs sorted by % ICs per IM > threshold');
print /home/julie/Manuscripts/Gamma/Frontiers/IMtemplateRMScutoff.jpg -djpeg

rmss=[];
for nxx = 1:length(subjlist)        
 nx = subjlist(nxx);
rmss= [rmss, reshape(rawrms{nx}(useims{nx},:),[1 length(useims{nx})*size(rawrms{nx},2)])];
end;
figure; plot(sort(rmss)); hold on;plot([get(gca,'xlim')],[.5 .5],'r-','linewidth',3);


%save /data/common1/emotion/HiFrqClustFacs.mat clustfacs templcell kptk freqs allbigs bigwts onebig mnspecs orivec modcorr origfreqs
%save /data/common1/emotion/AllClustFacs.mat clustfacs templcell kptk freqs allbigs bigwts onebig mnspecs orivec modcorr origfreqs
%save /data/common1/emotion/AllClustFacsWAVE.mat clustfacs templcell kptk freqs allbigs bigwts onebig mnspecs orivec modcorr origfreqs
%save /data/common1/emotion/AllClustFacsMuscle.mat clustfacs templcell kptk freqs allbigs bigwts onebig mnspecs orivec modcorr origfreqs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common1/emotion/AllClustFacs.mat    
load /data/common1/emotion/AllClustFacsWAVE.mat    
load /data/common1/emotion/AllClustFacsMuscle.mat    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster without sorting first:
% kptk -- [subject, IM, IC]
[ft,fm,fi,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(clustfacs,mnspecs,freqs,[3 128],kptk,[],34,1,4,'corr',freqscale);
keepclust = [1:25,27,29:32];
newtempls = zeros(0,length(freqs));
newmeans = zeros(0,length(freqs));
newidx=zeros(0,3);
for clss = 1:length(keepclust)
    cls = keepclust(clss);    
    newtempls(end+1:end+size(ft{cls},1),:) = ft{cls};
    newmeans(end+1:end+size(fm{cls},1),:) = fm{cls};
    newidx(end+1:end+size(fi{cls},1),:) = fi{cls};
end;
[finaltempls finalmeans finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(newtempls,newmeans,freqs,[3 128],newidx,[],20,1,4,'corr',freqscale);%  dotprod

save /data/common1/emotion/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first pass is to define delta, theta, alpha, beta or high freq:
% by finding freq of highest abs value of each template:
freqscale = 'quad';
[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,origfreqs,[3.5 120],[8.5 12.5 10],freqscale); % scale

[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,freqs,[0 128],[],freqscale); % no scale
[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,origfreqs,[3 125],[],freqscale); % no scale
%save /data/common1/emotion/SortedHiFrqClusts.mat gamaclust freqs
%save /data/common1/emotion/SortedClusts.mat deltaclust thetaclust alphaclust betaclust gamaclust freqs
%save /data/common1/emotion/SortedClustsWAVE.mat deltaclust thetaclust alphaclust betaclust gamaclust freqs origfreqs
%save /data/common1/emotion/SortedGammaClusts.mat gamaclust freqs
%save /data/common1/emotion/SortedMuscleClusts.mat deltaclust thetaclust alphaclust betaclust gamaclust freqs origfreqs

load /data/common1/emotion/SortedHiFrqClusts.mat gamaclust freqs
load /data/common1/emotion/SortedClustsWAVE.mat % regular decomp
load /data/common1/emotion/SortedClusts.mat % regular decomp
load /data/common1/emotion/SortedGammaClusts.mat %regular decomp, no freq shift
load /data/common1/emotion/SortedMuscleClusts.mat % muscle decomp (all shifted except gammaclust
freqs = origfreqs;
% cluster sorted templates:
% DELTA ******-----------------------------------------
%[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(deltaclust{1},deltaclust{2},freqs,[0 15],deltaclust{3},[],2,.5,3.5,'corr');
clear finaltempls finalidx finalmeans facvec
finaltempls{1} =deltaclust{1};finalidx{1} = deltaclust{3}; finalmeans{1} = deltaclust{2}; 
% eliminate by std cut:-----
% $$$ stdcut = 1.5; % 1.5 then 1
% $$$ zs = zscore(finaltempls{1});
% $$$ delmem = [];
% $$$ for mem = 1:size(finaltempls{1},1)
% $$$     if find(mean(abs(zs(mem,:))) > stdcut)
% $$$         delmem = [delmem mem];
% $$$     end;
% $$$ end;
% $$$ if ~isempty(delmem)
% $$$     finaltempls{1}(delmem,:) = [];
% $$$     finalmeans{1}(delmem,:) = [];
% $$$     finalidx{1}(delmem,:) = [];
% $$$ end;
%save /data/common1/emotion/DeltaClust.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/DeltaClustMusc.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/DeltaClustWave.mat finaltempls finalidx finalmeans freqs
     
% THETA   ******---------------------------
[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(thetaclust{1},thetaclust{2},freqs,[1 40],thetaclust{3},[],7,1,1.5,'corr',freqscale);

clear finaltempls finalidx finalmeans facvec
finaltempls{1} =thetaclust{1};finalidx{1} = thetaclust{3}; finalmeans{1} = thetaclust{2}; 
%save /data/common1/emotion/ThetaClust.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/ThetaClustMusc.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/ThetaClustWave.mat finaltempls finalidx finalmeans freqs
    
% ALPHA   ******-----------------------
[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(alphaclust{1},alphaclust{2},freqs,[1 40],alphaclust{3},[],6,1,1,'corr',freqscale);
[finaltempls finalmeans finalidx] = ClusterAlphaMods(alphaclust,freqs,[9 11],freqscale,1,2);
ph = textsc([savedat,' Alpha Clusters'],'title');set(ph,'fontsize',14);

%save /data/common1/emotion/AlphaClust.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/AlphaClustMusc.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/AlphaClustWave.mat finaltempls finalidx finalmeans freqs
[finaltempls finalmeans finalidx] = ClusterAlphaMods(alphaclust,freqs,[9.2 10.8],freqscale,1,4);
%save /data/common1/emotion/AlphaClustLongWave.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/AlphaClustWave.mat finaltempls finalidx finalmeans freqs
% BETA   ******--------------------------------
[ft,fm,fi,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(betaclust{1},betaclust{2},freqs,[1 40],betaclust{3},[],5,1,3.5,'corr',freqscale);%  dotprod ( put 3 and 5 together)
clear finaltempls finalidx finalmeans facvec
finaltempls{1} =[ft{3};ft{5};ft{2};ft{1}];finalidx{1}=[fi{3};fi{5};fi{2};fi{1}];finalmeans{1}=[fm{3};fm{5};fm{2};fm{1}]; 
finaltempls{2} =ft{4};finalidx{2} = fi{4}; finalmeans{2} = fm{4}; 
%save /data/common1/emotion/BetaClust.mat finaltempls finalidx finalmeans freqs
[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(betaclust{1},betaclust{2},freqs,[1 40],betaclust{3},[],5,1,1.5,'corr',freqscale);%  dotprod ( put 3 and 5 together)
finaltempls{1} =ft{1};finalidx{1} = fi{1}; finalmeans{1} = fm{1}; 
finaltempls{2} =ft{4};finalidx{2} = fi{4}; finalmeans{2} = fm{4}; 
finaltempls{3} =ft{2};finalidx{3} = fi{2}; finalmeans{3} = fm{2}; 
finaltempls{4} =ft{3};finalidx{4} = fi{3}; finalmeans{4} = fm{3}; 
%save /data/common1/emotion/BetaClustMusc.mat finaltempls finalidx finalmeans freqs
    
clear finaltempls finalidx finalmeans facvec
finaltempls{1} =betaclust{1};finalidx{1} = betaclust{3}; finalmeans{1} = betaclust{2}; 
%save /data/common1/emotion/BetaClustHighFrq.mat finaltempls finalidx finalmeans freqs


nclust = 20;perc = 5;freqtimes{1} = freqs; 
[IDX,newdata] = AffinityCluster(betaclust{1},nclust,perc,freqtimes,freqscale);



% GAMMA   ******-------------------------------------
% for brain decomp clustering, run clustcomodtempls from gamaclust sorted templates (non-shifted version). Then select non-noise templates from 100 preclusters. Then a simple sort by hi/low(1) and others(2) gives 2 broadband clusters.
[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(gamaclust{1},gamaclust{2},freqs,[2 128],gamaclust{3},[],50,2,2,'corr',freqscale);%  dotprod
keepclust = [1,10,11,12,16,30,31];
newtempls = zeros(0,size(gamaclust{1},2));
newmeans = zeros(0,size(gamaclust{2},2));
newidx=zeros(0,size(gamaclust{3},2));
for clss = 1:length(keepclust)
    cls = keepclust(clss);    
    newtempls(end+1:end+size(finaltempls{cls},1),:) = finaltempls{cls};
    newmeans(end+1:end+size(finalmeans{cls},1),:) = finalmeans{cls};
    newidx(end+1:end+size(finalidx{cls},1),:) = finalidx{cls};
end;
%[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(newtempls,newmeans,freqs,[2 128],newidx,[],length(keepclust),2,2.5,'corr',freqscale);%  dotprod
keepclust = [1,2,3];

clear finaltempls finalidx finalmeans facvec
finaltempls{1} =newtempls;finalidx{1} = newidx; finalmeans{1} = newmeans; 


% this is to cluster by hi/low and broadband by force
% $$$ fr1 = find(freqs<20); fr2 = find(freqs>45);
% $$$ relidx = zeros(0,size(kptk,2)); reltempls = zeros(0,size(newtempls,2)); relmeans = zeros(0,size(newmeans,2));
% $$$ relidx2 = zeros(0,size(kptk,2)); reltempls2 = zeros(0,size(newtempls,2)); relmeans2 = zeros(0,size(newmeans,2));
% $$$ for sp = 1:size(newidx,1)
% $$$     if mean(newtempls(sp,fr1))<0&mean(newtempls(sp,fr2))>0
% $$$         relidx(end+1,:) = newidx(sp,:);% hi/low cluster
% $$$         reltempls(end+1,:) = newtempls(sp,:);
% $$$         relmeans(end+1,:) = newmeans(sp,:);
% $$$     else % collect all others (not hi/low)
% $$$         relidx2(end+1,:) = newidx(sp,:); % regular broadband cluster
% $$$         reltempls2(end+1,:) = newtempls(sp,:);
% $$$         relmeans2(end+1,:) = newmeans(sp,:);        
% $$$     end;
% $$$ end;
% $$$ % eliminate by std cut:-----
% $$$ stdcut = 1.5;
% $$$ zs = zscore(reltempls2);
% $$$ delmem = [];
% $$$ for mem = 1:size(reltempls2,1)
% $$$     if find(mean(abs(zs(mem,:))) > stdcut)
% $$$         delmem = [delmem mem];
% $$$     end;
% $$$ end;
% $$$ if ~isempty(delmem)
% $$$     reltempls2(delmem,:) = [];
% $$$     relmeans2(delmem,:) = [];
% $$$     relidx2(delmem,:) = [];
% $$$ end;
% $$$ 
% $$$ clear finaltempls finalidx finalmeans facvec
% $$$ finaltempls{1} =reltempls;finalidx{1} = relidx; finalmeans{1} = relmeans; 
% $$$ finaltempls{2} =reltempls2;finalidx{2} = relidx2; finalmeans{2} = relmeans2; 


%save /data/common1/emotion/GammaClust.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/GammaClustWave.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/GammaClustMusc.mat finaltempls finalidx finalmeans freqs
%save /data/common1/emotion/GammaClustHighFrq.mat finaltempls finalidx finalmeans freqs
save /data/common1/emotion/GammaClust7-10-2009.mat finaltempls finalidx finalmeans freqs

%----------------------------------------------------------------------
% MUSCLE gamma diversion:--------------------------------
% go through gamaclust and separate out gdcomps, muscle and vent. frontal
freqs = origfreqs;
[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(gamaclust{1},gamaclust{2},freqs,[3.5 125],gamaclust{3},[],25,2,2.5,'corr',freqscale);%  dotprod
keepclust = [1,14,22,23];
newtempls = zeros(0,size(gamaclust{1},2));
newmeans = zeros(0,size(gamaclust{2},2));
newidx=zeros(0,size(gamaclust{3},2));
for clss = 1:length(keepclust)
    cls = keepclust(clss);    
    newtempls(end+1:end+size(finaltempls{cls},1),:) = finaltempls{cls};
    newmeans(end+1:end+size(finalmeans{cls},1),:) = finalmeans{cls};
    newidx(end+1:end+size(finalidx{cls},1),:) = finalidx{cls};
end;
clear finaltempls finalidx finalmeans facvec
finaltempls{1} =newtempls;finalidx{1} = newidx; finalmeans{1} = newmeans; 
%[finaltempls,finalmeans,finalidx,clustorigmeans,clustorigtempls,freqs,ndels] = ClustCoModTempls(newtempls,newmeans,freqs,[2 128],newidx,[],2,1,2.5,'corr','quad');%  dotprod

clustidx = 1; clear ft fm fi
for cls = 1:length(finaltempls)
    clear gamaclust
    % All gamma templs after first triage:
    gamaclust{1} = finaltempls{cls};
    gamaclust{2} = finalmeans{cls};
    gamaclust{3} = finalidx{cls};
    
    %% Search each dataset for brain/muscle/vf comps:
    clear deltaclust thetaclust alphaclust betaclust
    gdtempls = zeros(0,size(gamaclust{1},2));
    gdmeans = zeros(0,size(gamaclust{2},2));
    gdidx=zeros(0,size(gamaclust{3},2));
    mstempls = zeros(0,size(gamaclust{1},2));
    msmeans = zeros(0,size(gamaclust{2},2));
    msidx=zeros(0,size(gamaclust{3},2));
    vftempls = zeros(0,size(gamaclust{1},2));
    vfmeans = zeros(0,size(gamaclust{2},2));
    vfidx=zeros(0,size(gamaclust{3},2));
    EEG = pop_loadset('sources.set' ,newpaths{1});  nxx=1; 
    for sp = 1:size(gamaclust{3},1)
        nx = gamaclust{3}(sp,1);
        if nx ~= nxx
            EEG = pop_loadset('sources.set' ,newpaths{nx}); nxx = nx;
        end;
        if ~isempty(find(ismember(EEG.gdcomps, gamaclust{3}(sp,3))))& isempty(find(ismember(EEG.ventfrontal, gamaclust{3}(sp,3))))
            gdtempls(end+1,:) = gamaclust{1}(sp,:);
            gdmeans(end+1,:) = gamaclust{2}(sp,:);
            gdidx(end+1,:) = gamaclust{3}(sp,:);
        elseif ~isempty(find(ismember(EEG.muscle, gamaclust{3}(sp,3))))
            mstempls(end+1,:) = gamaclust{1}(sp,:);
            msmeans(end+1,:) = gamaclust{2}(sp,:);
            msidx(end+1,:) = gamaclust{3}(sp,:);
        elseif ~isempty(find(ismember(EEG.ventfrontal, gamaclust{3}(sp,3))))
            vftempls(end+1,:) = gamaclust{1}(sp,:);
            vfmeans(end+1,:) = gamaclust{2}(sp,:);
            vfidx(end+1,:) = gamaclust{3}(sp,:);
        end;
    end;
    ft{clustidx} = gdtempls; fm{clustidx} = gdmeans; fi{clustidx} = gdidx;
    ft{clustidx+1} = mstempls; fm{clustidx+1} = msmeans; fi{clustidx+1} = msidx;
    ft{clustidx+2} = vftempls; fm{clustidx+2} = vfmeans; fi{clustidx+2} = vfidx;
    clustidx = clustidx + 3;
end;

clear finaltempls finalmeans finalidx
% 1:3 = rapid rise, 4:6 = late rise, 7:9 = peaked (all = brain,ms,vf)
finaltempls = ft; 
finalmeans = fm; 
finalidx = fi;

% eliminate by std cut:-----
for cls = 1:length(finalidx)
    stdcut = 1.5; % 1.5 then 1
    zs = zscore(finaltempls{cls});
    delmem = [];
    for mem = 1:size(finaltempls{cls},1)
        if find(mean(abs(zs(mem,:))) > stdcut)
            delmem = [delmem mem];
        end;
    end;
    if ~isempty(delmem)
        finaltempls{cls}(delmem,:) = [];
        finalmeans{cls}(delmem,:) = [];
        finalidx{cls}(delmem,:) = [];
    end;
end;

%save /data/common1/emotion/GammaGdMsVf.mat finaltempls finalmeans finalidx freqs origfreqs

load /data/common1/emotion/GammaGdMsVf.mat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find N subjects who are in 1,2,3 alpha clusters:-----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load /data/common1/emotion/AlphaClustLongWave.mat finaltempls finalidx finalmeans freqs
%load /data/common1/emotion/AlphaClustWave.mat finaltempls finalidx finalmeans freqs
allNs = [];alleachs = [];
for nx = 1:length(gdcomps)
    neach = zeros(3,length(gdcomps{nx}));
    clustpart{nx} = zeros(1,length(gdcomps{nx}));
    for cls = 1:length(finalidx)
        oneclust = zeros(1,length(gdcomps{nx}));
        onesubj = finalidx{cls}(find(finalidx{cls}(:,1)==nx),:);
        for ic = 1:size(onesubj,1)            
            oneclust(ismember(gdcomps{nx},onesubj(ic,3))) = 1;
            neach(cls,ismember(gdcomps{nx},onesubj(ic,3))) = 1;
        end;
        clustpart{nx} = clustpart{nx}+ oneclust;
    end;
    alleachs=[alleachs,neach];
    allNs = [allNs,clustpart{nx}]; % collect all Ns
end;
sums = sum(alleachs,1);
for c = 1:2
    kind = find(sums == c);
    for cc = 1:3
        pieinfo(c,cc) = sum(alleachs(cc,kind));
    end;
end;
kind = find(sums == 2);% when two IMs
comb3 = [[1 1 0]',[0 1 1]',[1 0 1]'];
skind = alleachs(:,kind);
subpie = zeros(1,3);
for combs = 1:3
    for m = 1:size(skind,2)
        if skind(:,m) == comb3(:,combs)
            subpie(1,combs) = subpie(1,combs) + 1;
        end;
    end;
end;
figure; pie(subpie,[0 0 0],{'Low/Middle','Middle/High','Low/High'});

[hh,bins] = hist(allNs,4);
barcols = {'g','r',[.3 0 .7],[1 .6 0]};
figure;        ph = bar(pieinfo,'stacked');  hold on;
for b = 1:3:length(bins)
    ph = bar(b-1,hh(b),.7);hold on;
    set(ph,'facecolor',barcols{b});
end;
set(gca,'xticklabel',[]);set(gca,'xtick',[]);

str = ['print ',fullpaths{1}(1:end-5),savedat,'AlphaIMStats.eps -depsc -painters -adobe'];eval(str)
str = ['print ',fullpaths{1}(1:end-5),savedat,'AlphaIMSubPie.eps -depsc -painters -adobe'];eval(str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the clusters:*************************
clear facvec comods justcomps wtsmat1 jcwts denslist
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);
viewnum=[1,2,3,4];%length(viewnum) ;
zoom= 1;
figure;pl=1;row=6;col = 6;
for nx = 1:35
    try
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(ffs,allbigs,bigwts,orivec);
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot({finalidx{1}(find(finalidx{1}(:,1) == nx),:)},allbigs,bigwts,orivec);
clustcps{pl} = denslist{1};pl=pl+1;
    clust  =1;
    if pl > row*col
        figure; pl=1;
    end;
    sbplot(row,col,pl); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, denslist{1},[],[],[dipargs,'view',[0 0 1]],[],[],[],[]); camzoom(1.3);pl=pl+1;
    ph = title(['Subject ',int2str(nx)]);set(ph,'color','r');
     sbplot(row,col,pl); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, denslist{1},[],[],[dipargs,'view',[1 0 0]],[],[],[],[]); camzoom(1.3); pl=pl+1;
   %figure;[angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,1,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},[],[]); % next to last 1 plots solo IMs in black
    
end;
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/GammaIMdipoles1.jpg -djpeg
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/GammaIMdipoles2.jpg -djpeg

figure; PlotDipoleClusters('sources.set',fullpaths,gdcomps,clustcps,[1:length(clustcps)],row,col,1,'All Subj Gamma ICs',viewnum,[],'off');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/GammaIMdipolesAll.jpg -djpeg

%%%%%%%%%%**********************************************************************
% Plot scalp maps instead: (requires 'denslist' from loop above!)------------------
for clust = 1:length(denslist)
    PlotScalpMaps('sources.set',fullpaths,denslist{clust})
    textsc(['Cluster ',int2str(clust)],'title');
end;
%%%%%%%%%%**********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Plot dipole connections and relative template strength --------------
% plot one clust:
clust = 2; row = 2; col = 2; zoom = 1; viewnum=[1,2,3];
figure;[angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,1,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},[],[]); % next to last 1 plots solo IMs in black
textsc(['Cluster ',int2str(clust)],'title');
%%%------------------
row = 3%length(comods);
viewnum=[1,2,3];col = 3;%length(viewnum) ;
zoom= 1.1;
figure;pl = 1;
for clust = 1:length(comods)
figure;pl = 1;row=2;col=2;
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

str = ['print /home/julie/Manuscripts/Gamma/OMTdioples.tif -dtiff'];eval(str)

%%%DENSITY:------------------
row = 4; col = 3; 
sclusts = [1:3];pl=1;
for sc = 1:length(sclusts) 
    scls = sclusts(sc);
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,denslist{scls},[],[],[],{'mrislices',[63:-16:-17],'mriview','top','geom',[4,3]},'bred',[],'on');
    ph = textsc(['Emotion Alpha Clust ',int2str(scls)],'title');set(ph,'color','r')
str = ['print /home/julie/Meetings/CNS2009/EmoAlphaModDens',int2str(scls),'.jpg -djpeg'];eval(str)
end

save /home/julie/Meetings/CNS2009/EmoDens.mat dens

list_brain_areas( '/home/julie/Meetings/CNS2009/EmoDens.mat',.17,['/home/julie/Meetings/CNS2009/EmoClust',int2str(scls),'Dens.txt']);


[63:-8:-18]%-----------------------------------------------------------------------------------------------------
% plot spectral templates:***********
row=2;%round(sqrt(length(finaltempls))); 
col=2;%ceil(sqrt(length(finaltempls))); 
fr = find(freqs > 3.5 & freqs < 125);
figure; pl = 1; cls=2;
for cls = 1:length(finaltempls)
    sbplot(row,col,pl); pl = pl+1;
    [han,realx labelx] = quadplot(freqs(fr),finaltempls{cls}(:,fr)',2,[1 .5 0]); hold on;%[0 .75 .75][.2 1 .2]
    [han,realx labelx] = quadplot(freqs(fr),mean(finaltempls{cls}(:,fr),1)',3,'k');%[.16 .5 .3]
    %set(gca,'xlim',[realx(1) realx(end)]);
    set(gca,'ylim',[min(finaltempls{cls}(:)) max(finaltempls{cls}(:))]);
    set(gca,'ticklength',[.05 .05]);
    title(['Clust ',int2str(cls),' (',int2str(size(finaltempls{cls},1)),')']);      
    plot([realx(2) realx(2)],[get(gca,'ylim')],'g-');
    %plot([realx(4) realx(4)],[get(gca,'ylim')],'g-');
    plot([get(gca,'xlim')],[0 0],'k-');
end;
ph=textsc('Emotion Template clusters','title');set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
str = ['print ',fullpaths{1}(1:end-5),savedat,'AlphaTemplates.eps -depsc -painters -adobe'];eval(str)
%-----------------------------------------------------------------------------------------------------
%%%%%%%%%%**********************************************************************
% search all comods and see if any are between gdcomps and EEG.muscle-------
corrcut = .75; clear ngdcomps nmuscle nventfr gdmscomod gdvfcomod gdgdcomod comodics wtmat
for nx = 1:length(modcorr)
    s = load([newpaths{nx},savedat,'.mat']);  p1=1; p2=1; p3=1; 
    EEG = pop_loadset('sources.set' ,newpaths{nx});  
    for im = 1:length(modcorr{nx})
        cm = find(modcorr{nx}{im} > corrcut);
        if length(cm) > 1
            if length(find(ismember(EEG.gdcomps,s.complist(cm)))) > 0
                ngdcomps{nx}{im} = length(find(ismember(EEG.gdcomps,s.complist(cm))));
            else
                ngdcomps{nx}{im} = [];
            end;
            if length(find(ismember(EEG.muscle,s.complist(cm)))) > 0
                nmuscle{nx}{im} = length(find(ismember(EEG.muscle,s.complist(cm))));
            else
                nmuscle{nx}{im} = [];
            end;
            if length(find(ismember(EEG.ventfrontal,s.complist(cm)))) > 0
                nventfr{nx}{im} = length(find(ismember(EEG.ventfrontal,s.complist(cm)))); 
            else
                nventfr{nx}{im} = [];
            end;
            if ngdcomps{nx}{im} > 0&isempty(nmuscle{nx}{im}) &isempty(nventfr{nx}{im}) 
                gdgdcomod(nx,im) = 1;
                comodics{1}{nx}{p1} = s.complist(cm);
                justcmps{1}{nx} = [];
                justwts{1}{nx}=[];
                wtmat{1}{nx}{p1} =  ones(1,length(s.complist(cm))); p1 = p1+1;
            end;
            if ~isempty(ngdcomps{nx}{im}) & ~isempty(nmuscle{nx}{im}) 
                comodics{2}{nx}{p2} = s.complist(cm);
                justcmps{2}{nx} = [];
                justwts{2}{nx}=[];
                wtmat{2}{nx}{p2} =  ones(1,length(s.complist(cm))); p2 = p2+1;

                gdmscomod(nx,im) = 1;
            end;
            if ~isempty(ngdcomps{nx}{im}) & ~isempty(nventfr{nx}{im}) 
                gdvfcomod(nx,im) = 1;
                comodics{3}{nx}{p3} =  s.complist(cm); 
                justcmps{3}{nx} = [];
                justwts{3}{nx}=[];
                wtmat{3}{nx}{p3} =  ones(1,length(s.complist(cm))); p3 = p3+1;
            end;
        end;
    end;
end;
row = 3; col = 3; zoom = 1; viewnum=[1,2,3];
figure;clust = 3; place = 7;[angles] = PlotCoModasDipoles(comodics{clust},justcmps{clust},newpaths,'sources.set',row,col,place,zoom,0,viewnum,wtmat{clust},justwts{clust},1,[]); % next to last 1 plots solo IMs in black
str = ['print /data/common2/emotion/Figs/MuscleDecompComod.tif -dtiff'];eval(str)
str = ['print /data/common2/emotion/Figs/MuscleDecompBMcomod.tif -dtiff'];eval(str)
str = ['print /data/common2/emotion/Figs/MuscleDecompBVcomod.tif -dtiff'];eval(str)
lims=[-.5 1.4];
figure; 
sbplot(2,2,1);imagesc(gdgdcomod,lims);caxis([-2.5 1.5]);set(gca,'ticklength',[.04 .04]);
sbplot(2,2,2); imagesc(gdmscomod*-1,lims); caxis([-1.2 0.35]);set(gca,'ticklength',[.04 .04]);
sbplot(2,2,3);imagesc(gdvfcomod,lims);caxis([-.25 .8]);set(gca,'ticklength',[.04 .04]);
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

str = ['print /data/common2/emotion/Figs/MuscleDecompCoMods.eps -depsc -painters -adobe'];eval(str)

%%%----------------------------------------------------------------
%%%----------------------------------------------------------------

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Calculate correlations of wts/time between mod clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate trial to trial spectral modulation across components: 

shuffnum = 250; 
for nx = 1:length(fullpaths)
    [corr, bootstats] = CorrCoMod(fullpaths{nx},savedat,[1:15],shuffnum);
    subjcorr(:,:,nx) = corr;
end;
subjboots = bootstats;
%save /data/common4/emotion/IMCorrelations.mat subjcorr subjboots shuffnum
load /data/common4/emotion/IMCorrelations.mat subjcorr subjboots shuffnum

load /data/common4/emotion/AllCoModAlpha.mat
clslist1 = facvec;
load /data/common4/emotion/AllCoModBeta.mat
clslist1 = facvec; clslist1{end} = [];
load /data/common4/emotion/AllCoModGama.mat
clslist2 = facvec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear corrcell bootcellhi bootcelllo
for cls1 = 1:length(clslist1)
  for cls2 = 1:length(clslist2)
    onecorr = zeros(1,0);hiboot = zeros(1,0); loboot=zeros(1,0);
    for nx = 1:length(clslist1{cls1})
      snglims1 = unique(clslist1{cls1}{nx});
      snglims2 = unique(clslist2{cls2}{nx});
      if ~isempty(snglims1) & ~isempty(snglims2)
        for im1 = 1:length(snglims1)
          for im2 = 1:length(snglims2)
            onecorr(1,end+1) = corr{nx}(snglims1(im1),snglims2(im2));
            hiboot(1,end+1) = subjboots{nx}(snglims1(im1),snglims2(im2),1);
            loboot(1,end+1) = subjboots{nx}(snglims1(im1),snglims2(im2),2);
          end;
        end;
      end;
    end;
    corrcell{cls1,cls2} = onecorr;
    bootcellhi{cls1,cls2} = hiboot;
    bootcelllo{cls1,cls2} = loboot;
  end;
end;

figure; pl = 1; row = 3; col = 3;
for cls1 = 1:size(corrcell,1)
    for cls2 = 1:size(corrcell,2)
        sbplot(row,col,pl)
        [srtcorr srtidx] = sort(corrcell{cls1,cls2});
        plot(srtcorr,'g-');hold on;
        %plot(corrcell{cls1,cls2},'g-');
        plot(bootcellhi{cls1,cls2}(srtidx),'r-'); hold on;
        plot(bootcelllo{cls1,cls2}(srtidx),'b-');
        set(gca,'xlim',[0 length(srtcorr)+1]);
        title(['Alpha ',int2str(cls1),' vs Gama ',int2str(cls2)]);
        pl = pl+1;
    end;
end;
axcopy
