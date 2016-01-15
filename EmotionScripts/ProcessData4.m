% this is for after the emotion space decomp to discover what 
% median-weighted IM templates look like for each emotion.

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
savedat = 'SpecCoMod'; fullpaths = newpaths;
% see bottom of script for other subject lists
subjlist = [2:21,23:31,33:35];% all but 1,22,32

% for regular gdcomps decomposition:--------------------
load /data/common1/emotion/EmoWeights.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacs.mat    
load /data/common1/emotion/SpectralDecompStuff.mat 
strs = {
    'load /data/common1/emotion/ThetaClust.mat',
    'load /data/common1/emotion/AlphaClust.mat',
    'load /data/common1/emotion/BetaClust.mat',
    'load /data/common1/emotion/GammaClust.mat'};
ttl = {'Theta','AlphaLow','AlphaPeak','AlphaHigh','BetaLow','BetaHigh','Gamma'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for muscle decomposition:-------------------
load /data/common1/emotion/EmoWeightsMuscle.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacsMuscle.mat    
load /data/common1/emotion/SpectralDecompStuffMusc.mat 
strs = {'load /data/common1/emotion/GammaGdMsVf.mat'};
ttl = {'Brain','Muscle','OMT'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for Gamma cluster at least
 decomp = 7;
pl = 8;
 [fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,oricell,4,subjlist,useimssub,'pdist',0,['Gamma clusters'],0,'elim','noperm');
decomp = 8; pl = 9;
for subcls = 1:length(keeptrack{decomps})
    oricell = cell(1,max(subjlist)); 
    useimssub = cell(1,max(subjlist)); 
    for nx = 1:max(keeptrack{decomp}{subcls}(:,1))
        onesubj = keeptrack{decomp}{subcls}(find(keeptrack{decomp}{subcls}(:,1) == nx),2);
        if ~isempty(onesubj)
            pics{pl}{nx} = keeptrack{decomp}{subcls}(find(keeptrack{decomp}{subcls}(:,1) == nx),:);
            useimssub{nx} = unique(abs(onesubj))';
            for im = 1:length(useimssub{nx})
                idx = find(allbigs{nx}{useimssub{nx}(im)} == onebig{nx}{useimssub{nx}(im)});
                oricell{nx}(im) = orivec{nx}{useimssub{nx}(im)}(idx);
            end;
        end;
    end;
    [fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useimssub,'mds',0,['subcls ',int2str(subcls)],0,'elim','noperm'); pl = pl+1;
end;

decomp = 8; clear emotempls
 figure; row = 2; col = 3;
 for subcls = 1:length(keeptrack{decomp})
     sbplot(row,col,subcls);
     for e = 1:size(collmeans{decomp},1)
         alltempls = []; 
         for im = 1:size(keeptrack{decomp}{subcls},1)
             nx = keeptrack{decomp}{subcls}(im,1); 
             nxim = keeptrack{decomp}{subcls}(im,2); 
             alltempls = [alltempls;templcell{nx}(nxim,:)*fullwts{decomp}{subcls}(e,im)]; %weight by median wt
         end;
         emotempls{e} = alltempls;
         ph = logplot(origfreqs,mean(emotempls{e},1),1.5,cols(e,:)); hold on;
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-');
 end;
 
%---------------------------------------------
% plot dipoles color-coded by weights for a dimension

clustnames = {'Theta','AlphaLo','AlphaPeak','AlphaHi','BetaLo','BetaHi','Gamma'};
clustnames = {'Brain','Muscle','OMT'}; % muscle
dens = 'on'; % 'neg' or 'pos' weightings, [] to plot color-coded dipoles

cut = 0; % * std
clustinfo = 'load /data/common1/emotion/AllClustFacs.mat';    
realclusts = [1,2,2,2,3,3,4]; subclusts = [1,1,2,3,1,2,1]; % brain only

%clustinfo = 'load /data/common1/emotion/AllClustFacsMuscle.mat';  % muscle   
%realclusts = [1,1,1]; subclusts = [1,2,3]; % muscle

decomp 7;dim=1;
specclust = strs{realclusts(decomp)};
subclust = subclusts(decomp);
for e = 1:15
    EmoDimDipoles(subjlist,fullpaths,fullwts{decomp}(e,:),fullmd{decomp}(:,dim),keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens);        

    ph=textsc([clustnames{decomp},'-',emos{e},'-Dipole weightings -- Dim ',int2str(dim)],'title'); set(ph,'color','r');
    %ph=textsc([clustnames{decomp},' ',emos{e},' Dipole weightings -- Dim ',int2str(dim),' Masked by permuting weights only'],'title'); set(ph,'color','r');
 set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
    str = ['print /home/julie/Manuscripts/Gamma/',clustnames{decomp},'-',emos{e},'-EmoMeanWtdDens.tif -dtiff']; eval(str);
end;
