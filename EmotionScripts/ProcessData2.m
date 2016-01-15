%%  follows from PreProcess.m to run emotion analysis on spectral data

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
%savedat = 'SpecCoModWave'; fullpaths = newpaths;
savedat = 'SpecCoModMoreFreqs'; fullpaths = newpaths;
%savedat = 'SpecCoModWaveTest'; fullpaths = newpaths;

% see bottom of script for other subject lists
subjlist = [2:21,23:31,33:35];% all but 1,22,3
2
% for regular gdcomps decomposition:--------------------
load /data/common1/emotion/EmoWeights.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacs.mat    
strs = {
    'load /data/common1/emotion/DeltaClust.mat',
    'load /data/common1/emotion/ThetaClust.mat',
    'load /data/common1/emotion/AlphaClust.mat',
    'load /data/common1/emotion/BetaClust.mat',
    'load /data/common1/emotion/GammaClust.mat'};
basettls = {'Delta','Theta','Alpha','Beta','Gamma'};ttl = cell(1,0);
for s = 1:length(strs)
    eval(strs{s})
    for cls = 1:length(finalidx)
        ttl{end+1} = [basettls{s},int2str(cls)];
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for muscle decomposition:-------------------
load /data/common1/emotion/EmoWeightsMuscle.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacsMuscle.mat    
strs = {'load /data/common1/emotion/GammaGdMsVf.mat'};
ttl = {'Brain','Muscle','OMT'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear percnt
for nx = 1:35
    nims = [];
   for c = 1:length(gdcomps{nx})
        ims = find(finalidx{1}(:,1) == nx & finalidx{1}(:,3) == gdcomps{nx}(c));
        if ~isempty(ims)
            if length(ims) == 2
            nims = [nims length(ims)];
            end;
        end;
   end;
   percnt(nx) = round(100*length(nims)/length(gdcomps{nx}));
end;

for nx = 1:length(gdcomps)
    for ic = 1:length(gdcomps{nx})
        icc=gdcomps{nx}(ic);
        if ~isempty(find(ismember(icc,finalidx{1}(find(finalidx{1}(:,1)==nx&finalidx{1}(:,3)==icc),:,:))))
            mems{nx}(ic) = 1;
        else
            mems{nx}(ic) = 0;
        end;
    end;
end;
counts = [];
for nx = 1:length(gdcomps)
    if ismember(nx, subjlist)
    counts = [counts mems{nx}];
    end
end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find lists of clustered IMs for Emo Space decomp
allidx = []; whichcls = []; plotidx = cell(1,5);
for s = 1:length(strs)
    eval(strs{s})
    for cls = 1:length(finalidx)
        plotidx{s} = finalidx{cls};
        allidx = [allidx;finalidx{cls}];
        whichcls = [whichcls;[repmat(s,[size(finalidx{cls},1) 1]) repmat(cls,[size(finalidx{cls},1) 1])]];
    end;
end;
% find indices of all clustered IMs-----------------------
onesubjics = cell(1,max(subjlist)); 
whichcell = cell(1,max(subjlist)); 
oricell = cell(1,max(subjlist)); 
useims = cell(1,max(subjlist)); 
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    onesubjics{nx} = allidx(find(allidx(:,1) == nx),:);% all clustered templates
    whichcell{nx} = whichcls(find(allidx(:,1) == nx),:);% corresponding clust/subclust
    onesubjims = allidx(find(allidx(:,1) == nx),2);
    useims{nx} = unique(abs(onesubjims)');% list of clustered IMs for each subj
    for im = 1:length(useims{nx})
        idx = find(allbigs{nx}{useims{nx}(im)} == onebig{nx}{useims{nx}(im)});
        oricell{nx}(im) = orivec{nx}{useims{nx}(im)}(idx);
    end;
end;
pl = length(fullmd)+1;
[fullmd{pl},fullwts{pl},keeptrack{pl},collmeans,nsteps] = EmoSpace(emomeans,oricell,10,subjlist,useims,'pdist',0,['All IMs'],0,'elim','noperm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=2; 
[outmat emovec] = PredictEmoConstruct(savedat,fullpaths{nx},[1:20],'valence',{'compassion'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[emoorders] = FindEmoOrder(fullpaths,emos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put together emodeciles for only one cluster at a time:---
    eval(strs{s})
    bbdeciles = cell(1,max(subjlist)); 
    useimssub = cell(1,max(subjlist)); 
        cls=1;
        for nx = 1:max(finalidx{cls}(:,1))
            onesubj = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
            if ~isempty(onesubj)
                bbdeciles{nx} = zeros(0,size(emodeciles{nx},2),size(emodeciles{nx},3));
                useimssub{nx} = unique(abs(onesubj))';
                for im = 1:length(useimssub{nx})
                    bbdeciles{nx}(end+1,:,:) = emodeciles{nx}(useimssub{nx}(im),:,:);
                end;
            end;       
        end;
 %save /data/common1/emotion/BBdeciles.mat bbdeciles
%load /data/common1/emotion/BBdeciles.mat bbdeciles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makes a IM x emos (x quadrants) matrix/subj
[emodensity] = WtsDensity(savedat,fullpaths,gdcomps,[1:35]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makes a IM x emos (x decile) matrix/subj
nsamples = []; % only look at first x samples
[emomeans,emodeciles,subjpoints] = EmoWeights(savedat,fullpaths,gdcomps,[1:35],nsamples);
%save /data/common1/emotion/EmoWeights.mat emomeans emodeciles
%save /data/common1/emotion/EmoWeightsMuscle.mat emomeans emodeciles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selemos = {'anger','frustration'};
PairEmos(emomeans,subjlist,useims,selemos,emos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useims,'mds',0,['All IMs'],0,'elim','noperm');      
for nx = 1:max(keeptrack(:,1))
    onesubj = find(keeptrack(:,1) == nx);
    if ~isempty(onesubj)
        subjmeds{nx} = collmeans(:,onesubj);
    end;
end;

   
emo2 = [emo2 ' prebase', ' postbase'];

pl = 1; clear ktsubj fullmd fullwts keeptrack nsteps collmeans pics useimssub templs
% keep each subcluster separate:--------------
for s = 1:length(strs)
    eval(strs{s})
    useimssub = cell(1,35);  k = [];
    for cls = 1:length(finalidx)
        oricell = cell(1,max(subjlist)); 
        useimssub = cell(1,max(subjlist)); 
        for nx = 1:max(finalidx{cls}(:,1))
            onesubj = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
            if ~isempty(onesubj)
                pics{pl}{nx} = finalidx{cls}(find(finalidx{cls}(:,1) == nx),:);
                templs{pl}{nx} = finaltempls{cls}(find(finalidx{cls}(:,1) == nx),:);
                useimssub{nx} = unique(abs(onesubj))';
                for im = 1:length(useimssub{nx})
                    idx = find(allbigs{nx}{useimssub{nx}(im)} == onebig{nx}{useimssub{nx}(im)});
                    oricell{nx}(im) = orivec{nx}{useimssub{nx}(im)}(idx);
                end;
                k=[k; repmat(nx,[length(useimssub{nx}) 1]) useimssub{nx}'];
            end;
        end;
        ktsubj{pl} = k;
        [fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useimssub,'mds',0,['gamma'],0,'keep','noperm');      
        %[fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useimssub,'mds',0,[ttl{pl},'; Sub-Clust ',int2str(cls)],0,'elim','noperm');      
        %set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        %str = ['print /home/julie/Manuscripts/Emotion/figures/',ttl{pl},'EmoSpaceTemplates.jpg -djpeg']; eval(str)    
        pl = pl+1;
    end;
end;
close all
%%%-----------
%save /data/common1/emotion/SpectralDecompStuff.mat collmeans fullmd fullwts ktsubj keeptrack templs pics subjlist ttl strs 
%save /data/common1/emotion/SpectralDecompStuffMusc.mat collmeans fullmd fullwts ktsubj keeptrack templs pics subjlist ttl strs 
load /data/common1/emotion/SpectralDecompStuff.mat 
load /data/common1/emotion/SpectralDecompStuffMusc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test emotion space by jack-knife:--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=1;
for nd = .05:.1:.9
    [keeptemplates{r},keepwts{r}] = EmoSpaceStats(emomeans,oricell,[],subjlist,useimssub,'mds',0,'elim',500,nd);r=r+1;
    figure;
    for e = 1:size(keeptemplates,1)
        if length(dims)>2
            ph=plot3(squeeze(keeptemplates(e,dims(1),:)),squeeze(keeptemplates(e,dims(2),:)),squeeze(keeptemplates(e,dims(3),:)),'.');hold on;
        else            
            ph=plot(squeeze(keeptemplates(e,dims(1),:)),squeeze(keeptemplates(e,dims(2),:)),'.');hold on;
        end;
        set(ph,'markersize',5);set(ph,'color',cols(e,:));
    end;
    for e = 1:size(keeptemplates,1)
        if length(dims)>2
            ph=plot3(mean(keeptemplates(e,dims(1),:),3),mean(keeptemplates(e,dims(2),:),3),mean(keeptemplates(e,dims(3),:),3),'k.');hold on;
        else            
            ph = plot(mean(keeptemplates(e,dims(1),:),3),mean(keeptemplates(e,dims(2),:),3),'k.');
        end;
        set(ph,'markersize',40);
        if length(dims)>2
            ph = text(mean(keeptemplates(e,dims(1),:),3),mean(keeptemplates(e,dims(2),:),3),mean(keeptemplates(e,dims(3),:),3),emo2{e});
        else            
        ph = text(mean(keeptemplates(e,dims(1),:),3),mean(keeptemplates(e,dims(2),:),3),emo2{e});
        end;
        set(ph,'fontsize',20); set(ph,'color',cols(e,:));
    end;
    xlabel(['Dim ',int2str(dims(1))]);ylabel(['Dim ',int2str(dims(2))]);
    if length(dims) > 2
        zlabel(['Dim ',int2str(dims(3))]);
        zl = get(gca,'zlim');
        for e = 1:size(fullmd{clust},1)
            pl =plot3([mean(keeptemplates(e,dims(1),:),3) mean(keeptemplates(e,dims(1),:),3)],[mean(keeptemplates(e,dims(2),:),3) mean(keeptemplates(e,dims(2),:),3)],[zl(1)  mean(keeptemplates(e,dims(3),:),3)]);
            set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
        end;
        set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
    else        
        plot([0 0],[get(gca,'ylim')],'k-'); plot([get(gca,'xlim')],[0 0],'k-');
    end;
    title(['MD scaling of ',num2str(1-nd),'% of the data']);
    %str = ['print /home/julie/Manuscripts/Gamma/JackKnifeMDS',num2str(nd),'.jpg -djpeg'];eval(str);
    %str = ['print /home/julie/Manuscripts/Gamma/JackKnifeMDS',num2str(nd),'.eps -depsc -adobe -painters'];eval(str);    
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=1; clear rrs keeptemplates keepwts
dims = [1,2,3];clust = 11;
for nd = .05:.1:.9
    [keeptemplates{r},keepwts{r}] = EmoSpaceStats(emomeans,oricell,[],subjlist,useimssub,'mds',0,'elim',500,nd);
    % find correlation of each iteration with full decomp:
    for i = 1:size(keeptemplates{r},3)
        [rrs(r,i)] = corr2(keeptemplates{r}(:,dims,i),fullmd{clust}(:,dims));
    end;r=r+1;
end;
n=10;
figure; hist(keepwts{n},50)
markers = {'-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.'};
markcols = jet(size(rrs,1));
figure;
for nd=1:size(rrs,1)
    ph = plot(sort(abs(rrs(nd,:))),markers{nd},'linewidth',2);hold on;
    set(ph,'color',markcols(nd,:));
end;
legend({'5%','15%','25%','35%','45%','55%','65%','75%','85%','95%'},'location','southeast');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sig diff in weights between emos? :--------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[emomeans,emodeciles,subjpoints] = EmoWeights(savedat,fullpaths,gdcomps,[1:35]);% 
 pval = .00001; clear sallPs sallFs
for s = 1:length(strs)
    eval(strs{s})
    useimssub = cell(1,35);  k = []; clear bootcorr 
     clear allFs allPs
    for cls = 1:length(finalidx)
        keepF = []; keepP = [];emopoints = [];
        for nx = 1:max(finalidx{cls}(:,1))
            onesubj = unique(abs(finalidx{cls}(find(finalidx{cls}(:,1) == nx),2)));            
            for im = 1:length(onesubj)
                dim = onesubj(im);
                [P,atab,STATS] = anova1(subjpoints{nx}{dim},emos,'off');close; close;
                keepF = [keepF,atab{2,5}];
                keepP = [keepP,atab{2,6}];
                emopoints = [emopoints;subjpoints{nx}{dim}];
            end;            
        end;
        [P,atab,STATS] = anova1(emopoints,emos,'off');
        multcompare(STATS,'alpha',.01);
        clssig(1,cls) = atab{2,6};
        saveemos{cls} = emopoints;
        allFs{cls} = keepF;
        allPs{cls} = keepP;
    end;
    sallFs{s} = allFs;
    sallPs{s} = allPs;
end;
 percsig = cell(1,5);
for s = 1:length(strs)
    for cls = 1:length(sallPs{s})
        percsig{s} = [percsig{s} length(find(sallPs{s}{cls} < pval))/length(sallPs{s}{cls})];    
    end;
end;
figure;cols = jet(size(emopoints,2));
for e = 1:size(emopoints,2)
    ph = plot3(median(saveemos{1}(:,e)),median(saveemos{2}(:,e)),median(saveemos{3}(:,e)),'.','markersize',25);hold on;   set(ph,'color',cols(e,:));
    ph = text(median(saveemos{1}(:,e)),median(saveemos{2}(:,e)),median(saveemos{3}(:,e)),emo2{e});
    set(ph,'color',cols(e,:));set(ph,'fontsize',20);
end;
zl = get(gca,'zlim');
for e = 1:size(emopoints,2)
    pl =plot3([median(saveemos{1}(:,e)) median(saveemos{1}(:,e))],[median(saveemos{2}(:,e)) median(saveemos{2}(:,e))],[zl(1)  median(saveemos{3}(:,e))]);
    set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel('Low Alpha'); ylabel('Mid Alpha'); zlabel('High Alpha')
        [P,atab,STATS] = anova1(saveemos{2},emos,'off');
        figure; multcompare(STATS)
str = ['print /data/common1/emotion/AlphaMidANOVA.eps -depsc -adobe -painters']; eval(str)

% make a 3D plot of progressive correlations with less and less data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine subclusters instead:--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttl = {'Theta','Alpha','Beta','Gamma','All'};
pl = 1; clear ktsubj fullmd fullwts keeptrack nsteps templs pics

    useims = cell(1,max(subjlist)); k = [];onesubjics = [];
    oricell = cell(1,max(subjlist)); 
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        ix = [];tts=[]; ui=[];
        for cls = 1:length(finalidx)
            tts = [tts;finaltempls{cls}(find(finalidx{cls}(:,1) == nx),:)];
            onesubjics=[onesubjics; finalidx{cls}(find(finalidx{cls}(:,1) == nx),:)];
            onesubjims = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
            ui = [ui, unique(abs(onesubjims))'];% unique b/c only IMs
                for im = 1:length(onesubjims)
                    idx = find(allbigs{nx}{onesubjims(im)} == onebig{nx}{onesubjims(im)});
                    oricell{nx} = [oricell{nx} orivec{nx}{onesubjims(im)}(idx)];
                end;
                % more stringent:-----
            %delu = [];
            %sd = mean(std(emomeans{nx}'))+2*std(std(emomeans{nx}'));
            %for im = 1:length(ui)
            %    if std(emomeans{nx}(ui(im),:)) < sd
            %        delu = [delu im];
            %    end;                
            %end;
            %tts(delu,:) = [];            
            %onesubjics(delu,:) = [];
            %ui(delu) = [];
            
        end;
        ui = unique(ui);
        useimssub{nx} = ui;
        pics{pl}{nx} = onesubjics;
        templs{pl}{nx} = tts;        
        k=[k; [repmat(nx,[length(ui) 1]) ui']];
    end;
    ktsubj{pl} = k;
    if length(k) > 15        
        [fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useimssub,'mds',0,['All clusts'],0,'elim','noperm');      
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

% decomp with all clustered IMs (optional):-----------------------
[fullmd{pl},fullwts{pl},keeptrack{pl},collmeans{pl},nsteps] = EmoSpace(emomeans,[],subjlist,useims,'ica',0,['All Clustered IMs'],0);
pics{pl} = onesubjics; % if adding a decomp with all ims

%save /data/common1/emotion/SpectralDecompCombine.mat fullmd fullwts ktsubj keeptrack templs pics subjlist ttl strs flipvec whichcell 
pl = 1; clear ktsubj fullmd fullwts keeptrack nsteps collmeans pics useimssub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%collect permuted emotion space decomps for bootstrap limits:--------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shuffnum = 500;
%w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
w=load('/data/common1/emotion/BehavRatings/ZscoreTally.mat');
useemos = [1:7,9:15];% take out compassion for correlation
s=5;
%for s = 1:length(strs)
    eval(strs{s})
    useimssub = cell(1,35);  k = []; clear bootcorr 
    for cls = 1:length(finalidx)
        clear oricell useimssub
        for nx = 1:max(finalidx{cls}(:,1))
            onesubj = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
            useimssub{nx} = unique(abs(onesubj))';
            for im = 1:length(useimssub{nx})
                idx = find(allbigs{nx}{useimssub{nx}(im)} == onebig{nx}{useimssub{nx}(im)});
                oricell{nx}(im) = orivec{nx}{useimssub{nx}(im)}(idx);
            end;
        end;
        for b = 1:shuffnum
            [bootmd{b},bootwts{b},kt,collmeans,nsteps] = EmoSpace(emomeans,oricell,[],subjlist,useimssub,'mds',0,[''],0,'elim','perm');
            [bb,bint,r,rint,stats] = regress(w.emoval(useemos)',[bootmd{b}(useemos,1),bootmd{b}(useemos,2),ones(length(useemos),1)]);
            rcorr(b) = sqrt(stats(1));
            bootcorr{cls}(1,b) = corr2(w.emoval(useemos),bootmd{b}(useemos,1)');
        end;
    end;
%end;
figure; hist(bootcorr{cls},100); xlabel('Correlation coeff with subj ratings');
title('Brain Gamma IM decomps with permuted emotions');
str = ['print /home/julie/Manuscripts/Gamma/BootVal-GammaDim1CorrDist.jpg -djpeg']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/BootVal-GammaDim1CorrDist.eps -depsc -painters -adobe']; eval(str)
sortdist = sort(bootcorr{1});
sortdist(.01*shuffnum)
sortdist(end-.01*shuffnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot templates vs subjective ratings and calculate correlation
%w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
w=load('/data/common1/emotion/BehavRatings/ZscoreTally.mat');
erating = w.emoval;rw = 1;  % w.emoval or w.emoactiv
%erating = w.emoactiv; rw = 2; % w.emoval or w.emoactiv
clear cr
useemos = [1:7,9:15];% take out compassion for correlation
str = ['load ', fullpaths{1}(1:end-5),'AllSubjHeartInfo.mat hrates']; eval(str);
for clust = 1:length(fullmd)
    figure; row = 3; col = 3;
    % 2D correlation with dims 1 and 2:
    [b,bint,r,rint,stats] = regress(erating(useemos)',[fullmd{clust}(useemos,1),fullmd{clust}(useemos,2),ones(length(useemos),1)]);
    % 3D correlation with dims 1,2,3:
    [b3,bint3,r3,rint3,stats3] = regress(erating(useemos)',[fullmd{clust}(useemos,1),fullmd{clust}(useemos,2),fullmd{clust}(useemos,3),ones(length(useemos),1)]);
    for dim = 1:size(fullmd{clust},2)
        sbplot(row,col,dim)
        for e = 1:size(fullmd{clust},1)
            ph=plot(erating(e),fullmd{clust}(e,dim),'.');hold on;
            set(ph,'markersize',20);set(ph,'color',cols(e,:));
            ph = text(erating(e),fullmd{clust}(e,dim),emo2{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        end;
        plot([5 5],[get(gca,'ylim')],'k-');
        plot([get(gca,'xlim')],[0 0],'k-');
        cr(1,dim) = corr2(erating(useemos),fullmd{clust}(useemos,dim)');
        title(['Correlation ',num2str(cr(1,dim))]);
        xlabel('Subjective emotion ratings');
        ylabel([ttl{clust},' Dim ',int2str(dim),' Template']);
        % correlate with subj heart rate info
        for nx = 1:size(hrates,1)
            if ~isempty(find(hrates(nx,:)))
            hrcr(nx,dim) = corr2(hrates(nx,useemos),fullmd{clust}(useemos,dim)');
            end;
        end;  
    end;  
    [clustcorrs(rw,clust) corridx(rw,clust)] = max(abs(cr));
    textsc([ttl{clust},' correlations with subjective ratings'],'title');
end;
 close all
str = ['print /home/julie/Manuscripts/Gamma/Val-GammaDim1Correl.jpg -djpeg']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/Val-GammaDim1Correl.eps -depsc -painters -adobe']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/Val-BrainGammaDim1Correl.jpg -djpeg']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/Val-BrainGammaDim1Correl.eps -depsc -painters -adobe']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/Val-MuscleGammaDim1Correl.jpg -djpeg']; eval(str)
str = ['print /home/julie/Manuscripts/Gamma/Val-MuscleGammaDim1Correl.eps -depsc -painters -adobe']; eval(str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT Emotion space of specified dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clustnames = {'Delta','Theta1','Theta2','Theta3','AlphaLow','AlphaPeak','AlphaHigh','Beta1','Beta2','Beta3','Beta4','Gamma'};
%figure;
for clust = 1:length(fullmd)  % plot in 2D:----
         dims = corridx(:,clust)';         
         %figure;
         sbplot(3,4,clust)
         for e = 1:size(fullmd{clust},1)
             ph=plot(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),'.');hold on;
             set(ph,'markersize',10);set(ph,'color',cols(e,:));
             ph = text(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),emo2{e});
             set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
         end;
         xlabel(['Dim ',int2str(dims(1)),' (Valence: ',num2str(clustcorrs(1,clust)),')']);ylabel(['Dim ',int2str(dims(2)),' (Arousal: ',num2str(clustcorrs(2,clust)),')']);
         plot([0 0],[get(gca,'ylim')],'k-');         plot([get(gca,'xlim')],[0 0],'k-');
        title([clustnames{clust}]);
        %title([clustnames{clust},' Valence vs Arousal EmoSpace']);
         %str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{clust},'VvsAemospaceICA.tif -dtiff']; eval(str);
end;
textsc('ICA decomposition Valence vs Arousal dimensions','title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print ',fullpaths{1}(1:end-5),'results/AllVvsAemospace-ICA.eps -depsc']; eval(str);
str = ['print ',fullpaths{1}(1:end-5),'results/AllVvsAemospace-ICA.tif -dtiff']; eval(str);



clust = 12;% plot 3D emo space as spheres/cylinders:-------
dims = [1,2,3]; % 2
emo3 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sadness','  compassion','  love','  relief','  contentment','  awe','  happiness','  joy','  excitement'}; % for plotting purposes
% adjust color scale:----
w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
[xx yy] = sort(w.emoval);
emosval = emos{yy};
evals = round(w.emoval*100);
evals = evals-min(evals);evals = evals+1;% put all evals between 1 and 

cols = jet(max(evals));
cols(round(size(cols,1)/2) - (round(size(cols,1)/4)-1):round(size(cols,1)/2) + (round(size(cols,1)/4)-1),:)=[];
nvals = evals(1:7) - min(evals) +1;
pvals = evals(8:15) - min(evals(8:15)) +1;

evals = round(w.emoval*100);cols = jet(1000);  %cols(ceil(size(cols,1)/2)-9:ceil(size(cols,1)/2)+9,:) = [];
negs = cols(1:65:ceil(size(cols,1)/2)-109,:); poss = cols(ceil(size(cols,1)/2)+110:55:size(cols,1),:);
[vs ps] = sort(evals(8:15));
[vs ns] = sort(evals(1:7));
allord = [ns ps];
%------------------
[x,y,z]=sphere(20);x=x/30;y=y/30; z=z/30;
l=sqrt(x.*x+y.*y+z.*z);
normals = reshape([x./l y./l z./l],[21 21 3]);
figure; 
for e = 1:size(fullmd{clust},1)
    if e > 7
        cols = jet(max(pvals)*3); cols(1:max(pvals)*2,:) = [];
    colorarray = repmat(reshape(cols(pvals(e-7),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    %colorarray = repmat(reshape(poss(find(ps==(e-7)),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    else
        cols = jet(max(nvals)*3); 
    colorarray = repmat(reshape(cols(nvals(e),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    %colorarray = repmat(reshape(negs(allord(e),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    end;
    %colorarray = repmat(reshape(cols(e,:) , 1,1,3), [size(z,1) size(z,2) 1]);
    ph=surf(fullmd{clust}(e,dims(1))+x,fullmd{clust}(e,dims(2))+y,fullmd{clust}(e,dims(3))+z,colorarray, 'EdgeColor','none', 'VertexNormals', normals,'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);hold on;
end;
set(gca,'zlim',[-.7 .6]); 
set(gca,'fontsize',16);zl = get(gca,'zlim');
for e = 1:size(fullmd{clust},1)
    [cx,cy,cz] = cylinder([1 1],200);cx=cx/100;cy=cy/100; cz=cz/100;
    cx(1,:) = cx(1,:) + fullmd{clust}(e,dims(1));
    cy(1,:) = cy(1,:) + fullmd{clust}(e,dims(2));
    cz(1,:) = cz(1,:) + fullmd{clust}(e,dims(3));
    cx(2,:) = cx(2,:) + fullmd{clust}(e,dims(1));
    cy(2,:) = cy(2,:) + fullmd{clust}(e,dims(2));
    cz(2,:) = cz(2,:) + zl(1);
    if e > 7
        cols = jet(max(pvals)*3); cols(1:max(pvals)*2,:) = [];
    colorarray = repmat(reshape(cols(pvals(e-7),:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    %colorarray = repmat(reshape(poss(find(ps==(e-7)),:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    else
        cols = jet(max(nvals)*3); 
    colorarray = repmat(reshape(cols(nvals(e),:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    %colorarray = repmat(reshape(negs(allord(e),:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    end;
    %colorarray = repmat(reshape(cols(e,:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    ph=surf(cx,cy,cz,colorarray,'EdgeColor','none','backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.15);hold on;
end;
view(-19,32);
lighting phong; material shiny;
camlight right;camlight left; 
light;light;light;lightangle(-19,-32); 

%xlabel(['Dim ',int2str(dims(1))]);ylabel(['Dim ',int2str(dims(2))]);zlabel(['Dim ',int2str(dims(3))]);
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
%set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
set(gca,'xlim',[-.7 .6]); set(gca,'ylim',[-.8 .7]); set(gca,'zlim',[-.6 .4]);
set(gca,'xtick',[-.6:.2:.8]); set(gca,'ytick',[-1:.5:1]); set(gca,'ztick',[-1:.2:.6]);
%set(gca,'xticklabel',[{[],-.4,[],0,[],.4,[],.8}]);set(gca,'zticklabel',[{-.8,[],-.4,[],0,[],.4,[]}]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);set(gca,'zticklabel',[]);

print /home/julie/Manuscripts/Gamma/3DsphereEmoSpace.tif -dtiff
print /home/julie/Manuscripts/Gamma/3DsphereEmoSpace.eps -depsc -adobe
print /home/julie/Manuscripts/Gamma/3DsphereEmoSpaceStereo1.tif -dtiff
print /home/julie/Manuscripts/Gamma/3DsphereEmoSpaceStereo2.tif -dtiff

figure; clust = 12; % plot in 3D (traditional style):----
for e = 1:size(fullmd{clust},1)
    ph=plot3(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),fullmd{clust}(e,dims(3)),'.');hold on;
    set(ph,'markersize',30); set(ph,'color',cols(e,:));
    ph = text(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),fullmd{clust}(e,dims(3)),emo2{e});
    set(ph,'color',cols(e,:)); set(ph,'fontsize',40); 
end;
zl = get(gca,'zlim');
for e = 1:size(fullmd{clust},1)
    pl =plot3([fullmd{clust}(e,dims(1)) fullmd{clust}(e,dims(1))],[fullmd{clust}(e,dims(2)) fullmd{clust}(e,dims(2))],[zl(1)  fullmd{clust}(e,dims(3))]);
    set(pl,'color',cols(e,:)); set(pl,'linewidth',4)             
end;
xlabel(['Dim ',int2str(dims(1))]);ylabel(['Dim ',int2str(dims(2))]);zlabel(['Dim ',int2str(dims(3))]);
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
set(gca,'xlim',[-.7 .6]); set(gca,'ylim',[-.8 .7]); set(gca,'zlim',[-.6 .4]);

print /home/julie/Manuscripts/Gamma/3DsphereEmoSpaceNames.eps -depsc -adobe -painters

%%%%****************
clust = 12;% plot *2D* emo space as spheres/:-------
dims = [1,2]; % 2
emo3 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sadness','  compassion','  love','  relief','  contentment','  awe','  happiness','  joy','  excitement'}; % for plotting purposes
% adjust color scale:----
[xx yy] = sort(w.emoval);
emosval = emos{yy};
evals = round(w.emoval*100);evals = evals-min(evals);evals = evals+1;cols = jet(max(evals));cols(round(size(cols,1)/2) - (round(size(cols,1)/4)-1):round(size(cols,1)/2) + (round(size(cols,1)/4)-1),:)=[];

nvals = evals(1:7) - min(evals) +1;
pvals = evals(8:15) - min(evals(8:15)) +1;

cols = jet(1000);  %cols(ceil(size(cols,1)/2)-9:ceil(size(cols,1)/2)+9,:) = [];
negs = cols(1:65:ceil(size(cols,1)/2)-109,:); poss = cols(ceil(size(cols,1)/2)+110:55:size(cols,1),:);
[vs ps] = sort(evals(8:15));
[vs ns] = sort(evals(1:7),'descend');
allord = [ns ps];
%------------------
[x,y,z]=sphere(20);x=x/20;y=y/20; z=z/5;
l=sqrt(x.*x+y.*y+z.*z);
normals = reshape([x./l y./l z./l],[21 21 3]);
figure; colormap(jet(15));
for e = 1:size(fullmd{clust},1)
    if e > 7
        cols = jet(max(pvals)*3); cols(1:max(pvals)*2,:) = [];
    colorarray = repmat(reshape(cols(pvals(e-7),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    %colorarray = repmat(reshape(poss(find(ps==(e-7)),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    else
        cols = jet(max(nvals)*3); 
    colorarray = repmat(reshape(cols(nvals(e),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    %colorarray = repmat(reshape(negs(allord(e),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    end;
    %colorarray = repmat(reshape(cols(round(evals(e)/2),:) , 1,1,3), [size(z,1) size(z,2) 1]);
    ph=surf(fullmd{clust}(e,dims(1))+x,fullmd{clust}(e,dims(2))+y,z,colorarray, 'EdgeColor','none', 'VertexNormals', normals,'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);hold on;
    ph = text(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),emo3{e});
    %ph = text(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),emo3{e});
    if e > 7
        cols = jet(max(pvals)*3); cols(1:max(pvals)*2,:) = [];
    set(ph,'color',cols(pvals(e-7),:)); 
    else
        cols = jet(max(nvals)*3); 
    set(ph,'color',cols(nvals(e),:)); 
    end;set(ph,'fontsize',40); 
    %set(ph,'color',cols(e,:)); set(ph,'fontsize',40); 
end;
set(gca,'fontsize',16);
lighting phong; material shiny;
camlight right;camlight left; 
%xlabel(['Dim ',int2str(dims(1))]);ylabel(['Dim ',int2str(dims(2))]);
set(gca,'xgrid','on');  set(gca,'ygrid','on');
set(gca,'xlim',[-.65 .9]); set(gca,'ylim',[-.7 .7]); 
view(0,90)
ph = plot(b(1)*[-.5:.5]',b(2)*[-.5:.5]','k-');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); 
print /home/julie/Manuscripts/Gamma/2DsphereEmoSpace.tif -dtiff
print /home/julie/Manuscripts/Gamma/2DsphereEmoSpace2.tif -dtiff

figure; clust = 12; % plot in 3D (traditional style):----
for e = 1:size(fullmd{clust},1)
    ph=plot(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),'.');hold on;
    set(ph,'markersize',30); 
    if e > 7
    set(ph,'color',cols(pvals(e-7),:)); 
    else
    set(ph,'color',cols(nvals(e),:));
    end;
    ph = text(fullmd{clust}(e,dims(1)),fullmd{clust}(e,dims(2)),emo2{e});
    if e > 7
        cols = jet(max(pvals)*3); cols(1:max(pvals)*2,:) = [];
    set(ph,'color',cols(pvals(e-7),:)); 
    else
        cols = jet(max(nvals)*3); 
    set(ph,'color',cols(nvals(e),:)); 
    end;set(ph,'fontsize',40); 
     set(ph,'fontsize',40); 
end;
print /home/julie/Manuscripts/Gamma/2DsphereEmoSpaceNames.eps -depsc -adobe -painters

% Calculate linear discriminant based on first 2 dimensions:
 
[d p stats] = manova1(fullmd{clust}(:,dims),{'Negative','Negative','Negative','Negative','Negative','Negative','Negative','Negative','Positive','Positive','Positive','Positive','Positive','Positive','Positive'},.01);

[x,y,z]=sphere(20);x=x/4;y=y/4; z=z/4;
l=sqrt(x.*x+y.*y+z.*z);
normals = reshape([x./l y./l z./l],[21 21 3]);
figure; colormap(jet(15));
for e = 1:size(fullmd{clust},1)
    colorarray = repmat(reshape(cols(e,:) , 1,1,3), [size(z,1) size(z,2) 1]);
    ph=surf(stats.canon(e,1)+x,stats.canon(e,2)+y,stats.canon(e,1)+z,colorarray, 'EdgeColor','none', 'VertexNormals', normals, 'facecolor', 'interp');hold on;
    %ph = text(stats.canon(e,1),stats.canon(e,2),emo3{e});
    %set(ph,'color',cols(e,:)); set(ph,'fontsize',40); 
end;
set(gca,'fontsize',16);view(0,90)
lighting phong; material shiny;
camlight right;camlight left; 
 set(gca,'ylim',[-5 5]);
xlabel(['MDS ',int2str(dims(1))]);ylabel(['MDS ',int2str(dims(2))]);
set(gca,'xgrid','on');  set(gca,'ygrid','on');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/LdiscrimLine.tif -dtiff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------
% plot weighted templates:--------------------------------------
dim = 1; clust = 7; cut = 1;
cutval = cut*std(fullwts{clust}(dim,:)); % set cut value here****

figure;  wtedlist = cell(1,length(gdcomps));alltempls = [];
for x = 1:size(fullwts{clust},2)
    nx = keeptrack{clust}(x,1);
    im = keeptrack{clust}(x,2);
    if fullwts{clust}(dim,x) > cutval
        wtedlist{nx} = [wtedlist{nx} im];
        alltempls = [alltempls;emomeans{nx}(im,:)];
        for e = 1:size(fullmd{clust},1)
            ph=plot(e,emomeans{nx}(im,e),'.');hold on;
            set(ph,'markersize',12);set(ph,'color',cols(e,:));
        end;
    elseif fullwts{clust}(dim,x) < -cutval
        wtedlist{nx} = [wtedlist{nx} im];
        alltempls = [alltempls;emomeans{nx}(im,:)*-1];
        for e = 1:size(fullmd{clust},1)
            ph=plot(e,emomeans{nx}(im,e)*-1,'.');hold on;
            set(ph,'markersize',12);set(ph,'color',cols(e,:));
        end;        
    end;
end;
plot([get(gca,'xlim')],[0 0],'k-');
plot(mean(alltempls,1),'k-');
for e = 1:size(fullmd{clust},1)
    ph=plot(e,mean(alltempls(:,e),1),'.');hold on;
    set(ph,'markersize',30);set(ph,'color',cols(e,:));
end;        
title([ttl{clust},' weight cut off: ',num2str(cut),'*std(allweights)']);

str = ['print /home/julie/Manuscripts/Emotion/figures/',ttl{clust},'Dim',int2str(dim),'WtdTempls.eps -depsc -adobe -painters']; eval(str)

%---------------------------------------------
% find highest weight for each subject
dim =1; % same for all cls at the moment
clear clswts
for cls = 1:length(fullwts)
    clear subjwts 
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        subjinds = find(ktsubj{cls}(:,1) == nx);
        if ~isempty(subjinds)
            onesubj = fullwts{cls}(dim,subjinds);
            [x y] = max(abs(onesubj));
            subjwts(1,nx) = onesubj(y);
        end;
    end;
    clswts(cls,:) = subjwts;
end;
%---------------------------------------------
% plot dipoles color-coded by weights for a dimension
%---------------------------------------------
clustnames = ttl;
realclusts = [1,2,2,2,3,3,3,4,4,4,4,5]; subclusts = [1,1,2,3,1,2,3,1,2,3,4,1]; 
clustinfo = 'load /data/common1/emotion/AllClustFacs.mat';    
negs = [1:4,6:7]; poss = [9:10,12:15];
%clustnames = {'Brain','Muscle','OMT'}; % muscle
%realclusts = [1,1,1]; subclusts = [1,2,3]; % muscle
%clustinfo = 'load /data/common1/emotion/AllClustFacsMuscle.mat';  % muscle   

w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');useemos = [1:7,9:15];
dens = 'on'; % 'neg' or 'pos' weightings, [] to plot color-coded dipoles
cut = []; % * std
%---------------------------------------------
% WEIGHT BY MEDIAN IM WEIGHTS (pos vs neg emos grouped):
for decomp = length(fullwts):-1:1
    specclust = strs{realclusts(decomp)};
    subclust = subclusts(decomp);
    for ems = 1:2
        if ems == 1
            plotmeans = mean(collmeans{decomp}(poss,:),1);% pos emos
        else
            plotmeans = mean(collmeans{decomp}(negs,:),1);% neg emos
        end;
        EmoDimDipoles(subjlist,fullpaths,plotmeans,[],keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens); % weight by median wts
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
        if ems == 1
            ph=textsc([clustnames{decomp},' Dipole density wted by median IM wts -- avg of POSITIVE emos; Masked by permuting weights only'],'title'); set(ph,'color','r');
            str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'MedWtd-PosValence.tif -dtiff']; eval(str);
        else
            ph=textsc([clustnames{decomp},' Dipole density wted by median IM wts --  avg of NEGATIVE emos; Masked by permuting weights only'],'title'); set(ph,'color','r');
            str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'MedWtd-NegValence.tif -dtiff']; eval(str);
        end;
    end;
end;
%---------------------------------------------
% WEIGHT BY MEDIAN IM WEIGHTS:
for decomp = length(fullwts):-1:1
    specclust = strs{realclusts(decomp)};
    subclust = subclusts(decomp);
    for em = 1:15
        plotmeans = mean(collmeans{decomp}(em,:),1);
        EmoDimDipoles(subjlist,fullpaths,plotmeans,[],keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens); % weight by median wts
        ph=textsc([clustnames{decomp},' Dipole density wted by median IM wts -- ',emos{dim},'; Masked by permuting weights only'],'title'); set(ph,'color','r');
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
        str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'MedWtd-',emos{dim},'.tif -dtiff']; eval(str);
    end;
end;
%---------------------------------------------
% WEIGHT BY INDIVIDUAL CORRELATIONS WITH VALENCE:
clear mnspecs templcell templs kptk emodeciles emomeans clustfacs allbigs normals modcorr finaltempls finalmeans 
w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');useemos = [1:7,9:15];
for decomp = length(fullwts):-1:1
    specclust = strs{realclusts(decomp)};
    subclust = subclusts(decomp);
    for dim = 1:size(fullwts{decomp},1)% adds dim template when dipole plotting
        clear corr
        %for m = 1:size(collmeans{decomp},2)
        for m = 1:size(collmeans,2)
            %corr(1,m) = corr2(w.emoval(useemos),collmeans{decomp}(useemos,m)');
            corr(1,m) = corr2(w.emoval(useemos),collmeans(useemos,m)');
        end;
        %EmoDimDipoles(subjlist,fullpaths,corr,fullmd{decomp}(:,dim),keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens); % weight by correlation       

        EmoDimDipoles(subjlist,fullpaths,corr,fullmd(:,dim),keeptrack,pics,cut,specclust,subclust,clustinfo,dens); % weight by correlation       
        ph=textsc([clustnames{decomp},' Dipole density weighted by IM Corr w/Valence; Masked by permuting weights only'],'title'); set(ph,'color','r');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
        str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'ValenceCorrWtdDim',int2str(dim),'2.tif -dtiff']; eval(str);
        %str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'ValenceCorrWtdDim',int2str(dim),'Dipoles.tif -dtiff']; eval(str);
   end;
end;
%---------------------------------------------
% WEIGHT BY ICA EMO SPACE WEIGHTS:
for decomp = length(fullwts):-1:1
    specclust = strs{realclusts(decomp)};
    subclust = subclusts(decomp);
    for dim = 1:size(fullwts{decomp},1)% dims or emotions, depending on choice below.
        EmoDimDipoles(subjlist,fullpaths,fullwts{decomp}(dim,:),fullmd{decomp}(:,dim),keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens); % weight by correlation       
        ph=textsc([clustnames{decomp},' Dipole density weighted by ICA emospace weights; Masked by permuting weights only'],'title'); set(ph,'color','r');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
        str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'ICAWtdDensDim',int2str(dim),'.tif -dtiff']; eval(str);
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------
dipargs = {'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat','spherical','view',[0 0 1]};


%---------------------------------------------
%---------------------------------------------
% plot relevant IMs for a cluster/dimension:
clust = 1; % refering to strs
c = 2; % referring to cell array from EmoSpace decomps (all sub-clusters make this a different index)
cut = .8; % * std
savedat = 'SpecCoMod';
savedat = 'SpecCoModMuscle';  % includes muscle and inf frontal

cutval = cut*std(fullwts{c}(dim,:)); % set cut value here****
ifacs = find(abs(fullwts{c}(dim,:)) > cutval);
selwts = fullwts{c}(dim,ifacs);
subjfacs = keeptrack{c}(ifacs,:);       
   
for nx = 1:max(subjfacs(:,1))
    if ~isempty(find(subjfacs(:,1) == nx))
        ims = abs(subjfacs(find(subjfacs(:,1) == nx),2));
        for im = 1:length(ims)            
            SpecCoModPlot('sources.set',newpaths{nx},[],[],savedat,[3 128],'n',0,ims(im));
        end;
    end;
end;
    
%---------------------------------------------
%---------------------------------------------
clear subjfacs iclistp iclistn ctemplsp ctemplsn ilnxp ilnxn ptemplsp  ptemplsn
for c=1:13
    cutval = max(abs(fullwts{c}(dim,:))) -  max(abs(fullwts{c}(dim,:)))*.50;
    %cutval = mean(fullwts{c}(dim,:)) +  std(fullwts{c}(dim,:));
    ifacs = find(abs(fullwts{c}(dim,:)) > cutval);
    subjfacs{c} = ktsubj{c}(ifacs,:);
    multfac = fullwts{c}(dim,find(abs(fullwts{c}(dim,:)) > cutval))./abs(fullwts{c}(dim,find(abs(fullwts{c}(dim,:)) > cutval)));    
    subjfacs{c}(:,2) = subjfacs{c}(:,2).*multfac';% invert if neg weighting
    for nx = 1:35
        subidxs = subjfacs{c}(find(subjfacs{c}(:,1)==nx),:);
        ilp = []; iln = []; ptp = [];ptn = [];
        for im = 1:size(subidxs,1)
            assocics = pics{c}{nx}(find(abs(pics{c}{nx}(:,1))==abs(subidxs(im,2))),:);
            assoctempls = templs{c}{nx}(find(abs(pics{c}{nx}(:,1))==abs(subidxs(im,2))),:);
            for ic = 1:size(assocics,1)
                if  subidxs(im,2) < 0 % downward
                    iln = [iln assocics(ic,2)];
                    ptn = [ptn; assoctempls(ic,:)*-1];% reorient template
                elseif subidxs(im,2) > 0 % upward
                    ilp = [ilp assocics(ic,2)];
                    ptp = [ptp; assoctempls(ic,:)];
                end;
            end;
        end;
        ilnxp{nx} =  ilp;
        ilnxn{nx} =  iln;
        ptemplsp{nx} = ptp;
        ptemplsn{nx} = ptn;
    end;          
    iclistp{c} =  ilnxp;
    iclistn{c} =  ilnxn;
    ctemplsp{c} = ptemplsp;
    ctemplsn{c} = ptemplsn;
end;
% plot one clust: c is pre-defined
dipargs = {'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat','spherical'};

figure; row = 2; col = 2; viewnum = [1,2,3]; place = 1; 
sbplot(row,col,place); coltempls = [];
for nx = 1:length(ctemplsp{c})
    coltempls = [coltempls; ctemplsp{c}{nx}];
end;
ph = quadplot(freqs,coltempls,1.5,'r');place = place+1;
sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, iclistp{c},[],[dipargs,'view',[0 0 1]],[],{'r'},'off',[],[]);place = place +1;
sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, iclistp{c},[],[dipargs,'view',[1 0 0]],[],{'r'},'off',[],[]);place = place +1;
%PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistp,c,row,col,place,ttl{c},viewnum,[]);
set(gcf,'color','w');

figure; row = 2; col = 2; viewnum = [1,2]; place = 1; 
sbplot(row,col,place); coltempls = [];
for nx = 1:length(ctemplsp{c})
    coltempls = [coltempls; ctemplsn{c}{nx}];
end;
ph = quadplot(freqs,coltempls,1.5,'b');place = place+1;
PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistn,c,row,col,place,ttl{c},viewnum,[0 0 1]);
set(gcf,'color','w');

% plot the corresponding spectral templates
figure; row = 5; col = 6; viewnum = [1,2,3]; place = 1;pg=1;
for clust = 1:length(iclistp)
    if place > row*col
        ph = textsc(['Power increases-%max'],'title'); set(ph,'fontsize',14); set(ph,'color','r');
set(gcf,'color','w');set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        %str = ['print /data/common2/emotion/ClustFigs/Dim1PwrInc',int2str(pg),'.jpg -djpeg'];eval(str); %close
        figure; place = 1; pg = pg+1;
    end; 
    sbplot(row,col,place); coltempls = [];
    for nx = 1:length(ctemplsp{clust})
        coltempls = [coltempls; ctemplsp{clust}{nx}];
    end;
    ph = quadplot(freqs,coltempls,1.5,'r');place = place+1;
    PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistp,clust,row,col,place,ttl{clust},viewnum,[]);
    place = place+length(viewnum);
end;
ph = textsc(['Power increases-%max'],'title'); set(ph,'fontsize',14); set(ph,'color','r');set(gcf,'color','w');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%str = ['print /data/common2/emotion/ClustFigs/Dim1PwrInc',int2str(pg),'.jpg -djpeg'];eval(str); %close

figure; row = 5; col = 6; viewnum = [1,2]; place = 1;pg=1;% downward IMs/ICs
for clust = 1:length(iclistn)
    if place > row*col
        ph = textsc(['Power decreases-%max'],'title'); set(ph,'fontsize',14); set(ph,'color','r');set(gcf,'color','w');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        str = ['print /data/common2/emotion/ClustFigs/Dim1PwrDec',int2str(pg),'.jpg -djpeg'];eval(str); %close
        figure; place = 1; pg = pg+1;
    end;
    sbplot(row,col,place); coltempls = [];
    for nx = 1:length(ctemplsn{clust})
        coltempls = [coltempls; ctemplsn{clust}{nx}];
    end;
    ph = quadplot(freqs,coltempls,1.5,'b');place = place+1;
    PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistn,clust,row,col,place,ttl{clust},viewnum,[0 0 1]);
    place = place+length(viewnum);
end;
ph = textsc(['Power decreases-%max'],'title'); set(ph,'fontsize',14); set(ph,'color','r');set(gcf,'color','w');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /data/common2/emotion/ClustFigs/Dim1PwrDec',int2str(pg),'.jpg -djpeg'];eval(str); %close
%%%%% --- Are there different sets of subjects?-----
clear subjmat
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    for c = 1:length(iclistp)
        if ~isempty(iclistp{c}{nx})
            subjmat(nxx,c) = 1;
        else
            subjmat(nxx,c) = -1;
        end;
    end;
    for c = 1:length(iclistn)
        if ~isempty(iclistn{c}{nx})
            subjmat(nxx,c+length(iclistp)) = 1;
        else
            subjmat(nxx,c+length(iclistp)) = -1;
        end;
    end;
end;
[wts,sph,compvars,bias,signs,lrates,acts] = runica(subjmat','pca',rank(subjmat')-1,'extended',1,'stop',1e-7,'maxsteps',6000);
%%% can't see any obvious clusters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Individual Subject Emotion Spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%numdims = 2; % 3 accounts for between 77 and 94% of total variance
%5 accounts for between 91.0487% 97.9815% of total variance
%10 accounts for between 98.4593% and 99.8426% of total variance
% minimal pca on single-subj ICA, but only take ~3 dims from each subj. 
savedat = 'SpecCoModMoreFreqs'; fullpaths = newpaths;

clear data activations winv alldip kptk allbigs bigwts clustfacs mnspecs orivec
clustfacs = [];  mnspecs = []; kptk = [];
subjlist = [2:21,23:31,33:35];
freqscale = 'quad';
percmax =.5; % percent of max template to take as a 'comod'
 [clustfacs,templcell,mnspecs,kptk,allbigs,bigwts,onebig,orivec,modcorr,freqs,outcomods,rawrms,rawcorr] = CollectCoModTempls(savedat,newpaths,subjlist,percmax);
nsamples = []; % only look at first x samples
[emomeans,emodeciles,subjpoints] = EmoWeights(savedat,fullpaths,gdcomps,[1:35],nsamples);

% BUTTON press Groups:----------------------------
button = [2:12];  ttl = 'Repetitive button presses';% really subj 1 too
button = [21,23:26];  ttl = 'only when feeling it';
button = [13:20,27:31,33:35]; ttl = ' no button press';% (apart from the first one)
subjlist = button;
%---------------------------
% find lists of clustered IMs for Emo Space decomp
load /data/common1/emotion/GammaClust7-10-2009.mat finaltempls finalidx finalmeans freqs
allidx = []; whichcls = []; plotidx = cell(1,5);
for cls = 1:length(finalidx)
    allidx = [allidx;finalidx{cls}];
end;
% find indices of all clustered IMs-----------------------

%subjlist = [2:21,23:31,33:35];% all but 1,22,3
onesubjics = cell(1,max(subjlist)); 
oricell = cell(1,max(subjlist)); 
useims = cell(1,max(subjlist)); 
onlyims = [];
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    onesubjics{nx} = allidx(find(allidx(:,1) == nx),:);% all clustered templates
    onesubjims = allidx(find(allidx(:,1) == nx),2);
    useims{nx} = unique(abs(onesubjims)');% list of clustered IMs for each subj
    onlyims = [onlyims;[repmat(nx,[length(useims{nx}) 1]),useims{nx}']];
    for im = 1:length(useims{nx})
        idx = find(allbigs{nx}{useims{nx}(im)} == onebig{nx}{useims{nx}(im)});
        oricell{nx}(im) = orivec{nx}{useims{nx}(im)}(idx);
    end;
end;

%figure;    
[fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,subjlist,useims,'mds',0,'all subjs',0,'keep','noperm');      

w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');useemos = [1:7,9:15];
useemos = [1:7,9:15]; clear corrs% no compassion
for m = 1:size(collmeans,2)
    corrs(1,m) = corr2(w.emoval(useemos),collmeans(useemos,m)');
end;
sigcorr = find(abs(corrs)>.5); onlyims = onlyims(sigcorr,:); 
subjlist = unique(onlyims(:,1))'; % new subj list with only corr ims
onesubjics = cell(1,max(subjlist)); 
oricell = cell(1,max(subjlist)); 
useims = cell(1,max(subjlist)); 
allidx = onlyims; clear pics
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    onesubjics{nx} = allidx(find(allidx(:,1) == nx),:);% all clustered templates
    onesubjims = allidx(find(allidx(:,1) == nx),2);
    useims{nx} = unique(abs(onesubjims)');% list of clustered IMs for each subj
    for im = 1:length(useims{nx})
        idx = find(allbigs{nx}{useims{nx}(im)} == onebig{nx}{useims{nx}(im)});
        oricell{nx}(im) = orivec{nx}{useims{nx}(im)}(idx);
        pics{nx} = [repmat(nx,[length(allbigs{nx}{useims{nx}(im)}) 1]) repmat(useims{nx}(im),[length(allbigs{nx}{useims{nx}(im)}) 1]) allbigs{nx}{useims{nx}(im)}'];
    end;
end;
fidx{3} = pics;
ffs = cell(1,3);
for cls=1:3
    for nx = 1:length(fidx{cls})
        ffs{cls} = [ffs{cls}; fidx{cls}{nx}];
    end;
end;

figure;    [fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,subjlist,useims,'mds',0,['only val corr''d, ',int2str(length(subjlist)),' subjs'],0,'keep','noperm');      
set(gca,'fontsize',16);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/OnlyValCorrEmoSpace.eps -depsc
print /home/julie/Manuscripts/Gamma/Frontiers/OnlyValCorrWtdDips.jpg -djpeg


subjlist = button;
figure; pl = 1;
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    sbplot(3,2,pl);pl = pl+1;
    [fullmd,stresses(1,nx),keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,nx,useims,'mds',0,['Subj ',int2str(nx)],0,'keep','noperm');  hold on;     
end;
print /home/julie/Manuscripts/Gamma/Frontiers/OnlyValCorrIndivSubjs.eps -depsc
stresses(find(stresses==0)) = NaN;
bestest=find(stresses<.085)
worstest=find(stresses>=.085)
figure;[fullmd,stress,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,bestest,useims,'mds',0,['bestest'],0,'keep','noperm');  hold on;stress     
figure;[fullmd,stress,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,worstest,useims,'mds',0,['worstest'],0,'keep','noperm');  hold on;  stress   

%---------------------
subjlist = [2:21,23:31,33:35];% all but 1,22,3
subjlist = [2,5,6,9,11,14,16,17,18,19,20,21,23,25,26, 27,29,30,33,34,35]; % subjs with val-correlated IMs
buttgrps = {[13:20,27:31,33:35],[2:12,21,23:26]};% no butt., butt.

figure;[fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,subjlist,useims,'mds',0,'full mds',0,'keep','noperm');      
figure;[fullmdN,ostress1,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,buttgrps{1},useims,'mds',0,'',0,'keep','noperm');  % clf   
figure;[fullmdB,ostress2,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,buttgrps{2},useims,'mds',0,'',0,'keep','noperm'); 
cN = corr(fullmd,fullmdN);
cB = corr(fullmd,fullmdB);
if abs(cN(1,1)) < abs(cN(1,2))
    fullmdN(:,1:2) = fullmdN(:,2:-1:1);
    cN = corr(fullmd,fullmdN);
end;
if abs(cB(1,1)) < abs(cB(1,2))
    fullmdB(:,1:2) = fullmdB(:,2:-1:1);
    cB = corr(fullmd,fullmdB);
end;
if cN(1,1) < 0        
    fullmdN(:,1) = fullmdN(:,1)*-1;
end;
if cN(2,2) < 0        
    fullmdN(:,2) = fullmdN(:,2)*-1;
end;
if cB(1,1) < 0        
    fullmdB(:,1) = fullmdB(:,1)*-1;
end;
if cB(2,2) < 0        
    fullmdB(:,2) = fullmdB(:,2)*-1;
end;
cN = corr(fullmd,fullmdN); corrnobutt = [cN(1,1),cN(2,2)];
cB = corr(fullmd,fullmdB);corrbutt = [cB(1,1),cB(2,2)];
fullmdN = fullmdN.*repmat(negvec,[1 size(fullmdN,2)]); % reverse neg emos
fullmdB = fullmdB.*repmat(negvec,[1 size(fullmdB,2)]); % reverse neg emos
bd = fullmdN-fullmdB;
buttondiff = sum(sqrt(bd(:,1).^2 + bd(:,2).^2));

        
for e=1:15
fullthetas(1,e) = atan2(fullmd(e,2),fullmd(e,1));
end;
clear  distsone
figure;row=4; col=4;pl=1;

clear bootkeepmdists corr
for nsubj = 2:17
    %sbplot(row,col,pl);pl = pl+1;
    %cols = jet(15);randgrp1 = [];randgrp2 = []; 
    clear bootdists bootdiffs
    negvec = ones(15,1); negvec(1:8) = negvec(1:8)*-1;
    mdists = [];
    for b = 1:25%length(subjlist)
        %negvec = shuffle(negvec);
        %sbplot(row,col,b);
        randsubjs = randperm(length(subjlist));
        %randgrps{1} = subjlist(b);
        randgrps{1} = subjlist(randsubjs(1:nsubj));
        %randgrps{1} = subjlist(randsubjs(1:length(buttgrps{1})));
        %randgrps{2} = subjlist(randsubjs(length(buttgrps{1})+1:end));
        [fullmd1,stress1,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,randgrps{1},useims,'mds',0,'',0,'keep','noperm');      
        %[fullmd2,stress2,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,randgrps{2},useims,'mds',0,'',0,'keep','noperm');  
        if ~isempty(find(fullmd1)) & size(fullmd1,2) > 1
            c1 = corr(fullmd(:,1:2),fullmd1(:,1:2));
            %c2 = corr(fullmd(:,1:2),fullmd2(:,1:2));
            if abs(c1(1,1)) < abs(c1(1,2))
                fullmd1(:,1:2) = fullmd1(:,2:-1:1);
                c1 = corr(fullmd,fullmd1);
            end;
            if c1(1,1) < 0        
                fullmd1(:,1) = fullmd1(:,1)*-1;
            end;
            if c1(2,2) < 0        
                fullmd1(:,2) = fullmd1(:,2)*-1;
            end;
            %if abs(c2(1,1)) < abs(c2(1,2))
            %    fullmd2(:,1:2) = fullmd2(:,2:-1:1);
            %    c2 = corr(fullmd,fullmd2);
            %end;
            %if c2(1,1) < 0        
            %    fullmd2(:,1) = fullmd2(:,1)*-1;
            %end;
            %if c2(2,2) < 0        
            %    fullmd2(:,2) = fullmd2(:,2)*-1;
            %end;
            %c1 = corr(fullmd,fullmd1); randgrp1(b,:) = [c1(1,1),c1(2,2)];
            %c2 = corr(fullmd,fullmd2);randgrp2(b,:) = [c2(1,1),c2(2,2)];
            for rw = 1:size(fullmd1)
                indsums(b,rw) = sqrt(fullmd1(rw,1).^2 + fullmd1(rw,2).^2);
            end;            
            fullmd1 = fullmd1.*repmat(negvec,[1 size(fullmd1,2)]); % reverse neg emos
            fullmd1= sum(fullmd1,1);
            mdists = [mdists (sum(sqrt(fullmd1(:,1).^2 + fullmd1(:,2).^2)))/sum(indsums(b,:))];
            %fullmd2 = fullmd2.*repmat(negvec,[1 size(fullmd2,2)]); % reverse neg emos
            %bd = fullmd1-fullmd2;
            %bd = bd.*repmat(negvec,[1 size(bd,2)]); % reverse neg emos
            
            %bootdiffs(1,b) = sum(sqrt(bd(:,1).^2 + bd(:,2).^2));
            %for e=1:15
            %onerx = fullmd(e,1) - fullmd1(e,1);
            %onery = fullmd(e,2) - fullmd1(e,2);
            %bootdists(b,e) = sqrt(onerx^2 + onery^2);
            %distsone(b,e) = sqrt(onerx^2 + onery^2);
            %end;
            %for e = 1:size(fullmd,1)
            %ph=plot(fullmd1(e,1),fullmd1(e,2),'.');hold on;
            %set(ph,'markersize',5);set(ph,'color',cols(e,:)+(1-cols(e,:))*.75);
            %set(ph,'markersize',5);set(ph,'color',cols(e,:));
            %end;
            %for e = 1:size(fullmd,1)
            %    ph=plot(fullmd2(e,1),fullmd2(e,2),'.');hold on;
            %    set(ph,'markersize',5);set(ph,'color',cols(e,:)+(1-cols(e,:))*.75);
            %end;
        end;    
    end;
    bootkeepmdists{nsubj} = mdists;
    %set(gca,'xlim',[-1.2 1.2]);set(gca,'ylim',[-1.2 1.2]);
end;
bootdiffs = sort(bootdiffs);
mincut =  bootdiffs(round(.01*length(bootdiffs)));
maxcut = bootdiffs(end-round(.01*length(bootdiffs))+1);

figure; set(gca,'fontsize',16);
hist(bootdiffs,100); hold on;
plot([mincut mincut],[get(gca,'ylim')],'r-','linewidth',2);
plot([maxcut maxcut],[get(gca,'ylim')],'r-','linewidth',2);
plot([buttondiff buttondiff],[get(gca,'ylim')],'g-','linewidth',2);
xlabel('Weighted Vector Sum'); ylabel('Number of iterations');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/Butt-NoButtHist.jpg -djpeg


figure; dcols = cool(16);set(gca,'fontsize',16);
for nsubj = 2:17
    x =  sort(keepmdists{nsubj});
    plot(x,'-','linewidth',2,'color',dcols(nsubj-1,:)); hold on;
    x =  sort(bootkeepmdists{nsubj});
    plot([1:100],x([1:round(length(x)/100):length(x)]),':','linewidth',2,'color',dcols(nsubj-1,:)); hold on;
end;
set(gca,'xlim',[1 100]);
legend({'2 subjs','3 subjs','4 subjs','5 subjs','6 subjs','7 subjs','8 subjs','9 subjs','10 subjs','11 subjs','12 subjs','13 subjs','14 subjs','15 subjs','16 subjs','17 subjs'},'location','EastOutside');
xlabel('Sorted Permutation Iterations'); ylabel('Normalized Vector Sum (range, [0 1])');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeIncrSubjsMDScircCorrsNew.jpg -djpeg

%bootkeepmdists = keepmdists;

set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeIncrSubjsMDS.jpg -djpeg
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeIncrSubjsMDS.eps -depsc
distsone(find(distsone(:,1)==0),:) = [];
bootdists(find(bootdists(:,1)==0),:) = [];
figure; imagesc(bootdists);cbar; figure; plot(mean(bootdists,2));
figure; imagesc(distsone);cbar; figure; plot(mean(distsone,2));
figure; hist(distsone);
pval = .05; clear mask
mask(1,:) = bootdists(round(pval*size(bootdists,1)),:);
mask(2,:) = bootdists(end-(round(pval*size(bootdists,1))-1),:);
minmask = repmat(mask(1,:),[size(distsone,1) 1]);
maxmask = repmat(mask(2,:),[size(distsone,1) 1]);

sigmat = zeros(size(distsone));
sigmat(find(distsone < maxmask)) = 1;
figure; imagesc(sigmat);


set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeHalvesMDS.jpg -djpeg
figure; hist(mean(dists,2),50);

pval = .01;
x=sort(dists,1);
for e=1:15
    mask(1,e) = x(round(size(dists,1)*pval),e);
    mask(2,e) = x(end - (round(size(dists,1)*pval)-1),e);
end;
subjscore = zeros(size(distsone,1),size(distsone,2));
for nx = 1:length(subjlist)
    if distsone(nx,1) == 0 
        subjscore(nx,:) = NaN;
    else
       for e=1:15
           if distsone(nx,e) > mask(1,e)& distsone(nx,e) < mask(2,e)
               subjscore(nx,e) = 1;
           end;
       end;

    end;
end;
figure; plot(sort(sum(subjscore,2)'/15))

% Find vector sums:
x2 = fullmd(1,1);x1=0;
y2 = fullmd(1,2);y1=0;
thet = atan2(y2,x2);
distance = sqrt((x2 - x1)^2 + (y2-y1)^2);
        degs = (angles{2}*180)/pi;

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cols = jet(15);cols(10,:) = [1 .9 0];% for plotting purposes also
for e = 1:size(fullmd,1)
    %ph=plot(fullmdN(e,1),fullmdN(e,2),'.','markersize',1);set(ph,'color',cols(e,:));
    ph=text(fullmdN(e,1),fullmdN(e,2),'N');hold on;
    set(ph,'fontsize',22);set(ph,'fontweight','bold');
    set(ph,'color',cols(e,:));
end;
for e = 1:size(fullmd,1)
    %ph=plot(fullmdB(e,1),fullmdB(e,2),'.','markersize',1);set(ph,'color',cols(e,:));
    ph=text(fullmdB(e,1),fullmdB(e,2),'B');hold on;
    set(ph,'fontsize',22);set(ph,'fontweight','bold');
    set(ph,'color',cols(e,:));
end;
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeHalvesMDSsuperpose.jpg -djpeg


figure; hist([randgrp1;randgrp2],50); hold on;
legend({'Dim 1','Dim 2'},'location','northwest');
ylabel(['Number of iterations (out of ',int2str(b),')']);
xlabel('Correlation with full MDS');
plot([corrnobutt(1) corrnobutt(1)],[get(gca,'ylim')],'b-');
plot([corrnobutt(2) corrnobutt(2)],[get(gca,'ylim')],'r-');
plot([corrbutt(1) corrbutt(1)],[get(gca,'ylim')],'c-');
plot([corrbutt(2) corrbutt(2)],[get(gca,'ylim')],'m-');

figure; plot(sort([randgrp1(:,1);randgrp2(:,1)]),'b-','linewidth',4);hold on;
plot(sort([randgrp1(:,2);randgrp2(:,2)]),'r-','linewidth',4);
legend({'Dim 1','Dim 2'},'location','southeast');
plot([get(gca, 'xlim')],[corrnobutt(1) corrnobutt(1)],'g-','linewidth',3);
plot([get(gca,'xlim')],[corrnobutt(2) corrnobutt(2)],'g--','linewidth',3);
plot([get(gca,'xlim')],[corrbutt(1) corrbutt(1)],'c-','linewidth',3);
plot([get(gca,'xlim')],[corrbutt(2) corrbutt(2)],'c--','linewidth',3);
xlabel(['Iterations']);
ylabel('Correlation with full MDS');
print /home/julie/Manuscripts/Gamma/Frontiers/IterativeHalvesMDSsortvals.jpg -djpeg


plot([get(gca,'xlim')],[median(realcorr) median(realcorr)],'k-');
plot([get(gca,'xlim')],[origcorr origcorr],'r-');

median(realcorr) 0.6747
vs.
origcorr  0.6830

 
% How many subjects does it take to decrease stress to < .1 (.15 unacceptable)
figure; clear keepgrps
for n = 1:length(subjlist)
    clear nowgrp
    for b = 1:1000
        randsubjs = randperm(length(subjlist));
        nowgrp(b,:) = subjlist(randsubjs(1:n));
        [fullmd,stress1(n,b),keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,nowgrp(b,:),useims,'mds',0,'',0,'keep','noperm'); clf
    end;
    keepgrps{n} = nowgrp;
end;
save tmpstresstest.mat keepgrps stress1
stresstest = stress1;
stresstest(find(stress1==0)) = NaN;
clear nsubjs
for n=1:length(subjlist)
    [val idx] = sort(stresstest(n,:),'ascend');
    idx(find(isnan(val))) = [];
    val(find(isnan(val))) = [];
    bestsubjs =  keepgrps{n}(idx(1:10),:); % best stresses
    worstsubjs =  keepgrps{n}(idx(end-9:end),:);% worst stresses
    bestsubjs = reshape(bestsubjs,[1 size(bestsubjs,1)*size(bestsubjs,2)]);
    worstsubjs = reshape(worstsubjs,[1 size(worstsubjs,1)*size(worstsubjs,2)]);
    for s = 1:max(subjlist)
        nsubjs(1,s,n) = length(find(bestsubjs == s));
        nsubjs(2,s,n) = length(find(worstsubjs == s));
    end;
end;
x= squeeze(sum(nsubjs(1,:,:),3));
y=squeeze(sum(nsubjs(2,:,:),3));
z=x./y; [val idx] = sort(z,'descend');    idx(find(isnan(val))) = [];    val(find(isnan(val))) = [];

bestest = idx(1:5);
worstest = idx(end-4:end);
figure;[fullmd,stress,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,bestest,useims,'mds',0,['bestest'],0,'keep','noperm');  hold on;stress     
figure;[fullmd,stress,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,worstest,useims,'mds',0,['worstest'],0,'keep','noperm');  hold on;  stress   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KURTOSIS on intra-emotion IM weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /data/common1/emotion/GammaClust7-10-2009.mat finaltempls finalidx finalmeans freqs
cls = 1;
allskew = [];allkurt = []; allvar = [];allmedian= [];allmean= [];
figure; row = 10; col=10; pl=1;
for x = 1:size(finalidx{cls})
    nx = finalidx{cls}(x,1);
    im = abs(finalidx{cls}(x,2));
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;    winv = pinv(ws); 
    clear wts sph ws 
    speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  % templates   
    winv = specwts;  clear ekurt   eskew evar
    for e = 1:length(s.dstrials)
        evar(1,e) = var(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
        ekurt(1,e) = kurtosis(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
        eskew(1,e) = skewness(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
        emean(1,e) = mean(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
        emedian(1,e) = median(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
    end;
    if pl > row*col
        figure; pl=1;
    end;
    sbplot(row,col,pl); pl=pl+1;hist(winv(s.keeptrack(e,1):s.keeptrack(e,2),im));
    allskew = [allskew eskew];
    allkurt = [allkurt ekurt];
    allvar = [allvar evar];
    allmedian = [allmedian,emedian];
    allmean = [allmean,emean];
end;
figure; plot(allvar,allkurt,'b.');
figure; plot(allskew,allkurt,'b.');
figure; plot(allmean,allmedian,'b.');xlabel('mean'); ylabel('median');

sbplot(row,col,1)
hist(allskew); 
figure; hist(allkurt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSIFY emotions within-subject (see new script 'ClassifyEmos.m' instead)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cols = jet(15); cols= [cols;[.5 .5 .5]];

mtchpercents= zeros(0,15);   alldatmatchs= zeros(0,15); 
for nx = 1:35%length(fullpaths)
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;    winv = pinv(ws); 
    clear wts sph ws 
    speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  % templates   
    winv = specwts;  clear delpoints emeans  
    for e = 1:length(s.dstrials) % break up windows into emos...
        ndec = round(s.dstrials(e)*.05); % test data  = 5% of total for each emo
        epoints = sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e));
        rpoints = [round(length(epoints)/2-ndec/2):round(length(epoints)/2+ndec/2)];
        delpoints{e} = epoints(rpoints);% prediction windows
        epoints([rpoints(1)-1,rpoints,rpoints(end)+1]) = [];% take out buffer windows
        for dim = 1:size(winv,2) % for all (new)dims...
            emeans(dim,e) = mean(winv(epoints,dim));
        end;
    end; % makes a 15 dims X 15 emotions matrix
         % select only clustered IMs:--------------
% $$$         useims = [];subjims = cell(1,35);
% $$$         for ss = 1:length(strs)
% $$$             eval(strs{ss})    
% $$$             for cls = 1:length(finalidx)
% $$$                 onesubj = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
% $$$                 if ~isempty(onesubj)
% $$$                     useims = [useims,unique(abs(onesubj))'];
% $$$                 end;
% $$$             end;
% $$$         end;
% $$$         useims = unique(useims); 
    
    % OR ------------
    useims = [];grp = [];
    for e=1:length(s.dstrials)
        for t = s.keeptrack(e,1):s.keeptrack(e,2)
            grp{t} = emos{e};
        end;
    end; clear Fs
    for im = 1:size(winv,2)
        [P,table]=anovan(winv(:,im)',{grp},'display','off');
        Fs(1,im) = table{2,6};
    end;
    useims = find(Fs > mean(Fs));
    keepuseims{nx} = useims;
    nuseims(1,nx) = length(useims);
    %----------------
    svmdat = mean(emeans(useims,1:7),2)';
    svmdat = [svmdat;mean(emeans(useims,9:15),2)'];
    SVMStruct = svmtrain(svmdat,{'neg','pos'}); % no compassion
    if ~isempty(useims)
        clear imcombos colldata scores ntrials  percscores  mtchs svmgrp
        for emoidx = 1:length(s.dstrials)
            if s.dstrials(emoidx) > 130 % such that 1 min remains after test data removed
                clear mtch sgrp 
                for itr = 1:length(delpoints{emoidx})                    
                    nowtr = delpoints{emoidx}(itr); %use each point individually
                    onetrial = mean(winv(nowtr,useims),1);
                    nearemo = knnclassify(onetrial,emeans(useims,:)',emos);
                    if strcmp(emos{emoidx},nearemo)
                        mtch(itr) = 1; % how many correct with one time point?
                    else
                        mtch(itr) = 0;
                    end;
                    randtr = randperm(length(delpoints{emoidx}));
                    nowtr = delpoints{emoidx}(randtr(1:itr)); 
                    onetrial = mean(winv(nowtr,useims),1);
                    nearemo = knnclassify(onetrial,emeans(useims,:)',emos);
                    sgrp(1,itr) = svmclassify(SVMStruct,onetrial);
                end;
                mtchs(emoidx) = sum(mtch)/length(mtch); %
                svmgrp{emoidx} = sgrp{end};
                if strcmp(emos{emoidx},nearemo)
                    alldatmatch(emoidx) = 1; % 1/0 using all test data
                else
                    alldatmatch(emoidx) = 0;
                end;
            else
                mtchs(emoidx) = NaN;
                svmgrp{emoidx} = NaN;
                alldatmatch(emoidx) = NaN;
            end;
        end;
        mtchpercents(end+1,:) = mtchs*100; 
        alldatmatchs(end+1,:) = alldatmatch; 
        subjidxs(end+1) = nx;
    end;
end;

figure; 
boxplot(mtchpercents);hold on;
ph = plot(median(mtchpercents,1),'k.-','linewidth',3.5,'markersize',25); hold on;
ph = plot([get(gca,'xlim')],[50 50],'k--','linewidth',2.5); hold on;
xlabel('Emotions');ylabel('Percent matches');
title('Percent matches across subjects from 1-sec data tests');
set(gca,'xtick',[0:15]);
set(gca,'xticklabel',{'','anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'});
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

for e=1:size(alldatmatchs,2)
    robustemos(e) = 100*(length(find(alldatmatchs(:,e)))/size(alldatmatchs,1));
end;
figure; row = 4; col = 4;
ph = plot(robustemos,'.-','linewidth',1.5,'markersize',25); hold on;
set(gca,'xlim',[0 16]);set(gca,'ylim',[0 100]);
ph = plot([get(gca,'xlim')],[50 50],'k--','linewidth',2.5); hold on;
xlabel('Emotions');ylabel('Percent matches');
title('Percent of subjs with correct matches using 5% of data');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
set(gca,'xtick',[0:15]);
set(gca,'xticklabel',{'','anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'});

classfacs = [];classmeans=[];idxs=[];
for nx = 1:length(keepuseims)
     s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    dat=floatread([fullpaths{nx},savedat,'.fdt'],[s.pcs inf],[],0);        
    ws = wts*sph;    acts = ws*dat; 
   for im = 1:length(keepuseims{nx})
       cps = allbigs{nx}{keepuseims{nx}(im)};
       for c = 1:length(cps)
           rcp = find(s.complist == cps(c));
           classfacs = [classfacs;acts(keepuseims{nx}(im),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)];
           idxs = [idxs;[nx,keepuseims{nx}(im),cps(c)]];
           classmeans = [classmeans;s.meanpwr(rcp,:)];
       end;
   end;
end

[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(classfacs,idxs,classmeans,s.freqs,[3 128],[],'quad');
[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,s.freqs,[3 128],[],'quad');

 [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot([{thetaclust{3}},{alphaclust{3}},{betaclust{3}},{gamaclust{3}}],allbigs,bigwts,orivec);
          
row = length(comods);
viewnum=[1,2,3];col = 3;%length(viewnum) ;
zoom= 1.3;
figure;pl = 1;
for clust = 1:length(comods)
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

for cls = 1:length(denslist) 
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,denslist{cls},[],[],[],{'mrislices',[63:-16:-17],'mriview','top','geom',[4,3]},'bred',[],'on');
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,denslist{cls},[],gdcomps,[],{'mrislices',[63:-16:-17],'mriview','top','geom',[4,3]},'yred',[],'on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster emomedians with pdist
clustmeans = [];clustkeep = [];minmean = 1;
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    for t = 1:size(emomeans{nx},1)
        if max(abs(emomeans{nx}(t,:))) > minmean
            clustmeans = [clustmeans; emomeans{nx}(t,:)];
            clustkeep = [clustkeep; [nx t]];
        end;
    end;
    %clustkeep = [clustkeep; [repmat(nx,[size(emomeans{nx},1) 1]) [1:size(emomeans{nx},1)]']];
end;
nclust = 2;
alldist = pdist(clustmeans, 'correlation'); % euc better than seuc
links = linkage(alldist,'complete'); 
figure;[hnd,idx,perm]=  dendrogram(links,nclust);%close
figure;row = round(sqrt(nclust));col = ceil(sqrt(nclust));
row = 2; col = 2; cols = jet(15);pl = 1; clear ktsubj currwv
stdcut = .75; % 1.5 then 1
for cls = 1:nclust
    if pl == row*col+1
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        figure; pl = 1;  end;
    onecls = find(idx == cls);
    currwv{cls} = clustmeans(onecls,:);
    ktsubj{cls} = clustkeep(onecls,:);
% eliminate by std cut:-----
zs = zscore(currwv{cls});
delmem = [];
for mem = 1:size(currwv{cls},1)
    if find(median(abs(zs(mem,:))) > stdcut)
        delmem = [delmem mem];
    end;
end;
if ~isempty(delmem)
    currwv{cls}(delmem,:) = [];
    ktsubj{cls}(delmem,:) = [];
end;
    sbplot(row,col,pl);  pl = pl+1;
    %ph=plot(currwv{cls}','k-','linewidth',.5);hold on;    
    for e = 1:size(currwv{cls},2)
        ph=plot(e,currwv{cls}(:,e),'.');hold on;
        set(ph,'markersize',5);set(ph,'color',cols(e,:));
        %ph = text(e,mean(currwv{cls}(:,e),1),emo2{e});
        %set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-');
    ph=plot(mean(currwv{cls},1),'k-','linewidth',1.5);hold on;  
    set(gca,'xlim',[0 16]);
    title(['Cls ',int2str(cls)]);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

cls=1;plotlist = cell(1,35);templs = cell(1,35);
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*data;   
    clear wts sph ws allfacs alltemps onefac onemean kt 
    relims = ktsubj{cls}(find(ktsubj{cls}(:,1)==nx),2)';
    plotims{nx} = relims;
    for r=1:length(relims)
        clear wts sph ws allfacs alltemps onefac onemean kt 
        for rcp = 1:length(s.complist)
            alltemps(rcp,:) = activations(ims{nx}(1),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
        end;
        relcps = find(ismember(s.complist,allbigs{nx}{relims(r)}));
        templs{nx} = [templs{nx};alltemps(relcps,:)];
        plotlist{nx} = [plotlist{nx},allbigs{nx}{relims(r)}];
    end;
end;
figure;row = 1; col=2; 
[density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, plotlist,[],[dipargs,'view',[0 0 1]],[],{'r'},'off',[],[]);  view(60,30)

figure; row=6; col=6;pl=1;
for nxx = 1:length(subjlist)
    nx=subjlist(nxx);
    if ~isempty(templs{nx})
    if pl > row*col
        figure;pl=1;
    end;
    sbplot(row,col,pl);pl = pl+1;
    quadplot(s.freqs,templs{nx});title(['Subj ',int2str(nx)]);
    %SpecCoModPlot('sources.set',fullpaths{nx},[],[ims{nx}],savedat,[3 128],'n',0,[]);
    SpecCoModPlot('sources.set',fullpaths{nx},[],[],savedat,[3 128],'n',0,plotims{nx});
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useimssub = cell(1,max(subjlist)); 
oricell = cell(1,max(subjlist)); 
for nxx = 1:length(subjlist) % 3 and 33 need 12 dims, not rank(13)
    nx = subjlist(nxx);
    for clust = 5%1:length(strs)
        eval(strs{clust})
        for cls = 1:length(finalidx)
            onesubj = finalidx{cls}(find(finalidx{cls}(:,1) == nx),2);
            if ~isempty(onesubj)
                useimssub{nx} = [useimssub{nx} unique(abs(onesubj))'];  
            end;
        end;
        for im = 1:length(useimssub{nx})
            idx = find(allbigs{nx}{useimssub{nx}(im)} == onebig{nx}{useimssub{nx}(im)});
            oricell{nx}(im) = orivec{nx}{useimssub{nx}(im)}(idx);
        end;
    end;
end;

figure; row = 6;col = 6;
for nxx = 1:length(subjlist) % 3 and 33 need 12 dims, not rank(13)
    nx = subjlist(nxx);
    sbplot(row,col,nxx)
    [fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,oricell,2,nx,useimssub,'mds',0,['Subject ',int2str(nx)],1,'elim','noperm');
end;

% reorder fullmd and add deleted emos as mean            
% $$$             newtempls1{nx} = fullmd{nx}./repmat(std(fullmd{nx},0,1),[size(fullmd{nx},1) 1]);
% $$$             for e = 1:length(emos)
% $$$                 try
% $$$                     newtempls2{nx}(e,:) = newtempls1{nx}(find(strcmp(sems{nx},emos{e})),:);
% $$$                 catch % if that subject doesn't have that emo, replace with mean
% $$$                     newtempls2{nx}(e,:) = mean(newtempls1{nx},1);
% $$$                 end;   
% $$$             end;    
% $$$             close;
% $$$             clustwv(end+1:end+size(newtempls2{nx},2),:) = newtempls2{nx}';
% $$$             bfdim(1,nx) = input('Best fit dim: ');
% $$$             if bfdim(1,nx) < 0 | bfdim(1,nx) ~= 0
% $$$                 [val mnwt(1,nx)] = max(fullwts{nx}(abs(bfdim(1,nx)),:));% flip these to compare
% $$$                 [val mxwt(1,nx)] = min(fullwts{nx}(abs(bfdim(1,nx)),:));% across subjects 
% $$$             elseif bfdim(1,nx) ~= 0
% $$$                 [val mxwt(1,nx)] = max(fullwts{nx}(bfdim(1,nx),:));
% $$$                 [val mnwt(1,nx)] = min(fullwts{nx}(bfdim(1,nx),:));
% $$$             end;
      close;  
    kt = [kt; [repmat(nx,size(fullmd{nx},2),1) [1:size(fullmd{nx},2)]']];% subj,IM
    clustwv(end+1:end+size(fullmd{nx},2),:) = (fullmd{nx}./repmat(std(fullmd{nx},0,1),[size(fullmd{nx},1) 1]))'; % to skip all selection processes
end;
%%%%%%%%%%%%%%%
%%  Find dim with max(abs(diff in sad))
%%%%%%%%%%%%%%
em = 7; % sad
clear ims templs nxwts wts sph data activations 
for nxx = 1:length(subjlist)
    nx=subjlist(nxx);
    [x maxdim(1,nx)]=max(abs(fullmd{nx}(em,:)));
    %[x maxim] = max(abs(fullwts{nx}(maxdim(1,nx),:)));
    maxim = find(abs(fullwts{nx}(maxdim(1,nx),:)) > (mean(fullwts{nx}(maxdim(1,nx),:))+2*std(fullwts{nx}(maxdim(1,nx),:))));
    ims{nx} = keeptrack{nx}(maxim,2)';
    nxwts{nx} = fullwts{nx}(maxdim(1,nx),maxim);
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*data;   
    clear wts sph ws allfacs alltemps onefac onemean kt 
    for rcp = 1:length(s.complist)
        alltemps(rcp,:) = activations(ims{nx}(1),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
    end;
    relcps = find(ismember(s.complist,allbigs{nx}{ims{nx}(1)}));
    if nxwts{nx}(1) < 0
    templs{nx} = alltemps(relcps,:)*-1;
    plotlistn{nx} = allbigs{nx}{ims{nx}(1)};
    else
    templs{nx} = alltemps(relcps,:);
    plotlistp{nx} = allbigs{nx}{ims{nx}(1)};
    end;
end;
figure;row = 1; col=2; 
sbplot(row,col,1)
    [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, plotlistp,[],[dipargs,'view',[0 0 1]],[],{'r'},'off',[],[]);  view(60,30)
sbplot(row,col,2)
    [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, plotlistn,[],[dipargs,'view',[0 0 1]],[],{'b'},'off',[],[]);  view(60,30)
figure; row=6; col=6;pl=1;
for nxx = 1:length(subjlist)
    nx=subjlist(nxx);
    if pl > row*col
        figure;pl=1;
    end;
    sbplot(row,col,pl);pl = pl+1;
    quadplot(s.freqs,templs{nx});title(['Subj ',int2str(nx)]);
    %SpecCoModPlot('sources.set',fullpaths{nx},[],[ims{nx}],savedat,[3 128],'n',0,[]);
    %SpecCoModPlot('sources.set',fullpaths{nx},[],[],savedat,[3 128],'n',0,ims(1,nx));
end;


figure; cols = jet(15);cols(10,:) = [1 .9 0];row = 4; col=4;pl = 1;
for nx = 1:row*col
    %emo2 = sems{nx};
    %dd = pdist(clustwv', 'correlation') ;    
    dd = pdist(fullmd{nx}, 'euclidean') ;    
    [md,mwts] = cmdscale(dd);
    if pl > row*col
        set(gcf,'color','w');figure; pl=1;
    end;
    sbplot(row,col,pl);pl=pl+1;
    c1 = 1; c2 = 2; c3 = 3;
    for e = 1:size(md,1)
        ph=plot3(md(e,c1),md(e,c2),md(e,c3),'.');hold on;
        set(ph,'markersize',25);                set(ph,'color',cols(e,:));
        ph = text(md(e,c1),md(e,c2),md(e,c3),emo2{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
    end;
    zl = get(gca,'zlim');
    for e = 1:size(md,1)
        ph =plot3([md(e,c1) md(e,c1)],[md(e,c2) md(e,c2)],[zl(1)  md(e,c3)]);
        set(ph,'color',cols(e,:)); set(ph,'linewidth',2)             
    end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
    xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);zlabel(['Dim ',int2str(c3)]);
    title(['Subj ',int2str(nx)]);
end;
            
selemospc = [];clear mnallics mnreltemp mxreltemp mxallics
for nx = 1:length(mnwt)
    if mnwt(1,nx) ~= 0
        if bfdim(nx) < 0
        selemospc = [selemospc;fullmd{nx}(:,abs(bfdim(1,nx)))'/std(fullmd{nx}(:,abs(bfdim(1,nx))))*-1];
        else
        selemospc = [selemospc;fullmd{nx}(:,bfdim(1,nx))'/std(fullmd{nx}(:,bfdim(1,nx)))];
        end;
        mnallics{nx} = onebig{nx}{mnwt(1,nx)};
        rtidx = find(kptk(:,1)==nx&kptk(:,2)==mnwt(1,nx)&kptk(:,3)==mnallics{nx});
        mnreltemp(nx,:) = clustfacs(rtidx,:);
        mxallics{nx} = onebig{nx}{mxwt(1,nx)};% same for max
        rtidx = find(kptk(:,1)==nx&kptk(:,2)==mxwt(1,nx)&kptk(:,3)==mxallics{nx});
        mxreltemp(nx,:) = clustfacs(rtidx,:);
    end;
end
figure;
ph =plot([0 16],[0 0],'k-'); hold on;
for e = 1:size(md,1)
    ph=plot(e,selemospc(:,e),'.');hold on;
    set(ph,'markersize',25);                set(ph,'color',cols(e,:));
    ph = text(e,mean(selemospc(:,e)),emo2{e});set(ph,'rotation',45);
    set(ph,'color',cols(e,:)); set(ph,'fontsize',16); 
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

dipargs = {'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','on','projimg','off','coordformat','spherical'};
figure; pl=1;row = 4; col=6;
for nx = 1:size(mnreltemp,1)
    if mnwt(1,nx) ~= 0
        if pl > row*col
            set(gcf,'color','w');figure; pl=1;
        end;
        sbplot(row,col,pl); pl=pl+1;
        ph = quadplot(freqs,mnreltemp(nx,:)*-1,2,'b'); hold on;
        ph = quadplot(freqs,mxreltemp(nx,:),2,'r'); hold on;
        title(['Subj ',int2str(nx)]);
        sbplot(row,col,pl); pl=pl+1; clear onesubj
        onesubj{nx} = [mxallics{nx} mnallics{nx}];
    [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, onesubj,[],[dipargs,'view',[0 0 1]],[],{'r','b'},'off',[],[]);  view(60,30)
    end;
end;set(gcf,'color','w');
 figure;    [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, mnallics,[],[dipargs,'view',[0 0 1]],[],{'b'},'off',[],[]);
figure; quadplot(freqs,mnreltemp*-1)
  figure;    [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, mxallics,[],[dipargs,'view',[0 0 1]],[],{'r'},'off',[],[]);
figure; quadplot(freqs,mxreltemp)
 
sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[dipargs,'view',[0 0 1]],[],[],'off',[],[]);place = place +1;camzoom(1.1)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find angle between each subj pair emotion spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/jason/mbin/')
clear a
for nx1 = 1:35
    for nx2 = 1:35
        if ~isempty(fullmd{nx1}) & ~isempty(fullmd{nx2})
            a(nx1,nx2) = subspace(fullmd{nx1},fullmd{nx2});
            a(nx1,nx1) = -2;
        end;
    end;
end;
[mo,ord] = arr3(a);
%figure; imagesc(mo);
a(find(a == 0)) = 2; a(find(a == -2)) = 2; 
likesubjs = [];unrellist = [];others = [];
diffsubjs = [];subjmatch = zeros(0,2);
for nx = 1:size(a,1)
    [x y] = min(a(nx,:));
    if x < .75 % likeness cutoff         
        subjmatch(end+1,:) = [nx y];    
    end;
end;
for nx = 1:size(a,1)
    nmatch(1,nx) = length(find(subjmatch(:,2) == nx));
end;
[vals subjord] = sort(nmatch,'descend');
rellist = subjord(1); p = 1;clear corrsubjs
corrsubjs{1} = subjord(1);new= 0;
for nxx = 1:length(subjord)
    nx = subjord(nxx);
    %[x y] = min(a(nx,:));
    [x y] = find(a(nx,:)< .95);
    likesubjs{nx} = y;
    for pp = 1:length(corrsubjs)
        if ismember(nx,corrsubjs{pp})|~isempty(find(ismember(y,corrsubjs{pp})))
            corrsubjs{pp} = [corrsubjs{pp} nx y]; new= 0;break;
        else
            new = 1; 
        end;
    end;
    if new == 1
        p = p+1;
        corrsubjs{p} = [nx y];
    end;            
end;
% check for repeated subjs
repsub = [];
for p1 = 1:length(corrsubjs)
    for p2 = p1+1:length(corrsubjs)
        if ~isempty(find(ismember(corrsubjs{p1},corrsubjs{p2})))
            repsub = [repsub corrsubjs{p1}(ismember(corrsubjs{p1},corrsubjs{p2}))];
            p1del = ismember(corrsubjs{p1},corrsubjs{p2});
            p2del = ismember(corrsubjs{p2},corrsubjs{p1});
            corrsubjs{p1}(p1del) = [];
            corrsubjs{p2}(p2del) = [];
        end;
    end;
end;
for r = 1:length(repsub)
    nx = repsub(r);
    [x y] = min(a(nx,:));
    for pp = 1:length(corrsubjs)
        if ismember(y,corrsubjs{pp})
            corrsubjs{pp} = [corrsubjs{pp} nx]; break;            
        end;
    end;
end;
for pp = 1:length(corrsubjs)
    corrsubjs{pp} = unique(corrsubjs{pp});
end;
%---------------------------------------------- 
% ica on emo spaces:
c=[];k=[];
for n  = 1:length(corrsubjs{1})
    c=[c; clustwv(find(kt(:,1)==corrsubjs{1}(n)),:)];
    k = [k;kt(find(kt(:,1)==corrsubjs{1}(n)),:)];
end;
cols = jet(15);
[wts,sph,compvars,bias,signs,lrates,acts] = runica(c','pca',rank(clustwv'),'extended',1,'stop',1e-7,'maxsteps',6000);
winv = pinv(wts*sph); 
figure;
row = round(sqrt(size(winv,2)));col = ceil(sqrt(size(winv,2)));
for fd = 1:size(winv,2)
    sbplot(row,col,fd)
    ph = plot([0 size(winv,1)+1],[0 0],'k-'); hold on;
    for e = 1:size(winv,1)
        ph=plot(e,winv(e,fd),'.');hold on;
        set(ph,'markersize',15);set(ph,'color',cols(e,:));
        ph = text(e,winv(e,fd),emo2{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;          
end;         
figure; pl = 1; clear ktsubj currwv
for cls = 1:size(winv,2)
    if pl == row*col+1
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        figure; pl = 1;  end;
    onecls = find(abs(acts(cls,:)) > mean(acts(cls,:))+2*std(acts(cls,:)));
    currwv{cls} = clustwv(onecls,:);
    for w = 1:size(currwv{cls},1)
        [corr,indx,indy,corrs] = matcorr(currwv{cls}(w,:),winv(:,cls)');
        if corr < 0
            currwv{cls}(w,:) = currwv{cls}(w,:)*-1;
        end;
    end;
    for w = 1:size(currwv{cls},1)
        [corr,indx,indy,corrs] = matcorr(currwv{cls}(w,:),mean(currwv{cls},1));
        if corr < 0
            currwv{cls}(w,:) = currwv{cls}(w,:)*-1;
        end;
    end;
    ktsubj{cls} = kt(onecls,:); % keeps track of [subj,dim]
    sbplot(row,col,pl);  pl = pl+1;
    for e = 1:size(currwv{cls},2)
        ph=plot(e,currwv{cls}(:,e),'.');hold on;
        set(ph,'markersize',5);set(ph,'color',cols(e,:));
        %ph = text(e,mean(currwv{cls}(:,e),1),emo2{e});
        %set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-');
    ph=plot(mean(currwv{cls},1),'k-','linewidth',.5);hold on;  
    set(gca,'xlim',[0 16]);
    title(['Cls ',int2str(cls)]);
end;

% Cluster individual emo spaces:---------
% flip all vectors to double matrix:
x = clustwv*-1;
clustwv(end+1:end+size(x,1),:) = x;
y = kt; y(:,2) = y(:,2)*-1;
kt(end+1:end+size(x,1),:) = y;

nclust = 20;
alldist = pdist(clustwv, 'correlation'); % euc better than seuc
links = linkage(alldist,'complete'); 
figure;[hnd,idx,perm]=  dendrogram(links,nclust);%close
figure;row = round(sqrt(nclust));col = ceil(sqrt(nclust));
row = 4; col = 4; cols = jet(15);pl = 1; clear ktsubj currwv
%pclust = [1,2,3,4,6,7,8,9,13,16,17,18,19,21,26];
%pclust = [1:10,16:17,22,23,29];
pclust = [1:nclust]; % for no template flipping
for clss = 1:length(pclust)
    cls = pclust(clss);
    if pl == row*col+1
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        figure; pl = 1;  end;
    onecls = find(idx == cls);
    currwv{cls} = clustwv(onecls,:);
    ktsubj{cls} = kt(onecls,:);
    sbplot(row,col,pl);  pl = pl+1;
    ph=plot(currwv{cls}','k-','linewidth',.5);hold on;    
    for e = 1:size(currwv{cls},2)
        ph=plot(e,currwv{cls}(:,e),'.');hold on;
        set(ph,'markersize',5);set(ph,'color',cols(e,:));
        %ph = text(e,mean(currwv{cls}(:,e),1),emo2{e});
        %set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-');
    ph=plot(mean(currwv{cls},1),'k-','linewidth',.5);hold on;  
    set(gca,'xlim',[0 16]);
    title(['Cls ',int2str(cls)]);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%save /data/common2/emotion/ClusterDims.mat fullmd fullwts  currwv ktsubj clustwv kt idx
load /data/common2/emotion/ClusterDims.mat 
% prune clusters to be tig
for cls = 1:length(currwv)
    cutcorr = .05;
    while cutcorr < .4
        mntmp = median(currwv{cls},1); clear xc
        %figure;
        %for e = 1:size(mntmp,2)
        %    ph=plot(e,mntmp(:,e),'.');hold on;
        %    set(ph,'markersize',5);set(ph,'color',cols(e,:));
        %    ph = text(e,mntmp(:,e),emo2{e});
        %    set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
        %end;
        for m = 1:size(currwv{cls},1)
            xc(1,m) = (currwv{cls}(m,:)*mntmp')/(mntmp*mntmp');
        end;
        delsb = find(xc < cutcorr);
        cutcorr = cutcorr + .01;
        ktsubj{cls}(delsb,:)=[];
        currwv{cls}(delsb,:) = [];
    end;
end;
figure;
row = 4; col = 4; cols = jet(15);pl = 1; 
for clss = 1:length(pclust)
    cls = pclust(clss);
    if pl == row*col+1
        figure; pl = 1;  end;
    sbplot(row,col,pl);  pl = pl+1;
    for e = 1:size(currwv{cls},2)
        %ph=plot(e,currwv{cls}(:,e),'.');hold on;
        ph=plot(e,median(currwv{cls}(:,e),1),'.');hold on;
        set(ph,'markersize',15);set(ph,'color',cols(e,:));
        ph = text(e,median(currwv{cls}(:,e),1),emo2{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k-');
    title(['Cls ',int2str(cls)]);
end;
% plot clusters ordered by male/female of emo order:-----
% $$$ imgplot = ones(35,length(currwv))*-1;
% $$$ for cls = 1:length(ktsubj)
% $$$     imgplot(ktsubj{cls}(:,1),cls) = 1;
% $$$ end;
% $$$ currlist = alllist;
% $$$ currlist = malefem;
% $$$ sortlist = [];
% $$$ for sj = 1:length(currlist)
% $$$     sortlist = [sortlist currlist{sj}];
% $$$ end;
% $$$ figure; imagesc([1:length(currwv)],sortlist,imgplot(sortlist,:));  

% Collect IM indices for each cluster:--------------
clear clustims sides newkt
for cls = 1:length(ktsubj)
    kk = [];
    for nxx = 1:size(ktsubj{cls},1)
        nx = ktsubj{cls}(nxx,1);
        sbjdim = ktsubj{cls}(find(ktsubj{cls}(:,1)==nx),2);tpims=[];tpside = [];
        for dim = 1:length(sbjdim)
            cutval = mean(fullwts{nx}(sbjdim(dim),:)) + std(fullwts{nx}(sbjdim(dim),:));
            tpims = [tpims find(abs(fullwts{nx}(sbjdim(dim),:)) > cutval)];
            tpside = [tpside fullwts{nx}(sbjdim(dim),find(abs(fullwts{nx}(sbjdim(dim),:)) > cutval))./abs(fullwts{nx}(sbjdim(dim),find(abs(fullwts{nx}(sbjdim(dim),:)) > cutval)))];
        end;
        clustims{cls}{nx} =tpims;
        sides{cls}{nx} =tpside;
        kk = [kk;[repmat(nx,[length(tpims) 1]) tpims']];
    end;
    newkt{cls} = kk;
end;
% Plot results here:------------------------
% plot all templates from each member of a cluster
% plot mean spectra with back projections
frqlim = [3 125];
cls = 5;
for nx = 1:length(clustims{cls})    
    if ~isempty(clustims{cls}{nx})
        cps = [];
        for i = 1:length(clustims{cls}{nx})
            cps = [cps cell2mat(allbigs{nx}(clustims{cls}{nx}(i)))];
        end; cps = unique(cps);
        PlotSpecFacEnv('sources.set',savedat,fullpaths{nx},clustims{cls}{nx},cps,[],sides{cls}{nx},frqlim,1,.99,0);
    textsc(['Subject ',int2str(nx),'; Cls ',int2str(cls)],'title');
    end;
end;
          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Describe each cluster by the Spectral cluster representations:-------------
clear clustreps iclistp iclistn
for clust = 1:length(ktsubj)
    figure; row = 6; col = 2; 
    pl = 1; maxy = [];
    clustplus = [];clustminus = [];allfacs = [];
    clustsubjs = unique(ktsubj{clust}(:,1));
    for nxx = 1:length(clustsubjs)
        nx = clustsubjs(nxx); clear allidx2
        doi = ktsubj{clust}(find(ktsubj{clust}(:,1) == nx),2);% dims of interest
        for dimm = 1:length(doi)
            if pl > row*col
                textsc(['Cluster ',int2str(clust)],'title');        
                set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); close
                figure; pl=1;
            end;
            dim = abs(doi(dimm));                
            cutval = mean(fullwts{nx}(dim,:)) +  .5*std(fullwts{nx}(dim,:));% % of max
            %cutval = max(abs(fullwts{nx}(dim,:))) -  max(abs(fullwts{nx}(dim,:)))*.5;% % of max
            ifacs = find(abs(fullwts{nx}(dim,:)) > cutval);
            subjfacs = keeptrack{nx}(ifacs,:);
            multfac = fullwts{nx}(dim,find(abs(fullwts{nx}(dim,:)) > cutval))./abs(fullwts{nx}(dim,find(abs(fullwts{nx}(dim,:)) > cutval)));
            if doi(dimm) < 0
            subjfacs(:,2) = subjfacs(:,2).*(multfac*-1)';%account for inversion
            else
            subjfacs(:,2) = subjfacs(:,2).*multfac';% invert if neg weighting
            end;
            [clsmem,nims,allidx2,iclp,icln,tsn,tsp,maxy] = SortbyClusters2(subjfacs,row,col,pl,[]);
            allfacs = [allfacs; subjfacs]; % save for later
            title(['Subj ',int2str(nx),'; Dim ',int2str(dim)]);pl = pl+1;
            clustplus = [clustplus; clsmem(1,:)];
            clustminus = [clustminus; clsmem(2,:)];
        end;
        textsc(['Cluster ',int2str(clust)],'title');        
    end;close
    clustreps{clust}(1,:) = sum(clustplus,1);
    clustreps{clust}(2,:) = sum(clustminus,1);
    % run again to get iclist:
    figure; row=1; col = 1;pl=1;
    [clsmems{clust},nims,allidx2,iclistp{clust},iclistn{clust},templatesn{clust},templatesp{clust},maxy] = SortbyClusters2(allfacs,row,col,pl,[]);
            title(['Cluster ',int2str(clust)]);pl = pl+1;
    fprintf('\nCluster %s done.',int2str(clust));
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
% plot the output as the sum of ICs in each spectral category
figure; pl = 1;
for clust = 1:60
    if pl > 30
        figure; pl = 1;
    end;
    sbplot(5,6,pl);pl = pl+1;
    ph = plot([0 15],[0 0],'k-');hold on;
    ph = plot(clustreps{clust}(1,:),'r^');
    ph = plot(clustreps{clust}(2,:),'bv');
    ph = plot([2.5 2.5],[get(gca,'ylim')],'m-');
    ph = plot([6.5 6.5],[get(gca,'ylim')],'m-');
    ph = plot([9.5 9.5],[get(gca,'ylim')],'m-');
    ph = plot([12.5 12.5],[get(gca,'ylim')],'m-');    
    title(['Clust ',int2str(clust)]);
end;
%save /data/common2/emotion/ClusterDimICs.mat iclistp iclistn
load /data/common2/emotion/ClusterDimICs.mat iclistp iclistn

% plot the cluster dipoles:*************************

ttl = {'Delta1','Delta2','Theta1','Theta2','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Beta3','HiLow','BB','Peaked'};
viewnum=[1,2];
pclust = [1:15];
for clt = 1:length(iclistp) % emo space clusters
    clust = pclust(clt);
    figure;row = 5; col = 6;  place = 1;pg = 1;
    for c = 1:length(iclistp{clust})
        if clsmems{clust}(1,c) > 3
            if place > row*col
                ph = textsc(['Power increases for Cluster ',int2str(clust)],'title'); set(ph,'fontsize',14); set(ph,'color','c');set(gcf,'color','w');
                str = ['print /data/common2/emotion/ClustFigs/Clust',int2str(clust),ttl{c},'PwrInc',int2str(pg),'.jpg -djpeg'];eval(str); %close
                figure; place = 1;pg = pg+1;
            end;
            sbplot(row,col,place); coltempls = [];
            for nx = 1:length(templatesp{clust}{c})
                coltempls = [coltempls; templatesp{clust}{c}{nx}];
            end;
            ph = quadplot(freqs,coltempls,1.5,'r');place = place+1;
            PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistp{clust},c,row,col,place,ttl{c},viewnum,[]);
            place = place+length(viewnum);
        end;
    end;
    ph = textsc(['Power increases for Cluster ',int2str(clust)],'title'); set(ph,'fontsize',14); set(ph,'color','k');set(ph,'color','c');set(gcf,'color','w');
    str = ['print /data/common2/emotion/ClustFigs/Clust',int2str(clust),ttl{c},'PwrInc',int2str(pg),'.jpg -djpeg'];eval(str); 
    figure;row = 5; col = 6;  place = 1;pg = 1;
    for c = 1:length(iclistn{clust})
        if abs(clsmems{clust}(2,c)) > 3
            if place > row*col
                ph = textsc(['Power decreases for Cluster ',int2str(clust)],'title'); set(ph,'fontsize',14); set(ph,'color','c');set(gcf,'color','w');
                str = ['print /data/common2/emotion/ClustFigs/Clust',int2str(clust),ttl{c},'PwrDec',int2str(pg),'.jpg -djpeg'];eval(str); 
                figure; place = 1;pg = pg+1;
            end;
            sbplot(row,col,place); coltempls = [];
            for nx = 1:length(templatesn{clust}{c})
                coltempls = [coltempls; templatesn{clust}{c}{nx}];
            end;
            ph = quadplot(freqs,coltempls,1.5,'b');place = place+1;
            PlotDipoleClusters('sources.set',fullpaths,gdcomps,iclistn{clust},c,row,col,place,ttl{c},viewnum,[0 0 1]);
            place = place+length(viewnum);
        end;
    end;
    ph = textsc(['Power decreases for Cluster ',int2str(clust)],'title'); set(ph,'fontsize',14); set(ph,'color','k');set(ph,'color','c');set(gcf,'color','w');
    str = ['print /data/common2/emotion/ClustFigs/Clust',int2str(clust),ttl{c},'PwrDec',int2str(pg),'.jpg -djpeg'];eval(str); 
end;


%-------------------------------------------------------


%%%% correlate emo spaces pairwise, then plot by highest
for nx1 = 1:35
    for nx2=1:35
        [corr(nx1,nx2,:),indx(nx1,nx2,:),indy(nx1,nx2,:),corrs] = matcorr(fullmd{nx1}',fullmd{nx2}');          
    end;
end;
indx(find(corr == 1)) = 0;
indy(find(corr == 1)) = 0;
corr(find(corr == 1)) = 0;corr = abs(corr);
addpath('/home/jason/mbin/')
[mo,ord] = arr3(corr(:,:,1));% finds the order of subjects with the most in common with all others (according to the highest correlation)
% pick out mode of each subj highest correlated IM:
plotmds = zeros(15,0);plotcorrs = [];
for nx = 1:size(corr,1)
    mds(1,nx) = mode(indx(nx,:,1));
    corr(nx,find(indx(nx,:,1) == mds(1,nx)),1)
    plotmds{nx} = [fullmd(:,mds(1,nx),nx) fullmd(:,mds(1,nx),nx)];
    plotcorrs = [plotcorrs corr(nx,find(indx(nx,:,1) == mds(1,nx)),1)];
end;
figure; hist(plotcorrs,100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IM movie  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frqlim = [3 125]; 
datset = 'sources.set';
nx=2;
   
dim = 3;
weights = fullwts(dim,:);weights = weights/max(weights);
comp = [2,3,5,10]; im = useims{nx};  em = 10;
[M] = IMmovie(datset,savedat,fullpaths{nx},im,comp,frqlim,['Cp',int2str(comp),'EmSpc',int2str(dim),emos{em},'.avi'],[em],weights','off',{'r'},{emos{em}});
   
figure;movie(M)

SpecCoModPlot('sources.set',fullpaths{nx},[],[1:15],savedat,[3 128],'n',0,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% describe each dimension individually:
cols = jet(size(fullmd,1));
figure; row = 4; col = 4; pl = 1; maxy = [];
for dim = 1:size(fullwts,1)
    if pl > row*col
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        figure; pl=1;
    end;
    wtcut = max(abs(fullwts(dim,:))) -  max(abs(fullwts(dim,:)))*.30;
    ifacs{1,dim} = find(abs(fullwts(dim,:)) > wtcut);
    subjfacs{1,dim} = keeptrack(ifacs{1,dim},:);
    multfac = fullwts(dim,find(abs(fullwts(dim,:)) > wtcut))./abs(fullwts(dim,find(abs(fullwts(dim,:)) > wtcut)));
    subjfacs{1,dim}(:,2) = subjfacs{1,dim}(:,2).*multfac';
    sbplot(row,col,pl)
    plot([0 size(fullmd,1)+5],[0 0],'k-'); hold on;
    for e = 1:size(fullmd,1)
        ph=plot(e,fullmd(e,dim),'.');hold on;
        set(ph,'markersize',15);set(ph,'color',cols(e,:));
        ph = text(e,fullmd(e,dim),emos{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
    end;
    title(['Emo Space Dim ',int2str(dim)]); pl = pl+1;
    [clsmems{dim},nims,allidx{dim},allidx2{dim},maxy] = SortbyClusters2(subjfacs{dim},row,col,pl,[]);
    pl = pl+1;      
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- OR -- desribe each emotion individually:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find loading on each emotion. Common cutoff:
clear idims ifacs subjfacs allidx allidx2
cutoff = mean(std(fullmd(:,2:end),0,2)*1.5);% take out dim1 for deciles
wtcut = mean(std(fullwts,0,2)*2);
for em = 1:length(emos)
    idims{em} = find(abs(fullmd(em,:)) > cutoff);% important dims
    idims{em} = idims{em}.*(fullmd(em,idims{em})./abs(fullmd(em,idims{em})));
    % orient to be pos or neg
end;
for em = 1:length(emos)
    for f = 1:length(idims{em}) % number of dims emo is involved with
        ifacs{em,f} = find(abs(fullwts(abs(idims{em}(f)),:))>wtcut);
        mfact = fullwts(abs(idims{em}(f)),ifacs{em,f});
        mfact = mfact./abs(mfact);% retrieve actual orientation
        subjfacs{em}{f} = keeptrack(ifacs{em,f},:);
        if idims{em}(f) < 0 
            mfact = mfact*-1;% invert if emo is neg in emo space 
                             % because we want the case where the emo is
                             % positively weighted
        end;        
        subjfacs{em}{f}(:,2) = subjfacs{em}{f}(:,2) .*mfact';
        % record the factor sign if emo is pos in emo space dim
    end;
end;
figure; row = 4; col = 4; pl = 1; maxy = 195;
for em = 1:length(subjfacs)
    [clsmems,nims,allidx{em},allidx2{em},my] = SortbyClusters(subjfacs{em},row,col,pl,maxy);
    title(emos{em}); pl = pl+1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%  Visualizing results:----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%create new wts indicator column according to template orientation:
%%%****  Need to do this for most plotting applications below!***
for clust =1:5
    for c = 1:size(allidx2{em},2)
        for em = 1:length(emos)
            if ~isempty(allidx2{em}{clust,c})
                for mem = 1:size(allidx2{em}{clust,c},1)
                    currt = find(kptk(:,1) == allidx2{em}{clust,c}(mem,1)&...
                                 kptk(:,2) == abs(allidx2{em}{clust,c}(mem,2))&...
                                 kptk(:,3) == allidx2{em}{clust,c}(mem,3));
                    [bignum  idx] = max(abs(clustfacs(currt,:)));
                    % if pointed down, switch weighting
                    if clustfacs(currt,idx) < 0 & allidx2{em}{clust,c}(mem,2)>0
                        allidx2{em}{clust,c}(mem,4) = -1;% -1=pwr dec
                    elseif clustfacs(currt,idx) < 0 & allidx2{em}{clust,c}(mem,2)<0
                        allidx2{em}{clust,c}(mem,4) = 1;% 1=power inc
                    elseif clustfacs(currt,idx) > 0 & allidx2{em}{clust,c}(mem,2)<0
                        allidx2{em}{clust,c}(mem,4) = -1;% 1=power dec
                    elseif clustfacs(currt,idx) > 0 & allidx2{em}{clust,c}(mem,2)>0
                        allidx2{em}{clust,c}(mem,4) = 1;% 1=power inc
                    end;                    
                end;
            end;
        end;
    end;
end;
% plot the cluster dipoles:*************************
speccat = {'Delta','Theta','Alpha','Beta','Gamma'};
load /data/common4/emotion/AllClustFacsNEW.mat    
for clust = 1:5 % spectral categories (delta,theta,alpha...) 
    for c = 1:size(allidx2{em},2)  % clusters within categories
        figure;row = 4; col = 4; zoom = 1.1; viewnum=[1,2,3,4]; pl = 1;
        for pn = 1:2 % pos/neg factors
            for em = 1:length(allidx2) % emotions or dims
                if ~isempty(allidx2{em}{clust,c})
                    if pl > row*col
                        ph = textsc([speccat{clust},', Subclust ',int2str(c),' modulators'],'title'); set(ph,'fontsize',14);
                        figure; pl = 1;
                    end;
                    if pn == 1
                        clear facvec comods justcomps wtsmat1 jcwts denslist
                        ttl = [emos{em},';  Pwr incr'];
                        %str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',emos{em},'POS.jpg -djpeg'];
                    else
                        onegrp = allidx2{em}{clust,c}(find(allidx2{em}{clust,c}(:,4) == -1),:);
                        onegrp = abs(onegrp);
                        ttl = [emos{em},';  Pwr dec'];
                        %str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',emos{em},'NEG.jpg -djpeg'];
                    end;                
                    if ~isempty(onegrp)
                        clear facvec comods justcomps wtsmat1 jcwts denslist                
                        [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(gdcomps,fullpaths,{onegrp},allbigs,bigwts,orivec);
                        % plot one clust:        
                        [angles] = PlotCoModasDipoles(comods{1},justcomps{1},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{1},jcwts{1},[],[]); 
                        pl = pl+length(viewnum);ph=text(5,-115,150,ttl);
                        set(ph,'color','r');
                    end;
                end;
            end;
        end;
        ph = textsc([speccat{clust},', Subclust ',int2str(c),' modulators'],'title'); set(ph,'fontsize',14);
    end;
end;
ph = textsc([speccat{clust},', Subclust ',int2str(c),' modulators'],'title'); set(ph,'fontsize',14);
% plot the cluster dipole densities:*************************
for clust = 4:5 % spectral clusters
for em = 1:length(allidx2)% emotions or dims
    for c = 1:size(allidx2{em},2)
        if ~isempty(allidx2{em}{clust,c})
            for pn = 1:2 % pos/neg factors
                if pn == 1
                    onegrp = allidx2{em}{clust,c}(find(allidx2{em}{clust,c}(:,4) == 1),:);
                    %ttl = [emos{em},'; ',speccat{clust},' Cluster ',int2str(c),' Pwr incr Fac'];
                    %str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',emos{em},'POS.jpg -djpeg'];
                    ttl = ['Dim ',int2str(em),'; ',speccat{clust},' Cluster ',int2str(c),' Pwr incr Fac'];
                    str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',int2str(em),'POS.jpg -djpeg'];
                else
                    onegrp = allidx2{em}{clust,c}(find(allidx2{em}{clust,c}(:,4) == -1),:);
                    onegrp = abs(onegrp);
                    ttl = [emos{em},'; ',speccat{clust},' Cluster ',int2str(c),' Pwr dec Fac'];
                    %str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',emos{em},'NEG.jpg -djpeg'];
                    ttl = ['Dim ',int2str(em),'; ',speccat{clust},' Cluster ',int2str(c),' Pwr dec Fac'];
                    str = ['print ',fullpaths{1}(1:end-5),speccat{clust},'Cls',int2str(c),'-',int2str(em),'NEG.jpg -djpeg'];
                end;                
                if ~isempty(onegrp)
                    clear facvec comods justcomps wtsmat1 jcwts denslist                
                    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(gdcomps,fullpaths,{onegrp},allbigs,bigwts,orivec);
                    densargs ={'mrislices',[63:-6:-25],'mriview','top','geom',[4,4],'cmax',.1};
                    figure;[density] = PlotDipoles('sources.set', fullpaths,denslist{1},[],[],[],'alldistance',[],'off',densargs);
                    ph = textsc(ttl,'title'); set(ph,'fontsize',14);set(ph,'color','r')
                    eval(str);
                    close
                end;
            end;
        end;
    end;
end;
end;
% plot the associated modulations (with correct orienation):
speccat = {'Delta','Theta','Alpha','Beta','Gamma'};
load /data/common4/emotion/AllClustFacsNEW.mat    
for clust = 1:5
for c = 1:size(allidx2{1},2)
    figure; row = 4; col = 4; pl = 1;
    ttl = [speccat{clust},' Cluster, Sub-clust ',int2str(c)];
    for em = 1:length(allidx2)
        if pl > row*col
            ph = textsc(ttl,'title'); set(ph,'fontsize',14);set(ph,'color','r')
            figure; pl=1;
        end;
        plotfacsp = [];plotfacsn = [];
        if ~isempty(allidx2{em}{clust,c})
            for mem = 1:size(allidx2{em}{clust,c},1)
                currt = find(kptk(:,1) == allidx2{em}{clust,c}(mem,1)&...
                             kptk(:,2) == abs(allidx2{em}{clust,c}(mem,2))&...
                             kptk(:,3) == allidx2{em}{clust,c}(mem,3));
                [bignum  idx] = max(abs(clustfacs(currt,:)));
                 if clustfacs(currt,idx) > 0 & allidx2{em}{clust,c}(mem,4) == -1% pwr down
                    plotfacsn = [plotfacsn; clustfacs(currt,:)*-1]; 
                elseif clustfacs(currt,idx) < 0 & allidx2{em}{clust,c}(mem,4) == -1% pwr down
                    plotfacsn = [plotfacsn; clustfacs(currt,:)];                  
                elseif clustfacs(currt,idx) > 0&  allidx2{em}{clust,c}(mem,4) == 1% pwr up
                    plotfacsp = [plotfacsp; clustfacs(currt,:)];                    
                elseif clustfacs(currt,idx) < 0 & allidx2{em}{clust,c}(mem,4) == 1 % pwr up
                    plotfacsp = [plotfacsp; clustfacs(currt,:)*-1];                    
                end;
            end;
            sbplot(row,col,pl);pl = pl+1;hold on;
            if ~isempty(plotfacsp)
                [realx labelx handle] = quadplot(freqs,plotfacsp,1,[]);
                plot([get(gca,'xlim')],[0 0],'k-');
                %title([emos{em},' Pwr up']);
                title(['Dim ',int2str(em),' Pwr up']);
            end;
            sbplot(row,col,pl);pl = pl+1;hold on;
            if ~isempty(plotfacsn)
                [realx labelx handle] = quadplot(freqs,plotfacsn,1,[]);
                plot([get(gca,'xlim')],[0 0],'k-');
                %title([emos{em},' Pwr dwn']);
                title(['Dim ',int2str(em),' Pwr dwn']);
            end;
        end;
    end;
    ph = textsc(ttl,'title'); set(ph,'fontsize',14);set(ph,'color','r')
end;
end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for power assymetry between left and right frontal ics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 30; chans = [2 242];
freq = 10;
frqlim = [3 45];
datset = {'sad.set','happy.set'};
events = {{'sad',[]},{'happy',[]}};


EEG = pop_loadset('sources.set',fullpaths{nx},'all');  
[sources X Y Z XE YE ZE] = dipplot( EEG.dipfit.model , 'coordformat','spherical', 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.1/standard_BESA/avg152t1.mat', 'normlen', 'on','projlines','on');
%dipplot( EEG.dipfit.model(gdcomps{nx}) , 'coordformat','spherical', 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.1/standard_BESA/avg152t1.mat', 'normlen', 'on','gui','off','summary','on');
%  figure; dipplot( EEG.dipfit.model(gdcomps) , 'coordformat','spherical', 'mri', '/data/common/matlab/eeglab/plugins/dipfit2.1/standard_BESA/avg152t1.mat', 'normlen', 'on','gui','off','projimg','on');
for ic = 1:length(sources)
    x(1,ic) = sources(ic).talcoord(1,1);
    y(1,ic) = sources(ic).talcoord(1,2);
    rvs(1,ic) = sources(ic).rv;
end;    
Rics = find(x > 0 & y > -80 & rvs < .15);
Lics = find(x < 0 & y > -80 & rvs < .15);
%%% Plot the chosen ICs
 pop_topoplot(EEG,0, rtf , [fullpaths{nx}(end-4:end-1),' Right Hemisphere ICs'],[] ,0, 'electrodes', 'off', 'masksurf', 'on');
 pop_topoplot(EEG,0, ltf , [fullpaths{nx}(end-4:end-1),' Left Hemisphere ICs'],[] ,0, 'electrodes', 'off', 'masksurf', 'on');
 pop_topoplot(EEG,0, gdcomps, [fullpaths{nx}(end-4:end-1)],[] ,0, 'electrodes', 'off', 'masksurf', 'on');
%%%% find power for all chosen ICs
 [contrib,speccomp,cics,p,fbins] = LRicDiffs(datset,fullpaths{nx},events,Lics,Rics,frqlim,freq);
 

 pop_topoplot(EEG,0, [21,40,44,68,80,83] , [fullpaths{nx}(end-4:end-1),' Contributors to frontal 10 Hz'],[] ,0, 'electrodes', 'off', 'masksurf', 'on');
 pop_topoplot(EEG,0,combics, [fullpaths{nx}(end-4:end-1),' Contributors to frontal 10 Hz'],[] ,0, 'electrodes', 'off', 'masksurf', 'on');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjlist = [2:36]; % don't use tl81, events screwed up.
subjlist = [1:12]; % Repetitve button presses
subjlist = [21,23:26]; % only when 'feeling it' button presses (no mr72-2)
subjlist = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
subjlist = [13:20,36]; % no button press (apart from the first one)
subjlist = [1,2,4:6,8:12,14,17:21,23,25:30,31,33,34,35,36]; % all 'good' subjects (ones that said they got into it)
subjlist = [2:21,23:36];  % all but mr72-2
subjlist = [2:7,9:21,23:34,35,36];  % no mr72-2, tl81 or ar81
subjlist = [1,3:9,12,14,16,17,19,21,22,23,24,26,27,29,33,36]; % females
subjlist = [1,4:6,8,9,12,14,17,19,21,23,26,27,29,33,36]; % 'good' females
subjlist = [2,10,11,13,15,18,20,25,28,30,31,32,34,35]; % males
subjlist = [2,10,11,18,20,25,28,30,31,34]; % 'good' males

alllist{1} = [2,3,7,16,19,21,29]; % emo order 1
alllist{2} = [4,10,17,24,25,30,33]; % emo order 2
alllist{3} = [5,8,12,18,26,28,34]; % emo order 3
alllist{4} = [6,11,15,23,32,35]; % emo order 4
alllist{5} = [9,13,14,20,22,27,31]; % emo order 5
malefem{1} = [1,3:9,12,14,16,17,19,21,22,23,24,26,27,29,33]; % females
malefem{2} = [2,10,11,13,15,18,20,25,28,30,31,32,34,35]; % males
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
