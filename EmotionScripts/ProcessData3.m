% continued analysis of emotion study data-- CORRELATIONS******
% includes weights correlation and correlation with emotion
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
subjlist = [2:21,23:31,33:35];% all but 1,22,32

% for regular gdcomps decomposition:--------------------
load /data/common1/emotion/EmoWeights.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacs.mat    
%savedat = 'SpecCoMod'; fullpaths = newpaths;
savedat = 'SpecCoModMoreFreqs'; fullpaths = newpaths;
savedat = 'SpecCoModWave'; fullpaths = newpaths;
strs = {
    'load /data/common1/emotion/DeltaClust.mat',
    'load /data/common1/emotion/ThetaClust.mat',
    'load /data/common1/emotion/AlphaClust.mat',
    'load /data/common1/emotion/BetaClust.mat',
    'load /data/common1/emotion/GammaClust.mat'};
names = {'Delta','Theta','Alpha','Beta','Gamma'};




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for muscle decomposition:-------------------
load /data/common1/emotion/EmoWeightsMuscle.mat emomeans emodeciles
load /data/common1/emotion/AllClustFacsMuscle.mat    
savedat = 'SpecCoModMuscle'; fullpaths = newpaths;
strs = {'load /data/common1/emotion/GammaGdMsVf.mat'};eval(strs{1})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backproj gamma and take bi-coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=2;
ims = [2,4,7,10,54];
s = load([fullpaths{nx},savedat,'.mat']);
wts = floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);
sph = floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0);
data = floatread([fullpaths{nx},savedat,'.fdt'],[s.numrows s.numframes],[],0);
ws = wts*sph;winv = pinv(ws); acts = ws*data;
speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) inf],[],0);
specwts = speceig*winv; winv = specwts;     

backproj = winv(:,ims(im))*acts(ims(im),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate/mutual-info of IM weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allvec=cell(1,0);idxs = cell(1,0);
for s = 1:length(strs)
    eval(strs{s})
    clear facvec comods justcomps wtsmat1 jcwts denslist
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);  
    for cls = 1:2%length(facvec)
        allvec{end+1} = facvec{cls};
        idxs{end+1} = finalidx{cls};
    end;
end;
for nxx = 1:length(subjlist)
  nx= subjlist(nxx);
    clear corr bootstats
    %shuffnum = cell(1,0);
    %for f = 2:35
    %shuffnum{nx} = fullpaths{nx}; 
    %end;
    checkim = [];
    for cls = 1:length(allvec)
        checkim = [checkim allvec{cls}{nx}'];
    end;
    checkim = unique(checkim);
    if ~isempty(checkim)
        [corr, bootstats, MI, MIbootstats,varcorr,varstats] = CorrCoMod(fullpaths{nx},savedat,checkim,[],'off');
        %str =['save ',fullpaths{nx},'SpecCoModCorr.mat corr bootstats MI MIbootstats varcorr varstats allvec']; eval(str)
       %str =['save ',fullpaths{nx},'SpecCoModCorrEMOs.mat corr bootstats MI MIbootstats varcorr varstats allvec']; eval(str)
        %str =['save ',fullpaths{nx},'SpecCoModCorrMuscle.mat corr bootstats MI MIbootstats varcorr varstats allvec']; eval(str)
        %str =['save ',fullpaths{nx},'SpecCoModCorrEMOsMuscle.mat corr bootstats MI MIbootstats varcorr varstats allvec']; eval(str)
        str =['save ',fullpaths{nx},'SpecCoModCorrBrnMusc.mat corr bootstats MI MIbootstats varcorr varstats allvec']; eval(str)
    end;
    fprintf('\nSubject %s done.\n',int2str(nx));
end;
[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(corr,bootstats,subjlist,{'Brain','Muscle'},emos,allvec,'corr');

[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(corr,bootstats,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'corr');
[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(varcorr,varstats,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'varcorr');
[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(MI,MIbootstats,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'mi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = {'Theta','LowAlpha','PeakAlpha','HighAlpha','LowBeta','HighBeta','Broadband'};
%%%% pull out correlations WITHIN ICs for each cluster pair:
pl = 1; clear clustcorrs clustboots newlabels
for cls = 1:length(idxs)-1
    idx1 = idxs{cls};
    for cls2 = cls+1:length(idxs)
        idx2 = idxs{cls2};
        relcorrs = [];relboots = zeros(15,2,0);
        for nx = 1:length(gdcomps)
            try
              str =['load ',fullpaths{nx},'SpecCoModCorrEMOs.mat']; eval(str);
              ics = unique([idx1(find(idx1(:,1) == nx),3)',idx2(find(idx2(:,1) == nx),3)']);
            for ic = 1:length(ics)
                ims{1} = idx1(find(idx1(:,1) == nx & idx1(:,3) == ics(ic)),2); % keep orientation and deal in function
                ims{2} = idx2(find(idx2(:,1) == nx & idx2(:,3) == ics(ic)),2); 
                if ~isempty(ims{1}) & ~isempty(ims{2})
                    [rc,outboot,icothers] = IMcorrsbyIC(corr,bootstats,nx,ics(ic),ims,allbigs,'mi'); 
                    relcorrs = [relcorrs rc{1}];
                    relboots(:,:,end+1) = outboot{1};
                end;
            end;
            end;
        end;
        clustcorrs{pl} = relcorrs; 
        clustboots{pl} = mean(relboots,3); % mean of collected boots for clust comparison 
        newlabels{pl} = [names{cls},'/',names{cls2}];pl = pl+1;
    end;
end;
 
%%% Plot all emotions:----
figure; row = 3; col = 4; pl = 1;
cols = jet(size(clustcorrs{1},1));cols(10,:) = [1 .9 0];
for clss = 1:length(clustcorrs)
    if ~isempty(clustcorrs{clss})
        if pl > row*col
            textsc(ttl,'title');
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            figure; pl=1;
        end;
        sbplot(row,col,pl); pl = pl + 1;
        for ds = 1:size(clustcorrs{clss},1)
            ph = plot(ds,clustcorrs{clss}(ds,:),'.','markersize',4); hold on;
            set(ph,'color',cols(ds,:));                
            ph = plot(ds,mean(clustcorrs{clss}(ds,:)),'k.','markersize',10);
        end;
        ph = plot(mean(clustcorrs{clss},2),'k-');
        ph = plot(clustboots{clss},'r--');
        plot([get(gca,'xlim')],[0 0],'k-'); 
        set(gca,'xlim',[0 size(clustcorrs{clss},1)+1]);
        set(gca,'xticklabel',[]);
        title(newlabels{clss});
    end;
end;        
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
textsc(ttl,'title');
 
%%% Plot histograms with all emotions together:----
figure;  row = 3; col = 4; pl = 1;
for clss = 1:length(clustcorrs)
    if ~isempty(clustcorrs{clss})
        colmeans = mean(clustcorrs{clss},1);
        if pl > row*col
            textsc(ttl,'title');
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            figure; pl=1;
        end;
        sbplot(row,col,pl); pl = pl + 1;
        hist(colmeans,12); hold on;
        plot([0 0],[get(gca,'ylim')],'g-');
        plot([mean(colmeans) mean(colmeans)],[get(gca,'ylim')],'m-','linewidth',2);
        plot([mean(clustboots{clss}(:,1)) mean(clustboots{clss}(:,1))],[get(gca,'ylim')],'b-');
        plot([mean(clustboots{clss}(:,2)) mean(clustboots{clss}(:,2))],[get(gca,'ylim')],'r-');
        
        title(newlabels{clss});
    end;
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
textsc(ttl,'title');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find highest correlated IMs within subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjcorrs = [];keepims = []; nodata = 0; corr=[];
for nx = 1:35%length(fullpaths)
  try
    str =['load ',fullpaths{nx},'SpecCoModCorrBrnMusc.mat']; eval(str);
    %str =['load ',fullpaths{nx},'SpecCoModCorrEMOsMuscle.mat']; eval(str);
    %str =['load ',fullpaths{nx},'SpecCoModCorrMuscle.mat']; eval(str);
    %str =['load ',fullpaths{nx},'SpecCoModCorrEMOs.mat']; eval(str);
    %str =['load ',fullpaths{nx},'SpecCoModCorr.mat']; eval(str); 
    nodata = 0;  %corr=varcorr;bootstats=varstats;
  catch
    nodata = 1; corr = [];
  end;
  if nodata == 0
    for im = 1:size(corr,1)-1
      for imm = im+1:size(corr,2)
        if length(size(corr)) > 2 % split by emotion
          keepims = [keepims;[nx im imm squeeze(corr(im,imm,:))' squeeze(bootstats(im,imm,1,:))' squeeze(bootstats(im,imm,2,:))']];
        else % not split by emotion
          keepims = [keepims;[nx im imm corr(im,imm) bootstats(im,imm,1) bootstats(im,imm,2)]];
          end;
      end;
    end;
  end;
end;
% 1) find which frequencies are correlated
% 2) find where they are located.
justcomps = []; complist = []; templs = []; clustassign = zeros(size(keepims,1),3); nxlast = 0;
for k = 1:size(keepims,1)
    nx=keepims(k,1); im = keepims(k,2); imm = keepims(k,3);
    if nx > nxlast
        sets = 1;
    end;
    complist{1}{nx}{sets} = onebig{nx}{im};modwts{1}{nx}{sets} = 1;
    complist{2}{nx}{sets} = onebig{nx}{imm};modwts{2}{nx}{sets} = 1;
    justcomps{nx} = []; sets = sets +1;nxlast = nx;
    % collect templates:
    templs(1,:,k) = clustfacs(find(kptk(:,1) == nx & kptk(:,2) == im & kptk(:,3) == onebig{nx}{im}),:);
    templs(2,:,k) = clustfacs(find(kptk(:,1) == nx & kptk(:,2) == imm & kptk(:,3) == onebig{nx}{imm}),:);
    % check cluster assignment
    fullclust = 1;
    for s = 1:length(strs)
        eval(strs{s})
        for cls = 1:2%length(finalidx)
            clustnames{fullclust} = [names{cls}];
            %clustnames{fullclust} = [names{s},int2str(cls)];
            m = find(finalidx{cls}(:,1)==nx&abs(finalidx{cls}(:,2))==im&finalidx{cls}(:,3)==onebig{nx}{im});
            mm = find(finalidx{cls}(:,1)==nx&abs(finalidx{cls}(:,2))==imm&finalidx{cls}(:,3)==onebig{nx}{imm});
            if ~isempty(m)
                clustassign(k,1) = [fullclust];
            end;
            if ~isempty(mm)
                clustassign(k,2) = [fullclust];
                clustassign(k,3) = keepims(k,4);% keep correlation too
            end;
            if clustassign(k,1) ~=0 & clustassign(k,2) ~= 0
                clustpairs{k} = [clustnames{clustassign(k,1)},'-',clustnames{clustassign(k,2)}];
            end;
            fullclust = fullclust + 1;
        end;
    end;
    if k==1000 | k==2000|k==3000|k==4000|k==5000|k==6000|k==7000|k==8000|k==9000|k==10000
        fprintf('\nSamples 1 through %s completed.',int2str(k));
    end;
end;
save /data/common1/emotion/clustcorrvalsBrnMusc.mat clustassign keepims names
%save /data/common1/emotion/clustcorrvals.mat clustassign keepims
%save /data/common1/emotion/clustcorrvalsEMOs.mat clustassign keepims
%save /data/common1/emotion/clustcorrvalsVAR.mat clustassign keepims
%save /data/common1/emotion/clustcorrvalsEMOsVAR.mat clustassign keepims
%save /data/common1/emotion/clustcorrvalsMI.mat clustassign keepims
%save /data/common1/emotion/clustcorrvalsEMOsMI.mat clustassign keepims

% find cluster pairs when NOT separated by emotions:
fullclust = 1;
for s = 1:length(strs)
    eval(strs{s})
    for cls = 1:length(finalidx)
        clustnames{fullclust} = [names{s},int2str(cls)];
        fullclust = fullclust + 1;
    end;
end;
% or:
%clustnames = {'Brain','Muscle','OMT'}; % has to be the order of 'strs' and sub clusts

pairs = 1; clear pairnames pairboots paircorrs ncorrs  
for s = max(clustassign(:)):-1:1
    for ss = s:-1:1
        %ss=s;
        onegrp1 = find(clustassign(:,1)==s&clustassign(:,2)==ss)';
        onegrp2 = find(clustassign(:,2)==s&clustassign(:,1)==ss)';
        ncorrs(pairs) = length(onegrp1) + length(onegrp2);
        pairnames{pairs} = [clustnames{s},'-',clustnames{ss}];
        paircorrs{pairs} = keepims([onegrp1,onegrp2],4);
        pairboots{pairs} = [keepims([onegrp1,onegrp2],5),keepims([onegrp1,onegrp2],6)];
        meds(s,ss) = median(paircorrs{pairs});
        mns(s,ss) = mean(paircorrs{pairs});
        iqrs(s,ss) = iqr(paircorrs{pairs});
        pairs = pairs + 1;
    end;
end;
sigonly = 'off';     
markpairs = {'-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.'};
%markcols = hsv(length(paircorrs));% hsv for muscle
markcols = jet(length(paircorrs)+7);markcols(7:13,:)=[];
markcols = [[0 0 1];[0 1 0];[1 0 0]];markpairs = {'-','-','-'};
markcols = lines(length(paircorrs));markpairs = {'-.','--','-','-.','--','-','-.','--','-','-.','--','-','-.','--','-','-.','--','-',}; 
figure; newpairnames = cell(1,0);
for cp = 1:length(paircorrs) % 12 for gamma only
    if cp > 1 & ~strcmp(pairnames{cp}(1:4),pairnames{cp-1}(1:4))
        plot([get(gca,'xlim')],[0 0],'k-'); legend(newpairnames,'location','northwest');   
        xlabel('% of IM pairs'); ylabel('Correlation coefficient (r)');
        title('Correlations using all time windows');
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        figure;newpairnames = cell(1,0);
    end;
    if strcmp(sigonly,'on')
        sigs{cp} = find(paircorrs{cp} < pairboots{cp}(:,1) | paircorrs{cp} > pairboots{cp}(:,2));
        nsigs(cp) = length(sigs{cp});
        plotpairs=paircorrs{cp}(sigs{cp}); % only sig paircorrs
        %[hh bins] = hist(paircorrs{cp}(sigs{cp}),10);
    else
        plotpairs=paircorrs{cp}; % all paircorrs
        %[hh bins] = hist(paircorrs{cp},15);
        %[hh bins] = hist(paircorrs{cp},bins);
    end;
    if length(plotpairs) > 20
        [ploty idxy] = sort(plotpairs);
        plotx = [1:length(plotpairs)]; 
        ph = plot(plotx/max(plotx)*100,ploty,markpairs{cp},'linewidth',2); hold on;
        %hh = hh/ncorrs(cp); clear newh
        %for b = 1:length(hh)
        %    newh(b) = sum(hh(b:end))*100;
        %end;
        %saveh(cp,:) = newh;saveb(cp,:) = bins;
        %ph = plot(bins,newh,markpairs{cp},'linewidth',2); hold on;        
        set(ph,'color',markcols(cp,:));
        newpairnames{end+1} = pairnames{cp};
    end;
end;
legend(newpairnames,'location','northwest');
plot([get(gca,'xlim')],[0 0],'k-');
xlabel('% of IM pairs'); ylabel('Correlation coefficient (r)');
title('Correlations between broadband IMs using all time windows');
%xlabel('Cumulative Correlation'); ylabel('% of ''significant'' pairings');title('Correlations using all time windows');
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/Brain-MuscIMcorrs.jpg -djpeg
print /home/julie/Manuscripts/Gamma/IMcorrs.jpg -djpeg
print /home/julie/Manuscripts/Gamma/IMcorrs.eps -depsc -adobe -painters


% find cluster pairs when separated by emotions:
figure; row = 4; col = 4;
onepair = [];
for e = 1:15 % emotions
    pairs = 1; clear pairnames pairboots paircorrs ncorrs  
    for s = max(clustassign(:)):-1:2
        for ss = s-1:-1:1
            onegrp1 = find(clustassign(:,1)==s&clustassign(:,2)==ss)';
            onegrp2 = find(clustassign(:,2)==s&clustassign(:,1)==ss)';
            ncorrs(pairs) = length(onegrp1) + length(onegrp2);
            pairnames{pairs} = [clustnames{s},'-',clustnames{ss}];
            paircorrs{pairs} = keepims([onegrp1,onegrp2],e+3);
            onepair{pairs}(:,e) = paircorrs{pairs}; 
            if size(keepims,2) > 6
                pairboots{pairs} = [keepims([onegrp1,onegrp2],3+15+e),keepims([onegrp1,onegrp2],3+30+e)];
            else
                pairboots{pairs} = [keepims([onegrp1,onegrp2],3+1+e),keepims([onegrp1,onegrp2],3+2+e)];
            end;
            pairs = pairs + 1;
        end;
    end;
    sbplot(row,col,e)
    %figure; 
    sigonly = 'on';
    markpairs = {'-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.','--','-',':','-.'};
    %markcols = hsv(length(paircorrs));% hsv for muscle
    markcols = jet(length(paircorrs)+7);% 
    markcols(7:13,:)=[];
    for cp = 1:length(paircorrs)
        if strcmp(sigonly,'on')
            sigs{cp} = find(paircorrs{cp} < pairboots{cp}(:,1) | paircorrs{cp} > pairboots{cp}(:,2));
            nsigs(cp) = length(sigs{cp});
            [hh bins] = hist(paircorrs{cp}(sigs{cp}),10);
        else
            [hh bins] = hist(paircorrs{cp},20);
            %[hh bins] = hist(paircorrs{cp},bins);
        end;
        hh = hh/ncorrs(cp); clear newh
        for b = 1:length(hh)
            newh(b) = sum(hh(b:end))*100;
        end;
        %saveh(cp,:) = newh;
        %saveb(cp,:) = bins;
        ph = plot(bins,newh,markpairs{cp},'linewidth',2); hold on;

        %ph = plot(bins,hh,markpairs{cp},'linewidth',2); hold on;
        set(ph,'color',markcols(cp,:));
    end;
    set(gca,'ylim',[-1 101]); 
    plot([0 0],[get(gca,'ylim')],'k-','linewidth',2);
    plot([get(gca,'xlim')],[50 50],'k--');
    legend(pairnames);xlabel('Cumulative Correlation'); ylabel('% of possible pairings');title('Correlations using all time windows');
    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    title(emo2{e});
end;
print /home/julie/Meetings/Cosyne2008/ClustCorrs.jpg -djpeg
print /home/julie/Meetings/Cosyne2008/ClustCorrsMI.jpg -djpeg
print /home/julie/Meetings/Cosyne2008/ClustCorrsVAR.jpg -djpeg
print /home/julie/Meetings/Cosyne2008/ClustCorrsbyEMOsVAR.jpg -djpeg
print /home/julie/Meetings/Cosyne2008/ClustCorrsMuscle.jpg -djpeg

print /home/julie/Meetings/Cosyne2008/ClustCorrs.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/ClustCorrsMI.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/ClustCorrsVAR.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/ClustCorrsbyEMOsVAR.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/ClustCorrsMuscle.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/ClustCorrsbyEMOs.eps -depsc -painters -adobe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT correlation means across emotions and find sig diffs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; fg = figure(2); row=  5; col = 5; pois= [];
for pairs = 1:length(onepair)
    if ~isempty(onepair{pairs}) & size(onepair{pairs},1) > 5
    [P,atab,STATS] = anova1(onepair{pairs},emos,'off');
    figure(1); comp = multcompare(STATS,'alpha',.001,'ctype','bonferroni'); %close;
    comppairs = [];
    for cp = 1:size(comp,1)
        if length(find(comp(cp,[3,5]) == abs(comp(cp,[3,5])))) ~= 1 %(not straddling 0=sig)
            comppairs = [comppairs;[comp(cp,[1:2])]];
        end;
    end;
    sigemos = unique(comppairs);% allemos with sig diffs
    figure(fg); sbplot(row,col,pairs);
    ph = plot(mean(onepair{pairs},1),'k-');
    for e = 1:size(onepair{pairs},2)
        ph = plot(e,onepair{pairs}(:,e),'.','markersize',4); hold on;
        set(ph,'color',cols(e,:)); 
        ph = plot(e,mean(onepair{pairs}(:,e)),'.','markersize',14);
        set(ph,'color',cols(e,:)); 
        if ismember(e,sigemos)
            ph = plot(e,1,'*'); set(ph,'color',[.5 0 .8]);
            pois = [pois pairs];
        end;   
    end;
    yl = get(gca,'ylim'); set(gca,'ylim',[yl(1) 1.2]);
    set(gca,'xlim',[0 size(onepair{pairs},2)+1]);
    plot([get(gca,'xlim')],[0 0],'k-');set(gca,'xticklabel',[]);
    ph = title(pairnames{pairs});set(ph,'fontsize',8);
    end;
end; pois = unique(pois);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
textsc('Correlation coefficients betwee spectral clusters and across emotions; purple * is sig by ANOVA at p < .001','title');
print /home/julie/Meetings/Cosyne2008/CorrCoeffdiffsbyEMO.eps -depsc -painters -adobe
print /home/julie/Meetings/Cosyne2008/CorrCoeffdiffsbyEMO.jpg -djpeg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using correlations separated by emotion:
% now plot only the sig pairs up close
for prs = 1:length(pois)
    pairs = pois(prs);
    [P,atab,STATS] = anova1(onepair{pairs},emos,'off');
    figure(1); comp = multcompare(STATS,'alpha',.001,'ctype','bonferroni'); %close;
    comppairs = [];
    for cp = 1:size(comp,1)
        if length(find(comp(cp,[3,5]) == abs(comp(cp,[3,5])))) ~= 1 %(not straddling 0=sig)
            comppairs = [comppairs;[comp(cp,[1:2])]];
        end;
    end;
    sigemos = unique(comppairs);% allemos with sig diffs
    figure; set(gca,'fontsize',14);
    ph = plot(mean(onepair{pairs},1),'k-');hold on;  
    yl = get(gca,'ylim');set(gca,'ylim',[yl(1) yl(2)+abs(yl(2)*.15)]);
    yl = get(gca,'ylim'); 
    for e = 1:size(onepair{pairs},2)
        ph = plot(e,mean(onepair{pairs}(:,e)),'.','markersize',24);
        set(ph,'color',cols(e,:)); 
        if ismember(e,sigemos)
            ph = plot(e,yl(2)-abs(yl(2)*.05),'*'); set(ph,'color',[.5 0 .8]);
            set(ph,'markersize',14);
        end;   
    end;
    set(gca,'xlim',[0 size(onepair{pairs},2)+1]);
    plot([get(gca,'xlim')],[0 0],'k-');set(gca,'xticklabel',[]);
    ph = title(pairnames{pairs});set(ph,'fontsize',14);
    xlabel('Emotions'); ylabel('Mean Correlation');
    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    str = ['print /home/julie/Meetings/Cosyne2008/CorrCoeffdiffs',pairnames{pairs},'.eps -depsc -painters -adobe']; eval(str)
end;

% Plot gamma-related vs all else
bins = [-.6:.06:.9];gcorrs = [];rcorrs=[];
for pair = 1:6
    gcorrs = [gcorrs,onepair{pair}'];
end;
for pair = 7:length(onepair)
    rcorrs = [rcorrs,onepair{pair}'];
end;
[gh bins] = hist(gcorrs,bins);
[rh bins] = hist(rcorrs,bins);
gh = gh/length(gcorrs); clear newg
rh = rh/length(rcorrs); clear newr
for b = 1:length(gh)
    newg(b) = sum(gh(b:end))*100;
end;
for b = 1:length(rh)
    newr(b) = sum(rh(b:end))*100;
end;

% calculate bootstrap background
clear boot1 boot2 H P
allcorrs = [gcorrs,rcorrs];
for bt = 1:500
    randcorrs = allcorrs(randperm(length(allcorrs)));
    boot1(bt,:) = randcorrs(1:length(gcorrs));
    boot2(bt,:) = randcorrs(length(gcorrs)+1:end);
    [H(bt) P(bt)] = ttest2(boot1(bt,:),boot2(bt,:),.01);
end;
[reahH realP] = ttest2(gcorrs,rcorrs,.0000000001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using full time course for correlation:
% Plot gamma-related vs all else
bins = [-.6:.06:.9];gcorrs = [];rcorrs=[];
for pair = 2:12 % all gamma-related, not gamma-gamma
    if pair > 1 & pair < 13
    gcorrs = [gcorrs,paircorrs{pair}'];
    else
    end;
       rcorrs = [rcorrs,paircorrs{pair}'];
end;
[H P] = ttest(gcorrs); % is the gamma-other correlation <0?

[gh bins] = hist(gcorrs,bins);
gh = gh/length(gcorrs); clear newg

for b = 1:length(gh)
    newg(b) = sum(gh(b:end))*100;
end;

clear boot1 boot2 H P
allcorrs = [gcorrs,rcorrs];
for bt = 1:500
    randcorrs = allcorrs(randperm(length(allcorrs)));
    boot1(bt,:) = randcorrs(1:length(gcorrs));
    boot2(bt,:) = randcorrs(length(gcorrs)+1:end);
    [H(bt) P(bt)] = ttest2(boot1(bt,:),boot2(bt,:),.01);
end;
[reahH realP] = ttest2(gcorrs,rcorrs,.0000000001);


figure; 
for bt = 1:500
    [bt1 bins] = hist(boot1(bt,:),bins);
    [bt2 bins] = hist(boot2(bt,:),bins);
    bt1 = bt1/size(boot1,2);
    bt2 = bt2/size(boot2,2); 
    for b = 1:length(bt1)
        newbt1(b) = sum(bt1(b:end))*100;
    end;
    for b = 1:length(bt2)
        newbt2(b) = sum(bt2(b:end))*100;
    end;
    ph = plot(bins,newbt1,'b-','linewidth',1); hold on;
    set(ph,'color',[.6 .6 1]);
    ph = plot(bins,newbt2,'r-','linewidth',1); hold on;
    set(ph,'color',[1 .6 .6]);
end;

ph = plot(bins,newg,'b-','linewidth',4); hold on;
ph = plot(bins,newr,'r-','linewidth',4); hold on;
set(gca,'xlim',[-.6 1]); set(gca,'ylim',[-1 101]);
plot([0 0],[get(gca,'ylim')],'k-','linewidth',2);
plot([get(gca,'xlim')],[50 50],'k--');

print /home/julie/Meetings/Cosyne2008/MeanCorrs.eps -depsc -painters -adobe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find stats for all paircorrs

for cp = 1:length(paircorrs)
    sigs{cp} = find(paircorrs{cp} < pairboots{cp}(:,1) | paircorrs{cp} > pairboots{cp}(:,2));
    nsigs(cp) = length(sigs{cp});
    percsig(1,cp) = (nsigs(cp)/length(paircorrs{cp}))*100;
   % percsig(2,cp) = 100*mean(abs(paircorrs{cp}));
   percsig(2,cp) = 100*mean(abs(paircorrs{cp}(sigs{cp})));
end;
figure; bar(percsig','stacked'); hold on; colormap('winter')
set(gca,'xlim',[0 78]); set(gca,'xtick',[1:3:77]);set(gca,'xgrid','on');
 textsc('Emotion IM cluster time course correlations; blue=% sig corrs; green=mean(abs(sigcorrs));','title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

>> figure; ph=bar(percsig(:,1:12)'); hold on; colormap('winter') % just bb vs all others    
>> set(ph,'barwidth',1.25)
>> set(gca,'xlim',[0 13]);
>> set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
>> text(.8,-25,'BB-BB','rotation',90);  
>> text(1.8,-25,'BB-Beta4','rotation',90);
>> text(2.8,-25,'BB-Beta3','rotation',90);
>> text(3.8,-25,'BB-Beta2','rotation',90);
>> text(4.8,-25,'BB-Beta1','rotation',90);
>> text(5.8,-25,'BB-Alpha3','rotation',90);
>> text(6.8,-25,'BB-Alpha2','rotation',90);
>> text(7.8,-25,'BB-Alpha1','rotation',90);
>> text(8.8,-25,'BB-Theta3','rotation',90);
>> text(9.8,-25,'BB-Theta2','rotation',90);
>> text(10.8,-25,'BB-Theta1','rotation',90);
>> text(11.8,-25,'BB-Delta','rotation',90); 
>> print fig1.eps -depsc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try an 'emotion space' with pairing correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spcorrs = []; % use the 'EMOs' correlations
for pair = 1:6%length(onepair)
    spcorrs = [spcorrs,onepair{pair}'];
end;

w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
erating = w.emoactiv; % w.emoval or w.emoactiv
signames = cell(1,0); clear cr
for pair = 1:length(onepair)
    
    dd = pdist(onepair{pair}', 'euclidean') ;    
    [fullmd{pair},fullwts] = cmdscale(dd);

    figure; row = round(sqrt(size(fullmd{pair},2))); col = ceil(sqrt(size(fullmd{pair},2)));
    for fd = 1:size(fullmd{pair},2)
        sbplot(row,col,fd)
        for e = 1:size(fullmd{pair},1)
            ph=plot(e,fullmd{pair}(e,fd),'.');hold on;
            set(ph,'markersize',15);set(ph,'color',cols(e,:));
            ph = text(e,fullmd{pair}(e,fd),emo2{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
        end;  
        ph = plot([get(gca,'xlim')],[0 0],'k-'); hold on;
        set(gca,'xticklabel',[]);
        title(['Dim ',int2str(fd)]);
    end;
    textsc(pairnames{pair},'title');
    %set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    %str = ['print /home/julie/Manuscripts/Gamma/.eps -depsc -painters -adobe']; eval(str)
    figure; row = round(sqrt(size(fullmd{pair},2))); col = ceil(sqrt(size(fullmd{pair},2)));
    for dim = 1:size(fullmd{pair},2)
        sbplot(row,col,dim)
        for e = 1:size(fullmd{pair},1)
            ph=plot(erating(e),fullmd{pair}(e,dim),'.');hold on;
            set(ph,'markersize',20);set(ph,'color',cols(e,:));
            ph = text(erating(e),fullmd{pair}(e,dim),emo2{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        end;
        plot([5 5],[get(gca,'ylim')],'k-');
        plot([get(gca,'xlim')],[0 0],'k-');
        title(['Correlation ',num2str(cr(pair,dim))]);
        %xlabel('Subjective emotion ratings');
        %ylabel([' Dim ',int2str(dim),' Template']);
        cr(pair,dim) = corr2(erating,fullmd{pair}(:,dim)');
    end;  
    %textsc(pairnames{pair},'title');
    sigdims{pair} = find(cr(pair,:) > .55 | cr(pair,:) < -.55);
    actualsig{pair} = cr(pair,find(cr(pair,:) > .55 | cr(pair,:) < -.55));
    if ~isempty(sigdims{pair})
        signames{end+1} = pairnames{pair};
    end;
end;
% plot just the interesting dimensions:---------
figure; row =  round(sqrt(length(signames))); col = ceil(sqrt(length(signames)));pl = 1;
for pair = 1:length(sigdims)
    if ~isempty(sigdims{pair})
        for dm = 1:length(sigdims{pair})
            dim = sigdims{pair}(dm);            
            sbplot(row,col,pl); pl = pl+1;
            for e = 1:size(fullmd{pair},1)
                ph=plot(erating(e),fullmd{pair}(e,dim),'.');hold on;
                set(ph,'markersize',20);set(ph,'color',cols(e,:));
                ph = text(erating(e),fullmd{pair}(e,dim),emo2{e});
                set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
            end;
            plot([5 5],[get(gca,'ylim')],'k-');
            plot([get(gca,'xlim')],[0 0],'k-');
            title([pairnames{pair},', Dim ',int2str(dim),', Corr ',num2str(actualsig{pair})]);
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure; row = 2; col = 2; place = 1;viewnum = [1,2,3];zoom = 1;
%[angs,dists,newdips,realdips,countbil] = PlotCoModasDipoles(complist,justcomps,fullpaths,'sources.set',row,col,place,zoom,1,viewnum,modwts,[],1,[]);

figure; pl = 1; row = 8; col = 8;
for k = 1:size(templs,3)
    if pl > row * col
        figure; pl = 1;
    end;
    sbplot(row,col,pl); pl = pl+1;
    ph = logplot(origfreqs,templs(1,:,k),2,'r'); hold on;
    ph = logplot(origfreqs,templs(2,:,k),2,'b'); hold on;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate emo ratings (behavioral) with mean wts on gamma IMs
% this is basically unremarkable on the whole. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corrwith = 'Valence';%  'Valence' or 'Arousal'
clustnames = {'Delta','Theta1','Theta2','Theta3','AlphaLow','AlphaPeak','AlphaHigh','Beta1','Beta2','Beta3','Beta4','Gamma'};
pl = 1; clear corrs clustidxs keeptrack
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
        %[corrs{pl},clustidxs{pl},keeptrack{pl}] = EmoCorrs(emomeans,oricell,subjlist,useimssub,0,'elim',[2,3:5,6,8:12,15],corrwith);
        [corrs{pl},clustidxs{pl},keeptrack{pl}] = EmoCorrs(emomeans,oricell,subjlist,useimssub,0,'elim',[],corrwith);
        pl = pl+1;
    end;
    fprintf('.');
end;
% anger/happy
% anger/sad/grief vs happy/joy/love
cut = []; % * std
w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
clustinfo = 'load /data/common1/emotion/AllClustFacs.mat';    
realclusts = [1,2,2,2,3,3,3,4,4,4,4,5]; subclusts = [1,1,2,3,1,2,3,1,2,3,4,1]; 
dens = 'on'; % 'on' or [] to plot color-coded dipoles
for decomp = length(corrs):-1:1
    specclust = strs{realclusts(decomp)};
    subclust = subclusts(decomp);
    sqcorrs = corrs{decomp}.^2; sqcorrs=sqcorrs.*(corrs{decomp}./abs(corrs{decomp}));
        EmoDimDipoles(subjlist,fullpaths,sqcorrs,w.emoval,keeptrack{decomp},pics{decomp},cut,specclust,subclust,clustinfo,dens);        
    %ph=textsc([clustnames{decomp},' density wtd by correlation with ',corrwith,' -- all emos'],'title'); set(ph,'color','r');set(ph,'fontsize',14);
    ph=textsc([clustnames{decomp},' dipole density weighted by correlation with ',corrwith,' -- Anger/Sad/Happy/Joy'],'title'); set(ph,'color','r');set(ph,'fontsize',14);
 set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
    str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},'WtdDens',corrwith,'.tif -dtiff']; eval(str);
end;
str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},corrwith,'CorrWtd1.tif -dtiff']; eval(str);% all but compassion
str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},corrwith,'CorrWtd2.tif -dtiff']; eval(str);% all
str = ['print ',fullpaths{1}(1:end-5),'results/',clustnames{decomp},corrwith,'CorrWtd3.tif -dtiff']; eval(str);%Anger/Sad/Happy/Joy

figure;  % plot the corrletion distributions for all clusters
tmpcols = jet(length(corrs));
for x =1:12
    ph = plot([1:length(corrs{x})]/length(corrs{x}),sort(corrs{x}),'linewidth',2); 
    hold on;    set(ph,'color',tmpcols(x,:));
end; plot([get(gca,'xlim')],[0 0],'k-');
legend(clustnames,'location','northwest')
title(['Correlations with Emotion ',corrwith]);
 set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
str = ['print /home/julie/fig.eps -depsc']; eval(str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
regvec = repmat(zscore(w.emoval),[size(emomeans{1},2)]);; % emoval or emoactiv
clear corrwts
pval = .05;
for nx = 1:length(allbigs)
    if ~isempty(emomeans{nx})
        clear corrs
        for im = 1:length(allbigs{nx})
            emowts = zscore(emomeans{nx}(im,:));
            %[corrs(1,im),indx,indy,crs] = matcorr(regvec,emowts);
            [corrs(1,im)] = corr2(regvec,emowts);
            [b,BINT,R,RINT,STATS]  = regress(regvec',[ones(length(emowts),1) emowts']);
            %if STATS(3) > pval
            if abs(corrs(1,im)) < .3 % less than 50% correlation
                corrs(1,im) = 0;
            end;
            corrwts{nx}{im} = repmat(corrs(1,im),[1 length(allbigs{nx}{im})]);
        end;
    end;
end;
% plot:-----------------------
% maybe 6 sig corrs between valence and brain gamma, 8-10 with muscle gamma
specclust = 1; % referring to 'strs'
eval(strs{specclust})
[facvec,comods,wtsmat1,justcomps,jcwts,denslist,denswts] = Var4DipPlot(finalidx,allbigs,corrwts,orivec);%
figure; row = 2; col = 2;pl=1; 
for cls = 1:2%length(denslist) % subcluster
    sbplot(row,col,pl);[density,minmask,maxmask] = PlotDipoles('sources.set', newpaths,denslist{cls},denswts{cls},{'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat','spherical'},[],[],[],[],[],1);pl=pl+1;
    sbplot(row,col,pl);[density,minmask,maxmask] = PlotDipoles('sources.set', newpaths,denslist{cls},denswts{cls},{'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat','spherical','view',[1 0 0]},[],[],[],[],[],1);pl=pl+1;
end;
str =['print /home/julie/Manuscripts/Gamma/ValenceCorrelation.jpg -djpeg']; eval(str)
str =['print /home/julie/Manuscripts/Gamma/ArousalCorrelation.jpg -djpeg']; eval(str)

ph=textsc(['BRAIN Correlation with Valence (p < .05)'],'title');set(ph,'color','r');
ph=textsc(['MUSCLE Correlation with Valence (p < .05)'],'title');set(ph,'color','r');
ph=textsc(['BRAIN Correlation with Arousal (p < .05)'],'title');set(ph,'color','r');
ph=textsc(['MUSCLE Correlation with Arousal (p < .05)'],'title');set(ph,'color','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
