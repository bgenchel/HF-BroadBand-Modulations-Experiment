%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
%savedat = 'SpecCoModMoreFreqs'; fullpaths = newpaths;
savedat = 'SpecCoModNoOvrlap'; fullpaths = newpaths;
%savedat = 'SpecCoModNoFilt';  newpaths{22} = [newpaths{22}(1:22),'mr74/'];
%newdir = ['/data/projects/julieo/emotion/'];
%fullpaths = cell(1,length(newpaths));
%for nx =1:length(newpaths)
%  fullpaths{nx} = [newdir,newpaths{nx}(end-4:end)];
%end;

%savedat = 'SpecCoModMuscle';fullpaths = newpaths; % better classification
%savedat = 'SpecCoModWave';fullpaths = newpaths;
%savedat = 'SpecCoModWavePrePost';fullpaths = newpaths;
cols = jet(15); cols= [cols;[.5 .5 .5]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVM method
The function 'svmpredict' has three outputs. The first one,
predictd_label, is a vector of predicted labels. The second output,
accuracy, is a vector including accuracy (for classification), mean
squared error, and squared correlation coefficient (for regression).
The third is a matrix containing decision values or probability
estimates (if '-b 1' is specified). If k is the number of classes,
for decision values, each row includes results of predicting
k(k-1)/2 binary-class SVMs. For probabilities, each row contains k values
indicating the probability that the testing instance is in each class.
Note that the order of classes here is the same as 'Label' field
in the model structure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YP classification results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basedir = '/data/projects/YP/ResultsEmoWeightsNoFilt/';
basedir = '/data/projects/YP/ResultsEmoWeights1s/';
%basedir = 'C:\Users\julie\Documents\MatlabData\emotion\Classification\YP\EmoClass_1sec\';
%basedir = 'C:\Users\julie\Documents\MatlabData\emotion\Classification\YP\ResultsEmoWeightsNoFilt\';
%-----------------------------------
str = ['load ',basedir,'ConfusionTable.mat'];eval(str)
% ConfusionTable{subj,All/Val/Aro}.Strategy{Fscore/RandSel,1}{IM}(16,16)
%-----------------------------------
str = ['load ',basedir,'AccuracyMatrix.mat'];eval(str)
% TotalAccuracy{subj,All/Val/Aro}.Strategy{Fscore/RandSel,1}.All(IM,1)
% TotalAccuracy{subj,All/Val/Aro}.Strategy{Fscore/RandSel,1}.Single(IM,10,10)
%-----------------------------------
strs = {['load ',basedir,'Top_All_Fscore.mat'],['load ',basedir,'Top_Val_Fscore.mat'],['load ',basedir,'Top_Aro_Fscore.mat']};eval(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process confusion matrices
classes = {'All','Valence','Arousal'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max Acc/nfeat by subj graph
subjlist = [1:4,6:14,16:21,23:35];% elim subjs 5,15,22 is a repeat

%for classtype = 3:-1:1 % 1=all, 2=Val, 3=Aro
clear eacc subjacc meannumsampls
for sbj  =1:length(subjlist)
  eval(strs{1}) % load correct Table with Fscore order
  subj = subjlist(sbj);
  if ~isempty(TotalAccuracy{subj,1})
    r=find(Table(:,subj));r=r(end);
    q = ConfusionTable{subj,1}.Strategy{1,1}{r};% Fscore order, All
    for e = 1:size(q,1)-1
      eacc(subj,e) = 100* (q(e,e)/q(e,end));
    end;
    meannumsampls(1,subj) = mean(q(:,end));
    subjacc(1,subj) = TotalAccuracy{subj,1}.Strategy{1,1}.All(r,1); % only Fscore
    subjacc(2,subj) = r; % number of IMs for 'all'
    eval(strs{3})
    r=find(Table(:,subj));r=r(end); % 
    subjacc(3,subj) = TotalAccuracy{subj,3}.Strategy{1,1}.All(r,1); % Arousal
    subjacc(4,subj) = r; % number of IMs for 'Aro'
   eval(strs{2})
    r=find(Table(:,subj));r=r(end); 
    subjacc(5,subj) = TotalAccuracy{subj,2}.Strategy{1,1}.All(r,1); % Valence
    subjacc(6,subj) = r; % number of IMs for 'Val'
  end;
end;
mean(subjacc(1,find(subjacc(1,:))))
medacc = median(eacc,2);
%[aa bb] = sort(subjacc(1,:));% accuracy 'overall'
[a b] = sort(medacc); % mean or median accuracy across emotions
nonz = find(medacc(b)>0); b = b(nonz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1***************************************
figure; boxplot(eacc(b,:)', 'plotstyle' ,  'compact','colors','g');
set(gca,'xlim',[.5 35.5]);xlabel('Sorted subjects'); ylabel('Accuracy (%)'); hold on;
ph=plot([1:length(b)],subjacc(2,b),'bx','markersize',12); set(gca,'xlim',[.5 35.5]);
plot([1:length(b)],subjacc(1,b),'r^','linewidth',2,'markersize',10); hold on; % all
plot([1:length(b)],subjacc(3,b),'md','linewidth',2,'markersize',10); hold on; % arousal
plot([1:length(b)],subjacc(5,b),'cs','linewidth',2,'markersize',10); hold on;% valence
%plot([1:length(b)],meannumsampls(b),'cx','linewidth',2,'markersize',10); hold on;
xlabel('Sorted subjects');set(gca,'xlim',[.5 length(b)+.5]);set(gca,'ylim',[-1 101]);
textsc(classes{1},'title');
str = ['print ',basedir,'Accuracy',classes{1},'Medians.jpg -djpeg']; eval(str)
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
  str = ['print ',basedir,'Accuracy',classes{classtype},'Medians.eps -depsc -adobe -painters']; eval(str)
  %end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% include subjective ratings of emotion intensity, etc
% FIGURE ? (separate emos)***************************************
figure; boxplot(eacc,'colors',jet(15));
str = ['print C:\Users\julie\Documents\Manuscripts\EmotionClassification\Accuracy',classes{classtype},'byEmos.eps -depsc -adobe -painters']; eval(str)
% no statistical diffs by anova1

% FIGURE ?  -- confusion matrix
  eval(strs{classtype}) % load correct Table with Fscore order
figure; row=  6; col = 6;
confmean = zeros(15,15);
for subj  =1:35
      if ~isempty(TotalAccuracy{subj,1})
   x=TotalAccuracy{subj,classtype}.Strategy{ForRand,1}.Single;
   z = reshape(x(end,:,:),[1 size(x,2)*size(x,3)]);
   for im = 1:size(x,1)
      y = reshape(x(im,:,:),[1 size(x,2)*size(x,3)]);
      [h p(1,im)] = ttest2(y,z);
   end;
   r = find(p(2:end) > .01); r=r(1); % last consec under .01
    confsum = zeros(15,15);
   for e = 1:15
      allaccs(e,:) = 100* ConfusionTable{subj,classtype}.Strategy{ForRand,1}{r}(e,1:15)/ConfusionTable{subj,classtype}.Strategy{ForRand,1}{r}(e,end);
      corrmatch = allaccs(e,e); 
      allaccs(e,e) = 0; % zero out correct matches
      [v i]= max(allaccs(e,:));
      confmean(e,i) = confmean(e,i) + allaccs(e,i);
      if corrmatch < 50 % if correct < 50% of time
         confsum(e,i) = 1; % 1 = highest confusion
      end;
   end;
   subjconf(:,:,subj) = confsum;
   %sbplot(row,col,subj); imagesc(allaccs);cbar;   
end;
end;
confmean = confmean./35; % take mean
figure; imagesc(confmean); cbar;
figure; x=sum(subjconf,3); x(find(x<3)) = 0; imagesc(x); cbar;
% Emotions that were < 50% accurately classified within subject
% were only systematically confused for another emotion max 5 times
% a handful were mistook 3 times and 9 were confused 4 times.

%***************************************
% FIGURE 3  -- IM spectra/locations
%**********************************
s = load([fullpaths{1},savedat,'.mat']);
clear data activations winv alldip kptk allbigs bigwts clustfacs mnspecs orivec
clustfacs = [];  mnspecs = []; kptk = [];
%subjlist = [find(subjacc(1,:))];% elim subjs 5,15
%subjlist = [1:4,6:14,16:21,23:35];% elim subjs 5,15
freqscale = 'quad';
percmax =.5; % percent of max template to take as a 'comod'
 [clustfacs,templcell,mnspecs,kptk,allbigs,bigwts,onebig,orivec,modcorr,freqs,outcomods,rawrms,rawcorr] = CollectCoModTempls(savedat,fullpaths,subjlist,percmax);
[deltaclust, lthetaclust,hthetaclust,alphaclust,lbetaclust,hbetaclust,gamaclust,freqs,otherclust] = SortModTempls(clustfacs,kptk,mnspecs,s.freqs,[0 129],[],'quad',savedat,2.75);
fulltempls = {deltaclust{1},lthetaclust{1},hthetaclust{1},alphaclust{1},lbetaclust{1},hbetaclust{1},gamaclust{1},otherclust{1}}; 
fullmeans = {deltaclust{2},lthetaclust{2},hthetaclust{2},alphaclust{2},lbetaclust{2},hbetaclust{2},gamaclust{2},otherclust{2}}; 
fullidxs = {deltaclust{3},lthetaclust{3},hthetaclust{3},alphaclust{3},lbetaclust{3},hbetaclust{3},gamaclust{3},otherclust{3}}; 
%**********************************
load /data/common1/emotion/Classification/YuanPinEmoWeights1sec.mat eweights comment
str = ['load /data/common1/emotion/Classification/InitialParams1sec.mat AllCMD'];eval(str)
for nx = 1:length(eweights)    
    clear Fs
    for im = 1:size(eweights{nx},2)-1 % last is emo index
        [P,table]=anovan(eweights{nx}(:,im),{eweights{nx}(:,end)},'display','off');
        Fs(im,1) = table{2,6};
    end;
    allFs{nx} = Fs/std(Fs);
end;
%**********************************

eval(strs{1}) % load correct Table with Fscore order ('all')
clear subjIMs subjIMFs
for sbj  =1:length(subjlist)
  subj=  subjlist(sbj);
  if ~isempty(TotalAccuracy{subj,1})
    %x=TotalAccuracy{subj,classtype}.Strategy{ForRand,1}.Single;
    %z = reshape(x(end,:,:),[1 size(x,2)*size(x,3)]);
    %for im = 1:size(x,1)
    %  y = reshape(x(im,:,:),[1 size(x,2)*size(x,3)]);
    %  [h p(1,im)] = ttest2(y,z);
    %end;
    %r = find(p(2:end) > .01); r=r(1); % last consec under .01
    r=find(Table(:,subj));
    subjIMs{subj} = sort(Table(r,subj)'); % subject's best IMs
    subjIMFs{subj} = allFs{subj}(subjIMs{subj});
    for im=1:length(bigwts)
      bigwts{subj}{im} = repmat(allFs{subj}(im),[1 length(bigwts{subj}{im})]);
      orivec{subj}{im} = ones(1,length(orivec{subj}{im})); % no more flip
    end;
    %subjIMs{subj} = subjIMs{subj}(find(subjIMs{subj}));
  end;
end;
classfacs = [];classmeans=[];idxs=[];emopatterns = [];
for subj = 1:length(subjIMs)
  if ~isempty(subjIMs{subj})
    s = load([fullpaths{subj},savedat,'.mat']);  
    sph=floatread([fullpaths{subj},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{subj},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    dat=floatread([fullpaths{subj},savedat,'.fdt'],[s.pcs inf],[],0);        
    ws = wts*sph;    acts = ws*dat; 
    subjims = subjIMs{subj}(find(subjIMs{subj})); 
    for im = 1:length(subjims)
      cps = allbigs{subj}{subjims(im)};
      for c = 1:length(cps)
        %emopatterns = [emopatterns;subjemeans{subj}(keepuseims{subj}(im),:)];
        rcp = find(s.complist == cps(c));
        classfacs = [classfacs;acts(subjims(im),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)];
        idxs = [idxs;[subj,subjims(im),cps(c),subjIMFs{subj}(im)]];% last=F score
        classmeans = [classmeans;s.meanpwr(rcp,:)];
      end;       
    end;
  end
end
% make new sorted stuff with only templates that classify:
[deltaclust, lthetaclust,hthetaclust,alphaclust,lbetaclust,hbetaclust,gamaclust,freqs,otherclust] = SortModTempls(classfacs,idxs,classmeans,s.freqs,[0 129],[],'quad',savedat,2.75);
finaltempls = {deltaclust{1},lthetaclust{1},hthetaclust{1},alphaclust{1},lbetaclust{1},hbetaclust{1},gamaclust{1},otherclust{1}}; 
finalmeans = {deltaclust{2},lthetaclust{2},hthetaclust{2},alphaclust{2},lbetaclust{2},hbetaclust{2},gamaclust{2},otherclust{2}}; 
finalidxs = {deltaclust{3},lthetaclust{3},hthetaclust{3},alphaclust{3},lbetaclust{3},hbetaclust{3},gamaclust{3},otherclust{3}}; 
% Calculate % of total templates in each freq range

for cls = 1:length(finalidxs)
  perctmpls(1,cls) = size(finalidxs{cls},1)/size(fullidxs{cls},1);  
end;perctmpls=perctmpls*100;

% calculate median Fscores for each cluster:

for cls = 2:length(finalidxs)
  medFscores(1,cls) = median(finalidxs{cls}(:,4));  
end;
%%%%----------
row=2;col=3; 
fr = find(s.freqs > 0 & s.freqs <= 128);
figure; pl = 1; cls=2;
for cls = 2:7%length(finaltempls)
    sbplot(row,col,pl); pl = pl+1;
    [han,realx labelx] = quadplot(s.freqs(fr),finaltempls{cls}(:,fr)',2,[1 .5 0]); hold on;%[0 .75 .75][.2 1 .2]
    [han,realx labelx] = quadplot(s.freqs(fr),mean(finaltempls{cls}(:,fr),1)',3,'k');%[.16 .5 .3]
    set(gca,'ylim',[min(finaltempls{cls}(:)) max(finaltempls{cls}(:))]);
    set(gca,'ticklength',[.05 .05]);
    title(['Clust ',int2str(cls),' (N=',int2str(size(finaltempls{cls},1)),'; Median F=',num2str(medFscores(1,cls)),')']);      
    plot([realx(2) realx(2)],[get(gca,'ylim')],'g-');
    plot([get(gca,'xlim')],[0 0],'k-');
end;
textsc(['IM templates used for emotion classificaion'],'title');
str = ['print ',basedir,'FreqClustTmpls.eps -depsc -adobe'];eval(str)
str = ['print ',basedir,'FreqClustTmplsDelta.eps -depsc -adobe'];eval(str)
%--------------------------------------------------------
% Collect variables for dipole/density plotting:
%--------------------------------------------------------

[facvec,comods,wtsmat1,justcomps,jcwts,denslist,denswts] = Var4DipPlot(finalidxs,allbigs,bigwts,orivec);

%--------------------------------------------------------
% plot dipoles:

clustnames = {'delta','Lotheta','Hitheta','alpha','Lobeta','Hibeta','gamma','other'};
row = 3;%length(comods);
viewnum=[1,2,3];col = 3;%length(viewnum) ;
zoom= 1;
figure;pl = 1;
for clust = 2:4%length(comods)
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},[],[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print ',basedir,'FreqClustDipsThAlpha.jpg -djpeg'];eval(str)
str = ['print ',basedir,'FreqClustDipsBetaGamma.jpg -djpeg'];eval(str)

clear dlist
for cls = 1:length(denslist) 
  for nx = 1:length(denslist{cls}) % don't repeat ICs in a single     
    dlist{cls}{nx} = unique(denslist{cls}{nx});
  end;    
end;  
cmax=.0015;
for cls = 2:7%length(denslist) 
  if ~isempty(finalidxs{cls})
  figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[58:-20:-22],'mriview','top','geom',[2 3]},'bred',[],'on');% unmasked
   ph=textsc(clustnames{cls},'title'); set(ph,'color','r');
 str = ['print ',basedir,'FreqClustDens',clustnames{cls},'.jpg -djpeg'];eval(str); %close
  figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[-40:20:40],'mriview','side','geom',[5 1],'cmax',cmax},'bred',[],'on');% unmasked
  ph=textsc(clustnames{cls},'title'); set(ph,'color','r');
  str = ['print ',basedir,'FreqClustDens',clustnames{cls},'Side.jpg -djpeg'];eval(str); close
  % MASKED density plots:
  figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5 1],'cmax',cmax},'yred',[],'on');
  ph=textsc(clustnames{cls},'title'); set(ph,'color','r');
  str = ['print ',basedir,'FreqClustDens',clustnames{cls},'Masked.jpg -djpeg'];eval(str); close
  figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[-40:20:40],'mriview','side','geom',[5 1],'cmax',cmax},'yred',[],'on');
  ph=textsc(clustnames{cls},'title'); set(ph,'color','r');
  str = ['print ',basedir,'FreqClustDens',clustnames{cls},'MaskedSide.jpg -djpeg'];eval(str); close
  str = ['save ',basedir,'FreqClustDens',clustnames{cls},'Mask.mat dens minmask maxmask'];eval(str); close
end
end

% find brodmann areas of sig density
clscuts = [.002,.001,.001,.0002,.001];
for cls = 1:length(clustnames) 
  str = ['load ',basedir,'FreqClustDens',clustnames{cls},'Mask.mat dens minmask maxmask'];eval(str); 
  dens(find(dens>minmask&dens<maxmask)) = 0;
  save /home/julie/Manuscripts/Classification/EmoDens.mat dens
  
  list_brain_areas( '/home/julie/Manuscripts/Classification/EmoDens.mat',clscuts(cls),['/home/julie/Manuscripts/Classification/EmoClust',clustnames{cls},'Dens.txt']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load /data/common1/emotion/Classification/YuanPinEmoWeights.mat eweights comment
load /data/common1/emotion/Classification/YuanPinEmoWeights1sec.mat eweights comment
% last column for each subj is indicator variable (1:15)
addpath('/home/julie/MatlabScripts/libsvm-mat-2.9-1')
%addpath('C:\Users\julie\Documents\MatlabScripts\libsvm-mat-2.9-1')
testpercent = .1;
% (DataFormat: [NumOfSample x NumOfFeat] )
clear  AllFs AllDec_Values AllAccuracy AllClassLabels CMD AllSVMModel subjNdelpoints
%str = ['load /data/common1/emotion/Classification/InitialParams1sec.mat AllCMD'];eval(str)
%str = ['load /data/common1/emotion/Classification/InitialParams.mat AllCMD'];eval(str)
for nx = 1:length(eweights)    
    % Step 1: first do optimization based on whole data:
    %[C,Gamma]=SVMParaGridSearch(eweights{nx}(:,1:end-1),eweights{nx}(:,end));
    %CMD = ['-c ' num2str(C) ' -g ' num2str(Gamma)];
    %AllCMD{nx} = CMD;
    CMD = AllCMD{nx};
    clear Fs
    for im = 1:size(eweights{nx},2)-1 % last is emo index
        [P,table]=anovan(eweights{nx}(:,im),{eweights{nx}(:,end)},'display','off');
        Fs(im,1) = table{2,6};
    end;
    allFs{nx} = Fs/std(Fs);
    %[fsort idx] = sort(Fs,'descend'); 
    idx = randperm(size(eweights{nx},2)); % for random order
    clear PredictLabel Accuracy TLs MeanAcc 
    for r = 1:10 % random selections of test data start points        
      trainpoints = [];testpoints = [];
      for e = 1:length(unique(eweights{nx}(:,end))) % break up windows into emos...     
        ndec = round(length(find(eweights{nx}(:,end)==e))*testpercent);% # test points
        tpoints = find(eweights{nx}(:,end)==e)';
        randstart = randperm(length(tpoints) - ndec);randstart = randstart(1);
        rpoints = [randstart:randstart+ndec];% rel test windows
       %rpoints = [round(length(tpoints)/2-ndec/2):round(length(tpoints)/2+ndec/2)];% rel test windows
        testpoints = [testpoints tpoints(rpoints)];% abs test windows
        npoints(1,e) = length(rpoints);
        if randstart > 1 & testpoints(end)< tpoints(end)
           tpoints([rpoints(1)-1,rpoints,rpoints(end)+1]) = [];% take out buffer windows
        elseif testpoints(end)== tpoints(end)% against end edge
           tpoints([rpoints(1)-1,rpoints]) = [];% take out test and buffer windows
        elseif randstart == 1 % against beg end
           tpoints([rpoints,rpoints(end)+1]) = [];% take out test and buffer windows
        end;
        subjstarts(e,r,nx) = randstart;
        trainpoints = [trainpoints tpoints];
        subjNdelpoints{nx} = npoints; % number of test points per subj
      end;
      clear  PL Acc dv TL Models
      for n=1:length(idx) % add one IM in each time (highest to lowest F)
        % Step2: SVM Training Code
        SVMModel = svmtrain(eweights{nx}(trainpoints,end),eweights{nx}(trainpoints,idx(1:n)),CMD);
        % Step3: SVM Testing Code
        [PL(:,n), Acc(:,n), dv(:,:,n)] = svmpredict(eweights{nx}(testpoints,end), eweights{nx}(testpoints,idx(1:n)),SVMModel);
        TL(:,n) = eweights{nx}(testpoints,end);Models(r,n) = SVMModel;
      end;
      PredictLabel(:,:,r) = PL; Accuracy(:,:,r) = Acc; AllDec_Values{nx,r} = dv;
      MeanAcc(1,r) = mean(Acc(1,:)); MeanAcc(2,r) = median(Acc(1,:));
      TLs(:,:,r) = TL; 
    end;
    AllClassLabels{nx} = PredictLabel;% max accuracies
    AllTestLabels{nx} = TLs;
    AllAccuracy{nx} = Accuracy;
    AllMeanAcc{nx} = MeanAcc;
    AllModels{nx} = Models;
    AllFs{nx} = Fs; 
    str = ['save /data/common1/emotion/Classification/ClassVars1Sec10fold.mat subjstarts AllClassLabels AllTestLabels AllAccuracy AllMeanAcc AllDec_Values AllModels AllFs subjNdelpoints'];eval(str)
    str = ['save /data/common1/emotion/Classification/ClassVars1Sec10Rand.mat subjstarts AllClassLabels AllTestLabels AllAccuracy AllMeanAcc AllDec_Values AllModels AllFs subjNdelpoints'];eval(str)
    str = ['save /data/common1/emotion/Classification/ClassVars2Sec10fold.mat subjstarts AllClassLabels AllTestLabels AllAccuracy AllMeanAcc AllDec_Values AllModels AllFs subjNdelpoints'];eval(str)
    str = ['save /data/common1/emotion/Classification/ClassVars2Sec10Rand.mat subjstarts AllClassLabels AllTestLabels AllAccuracy AllMeanAcc AllDec_Values AllModels AllFs subjNdelpoints'];eval(str)
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['save /data/common1/emotion/Classification/InitialParams1sec.mat AllCMD'];eval(str)
%    str = ['save /data/common4/emotion/Fscores_weights.mat AllFs eweights allbigs finaltempls finalmeans finalidxs'];eval(str)
str = ['load C:\Users\julie\Documents\MatlabData\emotion\Fscores_weights.mat'];eval(str)
 str = ['load /data/common1/emotion/Classification/Fscores_weights.mat'];eval(str)

%-----------------------------------------------
% Plot F-scores for each subj
%---------------------------------------------------
for nx = 1:length(AllFs)
    AllFs{nx}(:,2:6) = 0;% initalize new columns to zero
end;
multcls = [];
for cls = 1:length(finalidxs)
    for m = 1:size(finalidxs{cls},1)
        nx = finalidxs{cls}(m,1);
        im = abs(finalidxs{cls}(m,2));
        AllFs{nx}(im,cls+1) = 1; % mark as used
    end;
end;
labels = {'Delta','Theta','Alpha','Beta','Gamma'};
figure; row = 6; col = 6;allvals = [];
for nx = 1:length(AllFs)    
    clsvals = [];
    for cls = 1:length(labels)
        clear cv
        cv(:,2) = AllFs{nx}(find(AllFs{nx}(:,cls+1)),1);
        cv(:,1) = cls;
        clsvals = [clsvals;cv];
    end;
    sbplot(row,col,nx)
    boxplot(clsvals(:,2),clsvals(:,1)); % Fscore vs cluster assignment
    %clsvals(:,2) = clsvals(:,2)/max(clsvals(:,2)); % normalize
    clsvals(:,2) = clsvals(:,2)/std(clsvals(:,2)); % normalize
    chooseF{nx} = clsvals(find(clsvals(:,2) > .75),:);
    allvals = [allvals;clsvals]; % collect across subjs
end;
% normalize by max F score and combine between subjs
figure; boxplot(allvals(:,2),allvals(:,1),'labels',labels); 
ylabel('Normalized F score');
str = ['print C:\Users\julie\Documents\MatlabData\emotion\Classification\FscoreVsFreq.eps -depsc -adobe -painters'];eval(str)
%---------------------------------------------------
% PlotAccuracies
%---------------------------------------------------
str = ['load C:\Users\julie\Documents\MatlabData\emotion\Classification\ClassVars1Sec10fold.mat'];eval(str)
str = ['load C:\Users\julie\Documents\MatlabData\emotion\Classification\ClassVars1Sec10Rand.mat'];eval(str)
str = ['load C:\Users\julie\Documents\MatlabData\emotion\Classification\ClassVars2Sec10fold.mat'];eval(str)
% LINUX ----
    str = ['load /data/common1/emotion/Classification/ClassVars1Sec10fold.mat'];eval(str)
    str = ['load /data/common1/emotion/Classification/ClassVars1Sec10Rand.mat'];eval(str)
    str = ['load /data/common1/emotion/Classification/ClassVars2Sec10fold.mat'];eval(str)
    str = ['load /data/common1/emotion/Classification/ClassVars2Sec10Rand.mat'];eval(str)

% Plot Accuracy for each subj
figure;row = 6; col = 6;collectaccX = [];collectaccN = [];
allsubjsX = zeros(0,0);allsubjsN = zeros(0,0);allsubjsMN = zeros(0,0);
minfeat = 100; clear medit allsubjsX allsubjsN allsubsMN
p = 1;
for nx = 1:length(AllAccuracy)
   clear mxval mnval
   minfeat = min(minfeat,size(AllAccuracy{nx},2));
   sbplot(row,col,nx); 
   for im = 1:size(AllAccuracy{nx},2) 
      [val mxval] = max(AllAccuracy{nx}(1,im,:));
      [val mnval] = min(AllAccuracy{nx}(1,im,:));
      val  = median(AllAccuracy{nx}(1,im,:));
      [val idx] = min(abs(AllAccuracy{nx}(1,im,:)-val));
      medit(nx,im) = AllAccuracy{nx}(1,im,idx);
      allsubjsX(nx,im) = AllAccuracy{nx}(1,im,mxval);
      allsubjsN(nx,im) = AllAccuracy{nx}(1,im,mnval);
      allsubjsMN(nx,im) = mean(AllAccuracy{nx}(1,im,:),3);
      plot(im,AllAccuracy{nx}(1,im,mxval),'r.'); hold on;
      plot([im im],[AllAccuracy{nx}(1,im,mxval)-AllAccuracy{nx}(2,im,mxval) AllAccuracy{nx}(1,im,mxval)+AllAccuracy{nx}(2,im,mxval)],'r-');
      plot(im,AllAccuracy{nx}(1,im,mnval),'b.'); hold on;
      plot([im im],[AllAccuracy{nx}(1,im,mnval)-AllAccuracy{nx}(2,im,mnval) AllAccuracy{nx}(1,im,mnval)+AllAccuracy{nx}(2,im,mnval)],'b-');
      plot([im im],mean(AllAccuracy{nx}(1,im,:),3),'g.'); % mean accuracy
   end;
    title(['Subj ',int2str(nx)]); set(gca,'xticklabel',[]);
   set(gca,'ylim',[0 101]);set(gca,'xlim',[0 im+1]);
end;
ylabel('Accuracy');xlabel('IMs -- F-sorted')
allsubjsX(:,minfeat+1:end) = []; % delete feat #'s above min feats
allsubjsN(:,minfeat+1:end) = []; % delete feat #'s above min feats
allsubjsMN(:,minfeat+1:end) = []; % delete feat #'s above min feats
textsc(str,'title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print C:\Users\julie\Documents\MatlabData\emotion\Classification\SubjAccuracies.eps -depsc -adobe -painters'];eval(str)
%str = ['print /home/julie/Manuscripts/Classification/SubjAccuracies1sec.eps -depsc -adobe -painters'];eval(str)
%str = ['print /home/julie/Manuscripts/Classification/SubjAccuracies2sec.eps -depsc -adobe -painters'];eval(str)
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%str = ['print /home/julie/Manuscripts/Classification/SubjAccuracies1sec.jpg -djpeg'];eval(str)
%str = ['print /home/julie/Manuscripts/Classification/SubjAccuracies2sec.jpg -djpeg'];eval(str)

% Plot mean accuracy across subjs
figure;mnsubjsX = mean(allsubjsX,1);stsubjsX=std(allsubjsX,1);
mnsubjsN = mean(allsubjsN,1);stsubjsN=std(allsubjsN,1);
mnsubjsMN = mean(allsubjsMN,1);stsubjsMN=std(allsubjsMN,1);
for f = 1:size(allsubjsX,2)
   stderrX = stsubjsX(f)/sqrt(size(allsubjsX,1));
   stderrN = stsubjsN(f)/sqrt(size(allsubjsN,1));
   stderrMN = stsubjsMN(f)/sqrt(size(allsubjsMN,1));
   plot([f f],[mnsubjsX(f)-stderrX mnsubjsX(f)+stderrX],'r-');hold on;
  plot([f f],[mnsubjsN(f)-stderrN mnsubjsN(f)+stderrN],'b-');hold on;
  plot([f f],[mnsubjsMN(f)-stderrMN mnsubjsMN(f)+stderrMN],'g-');hold on;
end;
plot(mnsubjsX,'r.-');plot(mnsubjsN,'b.-');plot(mnsubjsMN,'g.-');
set(gca,'ygrid','on');title('Mean/std of all subjs, all emotions')
xlabel('F-score sorted IMs'); ylabel('Accuracy');set(gca,'xlim',[0 f+1])
set(gca,'ylim',[0 101]);%textsc(str,'title');
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print C:\Users\julie\Documents\MatlabData\emotion\Classification\MeanAccuracies.eps -depsc -adobe -painters'];eval(str)
%str = ['print /home/julie/Manuscripts/Classification/MeanAccuracies1sec.jpg -djpeg'];eval(str)

 %str = ['print /home/julie/Manuscripts/Classification/MeanAccuracies2sec.jpg -djpeg'];eval(str)
%str = ['print /home/julie/Manuscripts/Classification/MeanAccuracies1sec.eps -depsc -adobe -painters'];eval(str)

 %str = ['print /home/julie/Manuscripts/Classification/MeanAccuracies2sec.eps -depsc -adobe -painters'];eval(str)
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate percent correct for each emotion (each subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get minfeat from all subj accuracies above
 clear corrperc meanacc subjbest nfeat
for nx = 1:35
   % find min optimal nfeat for all subjs
   tmpfeat = find(diff(movav(medit(nx,:),[],4)) < 0); 
   nfeat(1,nx) = tmpfeat(1); % first time slope is < 0
   subjbest(1,nx) = medit(nx,nfeat(1,nx));
   [corrperc{nx}] = FindEmoAccuracy(AllClassLabels{nx},AllTestLabels{nx});
   for e = 1:size(corrperc{nx},3)
      meanacc{nx}(:,e) = mean(corrperc{nx}(1:nfeat(1,nx),:,e),2)';% mean of iterations
   end;
end;

figure; row = 4; col = 4; cols = jet(length(meanacc));
for e = 1:size(corrperc{nx},3)
   sbplot(row,col,e);
  for nx = 1:35
       plot(meanacc{nx}(:,e)','color',cols(nx,:)); hold on;
   end;
end;
% plot 'best','mid' and 'worst' subjs as Acc by Emotion
figure; hist(subjbest,100);
bsubjs = find(subjbest>70); % best subjs, natural break
msubjs = find(subjbest>=50 & subjbest<=70); % mid subjs
wsubjs = find(subjbest<50); %worst subjs
allsubjs = {bsubjs,msubjs,wsubjs}; cols = {'r','g','b'};
figure; 
for s = 1:length(allsubjs)
   tmpsubjs = [];
   for nxx = 1:length(allsubjs{s})
      nx = allsubjs{s}(nxx);
      tmpsubjs = [tmpsubjs;meanacc{nx}(end,:)];
   end;
   plot(mean(tmpsubjs),'.-','markersize',20,'linewidth',2,'color',cols{s});
   hold on;
end;
legend({['Best subjects (N=',int2str(length(allsubjs{1})),')'],['Middle subjects (N=',int2str(length(allsubjs{2})),')'],['Worst subjects (N=',int2str(length(allsubjs{3})),')']},'location','southeast');
for s = 1:length(allsubjs)
   tmpsubjs = [];
   for nxx = 1:length(allsubjs{s})
      nx = allsubjs{s}(nxx);
      tmpsubjs = [tmpsubjs;meanacc{nx}(end,:)];
   end;
   for e = 1:size(tmpsubjs,2) % plot std
      stderr = std(tmpsubjs(:,e),1)/sqrt(length(allsubjs{s}));
      plot([e e],[mean(tmpsubjs(:,e),1)-stderr mean(tmpsubjs(:,e),1)+stderr],'linewidth',2,'color',cols{s});
   end;
end;
set(gca,'xlim',[0 e+1]);set(gca,'ylim',[0 101]);
ylabel('Accuracy (%)'); xlabel('Emotions');title('Mean +/- Std Err of the Mean');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print C:\Users\julie\Documents\MatlabData\emotion\Classification\AccAcrossEmos.eps -depsc -adobe -painters'];eval(str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCALIZE top feature IMs to brain ICs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classfacs = [];classmeans=[];idxs=[];emopatterns = [];
for nx = 1:length(keepuseims)
     s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    dat=floatread([fullpaths{nx},savedat,'.fdt'],[s.pcs inf],[],0);        
    ws = wts*sph;    acts = ws*dat; 
    [subjFs subjims] = sort(Fs{nx});
    subjims = subjims(1:nfeat(1,nx)); % limit to top features
   for im = 1:length(subjims)
       cps = allbigs{nx}{subjims(im)};
       for c = 1:length(cps)
           %emopatterns = [emopatterns;subjemeans{nx}(keepuseims{nx}(im),:)];
           rcp = find(s.complist == cps(c));
           classfacs = [classfacs;acts(subjims(im),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)];
           idxs = [idxs;[nx,subjims(im),cps(c)]];
           classmeans = [classmeans;s.meanpwr(rcp,:)];
       end;
   end;
end

[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(classfacs,idxs,classmeans,s.freqs,[3 128],[],'quad');
%[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,s.freqs,[3 128],[],'quad');
% plot templates:-----------------------------------------
finaltempls = {deltaclust{1},thetaclust{1},alphaclust{1},betaclust{1},gamaclust{1}};        
row=2;col=2; 
fr = find(s.freqs > 3 & s.freqs < 128);
figure; pl = 1; cls=2;
for cls = 1:length(finaltempls)
    sbplot(row,col,pl); pl = pl+1;
    [han,realx labelx] = quadplot(s.freqs(fr),finaltempls{cls}(:,fr)',2,[1 .5 0]); hold on;%[0 .75 .75][.2 1 .2]
    [han,realx labelx] = quadplot(s.freqs(fr),mean(finaltempls{cls}(:,fr),1)',3,'k');%[.16 .5 .3]
    set(gca,'ylim',[min(finaltempls{cls}(:)) max(finaltempls{cls}(:))]);
    set(gca,'ticklength',[.05 .05]);
    title(['Clust ',int2str(cls),' (',int2str(size(finaltempls{cls},1)),')']);      
    plot([realx(2) realx(2)],[get(gca,'ylim')],'g-');
    plot([get(gca,'xlim')],[0 0],'k-');
end;
%str = ['print FreqClustTmpls.eps -depsc -painters -adobe'];eval(str)
%--------------------------------------------------------
% Collect variables for dipole/density plotting:
%--------------------------------------------------------

[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot([{deltaclust{3}},{thetaclust{3}},{alphaclust{3}},{betaclust{3}},{gamaclust{3}}],allbigs,bigwts,orivec);
          
%--------------------------------------------------------
% plot dipoles:

clustnames = {'delta','theta','alpha','beta','gamma'};
row = length(comods);
viewnum=[1,2,3];col = 3;%length(viewnum) ;
zoom= 1.3;
figure;pl = 1;
for clust = 1:length(comods)
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/FreqClustDips.jpg -djpeg'];eval(str)

for cls = 1:length(denslist) 
    for nx = 1:length(denslist{cls}) % don't repeat ICs in a single 
         dlist{cls}{nx} = unique(denslist{cls}{nx});
    end;    
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5 1]},'bred',[],'on');
    %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[-55:25:55],'mriview','side','geom',[5 1]},'bred',[],'on');
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'Side.jpg -djpeg'];eval(str); close
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'.jpg -djpeg'];eval(str); close
    str = ['print /data/common4/emotion/Classification/FreqClustDens',clustnames{cls},'.jpg -djpeg'];eval(str); close
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5 1]},'yred',[],'on');
    %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[-55:25:55],'mriview','side','geom',[5 1]},'yred',[],'on');
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'MaskedSide.jpg -djpeg'];eval(str); close
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'Masked.jpg -djpeg'];eval(str); close
    str = ['print /data/common4/emotion/Classification/FreqClustDens',clustnames{cls},'Masked.jpg -djpeg'];eval(str); close
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSIFY emotions within-subject ( Old method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; row = 4; col = 4; pl = 1; pt = 1; clear testnims robustemos
clear perccorr pmatchs testnims
for thresh = [10000:-100:1100,1000:-10:10,9:-1:1]
    threshold = repmat(thresh,[1,length(subjlist)]);
    [mtchpercents,meantestmatchs,subjoneconf,fullconf,meanconf,subjemeans,keepuseims,nims] = ClassifyIMdata(savedat,fullpaths,subjlist,minsize,testpercent,threshold,method);
    for nx = 1:size(mtchpercents,1)
        tmpmatch = mtchpercents(nx,~isnan(mtchpercents(nx,:)));
        perccorr(nx,pt) = mean(tmpmatch(find(tmpmatch)));  % single epoch classification 
        tmpmatch2 = meantestmatchs(nx,~isnan(meantestmatchs(nx,:)));
        pmatchs(nx,pt) = 100*(length(find(tmpmatch2))/length(tmpmatch2));% 'test' data mean
    end;
    testnims(pt,:) = nims;    pt = pt+1;     
% $$$     if pl > row*col
% $$$         figure; pl = 1;
% $$$     end;
% $$$     sbplot(row,col,pl); pl = pl+1;
% $$$     for e=1:size(alldatmatchs,2)
% $$$         tmpmatch2 = meantestmatchs(nx,~isnan(meantestmatchs(nx,:)));
% $$$         pmatchs(nx,pt) = 100*(length(find(tmpmatch2))/length(tmpmatch2));% 'test' data mean
% $$$     end;

% $$$     boxplot(mtchpercents);hold on;
% $$$     ph = plot(robustemos,'.-','linewidth',2,'markersize',15); hold on;set(ph,'color',[0 1 0])
% $$$     ph = plot([get(gca,'xlim')],[50 50],'k--','linewidth',2.5); hold on;
% $$$     title([int2str(thresh),' top Fscores']);    
end;
% Use all trials to find best F for each subject:
allfs = [10000:-100:1100,1000:-10:10,9:-1:1];
for nx = 1:size(pmatchs,1)
    [val bestfs(1,nx)] = max(pmatchs(nx,:)); % find highest F that returns best results
end;
threshold = allfs(bestfs);
str = ['save /home/julie/Meetings/sfn2009/FindBestFs.mat pmatchs perccorr allfs bestfs threshold testnims'];eval(str)
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesBoxAllnIMs.eps -depsc -adobe -painters'];eval(str)
%---------------------------------------------------
str = ['load /home/julie/Meetings/sfn2009/FindBestFs.mat'];eval(str)

subjlist = [1:35]; ttl='All Subjs'; % all subjs
[mtchpercents,meantestmatchs,subjoneconf,fullconf,meanconf,subjemeans,keepuseims,nims,subjNdelpoints] = ClassifyIMdata(savedat,fullpaths,subjlist,minsize,testpercent,threshold,method);
for nx = 1:size(mtchpercents,1)
    tmpmatch2 = meantestmatchs(nx,~isnan(meantestmatchs(nx,:)));
    tmatchs(nx,1) = 100*(length(find(tmpmatch2))/length(tmpmatch2));% 'test' data mean
end;

[val i] = max(tmatchs,[],2);%figure; plot(val); % find highest F that returns best results

% how much is 5% of the data (on avg)?
npoints = []; % color magenta
for nx = 1:length(subjNdelpoints)
    npoints = [npoints, subjNdelpoints{nx}];
end;
figure; hist(npoints); xlabel('Number of spectral windows'); 
ylabel('Number of emotions across subjects');
title('Number of spectral windows in 5% of total data for all emotions and all subjects');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/HowLongIs5%.eps -depsc -adobe -painters'];eval(str)

% my best/worst choices by confusion matrices (ey eye):
%subjlist = [1,4,6,9,18,19,21,25,29]; % best classification subjects
%subjlist = [7,8,12,13,27,31,32]; %worst subjects
%subjlist = [2,3,5,10,11,14,15,16,17,20,22,23,24,26,28,30,33,34,35]; % the rest

%-----------
% LINE plots of 5% data avg percent classification accuracy
%-----------

figure; cols = {'r',[.6 0 1],'b','g'};
for x=1:4 % 4 subject groups
    if x == 1
        subjlist = find(val >= 99)'; 
        ttl=['Best Subjects (N=',int2str(length(subjlist)),')'];% best classification subjects
    elseif x == 2
        subjlist = find(val >= 80 & val < 99)'; 
        ttl=['Middle Subjects (N=',int2str(length(subjlist)),')'];% the rest
    elseif x == 3
        subjlist = find(val <= 79)';
        ttl=['Worst Subjects (N=',int2str(length(subjlist)),')'];%worst subjects
    elseif x == 4
        subjlist = [1:35]; ttl=['All Subjs (N=',int2str(length(subjlist)),')']; % all subjs
    end
    ttls{x} = ttl;
    tmpthresh = threshold(subjlist);
    [mtchpercents,meantestmatchs,subjoneconf,fullconf,meanconf,subjemeans,keepuseims,nims,subjNdelpoints] = ClassifyIMdata(savedat,fullpaths,subjlist,minsize,testpercent,tmpthresh,method);
    meantestmatchs = meantestmatchs(subjlist,:);
    mtchpercents =mtchpercents(subjlist,:);
    for e=1:size(meantestmatchs,2)
        tmpmatch2 = meantestmatchs(~isnan(meantestmatchs(:,e)),e);
        robustemos(e) = 100*(length(find(tmpmatch2))/length(tmpmatch2));
    end;
    ph = plot(robustemos,'.-','linewidth',2,'markersize',30); 
    hold on;set(ph,'color',cols{x})
end;
set(gca,'xtick',[0:15]);set(gca,'xlim',[0 16]);
set(gca,'xticklabel',{'','anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'});
ph = plot([get(gca,'xlim')],[50 50],'k--','linewidth',2.5); hold on;
legend(ttls,'location','southeast');
title('Classification accuracy across subjects for best/worst subjects');
ylabel('Percent correct matches'); set(gca,'ylim',[30 110]);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesLinesOnly.eps -depsc -adobe -painters'];eval(str)
%-----------
% Plot box plot for % matches for each emotion (2-sec and 5% tests)--------------
%-----------
figure; 
for e=1:size(meantestmatchs,2)
    tmpmatch2 = meantestmatchs(~isnan(meantestmatchs(:,e)),e);
    robustemos(e) = 100*(length(find(tmpmatch2))/length(tmpmatch2));
end;
boxplot(mtchpercents,'labelorientation' ,'inline');hold on;
ph = plot(robustemos,'.-','linewidth',2,'markersize',30); hold on;set(ph,'color',[0 1 0])
legend('5% data avg classification');
ph = plot([get(gca,'xlim')],[50 50],'k--','linewidth',2.5); hold on;
ylabel('Percent matches'); set(gca,'ylim',[-5 110]);
title(['Percent correct matches from single 2-sec epochs']);
title([ttl]);
set(gca,'xtick',[0:15]);
set(gca,'xticklabel',{'','anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'});
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesBox.eps -depsc -adobe -painters'];eval(str)
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesBoxBest.eps -depsc -adobe -painters'];eval(str)
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesBoxMid.eps -depsc -adobe -painters'];eval(str)
str = ['print /home/julie/Meetings/sfn2009/BestIMsMatchesBoxWorst.eps -depsc -adobe -painters'];eval(str)

%---------------------------------------------------
% Plot confusion matrix (2-sec tests)--------------
%---------------------------------------------------


figure; row = round(sqrt(length(subjoneconf))); col = round(sqrt(length(subjoneconf))); 
mnconf = zeros(15,15); histconf = [];
for nx = 1:length(subjoneconf)
    % zero out diagonal:
    for x = 1:size(subjoneconf{nx},1)
        subjoneconf{nx}(x,x) = 0;
    end;
    histconf = [histconf,subjoneconf{nx}];
    sbplot(row,col,nx);
    imagesc(subjoneconf{nx},[0 50]); 
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    title([int2str(nims(nx)),' IMs; F=',int2str(threshold(nx))]);
end;cbar;
textsc('Confusion matrices for all subjects','title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/IndivSubjConfMats.eps -depsc -adobe -painters'];eval(str)
% zero out diagonal:
for x = 1:size(meanconf,1)
    meanconf(x,x) = 0;
end;
figure; imagesc(meanconf,[0 10]);
title('Mean confusion matrix (all subjs; 5% data means)');cbar;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/MeanConfMats.eps -depsc -adobe -painters'];eval(str)

%---------------------------------------------------
% Plot F-scores vs corresponding number of IMs used 
%---------------------------------------------------
figure; plot(threshold(find(val >= 99)),nims(find(val >= 99)),'r*','markersize',10); hold on;
plot(threshold(find(val >= 80 & val < 99)),nims(find(val >= 80 & val < 99)),'go','markersize',10)
plot(threshold(find(val <= 79)),nims(find(val <= 79)),'bv','markersize',10);
legend('Best subjs','Middle subjs','Worst subjs'); set(gca,'xlim',[-50 max(threshold)+50]);
xlabel('F-score'); ylabel('Number of IMs');
textsc('F-score cut-off and number of IMs used for classification','title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/FscoreByNumIMs.eps -depsc -adobe -painters'];eval(str)


%-----------
% plot performance of 5% and 2-sec epochs as fxn of F score
%-----------
yax = [30 100];xax = [0 2000]; clear plotmean
figure;
for f = 1:size(perccorr,2)
    plotmean(1,f) = mean(perccorr(find(~isnan(perccorr(:,f))),f));
    ploterr(1,f) = plotmean(1,f)+ std(perccorr(find(~isnan(perccorr(:,f))),f));
    ploterr(2,f) = plotmean(1,f)- std(perccorr(find(~isnan(perccorr(:,f))),f));
end;
ph = plot(allfs,plotmean);hold on;
set(ph,'marker','.');set(ph,'markersize',25);
set(ph,'color','b'); clear plotmean
ph = plot(allfs,ploterr,'b--','linewidth',.5);
for f = 1:size(perccorr,2)
    plotmean(1,f) = mean(pmatchs(find(~isnan(pmatchs(:,f))),f));
    ploterr(1,f) = plotmean(1,f)+ std(perccorr(find(~isnan(perccorr(:,f))),f));
    ploterr(2,f) = plotmean(1,f)- std(perccorr(find(~isnan(perccorr(:,f))),f));
end;
ph = plot(allfs,plotmean); 
set(ph,'markersize',25);set(ph,'marker','.');
set(ph,'color','g');
legend('2-sec epoch','2-sec std','5% data avg','5% std');
ph = plot(allfs,ploterr,'g--','linewidth',.5);
set(gca,'ylim',yax);set(gca,'xlim',xax);
xlabel('F-score threshold'); 
ylabel('Correct classification (%)');
title('Mean percent correct matches across subjects and emotions')
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/MeanPercVsFscore.eps -depsc -adobe -painters'];eval(str)

%---------------------------------------
% COLLECT templates and info for classification IMs
%---------------------------------------

classfacs = [];classmeans=[];idxs=[];emopatterns = [];
for nx = 1:length(keepuseims)
     s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    dat=floatread([fullpaths{nx},savedat,'.fdt'],[s.pcs inf],[],0);        
    ws = wts*sph;    acts = ws*dat; 
   for im = 1:length(keepuseims{nx})
       cps = allbigs{nx}{keepuseims{nx}(im)};
       for c = 1:length(cps)
           emopatterns = [emopatterns;subjemeans{nx}(keepuseims{nx}(im),:)];
           rcp = find(s.complist == cps(c));
           classfacs = [classfacs;acts(keepuseims{nx}(im),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)];
           idxs = [idxs;[nx,keepuseims{nx}(im),cps(c)]];
           classmeans = [classmeans;s.meanpwr(rcp,:)];
       end;
   end;
end

[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(classfacs,idxs,classmeans,s.freqs,[3 128],[],'quad');
%[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(clustfacs,kptk,mnspecs,s.freqs,[3 128],[],'quad');
% plot templates:-----------------------------------------
finaltempls = {thetaclust{1},alphaclust{1},betaclust{1},gamaclust{1}};        
row=2;col=2; 
fr = find(s.freqs > 3 & s.freqs < 128);
figure; pl = 1; cls=2;
for cls = 1:length(finaltempls)
    sbplot(row,col,pl); pl = pl+1;
    [han,realx labelx] = quadplot(s.freqs(fr),finaltempls{cls}(:,fr)',2,[1 .5 0]); hold on;%[0 .75 .75][.2 1 .2]
    [han,realx labelx] = quadplot(s.freqs(fr),mean(finaltempls{cls}(:,fr),1)',3,'k');%[.16 .5 .3]
    set(gca,'ylim',[min(finaltempls{cls}(:)) max(finaltempls{cls}(:))]);
    set(gca,'ticklength',[.05 .05]);
    title(['Clust ',int2str(cls),' (',int2str(size(finaltempls{cls},1)),')']);      
    plot([realx(2) realx(2)],[get(gca,'ylim')],'g-');
    plot([get(gca,'xlim')],[0 0],'k-');
end;
str = ['print /home/julie/Meetings/sfn2009/FreqClustTmpls.eps -depsc -painters -adobe'];eval(str)
str = ['print /home/julie/Meetings/sfn2009/FreqClustTmplsAll.eps -depsc -painters -adobe'];eval(str)
%--------------------------------------------------------
% Collect variables for dipole/density plotting:
%--------------------------------------------------------

[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot([{deltaclust{3}},{thetaclust{3}},{alphaclust{3}},{betaclust{3}},{gamaclust{3}}],allbigs,bigwts,orivec);
          
%--------------------------------------------------------
% plot dipoles:

clustnames = {'delta','theta','alpha','beta','gamma'};
row = length(comods);
viewnum=[1,2,3];col = 3;%length(viewnum) ;
zoom= 1.3;
figure;pl = 1;
for clust = 1:length(comods)
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/FreqClustDips.jpg -djpeg'];eval(str)

for cls = 1:length(denslist) 
    for nx = 1:length(denslist{cls}) % don't repeat ICs in a single 
         dlist{cls}{nx} = unique(denslist{cls}{nx});
    end;    
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5 1]},'bred',[],'on');
    %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],[],[],{'mrislices',[-55:25:55],'mriview','side','geom',[5 1]},'bred',[],'on');
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'Side.jpg -djpeg'];eval(str); close
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'.jpg -djpeg'];eval(str); close
    str = ['print /data/common4/emotion/Classification/FreqClustDens',clustnames{cls},'.jpg -djpeg'];eval(str); close
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5 1]},'yred',[],'on');
    %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,dlist{cls},[],gdcomps,[],{'mrislices',[-55:25:55],'mriview','side','geom',[5 1]},'yred',[],'on');
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'MaskedSide.jpg -djpeg'];eval(str); close
    %str = ['print /home/julie/Meetings/sfn2009/FreqClustDens',clustnames{cls},'Masked.jpg -djpeg'];eval(str); close
    str = ['print /data/common4/emotion/Classification/FreqClustDens',clustnames{cls},'Masked.jpg -djpeg'];eval(str); close
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARE IMs pos or neg weighted for each classification?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clustnames = {'theta','alpha','beta','gamma'};
cols = jet(15); clim = .25; % max/min emomeans locked to 1
for f = 4:-1:1 % theta,alpha,beta,gamma
    for e = 1:length(emos)% 15 emos
        for nx = 1:length(denslist{f})% [nx im ic]
            wtcell{nx} = [];
            for im = 1:length(denslist{f}{nx})
                id = find(idxs(:,1)==nx&idxs(:,2)==facvec{f}{nx}(im)&idxs(:,3)==denslist{f}{nx}(im));
                wtcell{nx} = [wtcell{nx} subjemeans{nx}(facvec{f}{nx}(im),e)];
            end;
        end;
        %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,denslist{f},wtcell,denslist{f},[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5,1]},'yred',[],'on');
        figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,denslist{f},wtcell,[],{'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','view', [0 0 1]},[],{'r'},1,[58:-20:-22],[]);
        %textsc([clustnames{f},' influence on ',emos{e}],'title');
        %str = ['print /home/julie/Meetings/sfn2009/',clustnames{f},'on',emos{e},'Dens.jpg -djpeg'];eval(str); close 
        str = ['print /home/julie/Meetings/sfn2009/',clustnames{f},'on',emos{e},'Dips.jpg -djpeg'];eval(str); close all
    end;
end;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are there clusterable patterns of weightings across emotions?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emopatterns = [];newidxs = [];
for nx = 1:length(keepuseims)% [nx im ic]
    emopatterns = [emopatterns;subjemeans{nx}(keepuseims{nx},:)];
    newidxs = [newidxs; [repmat(nx,[length(keepuseims{nx}) 1]) keepuseims{nx}']];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
% First calculate correlation distances:---------------
alldist = pdist(emopatterns, 'correlation'); % context vectors

%--pdist clustering-------------------  
% optional: prune for highest question weighting before clustering
nclust =16;
links = linkage(alldist,'complete');
figure;[hnd,idx,perm]=  dendrogram(links,nclust);close 
labels = emo2;
cols= jet(length(labels));method = 'Correlation';
figure; row = 6;col =6;
 pl = 1;pg = 1;  clear allidx icidx tmpls patts idxfreq newidx
for cls = 1:max(idx)
    allidx = find(idx==cls);
    icidx{cls} = idxs(allidx,:);% nx, ic, dim
    templs{cls} = classfacs(allidx,:);
    clustmeans{cls} = classmeans(allidx,:);
    patts{cls} = emopatterns(allidx,1:length(labels)); 
[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(templs{cls},icidx{cls},clustmeans{cls} ,s.freqs,[3 128],[],'quad'); close
ntempls = [size(thetaclust{1},1),size(alphaclust{1},1),size(betaclust{1},1),size(gamaclust{1},1)];
newidx{cls} = [abs(thetaclust{3});abs(alphaclust{3});abs(betaclust{3});abs(gamaclust{3})];
idxfreq{cls} = [ones(1,size(thetaclust{1},1)),ones(1,size(alphaclust{1},1))*2,ones(1,size(betaclust{1},1))*3,ones(1,size(gamaclust{1},1))*4];
    if pl > row*col
        textsc(['Distance clusters -- ',method,' method'],'title');
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        %str = ['print /home/julie/Meetings/sfn2009/EmoPatterns',int2str(pg),'.eps -depsc -painters -adobe']; eval(str)
        set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        %str = ['print /home/julie/Meetings/sfn2009/EmoPatterns',int2str(pg),'.jpg -djpeg']; eval(str)
        figure;pl=1;pg = pg+1;
    end;    
    if ~isempty(templs{cls})
        sbplot(row,col,pl); pl=pl+1;
        ph = plot(s.freqs,templs{cls}); set(gca,'xlim',[s.freqs(1) s.freqs(end)]);hold on;
        title(['# tmpls: ',int2str(ntempls)]);    
        set(gca,'xticklabel',[]);
        sbplot(row,col,pl); pl=pl+1;
        ph = plot([1:length(labels)],mean(patts{cls},1),'k-','linewidth',1); hold on;
        for c = 1:size(patts{cls},2)
            ph = plot(c,patts{cls}(:,c)','k.');hold on;
            set(ph,'color',cols(c,:));set(ph,'markersize',10);         
            ph=text(c,max(patts{cls}(:,c)'),labels{c});;
            set(ph,'color',cols(c,:));set(ph,'fontsize',8);
            set(ph,'rotation',90);
        end;
        ph = plot([1:length(labels)],mean(patts{cls},1),'k-','linewidth',1); hold on;
        set(gca,'xlim',[0 length(labels)+1]);set(gca,'xticklabel',[]);
        ph = plot([get(gca,'xlim')],[0 0],'k-');set(gca,'ylim',[-1 1]);
        title(['Cls ',int2str(cls),'-- Ss ',int2str(length(unique(icidx{cls}(:,1))))]);    
    end;
end;
textsc(['Distance clusters -- ',method,' method'],'title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/EmoPatterns',int2str(pg),'.eps -depsc -painters -adobe']; eval(str)
set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/EmoPatterns',int2str(pg),'.jpg -djpeg']; eval(str)

clustcps = cell(1,length(newidx));wtcell = cell(1,length(newidx));
for cls = 1:length(newidx)
    for nx = 1:length(gdcomps)
        clustcps{cls}{nx} = [];wtcell{cls}{nx} = [];
        if ~isempty(find(newidx{cls}(:,1) == nx))
            imidx = find(newidx{cls}(:,1)==nx);
            subjidx = newidx{cls}(imidx,2);
            for im = 1:length(subjidx)
                id = find(idxs(:,1) == nx & idxs(:,2) == subjidx(im));
                clustcps{cls}{nx} = [clustcps{cls}{nx} idxs(id,3)'];
                wtcell{cls}{nx} = [wtcell{cls}{nx} repmat(idxfreq{cls}(imidx(im)),[1 length(id)])];
            end;
        end;
    end;
end;
clustvec = [1:length(clustcps)]; row = 2; col = 2; place = 1;
figure; PlotDipoleClusters('sources.set',fullpaths,gdcomps,clustcps,clustvec,row,col,place,[method,' Emo Pattern Clusts '],[1,2,3,4],[],'off');
for cls = 1:length(clustcps) % this has repeated dipoles, so looks like all gamma (last)
    figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,clustcps{cls},wtcell{cls},[],{'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','view', [0 0 1]},[],[],4,[58:-20:-22],[]);
    ph=textsc(['Cluster ',int2str(cls)],'title');  set(ph,'color','r')
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    str = ['print /home/julie/Meetings/sfn2009/EmoPatternDipsCls',int2str(cls),'.jpg -djpeg']; eval(str)
    %figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,clustcps{cls},wtcell{cls},[],{'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','view', [0 0 1]},[],[],[],[],[]);
end;
for cls = 1:length(clustcps)
figure; [dens,minmask,maxmask] = PlotDipoles('sources.set', fullpaths,clustcps{cls},[],[],[],{'mrislices',[58:-20:-22],'mriview','top','geom',[5,1]},'yred',[],'on');
ph=textsc(['Cluster ',int2str(cls)],'title');set(ph,'color','r')
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Were the best classified also the best rated emotions?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not correlated and separating into pos/neg zscore also doesn't show consistency.


r = load('/data/common1/emotion/SubjRatingsMatrix.mat');% Not in valence order

for e=1:size(meantestmatchs,2)
    tmpmatch2 = meantestmatchs(~isnan(meantestmatchs(:,e)),e);
    robustemos(e) = 100*(length(find(tmpmatch2))/length(tmpmatch2));
    tt = r.fullrats(subjlist,e);
    mnrats(e) = median(tt(~isnan(tt)));
    stdrats(e) = std(tt(~isnan(tt)));
end;

figure; [ax h1 h2] = plotyy([1:15],mnrats,[1:15],robustemos);hold on; 
set(h2,'color','g');set(h1,'color',[.5 0 1]); set(ax(2),'ticklength',[0 0]);
set(ax(2),'ycolor','g');set(ax(1),'ycolor',[.5 0 1]);set(ax(2),'ylim',[60 110]);
set(ax,'xlim',[0 16]); set(ax,'xtick',[1:15]);
set(ax(1),'ylim',[.5 9.5]);set(ax(1),'ytick',[0:10]);set(ax(1),'yticklabel',[0:10]);
set(h2,'linewidth',2);
set(h2,'marker','.'); set(h2,'markersize',25);
boxplot(r.fullrats,'colors',[.5 0 1]);set(gca,'xtick',[1:15]);set(gca,'xlim',[0 16]);
set(gca,'ticklength',[0 0]);
ylabel('Intensity rating (z-score)');xlabel('Emotions');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/IntensityVsAccuracyBoxLine.eps -depsc'];eval(str); 


for nx = 1:size(r.zrats,1)
    %goodemos=  find(r.zrats(nx,:) > 0);
    %bademos = find(r.zrats(nx,:) < 0);
    goodemos=  find(r.fullrats(nx,:) >7); 
    bademos = find(r.fullrats(nx,:) <=7);  
    gdmean(1,nx)  = 100*(length(find(meantestmatchs(nx,goodemos) == 1))/length(goodemos));
    badmean(1,nx) = 100*(length(find(meantestmatchs(nx,bademos) == 1))/length(bademos)); % percent correct matches
end;
figure; plot(gdmean,'r.-'); hold on; plot(badmean,'b.-');
ylabel('Percent of correct matches'); xlabel('Subjects');

figure; hist([gdmean-badmean]'); ylabel('Number of subjects'); % color orange
xlabel('% accuracy difference (high - low intensity emotions)');
title('Accuracy of emotions rated > 7 (high) and <= 7 (low)');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/IntensityVsAccuracyHist.eps -depsc'];eval(str); 


histvals = zeros(2,0)
for rr = min(round(r.fullrats)):max(round(r.fullrats))
    [x y] = find(r.fullrats == rr);
    histvals(1,end+1) = length(find(meantestmatchs(x,y) == 1));
    histvals(2,end) = length(find(meantestmatchs(x,y) == 0));
end;
figure; ph=bar([min(round(r.fullrats)):max(round(r.fullrats))],histvals');
legend('Correct','Incorrect')
ylabel('Number of classifications'); xlabel('Intensity ratings');
title('Classification accuracy relative to self-reported intensity')
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/AccuracyWinIntensityRating.eps -depsc'];eval(str); 
%------------------------------------------------------------------
% of the zeros for each subj, what are the mean ratings (compared to ratings for 1's)?
%------------------------------------------------------------------
whatrat = r.fullrats; % zrats or fullrats
for nx = 1:size(meantestmatchs,1)
    nomatch = find(meantestmatchs(nx,:)==0); % incorrect classification
    bdrats = whatrat(nx,nomatch);
    nomatchrat(1,nx) = mean(bdrats(~isnan(bdrats))); % find corresponding emo ratings
    nomatchstd(1,nx) = std(bdrats(~isnan(bdrats))); % find corresponding emo ratings
    match= find(meantestmatchs(nx,:)==1); % incorrect classification
    gdrats = whatrat(nx,match);
    matchrat(1,nx) = mean(gdrats(~isnan(gdrats))); % find corresponding emo ratings
    matchstd(1,nx) = std(gdrats(~isnan(gdrats))); % find corresponding emo ratings
end;
diffmat = matchrat-nomatchrat;diffmat=diffmat(~isnan(diffmat));

figure; 
ph = bar(find(diffmat>0),diffmat(find(diffmat>0))); hold on;
set(ph,'facecolor','b');
ph = bar(find(diffmat<0),diffmat(find(diffmat<0)));
set(ph,'facecolor','r');
set(gca,'xtick',[1:length(diffmat)+1]); set(gca,'xticklabel',[]);
set(gca,'xlim',[0 length(diffmat)+1]);
ylabel('Mean intensity ratings for correct - incorrect classifications');
xlabel('Subjects with correct and incorrect classifications');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/Corr-IncorrIntenseRating.eps -depsc'];eval(str); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices for Yuan-Pin Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate 1-sec freq measures for all gdcomps and chans using same transform as IM decomp:
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab file loads all subject info needed
fullpaths = newpaths;
subjlist = [1:4,6:14,16:21,23:35];% elim subjs 5,15,22 is a repeat
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
savedat = 'ChanPwrDecomp'; frqlim = [3 125];overlap = 1;nfreqs = 370;freqscale = 'quad';wsize = 1;chancomp= 'chan';
savedat = 'SpecCoModNoOvrlap'; frqlim = [3 125];overlap = 1;nfreqs = 370;freqscale = 'quad';wsize = 1; chancomp= 'comp';

for nxx = length(subjlist):-1:1
  nx=subjlist(nxx);
  CalcFFTonly(datset,fullpaths{nx},[],savedat,frqlim,freqscale,overlap,nfreqs,wsize,chancomp);
end;
  
% Compile data into subj cell array
s = load([fullpaths{1},savedat,'.mat']);  
specranges = {[find(s.freqs >3 & s.freqs < 8)],[find(s.freqs >=  8 & s.freqs <= 13.5)],[find(s.freqs > 13.5 & s.freqs <= 20)],[find(s.freqs > 20 & s.freqs <= 35)],[find(s.freqs > 35)]};
  
clear eweights
for nxx = 1:length(subjlist)
  nx=subjlist(nxx);
  
  s = load([fullpaths{nx},savedat,'.mat']);  
  dat = floatread([fullpaths{nx},savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);
  allCs = [];
  for c = 1:length(s.complist)
    onec = dat(:,(c-1)*length(s.freqs)+1:c*length(s.freqs));
    clear specest
    for sp = 1:length(specranges)
      specest(:,sp) = mean(onec(:,specranges{sp}'),2);
    end;
    allCs = [allCs,specest];
  end;
  eidx = [];
  for e = 1:length(s.dstrials) % break up windows into emos...
    eidx = [eidx; ones(length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))),1)*e];
  end; 
  eweights{nx} = [allCs,eidx];
end;

comment = 'The variable ''eweights'' is a cell array with all 35 subjects. Each cell (subject) contains a matrix with rows equal to the number of spectral windows (0% overlap of 1-sec windows) by number of ICs*number of freq bins.The FIVE Freq bins are: theta(3-8 Hz), Alpha (8-13.5 Hz), Low Beta(13.5-20 Hz), Hi Beta (20-35 Hz), and Gamma (>35 Hz). The columns are arranged as the five freq bands for each IC, followed by the 5 freq bands (in order) for the next IC. Each spectral estimate was obtained by FFT (pwelch method). The last column of the matrix is an indicator matrix with numbers between 1 and 15 indicating which of the 15 emotions that spectral window (row) was taken from.'; 
save /data/common1/emotion/YuanPinICSpecPwr1Sec.mat eweights comment

comment = 'The variable ''eweights'' is a cell array with all 35 subjects. Each cell (subject) contains a matrix with rows equal to the number of spectral windows (0% overlap of 1-sec windows) by number of Chans*number of freq bins.The FIVE Freq bins are: theta(3-8 Hz), Alpha (8-13.5 Hz), Low Beta(13.5-20 Hz), Hi Beta (20-35 Hz), and Gamma (>35 Hz). The columns are arranged as the five freq bands for each Chan, followed by the 5 freq bands (in order) for the next Chan. Each spectral estimate was obtained by FFT (pwelch method). The last column of the matrix is an indicator matrix with numbers between 1 and 15 indicating which of the 15 emotions that spectral window (row) was taken from.'; 
save /data/common1/emotion/YuanPinChanSpecPwr1Sec.mat eweights comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedat ='SpecCoModNoOvrlap'; % for 1-sec, no overlap
%savedat = 'SpecCoModMoreFreqs';% for 2-sec, 50% overlap

for nx = 1:35%length(fullpaths)
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;    winv = inv(ws); 
    clear wts sph ws 
    speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  % templates   
    winv = specwts;  clear delpoints emeans  
    eidx = [];
    for e = 1:length(s.dstrials) % break up windows into emos...
        eidx = [eidx; ones(length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))),1)*e];
    end; 
    eweights{nx} = [winv,eidx];
end;
comment = 'The variable ''eweights'' is a cell array with all 35 subjects. Each cell (subject) contains a matrix with rows equal to the number of spectral windows (50% overlap of 2-sec windows) by number of spectral modulators as columns. Each spectral estimate was obtained by FFT (pwelch method) and then decomposed by PCA/ICA into independent spectral modulators. The last column of the matrix is an indicator matrix with numbers between 1 and 15 indicating which of the 15 emotions that spectral window (row) was taken from.'; 
save /data/common1/emotion/YuanPinEmoWeights.mat eweights comment
comment = 'The variable ''eweights'' is a cell array with all 35 subjects. Each cell (subject) contains a matrix with rows equal to the number of spectral windows (0% overlap of 1-sec windows) by number of spectral modulators as columns. Each spectral estimate was obtained by FFT (pwelch method) and then decomposed by PCA/ICA into independent spectral modulators. The last column of the matrix is an indicator matrix with numbers between 1 and 15 indicating which of the 15 emotions that spectral window (row) was taken from.'; 
save /data/common1/emotion/YuanPinEmoWeights1sec.mat eweights comment
