% PCAs emotion scores from ica clusters within and between subject
% input variables and nxlists from PlotEmoClusters.m
% input scores from DescriptStats.m

button = [2:36]; % don't use tl81, events screwed up.
button = [1:12]; % Repetitve button presses
button = [21,23:26]; % only when 'feeling it' button presses (no mr72-2)
button = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
button = [13:20,36]; % no button press (apart from the first one)
button = [1,2,4:6,8:12,14,17:21,23,25:30,31,33,34,35,36]; % all 'good' subjects (ones that said they got into it)
button = [2:21,23:36];  % all but mr72-2
button = [2:7,9:21,23:34,35,36];  % no mr72-2, tl81 or ar81
button = [1,3:9,12,14,16,17,19,21,22,23,24,26,27,29,33,36]; % females
button = [1,4:6,8,9,12,14,17,19,21,23,26,27,29,33,36]; % 'good' females
button = [2,10,11,13,15,18,20,25,28,30,31,32,34,35]; % males
button = [2,10,11,18,20,25,28,30,31,34]; % 'good' males

button = [2,3,7,16,19,21,29]; % emo order 1
button = [4,10,17,24,25,30,33]; % emo order 2
button = [5,8,12,18,26,28,34]; % emo order 3
button = [6,11,15,23,32,35]; % emo order 4
button = [9,13,14,20,22,27,31]; % emo order 5

load /data/common4/emotion/PCASpacegdsubjs.mat ws winv activations button allpcs

savedat = 'SpecCoMod';
savedat = 'SpecCoMod-2';
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
str = ['load /data/common4/emotion/GoodComps.mat ']; eval(str);
subjlist = button;
%%%%%%%%%%%%%%
eeglab
subjlist = button;
%[emomeans] = EmoSpacebyICA(savedat,fullpaths,gdcomps,subjlist); % antiquated
[emomeans,emodeciles] = EmoWeights(savedat,fullpaths,gdcomps,subjlist);
%save /data/common4/emotion/GroupEmoFacAnal.mat emomeans    
load /data/common4/emotion/GroupEmoFacAnal.mat 
 
%%% Option 1:  ***********
%%%%%% Cluster emomeans by diff methods and cut down by emo vectors that are correlated with others %%%%
plotall = 0; % 0: plot only clusters, 1: plot factors/dipoles for all clusters
nclusts = 3; % decimation factor for determining number of clusters.: num input/decfac
[clustmeans, clusttrack,collmeans, keeptrack, polarvec] = CorrEmoMeans(savedat,gdcomps,emos,emomeans,fullpaths,.5,nclusts,plotall,'p');

PlotCorrClust(gdcomps,fullpaths,clustmeans,clusttrack,emomeans,polarvec)

%%% Option 2:  ***********
%%%%% multi-dimensional scaling

subjlist = button;
numdims = 3; mvon = 0; % don't make a movie
corrcut = .5; % corr coeff requirement to be added to MD scaling
[fullwts] = EmoSpacebyMDscale(emodeciles,numdims,subjlist,corrcut,mvon,'Good Males');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cols = jet(15);cols(10,:) = [.9 .9 0];
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure;  row = 6; col = 6; pl=1;
 for nx=1:35
    fullmat = emomeans{nx}';
    dd = pdist(collmeans', 'correlation') ;    
    [md,STRESS,DISPARITIES] = mdscale(dd,2);
    winv = md;
    % for one point per emotion (mean)
    %for wv = 1:size(winv,2)
        sbplot(row,col,pl)
    %    for e = 1:size(winv,1)
    %        ph=plot(e,winv(e,wv),'.');hold on;
    %        set(ph,'markersize',20);set(ph,'color',cols(e,:));
    %        %ph = text(e,winv(e,wv),emo2{e});
    %        %set(ph,'color',cols(e,:));         
    %    end;set(gca,'xlim',[1 16]);  pl = pl+1;
c1 = 1; c2 = 2;
  % for one point per emotion (mean)
    for e = 1:size(winv,1)
    ph=plot(winv(e,c1),winv(e,c2),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    ph = text(winv(e,c1),winv(e,c2),emo2{e});
    set(ph,'color',cols(e,:)); 
    end;
        title(['Sbj ',int2str(nx)]); set(gca,'xticklabel',[]);pl = pl+1;
    %end;
    axcopy
end;
     set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
   textsc(['Multidimensional Scaling of Individual Subjects'],'title');
print /data/common4/emotion/MDscaleEmoSpace.eps -depsc2 -adobecset -painter 

%%%%%%%%%%%%%%%%%%
% this one is to see if spectral cluster factors basically correlate with emotions
PlotSpecClustEmos(savedat,fullpaths,emomeans); % plots emo weights for all 8 spectral clusters
%%%%%%%%%%%%%%%%%%

 subjlist = [2:3];
 PlotEnvEmoSpace(savedat,fullpaths,origwinvs,spaceacts,savewinvs,subjlist); % plots envelopes of dim back-proj
 EmoSpacebyCorr(savedat,fullpaths,subjlist); %maybe a similar thing to second ICA

[faclist, allfacs,clustlist,specs] = Facs2Clusts(gdcomps,savedat,fullpaths,emomeans,whichfacs);
%save /data/common4/emotion/IndivEmoInfo.mat clustlist specs allfacs faclist 
load /data/common4/emotion/IndivEmoInfo.mat clustlist specs allfacs faclist 
emlist = [1:15];
PlotIndivEmoSpace(clustlist,specs,gdcomps,fullpaths,emlist,'/data/common4/emotion/');

gdsubjs = [2,4,7,11,14,16,  12,22,27,28,29,6,10,33];% valence,arousal,singletons
cutline = 2;
for nx = 1:length(gdsubjs)
PlotEmoSpace(emomeans,maxfacs,gdsubjs(nx),[],cutline);
end;
PlotEmoSpace(emomeans,maxfacs,2,[1,2]);
%%%%%%%%%%%%
% for figure purposes:
 PlotEnvEmoSpace(savedat,fullpaths{nx},origwinvs,spaceacts,savewinvs,subjlist);
row = 4; col = 4; place = 1; 
 [emoorders] = FindEmoOrder(fullpaths,emos);
 for nx = 1:35    
    %PlotSubjBackProjs(fullpaths{nx},savedat,emos,1,emos); % valence order
    [newpl] = PlotSubjBackProjs(fullpaths{nx},savedat,emos,[1:15],1,emoorders{nx},row,col,place,1); % presentation order
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    %str  = ['print /data/common4/emotion/BackProjWts',int2str(nx),'.jpg -djpeg']; eval(str)
    %close
end


%%%%%%%%%%%%
%%%%%%%%%%%%
emopair = [7 13]; % sad vs happy
[diffims,P] = CheckEmoShifts(savedat,fullpaths,subjlist,emopair,.00001);
PlotSpecCoModAcrSubj('sources1.set',fullpaths,savedat,gdcomps,diffims,[3 45],12,'IMs with Sig Diff Happy/Sad Weights');
PlotSpecCoModClusts(fullpaths,savedat,gdcomps,diffims,figinfo);
%%%%%%%%%%%%
%%%%%%%%%%%%

    pcamatall = zeros(15,0);clear mnemodiff forstats
alldat = zeros(45,0);moredat = zeros(0,300); clear allpcs
for nxs = 1:length(button)
    nx = button(nxs);
    str = ['load ',fullpaths{nx},savedat,'.mat '];eval(str);    
    allpcs(nx) = pcs;
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
    ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws allfacs alltemps
    
    % find RMS of each factor activation and multiply with winv (doesn't actually change results,
    %for f = 1:size(winv,2)                                      at least when you normalize by 
    %    rms = sqrt(mean(activations(f,:).^2));                  rms to make 'pcanorm')
    %    winv(:,f) = winv(:,f)*rms;
    %end;
    %alltps = zeros(45,0);
    %percscores = zeros(300,0);mnemo = [];
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; clear allemos
        allemos = zeros(0,1);pscores = zeros(0,1);
        for e = 1:length(emos)  % start with 2 for2 straight nums (not diffs)
            clear newmat newmat2
            tempmat = tpwts(sum(dstrials(1:e-1))+1:sum(dstrials(1:e))); %tempmat = sort(tempmat);
            mnemo(e,tp) = mean(tempmat);%instead of diff
            %pl=1;      % for inputing more than just mean           
            %for pk = .25:.25:.75
            %    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
            %end;
            %allemos(end+1:end+3,1) = newmat'; % wts by quartiles
            %pl=1;      % for inputing more than just mean           
            %for pk = .01:.05:.99
            %    newmat2(1,pl) = tempmat(ceil(size(tempmat,2)*pk));pl = pl+1;
            %end;
            %pscores(end+1:end+size(newmat2,2),1) = newmat2'; % wts by fractions
        end; 
        %alltps(:,end+1) = allemos; % add quartiles to a template matrix
        %percscores(:,end+1) = pscores;% add fractions to a template matrix
    end;
    %moredat(end+1:end+size(percscores,2),:) = percscores';% add quartiles to the grand matrix
    %alldat(:,end+1:end+size(alltps,2)) = alltps;% add fractions to the grand matrix
    pcamatall(:,end+1:end+size(mnemo,2)) = mnemo; 
    fprintf('\n One More SUBJECT Done: %i',nx);clear sph winv wts ws activations icamatall
end;
load /data/common4/emotion/FacMatAllFemale.mat pcamatall button allpcs ttl
load /data/common4/emotion/FacMatAllGoodFemale.mat pcamatall button allpcs ttl
load /data/common4/emotion/FacMatAllMale.mat pcamatall button allpcs ttl
load /data/common4/emotion/FacMatAllGoodMale.mat pcamatall button allpcs ttl
load /data/common4/emotion/FacMatAllSubj.mat pcamatall button allpcs ttl
load /data/common4/emotion/FacMatAllGoodSubj.mat pcamatall button allpcs ttl
%%%%%%%%%%%%%%%%%
EmoSpacebyContext(savedat,fullpaths,button);
EmoSpacebyCorr(savedat,fullpaths,[2])
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[allsubjmat] = DescriptMat(savedat,fullpaths,button);

%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize within factor
clear pcanorm
for f = 1:size(pcamatall,2)
    pcanorm(:,f) = pcamatall(:,f)/sqrt(mean(pcamatall(:,f).^2));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%
% run individual emotion space analysis, plus cluster assignment
eeglab
subjlist = zeros(1,0);
for f = 1:length(button)
    subjlist(end+1:end+75) = repmat(button(f),[1,75]);
end;

%%%%%%%%%%%%%%%%%%%%%
% run pca instead
%%[pc,eigvec,sv] = runpca(pcanorm,1);
activations = pc; % now go to cut subejcts part
%%%%%%%%%%%%%%%%%%%%%
%%%%% multi-dimensional scaling
cols = jet(15);cols(10,:) = [.9 .9 0];
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure;  row = 6; col = 6; pl=1;
 for nx=25:35
    fullmat = emomeans{nx}';
    dd = pdist(collmeans', 'correlation') ;
    
    [md,STRESS,DISPARITIES] = mdscale(dd,3);
    winv = md;
    % for one point per emotion (mean)
    for wv = 1:size(winv,2)
        sbplot(row,col,pl)
        for e = 1:size(winv,1)
            ph=plot(e,winv(e,wv),'.');hold on;
            set(ph,'markersize',20);set(ph,'color',cols(e,:));
            %ph = text(e,winv(e,wv),emo2{e});
            %set(ph,'color',cols(e,:));         
        end;set(gca,'xlim',[1 16]);  pl = pl+1;
        title(['Sb ',int2str(nx),' Dim ',int2str(wv)]); set(gca,'xticklabel',[]);
    end;axcopy
end;
     set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
   textsc(['Multidimensional Scaling of Individual Subjects'],'title');
print /data/common4/emotion/MDscaleEmoSpace.eps -depsc2 -adobecset -painter 
%%%%%%%%%%%%%%%%%%%%%

% cluster emoscores by kmeans:
ids = kmeans(pcamatall',15,'replicates',5);
figure;for k = 1:15
    subplot(4,4,k)
    clustemos = pcamatall(:,find(ids==k));
    clustemos = mean(clustemos,2);
    ph=plot(clustemos,'.');
    set(ph,'markersize',10);  
end;
%  save /data/common4/emotion/PCASpacegdsubjs.mat ws winv activations button allpcs
% dimensions:  2, 5, 8 !!!
% print /data/common4/emotion/3DEmoSpace.eps -depsc

% remove subjects without highly weighted factors for specific dimensions
dims = [1,2,3,4,5];cutwt = 2; 
[badsubj] = RmSubjsbyWts(activations,subjlist,dims,cutwt)
%%%%%%%%%%%%%%%%%
 [elimsubj] = ElimSubjbyBackProj(pcanorm,subjlist)
button(find(ismember(button,elimsubj))) = [];
pcanorm2 = pcanorm;
pcanorm2(:,find(ismember(subjlist,elimsubj))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
 [elimsubj2] = ElimSubjbyBackProj(pcanorm2,subjlist)
button(find(ismember(button,elimsubj2))) = [];
pcanorm3 = pcanorm2;
pcanorm3(:,find(ismember(subjlist,elimsubj2))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
 [elimsubj3] = ElimSubjbyBackProj(pcanorm3,subjlist)
button(find(ismember(button,elimsubj3))) = [];
pcanorm4 = pcanorm3;
pcanorm4(:,find(ismember(subjlist,elimsubj3))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
 [elimsubj4] = ElimSubjbyBackProj(pcanorm4,subjlist)
button(find(ismember(button,elimsubj4))) = [];
pcanorm5 = pcanorm4;
pcanorm5(:,find(ismember(subjlist,elimsubj4))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
%%%%%%%%%%%%%%%%%

[suggestrm1,subjprob] = ElimSubjbyWinv(pcanorm,subjlist,5000);

button(find(ismember(button,suggestrm1))) = [];
pcanorm2 = pcanorm;
pcanorm2(:,find(ismember(subjlist,suggestrm1))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
[suggestrm2,subjprob] = ElimSubjbyWinv(pcanorm2,subjlist,5000);

button(find(ismember(button,suggestrm2))) = [];
pcanorm3 = pcanorm2;
pcanorm3(:,find(ismember(subjlist,suggestrm2))) = [];
subjlist = zeros(1,0); % redefine
for f = 1:length(button)
    subjlist(end+1:end+15) = repmat(button(f),[1,15]);
end;
[suggestrm3,subjprob] = ElimSubjbyWinv(pcanorm3,subjlist,5000);
button(find(ismember(button,suggestrm3))) = [];
pcanorm4 = pcanorm3;
pcanorm4(:,find(ismember(subjlist,suggestrm3))) = [];

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
  load /data/common4/emotion/PCASpacegdsubjs.mat ws winv activations button allpcs

  % load /data/common4/emotion/PCASpace.mat ws winv activations
%%%%%%%%%%%%%%%%%%
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(fullmat,'pca',4,'stop',1e-7,'maxsteps',1500);
 ws = weights*sphere; winv = pinv(ws);
%%%%%%%%%%%%%%%%%%
% plot each dimension separately with color-coded emotions
cols = jet(15);cols(10,:) = [.9 .9 0];
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
    title(['Dim ',int2str(wv)]); set(gca,'xticklabel',[]);
end;axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
textsc(ttl,'title');
%%%%%%%%%%%
%%  Plot 2 Dims vs each other:
c1 = 1; c2 = 2;
figure;  % for one point per emotion (mean)
    for e = 1:size(winv,1)
    ph=plot(winv(e,c1),winv(e,c2),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    ph = text(winv(e,c1),winv(e,c2),emo2{e});
    set(ph,'color',cols(e,:)); 
    end;
xlabel(['Winv ',int2str(c1)]);ylabel(['Winv ',int2str(c2)]);

set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
textsc(ttl,'title');

%%%%%%%%%%%
%%  Plot 3 Dims vs each other:
cols = jet(15);cols(10,:) = [.9 .9 0];
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure; % just 3  dims vs each other
c1 = 2; c2 = 3; c3 = 1;
for e = 1:size(winv,1)
    ph=plot3(winv(e,c1),winv(e,c2),winv(e,c3),'.');hold on;
    set(ph,'markersize',25);                set(ph,'color',cols(e,:));
    ph = text(winv(e,c1),winv(e,c2),winv(e,c3),emo2{e});
    set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
end;
zl = get(gca,'zlim');
for e = 1:size(winv,1)
    pl =plot3([winv(e,c1) winv(e,c1)],[winv(e,c2) winv(e,c2)],[zl(1)  winv(e,c3)]);
    set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['Winv ',int2str(c1)]);ylabel(['Winv ',int2str(c2)]);zlabel(['Winv ',int2str(c3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  find contributing factors for good dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dim = dimension of data to plot, dir = positive or negative influence on the dimension
eeglab
dim = 5;  dir =  1;  clear subjwts maxwt idx emowts
emowts = winv(:,dim)*dir;
if dir == 1
cutoff = median(activations(dim,:)) + std(activations(dim,:));
elseif dir == -1
    cutoff = median(activations(dim,:)) - std(activations(dim,:));
end;    
for nxs = 1:length(button)
    nx=button(nxs);
    subjwts = activations(dim,sum(allpcs(1:nx-1))+1:sum(allpcs(1:nx)));    
    if dir == 1
        idx{nx} = find(subjwts> cutoff );
        %idx{nx} = find((subjwts> mean(subjwts)+ std(subjwts)) & subjwts > 0 );
        %[val idx{nx}] = max(subjwts);
        %if val < 0
        %    idx{nx} = [];
        %else      
            maxwt{nx} = subjwts(idx{nx});            
            %emowts{nx} = pcamatall(:,idx{nx});
        %end;        
    else
        idx{nx} = find(subjwts < cutoff );
        %idx{nx} = find((subjwts < mean(subjwts)- std(subjwts)) & subjwts < 0 );
        %[val idx{nx}] = min(subjwts);   
        %if val > 0
        %    idx{nx} = [];
        %else      
            maxwt{nx} = subjwts(idx{nx});            
            %emowts{nx} = pcamatall(:,idx{nx});
        %end;        
    end; 
end;
forsum= zeros(1,0);
for x=1:length(maxwt)
    forsum(1,end+1) = length(maxwt{x});
end;
numfacs = sum(forsum);

%PlotSpecCoModAcrSubj(fullpaths,gdcomps,idx,freqs,12,'Factors separating emotion');
% find highest var comp and overplot/average...
clear data activations alldip kptk allspecs allcmpscales subjspecs allbesa
pl = 1;  allwts = zeros(1,0); new=1;
for  nxs = 1:length(button)
    nx=button(nxs); ALLEEG=[];EEG=[];
    if ~isempty(idx{nx})       
        EEG = pop_loadset('sources.set', fullpaths{nx});    
        if isfield(EEG.dipfit.model,'diffmap')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
        end;
        if isfield(EEG.dipfit.model,'active')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
        end;
        if isfield(EEG.dipfit.model,'select')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
        end;
        str = ['load ',fullpaths{nx},'SpecCoModStuff.mat numtrials numframes freqs keeptrack rmepochs dstrials pcs savedat comment'];eval(str);  
        sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
        data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
        ws = wts*sph;    activations = ws*data;  winv = pinv(ws); clear wts sph ws data
        backprojdat = winv*activations ;
         clear cmpscale
        for fc = 1:length(idx{nx}) 
            dipsources = EEG.dipfit.model(gdcomps{nx}(1));
            x=[1:size(activations,1)];x(idx{nx}(fc)) = [];
            acts = activations; acts(x,:) = 0;
            backproj = winv*acts ;
            allspecs = zeros(length(gdcomps{nx}),length(freqs),15);
            for cmp = 1:length(gdcomps{nx})
                dipsources(1,cmp) = EEG.dipfit.model(gdcomps{nx}(cmp));
                oneorig = backprojdat(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                onebkprj = backproj(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                cmpscale(cmp,fc) = var(onebkprj)/var(oneorig);
                for em = 1:length(dstrials)
                    oneemo = onebkprj(sum(dstrials(1:em-1))+1:sum(dstrials(1:em)),:);%by emo
                    allspecs(cmp,:,em) = mean(oneemo,1);
                end;
            end;                
            allwts(1,end+1:end+size(cmpscale,1)) =  cmpscale(:,fc)';
            if new == 1
                allbesa = dipsources;new = 0;
            else
                allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            end;
            subjspecs{nx}{fc} = allspecs;
        end;
        allcmpscales{nx} = cmpscale;
        %for t = 1:length(idx{nx})
        %    for rcp = 1:length(gdcomps{nx})
        %        alltemps(rcp,:) = activations(idx{nx}(t),length(freqs)*(rcp-1)+1:length(freqs)*rcp);  
        %    end;
        %    x=var(alltemps'); 
        %    y = mean(x);
        %    seltemps = alltemps(find(x > y+std(x)),:);
        %    bigcomp = zeros(1,0);
        %    for tp = 1:length(gdcomps{nx})
        %        [corr,indx,indy,corrs] = matcorr(mean(seltemps,1),alltemps(tp,:));
        %        if corr > .96
        %            bigcomp(1,end+1) = gdcomps{nx}(tp);clear dip
        %            kptk(pl,:) = [nx,idx{nx}(t),gdcomps{nx}(tp)];   pl = pl+1;   
        %        end;
        %    end;
        %    if ~isempty(bigcomp)                
        %        alldips{end+1} = EEG.dipfit.model(bigcomp);
        %        allfacs(end+1,:) = mean(seltemps,1);
        %        allwts(1,end+1) =  maxwt{nx}(t);
        %        factrack(end+1,:) = [nx,idx{nx}(t)]; 
    end;
    fprintf('.');
end;

% plot only the largest cmp from each factor
plotspecs = zeros(0,length(freqs),15); allspecs = zeros(0,length(freqs));kptk = zeros(0,2); pl = 1;
for nx=1:length(allcmpscales)
    if ~isempty(allcmpscales{nx})
        for fc = 1:length(subjspecs{nx})
            for em = 1:size(subjspecs{nx}{fc},3)
                plotspecs(pl,:,em) = subjspecs{nx}{fc}(find(max(allcmpscales{nx}(:,fc))),:,em);
            end;
                allspecs(end+1,:) = mean(subjspecs{nx}{fc}(find(max(allcmpscales{nx}(:,fc))),:,:),3);
                kptk(end+1,:) = [nx,fc];
            pl = pl+1;
        end;
    end;
end;

figure; row = 4; col = 4; cols = jet(15);cols(10,:) = [.9 .9 0];
cutwt = sort(allwts);cutwt = cutwt(round(length(allwts)*.25));
for em = 1:length(emos)
    colspecs = zeros(0,99);
    %for fc = 1:size(plotspecs,1)
        subplot(row,col,em)
        ph = plot(freqs,plotspecs(:,:,em));hold on;
        %colspecs(end+1,:) = plotspecs(:,:,em);
        %ph = plot(freqs,plotspecs(fc,:)*emowts(emo,1));hold on;
        set(ph,'color',cols(em,:));
        %colspecs(end+1,:) = plotspecs(fc,:)*emowts(emo,1);
    %end;
    ph = plot(freqs,mean(plotspecs(:,:,em),1),'k-'); set(ph,'linewidth',2);
    title(emos{em});
    %set(gca,'ylim',[-1.5 1.5]);
end;axcopy
if dir == 1
    textsc(['Positively weighted template factors for each emotion for Dimension ',int2str(dim),'; showing ',int2str(size(colspecs,1)),' Facs'],'title');
else
    textsc(['Negatively weighted template factors for each emotion for Dimension ',int2str(dim),'; showing ',int2str(size(colspecs,1)),' Facs'],'title');
end;    

% Then do a density plot weighting for variance in each factor
 pl = 1; allbesa=[];
for xx = 1:length(alldips)
    if isfield(alldips{xx},'diffmap')
    alldips{xx} = rmfield(alldips{xx},'diffmap');      
    end;
    if isfield(alldips{xx},'active')
    alldips{xx} = rmfield(alldips{xx},'active');      
    end;
    if isfield(alldips{xx},'select')
    alldips{xx} = rmfield(alldips{xx},'select');      
    end;
    allbesa = [allbesa alldips{xx}];
    for hf = 1:length(alldips{xx})    
        wts(1,pl) = allwts(xx); pl = pl+1;
    end;
end;
optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','dipolelength',0,'spheres','on'};
figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',allwts); 
ph = textsc(['Density of dipoles weighted for contribution to Dim ',num2str(dim), '; ',int2str(numfacs),' Facs'],'title');
    set(ph,'color','r');set(ph,'fontsize',14);

    % plot the factor templates corresponding
figure; row = 4; col = 4; cols = jet(15);cols(10,:) = [.9 .9 0];
cutwt = sort(allwts);cutwt = cutwt(round(length(allwts)*.25));
for emo = 1:length(emos)
    colspecs = zeros(0,99);
    for fc = 1:size(allfacs,1)
            if allwts(fc) >= cutwt
                subplot(row,col,emo)
                ph = plot(freqs,allfacs(fc,:)*emowts(emo,1));hold on;
                set(ph,'color',cols(emo,:));
                colspecs(end+1,:) = allfacs(fc,:)*emowts(emo,1);
            end;            
    end;
    ph = plot(freqs,mean(colspecs,1),'k-'); set(ph,'linewidth',2);
    title(emos{emo});
    set(gca,'ylim',[-1.5 1.5]);
end;axcopy
if dir == 1
textsc(['Positively weighted template factors for each emotion for Dimension ',int2str(dim),'; showing ',int2str(size(colspecs,1)),' Facs'],'title');
else
textsc(['Negatively weighted template factors for each emotion for Dimension ',int2str(dim),'; showing ',int2str(size(colspecs,1)),' Facs'],'title');
end;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cluster factor spectra by kmeans   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cfr = find(freqs > 6 & freqs < 34.5);
clear normfacs 
 for x=1:size(allspecs,1)
    normfacs(x,:) = allspecs(x,:)/std(allspecs(x,cfr));
end;
pcnum = 5; 
[pc,eigvec,sv] = runpca(normfacs(:,cfr)',pcnum);

reps=10;optk = 4;
[kout, C,sumd, allsums] = kmeans(pc',optk,'replicates',reps); 
clear clustidxs clustspecs  outliers facvec
for cl = 1:optk
    tmpvec = cell(1,0);
    oneclust = find(kout == cl);
    oneclust';
    fout = allsums(oneclust)/std(allsums(oneclust));
    outliers{cl} = find(fout > 4)';
    kout(oneclust(outliers{cl})) = 0;
    oneclust = find(kout == cl);
    clustspecs{cl} = allspecs(oneclust,:);
    %clustspecs{cl} = normfacs(oneclust,:);
    for nx = 1:length(gdcomps)
        if ~isempty(find(kptk(oneclust,1) == nx))
            tmpvec{nx} = kptk(oneclust(kptk(oneclust,1) == nx),2)';
        end;
    end;
    facvec{cl} = tmpvec;
end;
outliers   % repeat until outliers minimal
figure; row = 2; col =2; pl = 0;
for clust =1:length(clustspecs)
    if ~isempty(clustspecs{clust})
    pl = pl+1;
    subplot(row,col,pl)
    %plot(freqs(cfr),clustspecs{clust}(:,cfr));hold on;
    %ph = plot(freqs(cfr),mean(clustspecs{clust}(:,cfr),1),'k-');set(ph,'linewidth',1.5);
    plot(freqs,clustspecs{clust});hold on;
    ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',1.5);
    %set(gca,'ylim',[-5 15]); set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
    title(['Cluster ',int2str(clust)]);
    end;
end;
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
 
for clust = 1:4
    pl = 1; allbesa=[]; newfacs = zeros(0,99);clear wts
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            subjidx = find(factrack(:,1) == nx);
            for fc = 1:length(facvec{clust}{nx})
                facidx = find(factrack(subjidx,2) == facvec{clust}{nx}(fc));
                realidx = subjidx(facidx);
                
                if isfield(alldips{realidx},'diffmap')
                    alldips{realidx} = rmfield(alldips{realidx},'diffmap');      
                end;
                if isfield(alldips{realidx},'active')
                    alldips{realidx} = rmfield(alldips{realidx},'active');      
                end;
                if isfield(alldips{realidx},'select')
                    alldips{realidx} = rmfield(alldips{realidx},'select');      
                end;
                allbesa = [allbesa alldips{realidx}];
                for hf = 1:length(alldips{realidx})    
                    wts(1,pl) = allwts(realidx); pl = pl+1;
                end;
                newfacs(end+1,:) = allfacs(realidx,:);
            end;
        end;    
    end;
    optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','dipolelength',0,'spheres','on'};
    figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',wts,'plotargs',{'cmax',.06,'cmap',jet}); 
    ph = textsc(['Density of dipoles weighted for contribution to Dim ',num2str(dim), '; Clust ',int2str(clust)],'title');
    set(ph,'color','r');set(ph,'fontsize',14);
    % plot the factor templates corresponding
    figure; row = 4; col = 4; cols = jet(15);cols(10,:) = [.9 .9 0];
    cutwt = sort(wts);cutwt = cutwt(round(length(wts)*.25));
    for emo = 1:length(emos)
        colspecs = zeros(0,99);
        for fc = 1:size(newfacs,1)
            if wts(fc) >= cutwt
                subplot(row,col,emo)
                ph = plot(freqs,newfacs(fc,:)*emowts(emo,1));hold on;
                set(ph,'color',cols(emo,:));
                colspecs(end+1,:) = newfacs(fc,:)*emowts(emo,1);
            end;            
        end;
        ph = plot(freqs,mean(colspecs,1),'k-'); set(ph,'linewidth',2);
        title(emos{emo});
        set(gca,'ylim',[-1 1]);
    end;axcopy
    textsc(['Weighted template factors for each emotion for Dimension ',int2str(dim),'; Cluster ',int2str(clust),'; showing ',int2str(size(colspecs,1)),' Facs'],'title');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteratively eliminate subjects who are not highly weighted in dimensions of interest2% need allpcs from loop at top
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutsubjs = [];gdsubjs = [];

dim =5; clear subjwts maxwt
for nxs = 1:length(button)
    nx=button(nxs);
    subjwts{nx} = activations(dim,sum(allpcs(1:nx-1))+1:sum(allpcs(1:nx)));
    maxwt(1,nx) = max(abs(subjwts{nx}));
end;
cutwt = maxwt(find(maxwt));
tocut = sort(cutwt);   
for x = 1:length(cutwt)-1
    ints(1,x) = tocut(x+1) - tocut(x);
end;
cuts= [];
for x = 1:length(ints)
    if ~isempty(find(ints(x) > median(ints)))
        cuts(end+1) = x;
    else
        break
    end;
end;
if isempty(cuts)
    cutval = tocut(2);    
else
   cutval = tocut(cuts(end)); 
end;
cutsubjs = [cutsubjs find(maxwt<=cutval)] ; % includes subjects who were never a part of analysis


gdsubjs = button;
gdsubjs(find(ismember(gdsubjs,cutsubjs))) = [];

button = unique(gdsubjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 'activations', or factor weights of selected dimensions, coded by subject
dims = [1,4];
figure; scols = jet(length(button));
for nx = 1:length(button)    
    goto2 = {pcadims{button(1:nx)}};    goto2 = cell2mat(goto2);clear subjwts
    for dim = 1:length(dims)
        if nx==1
            subjwts(dim,:) = abs(activations(dims(dim),1:sum(goto2)));
        else
            goto1 = {pcadims{button(1:(nx-1))}};goto1 = cell2mat(goto1);
            subjwts(dim,:) = abs(activations(dims(dim),sum(goto1)+1:sum(goto2)));        
        end;
    end;
    ph = plot(subjwts(1,:),subjwts(2,:),'k.'); hold on; set(ph,'markersize',15);set(ph,'color',scols(nx,:));
end;
% cluster factor wts based on all 6 dims
kmwts = zeros(0,size(winv,2));
for nx = 1:length(button) 
    goto2 = {pcadims{button(1:nx)}};    goto2 = cell2mat(goto2);clear subjwts
    if nx == 1
        kmwts(end+1:end+pcadims{button(nx)},:) = abs(activations(:,1:sum(goto2)))';
    else
        goto1 = {pcadims{button(1:(nx-1))}};goto1 = cell2mat(goto1);
        kmwts(end+1:end+pcadims{button(nx)},:) = abs(activations(:,sum(goto1)+1:sum(goto2)))';
    end;
end;

[ids, C, SUMD, D] = kmeans(kmwts,15,'replicates',5);
figure;for k = 1:15
    subplot(4,4,k)
clustemos = kmwts(find(ids==k),:);
clustemos = mean(clustemos,1);
 ph=plot(clustemos,'-');
     set(ph,'markersize',10);  
end;

% run pca instead:
 [pc,eigvec,sv] = runpca(pcanorm',3);
cols = jet(15);
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure; c1 = 1;  c2 = 2;  c3 = 3;
for e = 1:size(pc,2)
    ph =plot3(pc(c1,e),pc(c2,e),pc(c3,e),'.');hold on;
    set(ph,'markersize',20); set(ph,'color',cols(e,:));
    ph = text(pc(c1,e),pc(c2,e),pc(c3,e),emo2{e});
    set(ph,'color',cols(e,:));     
end;
zl = get(gca,'zlim');
for e = 1:size(pc,2)
    pl =plot3([pc(c1,e) pc(c1,e)],[pc(c2,e) pc(c2,e)],[zl(1)  pc(c3,e)]);
    set(pl,'color',cols(e,:))             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['Dim ',int2str(c1)]); ylabel(['Dim ',int2str(c2)]); zlabel(['Dim ',int2str(c3)]); 

figure;  % for one point per emotion (mean)
for wv = 1:size(pc,1)
    subplot(round(sqrt(size(pc,1))),ceil(sqrt(size(pc,1))),wv)
    for e = 1:size(pc,2)
    ph=plot(e,pc(wv,e),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    ph = text(e,pc(wv,e),emo2{e});
    set(ph,'color',cols(e,:)); 
    end;
end;

%%%%%%%%%%%%%%%%
% plot for three points per emotion:
emo3 = {' anger',' anger',' anger',' frustration',' frustration',' frustration',' jealousy',' jealousy',' jealousy',' fear' ,' fear' ,' fear' ,' disgust',' disgust',' disgust',' grief',' grief',' grief',' sad',' sad',' sad',' compassion',' compassion',' compassion',' love',' love',' love',' relief', ' relief',' relief',' content',' content',' content',' awe',' awe',' awe',' happy',' happy',' happy',' joy',' joy',' joy',' excited',' excited',' excited'};
cols(1:3:43,:) = jet(15);cols(2:3:44,:) = jet(15);cols(3:3:45,:) = jet(15);
row=4; col = 5;            
figure;pp = 1;cb=1; bc = 2; cbc = 2;            
for cb = 1:size(winv,2)-2
    for bc = cb+1:size(winv,2)-1
        for cbc = bc+1:size(winv,2)
            sbplot(row,col,pp)
            for e = 1:3:size(winv,1)-2
                ph =plot3(winv(e,cb),winv(e,bc),winv(e,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e,cb),winv(e,bc),winv(e,cbc),emo3{e});
                set(ph,'color',cols(e,:)); 
                ph =plot3(winv(e+1,cb),winv(e+1,bc),winv(e+1,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e+1,cb),winv(e+1,bc),winv(e+1,cbc),emo3{e});
                set(ph,'color',cols(e,:)); 
                ph =plot3(winv(e+2,cb),winv(e+2,bc),winv(e+2,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e+2,cb),winv(e+2,bc),winv(e+2,cbc),emo3{e});
                set(ph,'color',cols(e,:)); 
            end;
            mx =  max(winv(:,cb));mnx = min(winv(:,cb));
            my =  max(winv(:,bc));mny = min(winv(:,bc));
            mz =  max(winv(:,cbc));mnz = min(winv(:,cbc));
            set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
            zl = get(gca,'zlim');
            for e = 1:size(winv,1)-2
                pl =plot3([winv(e,cb) winv(e,cb)],[winv(e,bc) winv(e,bc)],[zl(1)  winv(e,cbc)]);
                set(pl,'color',cols(e,:))             
                pl =plot3([winv(e+1,cb) winv(e+1,cb)],[winv(e+1,bc) winv(e+1,bc)],[zl(1)  winv(e+1,cbc)]);
                set(pl,'color',cols(e,:))             
                pl =plot3([winv(e+2,cb) winv(e+2,cb)],[winv(e+2,bc) winv(e+2,bc)],[zl(1)  winv(e+2,cbc)]);
                set(pl,'color',cols(e,:))             
            end;
            set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             %set(gca,'xticklabel',[]);              set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);
            
            xlabel(['Winv ',int2str(cb)]); ylabel(['Winv ',int2str(bc)]); zlabel(['Winv ',int2str(cbc)]); 
            pp=pp+1;
        end;
    end;
end;axcopy;set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%5
% run pca  on the winv:
 [pc,eigvec,sv] = runpca(winv,3);
 figure;
 for e = 1:3:size(eigvec,1)
     ph=plot(eigvec(e,1),eigvec(e,2),'.');hold on;
     set(ph,'markersize',20);   set(ph,'color',cols(e,:));
     ph = text(eigvec(e,1),eigvec(e,2),emo3{e}); set(ph,'color',cols(e,:)); 
     ph=plot(eigvec(e+1,1),eigvec(e+1,2),'.');hold on;
     set(ph,'markersize',20);   set(ph,'color',cols(e,:));
     ph = text(eigvec(e+1,1),eigvec(e+1,2),emo3{e}); set(ph,'color',cols(e,:)); 
     ph=plot(eigvec(e+2,1),eigvec(e+2,2),'.');hold on;
     set(ph,'markersize',20);   set(ph,'color',cols(e,:));
     ph = text(eigvec(e+2,1),eigvec(e+2,2),emo3{e}); set(ph,'color',cols(e,:)); 
 end;
 
     
% plot each dimension alone
cols(1:3:43,:) = jet(15);cols(2:3:44,:) = jet(15);cols(3:3:45,:) = jet(15);
figure;  % for three points per emotion
for wv = 1:size(winv,2)
    subplot(3,4,wv)
    for e = 1:3:size(winv,1)-2
    ph=plot(e,winv(e,wv),'.');hold on;
                set(ph,'markersize',20);                set(ph,'color',cols(e,:));
    ph=plot(e,winv(e+1,wv),'.');hold on;
                set(ph,'markersize',20);                set(ph,'color',cols(e,:));
    ph=plot(e,winv(e+2,wv),'.');hold on;
                set(ph,'markersize',20);                set(ph,'color',cols(e,:));    
    end;
end;
cols = jet(15);
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure;  % for one point per emotion (mean)
for wv = 1:size(winv,2)
    sbplot(round(sqrt(size(winv,2))),ceil(sqrt(size(winv,2))),wv)
    for e = 1:size(winv,1)
    ph=plot(e,winv(e,wv),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    ph = text(e,winv(e,wv),emo2{e});
    set(ph,'color',cols(e,:)); 
    end;
end;
%%%%%%%%%%%
figure; % just 2 dims vs each other
c1 = 5; c2 = 2; cols = jet(15);
for e = 1:size(winv,1)
    ph=plot(winv(e,c1),winv(e,c2),'.');hold on;
    set(ph,'markersize',20);                set(ph,'color',cols(e,:));
    ph = text(winv(e,c1),winv(e,c2),emo2{e});
    set(ph,'color',cols(e,:)); 
end;
xlabel(['Winv ',int2str(c1)]);ylabel(['Winv ',int2str(c2)]);
%%%%%%%%%%%
cols = jet(15);
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure; % just 3  dims vs each other
c1 = 2; c2 = 5; c3 = 8;
for e = 1:size(winv,1)
    ph=plot3(winv(e,c1),winv(e,c2),winv(e,c3),'.');hold on;
    set(ph,'markersize',25);                set(ph,'color',cols(e,:));
    ph = text(winv(e,c1),winv(e,c2),winv(e,c3),emo2{e});
    set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
end;
set(gca,'zlim',[-.14 .12]);
zl = get(gca,'zlim');
for e = 1:size(winv,1)
    pl =plot3([winv(e,c1) winv(e,c1)],[winv(e,c2) winv(e,c2)],[zl(1)  winv(e,c3)]);
    set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['Winv ',int2str(c1)]);ylabel(['Winv ',int2str(c2)]);zlabel(['Winv ',int2str(c3)]);
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% plots dimensions vs each other from ICA on just mean emo scores
 
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
cols = jet(15);
figure;  row = 4; col = 5;pp = 1;
for cb = 1:size(winv,2)-2
    for bc = cb+1:size(winv,2)-1
        for cbc = bc+1:size(winv,2)
            subplot(row,col,pp)
            for e = 1:size(winv,1)
                ph =plot3(winv(e,cb),winv(e,bc),winv(e,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e,cb),winv(e,bc),winv(e,cbc),emo2{e});
                set(ph,'color',cols(e,:)); 
            end;
            mx =  max(winv(:,cb));mnx = min(winv(:,cb));
            my =  max(winv(:,bc));mny = min(winv(:,bc));
            mz =  max(winv(:,cbc));mnz = min(winv(:,cbc));
            set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
            zl = get(gca,'zlim');
            for e = 1:size(winv,1)
                pl =plot3([winv(e,cb) winv(e,cb)],[winv(e,bc) winv(e,bc)],[zl(1)  winv(e,cbc)]);
                set(pl,'color',cols(e,:))             
            end;
            set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             %set(gca,'xticklabel',[]);              set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);
            
            xlabel(['Winv ',int2str(cb)]); ylabel(['Winv ',int2str(bc)]); zlabel(['Winv ',int2str(cbc)]); 
            pp = pp+1;
        end;
    end;
end;axcopy
set(gcf,'color','w');

clustact = activations;
clustwinv = winv;
%save /data/common2/emotion/clusters/GrandICAonScores.mat clustact clustwinv  pcamatall
%load /data/common2/emotion/clusters/GrandICAonScores.mat clustact clustwinv  pcamatall


emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excited'};
cols = jet(15);
selfacts = [1:6];
figure;  row = 6; col = 6;pp = 1;
for cb = 1:length(selfacts)
    for bc = cb+1:length(selfacts)
        for cbc = bc+1:length(selfacts)
            subplot(row,col,pp)
            for e = 1:size(clustwinv,1)
                ph =plot3(clustwinv(e,selfacts(cb)),clustwinv(e,selfacts(bc)),clustwinv(e,selfacts(cbc)),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(clustwinv(e,selfacts(cb)),clustwinv(e,selfacts(bc)),clustwinv(e,selfacts(cbc)),emo2{e});
                set(ph,'color',cols(e,:)); 
            end;
            mx =  max(clustwinv(:,selfacts(cb)));mnx = min(clustwinv(:,selfacts(cb)));
            my =  max(clustwinv(:,selfacts(bc)));mny = min(clustwinv(:,selfacts(bc)));
            mz =  max(clustwinv(:,selfacts(cbc)));mnz = min(clustwinv(:,selfacts(cbc)));
            set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
            zl = get(gca,'zlim');
            for e = 1:size(clustwinv,1)
                pl =plot3([clustwinv(e,selfacts(cb)) clustwinv(e,selfacts(cb))],[clustwinv(e,selfacts(bc)) clustwinv(e,selfacts(bc))],[zl(1)  clustwinv(e,selfacts(cbc))]);
                set(pl,'color',cols(e,:))             
            end;
            set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
            % set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);
            
            xlabel(['Winv ',int2str(selfacts(cb))]); ylabel(['Winv ',int2str(selfacts(bc))]); zlabel(['Winv ',int2str(selfacts(cbc))]); 
            pp = pp+1;
        end;
    end;
end;axcopy
set(gcf,'color','w');
figure;plot3(clustact(3,:),clustact(5,:),clustact(6,:),'.k');hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot ICA activations (subj*factors) color coded by spectral clusters:
load /data/common2/emotion/clusters/KmeanClustersGdSubjs.mat
load /data/common2/emotion/clusters/KmeanClustersGdSubjs10.mat

figure;plot(activations(3,:),clustact(3,:),'.r');hold on;
 plot(activations(,:)/std(clustact(5,:)),'.b');hold on;       
 plot(activations(3,:)/std(clustact(6,:)),'.g');hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   find factor contributions
PCAdim = 3;
normact = clustact(PCAdim,:)/std(clustact(PCAdim,:));

sigfacts = find(normact > 1 | normact < -1); % chooses factors about z score of 1
subjfactors = cell(1,length(gdcomps));  
for nx = 1:length(button)
    subjfactors{button(nx)} = sigfacts(find(sigfacts>(nx-1)*10+1 & sigfacts<nx*10))-(nx-1)*10;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now plot dipoles, spectral templates and emotion weightings for only clustered factors
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
qt = 0;wt = 1; clear dipcols  mnemodiff subjactcell actcell dipsource subjwts subjactcell allnames allidx
for nx = 1:length(gdcomps)
    if ~isempty(subjfactors{nx})
        cpcols{nx} = hsv(length(gdcomps{nx})); clear actcell factwts factname actidx
        EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
        tw = 1;   
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{button(nx)}],[pcadims{button(nx)} subjdims{button(nx)}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
        for e = 2:length(numtrials{nx})+1
            emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
        end;
        if wt == 1 
            dipsource{nx,1}(1) = EEG.dipfit.model(1);
        end;
    end; 
    for tp = 1:length(subjfactors{nx})
        for cp = 1:length(nxlists{nx}{subjfactors{nx}(tp)})
            rcp = find(nxlists{nx}{subjfactors{nx}(tp)}(cp) == gdcomps{nx});
            dipsource{nx,tw}(cp) = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));
            dipcols{nx,tw}{cp} = cpcols{nx}(rcp,:);
            actidx{tw}(cp) = rcp;
        end;        
        
        tpwts =  winv(:,subjfactors{nx}(tp))'; 
        na = activations(subjfactors{nx}(tp),:);  % had flipped based on mean tpwts 
        actcell{tw} = na;
        factname{tw} = subjfactors{nx}(tp);
        factwts{tw} = tpwts;tw = tw+1; wt = wt+1;
        if wt == 11
            qt = 1;
            break;break
        end;    
    end;
    if ~isempty(subjfactors{nx})
        subjactcell{nx} = actcell;
        subjwts{nx} = factwts;
        allnames{nx} = factname;
        allidx{nx} = actidx;
    end;
    if qt == 1
        break;break;
    end;    
    ALLEEG=[];EEG=[];
end; 
tl = ['Subject Factors contributing to Dimension ',int2str(PCAdim),' of ICA clustering on EmoScores'];
row = 5; col = 6;     
figure; op = 7%28;
for nx = size(dipsource,1):-1:1
    for tw = size(dipsource,2):-1:1
        if ~isempty(dipsource{nx,tw})
        subplot(row,col,op);       
        dipplot(dipsource{nx,tw},'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',dipcols{nx,tw},'projlines','on','projimg','on'); op = op-3; view(60,20); 
        end;
    end;
end;
p=2;
for nx = 1:length(subjactcell)
    actcell = subjactcell{nx};
    if ~isempty( subjactcell{nx})
    for tw = 1:length(actcell)
        if ~isempty(actcell{tw})
            subplot(row,col,p)        
            mxcp = 0;         mncp = 0;  
            for cp = 1:length(allidx{nx}{tw})
                rcp = allidx{nx}{tw}(cp);
                %rcp = find(nxlists{nx}{tw}(cp) == gdcomps{nx});
                ph=plot(freqs,actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp));hold on; 
                set(ph,'color',cpcols{nx}(rcp,:)); set(ph,'linewidth',1.5); 
                mxcp1 = max(actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                mncp1 = min(actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                if mxcp1 > mxcp
                    mxcp = mxcp1;
                end;
                if mncp1< mncp
                    mncp = mncp1;
                end;        
            end;p = p+1;
            set(gca,'xlim',[0 40]);     set(gca,'ylim',[mncp+mncp*.01 mxcp+mxcp*.01]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on'); set(gca,'fontsize',10);    
            set(gca,'ycolor','c');
            subplot(row,col,p);             
            clear basemat basemat1
            % need to calculate new emomap for each subject
            for e = 2:length(numtrials{nx})+1
                emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
            end;          
            emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
            for e = 1:length(numtrials{nx})-1
                e=emoorder(e);clear newmat
                tempmat = subjwts{nx}{tw}(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
                if e == 1
                    pl=1; for pk = .1:.1:.9
                        basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat1,basemat1-basemat1,'.');hold on;
                    set(ph,'color',[.5 .5 .5]);        set(ph,'markersize',16);
                    pl=1; 
                    for pk = .03:.01:.98
                        basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat,basemat-basemat);hold on;
                    handvec(1,e) = ph;
                    set(ph,'color',[.5 .5 .5]);        set(ph,'linewidth',1.5); 
                    mx = max(basemat-basemat); mn = min(basemat-basemat);
                else
                    pl=1;
                    for pk = .1:.1:.9
                        newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat1,newmat1-basemat1,'.');hold on;
                    set(ph,'color',cols(e,:));        set(ph,'markersize',16);pl=1;
                    for pk = .03:.01:.98
                        newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat,newmat-basemat);hold on;
                    handvec(1,e) = ph;
                    set(ph,'color',cols(e,:));         set(ph,'linewidth',1.5); 
                    mx1 = max(newmat-basemat); mn1 = min(newmat-basemat);
                    if mx1 > mx
                        mx = mx1;
                    end;
                    if mn1 < mn
                        mn = mn1;
                    end;                
                end; 
            end;
            set(gca,'xlim',[basemat(1)-basemat(1)*.02 basemat(end)+basemat(end)*.02]);
            set(gca,'box','off'); set(gca,'fontsize',10); 
            set(gca,'ylim',[mn+mn*.01 mx+mx*.01]);   p = p+2;
            title(['Subj ',int2str(nx),'; Fact ',int2str(allnames{nx}{tw})]);
        end; set(gcf,'color','w');    
    end;
    end;    
end;
ph = textsc(tl,'title'); set(ph,'color','r');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

print  -dpsc2 -Pcoloring 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   find factor contributions
emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excited'};
load /data/common2/emotion/clusters/GrandICAonScores.mat clustact clustwinv  pcamatall
PCAdim = 5;
normact = clustact(PCAdim,:)/std(clustact(PCAdim,:));

sigfacts = find(normact > 1 | normact < -1); % chooses factors about z score of 1
subjfactors = cell(1,length(gdcomps));  
for nx = 1:length(button)
    subjfactors{button(nx)} = sigfacts(find(sigfacts>(nx-1)*10+1 & sigfacts<nx*10))-(nx-1)*10;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;pl = 1; fact = 1;page = 1; bd = 1;gd = 1; pd = 1; dp=1;clear emogood emobad  dplistgd  dplistbd
for nx = 1:length(subjfactors)
    if ~isempty(subjfactors{nx})
        %if fact == 14
        %    set(gcf,'PaperOrientation','landscape');
        %    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        %    set(gcf,'color','w'); axcopy
        %    textsc(['From Dimension ',int2str(PCAdim),'. Scaled contributions to each emotion (Clustact*Clustwinv values * spectra). Page: ',int2str(page)],'title');
            %print  -dpsc2 -Pcoloring 
       %     figure; pl = 1;fact = 1;page = page+1;
       % end;       
        cpcols{nx} = hsv(length(gdcomps{nx})); clear actcell factwts factname actidx
        EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
        sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
        wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
        ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
        for e = 2:length(numtrials{nx})+1
            emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
        end;
        %if nx == 1 
       %     dipsource{nx,1}(1) = EEG.dipfit.model(1);
        %end;
        for tp = 1:length(subjfactors{nx})
            for cp = 1:length(nxlists{nx}{subjfactors{nx}(tp)})
                rcp = find(nxlists{nx}{subjfactors{nx}(tp)}(cp) == gdcomps{nx});
                dipsource{nx,tp}(cp) = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));
                dipcols{nx,tp}{cp} = cpcols{nx}(rcp,:);
                actidx{tp}(cp) = rcp;
            end;      
        end;    
        for tp = 1:length(subjfactors{nx})
        %if fact == 14
        %    set(gcf,'PaperOrientation','landscape');
        %    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
       %     set(gcf,'color','w'); axcopy
        %    textsc(['From Dimension ',int2str(PCAdim),'. Scaled contributions to each emotion (Clustact*Clustwinv values * spectra). Page: ',int2str(page)],'title');
            %print  -dpsc2 -Pcoloring 
        %    figure; pl = 1;fact = 1;page = page+1;
        %end;       
            % clear yl nyl
            %subplot(13,16,pl)
            %dipplot(dipsource{nx,tp},'image','mri','gui','off','normlen','on','dipolesize',42,'dipolelength',0,'spheres','on','color',dipcols{nx,tp},'projlines','on','projimg','on'); 
            %camzoom(.95);view(60,20);pl=pl+1;
            %title(['Subject: ',int2str(nx),' Factor: ',int2str(subjfactors{nx}(tp))]);
            %for em = 1:size(clustwinv,1)
                %subplot(13,16,pl);set(gca,'fontsize',7);
                for cp = 1:length(nxlists{nx}{subjfactors{nx}(tp)})
                    rcp = find(nxlists{nx}{subjfactors{nx}(tp)}(cp) == gdcomps{nx});
                    emoact = activations(subjfactors{nx}(tp),(length(freqs)*(rcp-1)+1:length(freqs)*rcp));
                    %emoact = emoact*(clustwinv(em,PCAdim)*clustact(PCAdim,((find(button==nx))-1)*15+subjfactors{nx}(tp)));
                    emoact = emoact*clustact(PCAdim,((find(button==nx))-1)*15+subjfactors{nx}(tp))*-1;
                    %if ismember(em,[8,9,10,12,14,15]) 
                    %    emogood(gd,:) = emoact;gd = gd+1;
             %dplistgd{dp} = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));dp= dp+1;
                    %elseif ismember(em,[2,3,5])
                        emobad(bd,:) = emoact; bd = bd+1;
             dplistbd{pd} = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));pd= pd+1;
             
             
                    %end;                    
                    %ph = plot(freqs,emoact,'k-'); hold on;
                    %set(ph,'color',cpcols{nx}(rcp,:));set(ph,'linewidth',.75);                    
                    %if cp == 1
                    %    yyl = [min(emoact) max(emoact)];
                    %end;                    
                    %tmpyl = [min(emoact) max(emoact)];
                    %if tmpyl(1) < yyl(1)
                    %    yyl(1) = tmpyl(1);
                    %end;
                    %if tmpyl(2) > yyl(2)
                    %    yyl(2) = tmpyl(2);
                    %end;                    
                    %set(gca,'ylim',[yyl(1) yyl(2)]);
                end;
                %if em == 1
                %    yl = get(gca,'ylim'); 
                %else
                %    nyl= get(gca,'ylim'); 
                %    if nyl(1)< yl(1)
                %        yl(1) = nyl(1);
                %    end;
                %    if nyl(2)> yl(2)
                %        yl(2) = nyl(2);
                %    end;
                %end;
                %set(gca,'xtick',[10:20:50]);
                %if fact == 1
                %title(emo2{em});
                %end;
                %set(gca,'xgrid','on');set(gca,'xlim',[0 50]);pl=pl+1;
            %end;
            %for sb = pl-15:pl-1
            %    subplot(13,16,sb); set(gca,'ylim',[yl(1) yl(2)]);
            %end; 
            %fact = fact+1;
        end; % to tp
        ALLEEG=[];EEG=[];
    end; % to if subjfactors ~ isempty
    fprintf('\n One More SUBJECT Done: %i',nx);
end;
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color','w');axcopy
textsc(['From Dimension ',int2str(PCAdim),'. Scaled contributions to each emotion (Clustact*Clustwinv values * spectra). Page: ',int2str(page)],'title'); 
print  -dpsc2 -Pcoloring 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k means cluster good and bad spectra:
reps = 10;  % currently good is comp,love,relief,awe,joy,excite; bad = frust, jeal.,disgust
fr = find(freqs> 0 & freqs < 30);
egood = emogood;
for dip = 1:size(emogood,1)
egood(dip,:) = egood(dip,:)/max(abs(egood(dip,:)));
end;

[optkgood,centr,clst,Cg] =Kmeangrp(egood,22,reps,1); % find optimal cluster number, don't plot
[koutgood, C,sumd, allsums] = kmeans(egood,optkgood,'replicates',reps); % 10 clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster by dipoles instead %%%%%%
locspec = egood;locspec(:,20:22) = zeros(72,3);
for dip = 1:length(dplistgd)
    subjwinv = dplistgd{dip}.posxyz;
    if size(subjwinv,1) > 1
        subjwinv(1,end+1:end+size(subjwinv,2)) = subjwinv(2,:); subjwinv(2,:) = [];
    else
        subjwinv(1,end+1:end+3) = zeros(1,3);
    end;   
    locspec(dip,:) = locspec(dip,:)/max(abs(locspec(dip,:)));
    %subjwinv = subjwinv(1,1:3)/max(abs(subjwinv(1,1:3)));
    locspec(dip,20:22) = subjwinv(1,1:3);
end;
[optkgood,centr,clst,Cg] =Kmeangrp(allloc,6,reps,1); % find optimal cluster number, don't plot
[koutgood, C,sumd, allsums] = kmeans(locspec,optkgood,'replicates',25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for cl = 1:optkgood
    oneclust = find(koutgood == cl);
    fout = allsums(oneclust)/std(allsums(oneclust));
    outliers = find(fout > 6);
    koutgood(oneclust(outliers)) = 0;
end;
figure; col = 6;  pl = 1;  clustcol = jet(optkgood);row = 3;%round(optkgood/3);
for cl = 1:optkgood
    oneclust = find(koutgood == cl);
    if ~isempty(oneclust)
        plotdip = dplistgd{oneclust(1)};
        for w=  2:length(oneclust)
            plotdip(end+1) = dplistgd{oneclust(w)};
        end;
        subplot(row,col,pl)
        dipplot(plotdip,'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','projlines','on','projimg','on','color',{clustcol(7-cl,:)}); %camzoom(1.3)
        view(60,20); pl = pl+1;
        subplot(row,col,pl);
        x=mean(egood(oneclust,:),1); ph=plot(freqs,x);hold on;set(ph,'color',clustcol(7-cl,:));pl = pl+1; 
        set(ph,'linewidth',1.75);set(gca,'xtick',[10:10:50]);set(gca,'xlim',[1 50]); set(gca,'xgrid','on');
        title(['Cluster ',int2str(cl)]);
        subplot(row,col,pl); plall = hsv(length(oneclust));
        for w = 1:length(oneclust)
            ph=plot(freqs,egood(oneclust(w),:)');hold on; 
            set(ph,'color',plall(w,:)); set(ph,'linewidth',1);
        end;pl = pl+1;    set(gca,'xtick',[10:10:50]);    
        set(gca,'xlim',[1 50]); set(gca,'xgrid','on');    title(['Cluster ',int2str(cl)]);
    end;
end;axcopy
ph = textsc('Positive Emotion Spectral Change Clusters','title');set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color','w');axcopy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ebad = emobad;
for dip = 1:size(emobad,1)
ebad(dip,:) = ebad(dip,:)/max(abs(ebad(dip,:)));
end;
[optkbad,centr,clst,Cg] =Kmeangrp(emobad(:,fr),15,reps,1); % find optimal cluster number, don't plot
[koutbad, C,sumd, allsums] = kmeans(emobad(:,fr),optkbad,'replicates',reps); % 6 clusters
for cl = 1:optkbad
    oneclust = find(koutbad == cl);
    fout = allsums(oneclust)/std(allsums(oneclust));
    outliers = find(fout > 5);
    koutbad(oneclust(outliers)) = 0;
end;
figure; col = 6; pl = 1;  clustcol = jet(optkbad);row = 3;%optkbad/3; 
for cl = 1:optkbad
    oneclust = find(koutgood == cl);
    plotdip = dplistbd{oneclust(1)};
    for w=  2:length(oneclust)
        plotdip(end+1) = dplistbd{oneclust(w)};
    end;
    subplot(row,col,pl)
    dipplot(plotdip,'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','projlines','on','projimg','on','color',{clustcol(7-cl,:)}); %camzoom(1.3)
    view(60,20); pl = pl+1;
    subplot(row,col,pl); 
    x=mean(ebad(oneclust,:),1); ph=plot(freqs,x);hold on;set(ph,'color',clustcol(7-cl,:));pl = pl+1;   
    set(ph,'linewidth',1.75);set(gca,'xtick',[10:10:50]);set(gca,'xlim',[1 50]); set(gca,'xgrid','on');
    title(['Cluster ',int2str(cl)]);
    subplot(row,col,pl); plall = hsv(length(oneclust));
    for w = 1:length(oneclust)
        ph=plot(freqs,ebad(oneclust(w),:)');hold on; 
        set(ph,'color',plall(w,:)); set(ph,'linewidth',1);
    end;pl = pl+1;    set(gca,'xtick',[10:10:50]);
    set(gca,'xlim',[1 50]); set(gca,'xgrid','on');
    title(['Cluster ',int2str(cl)]);
end;axcopy
ph = textsc('Negative Emotion Spectral Change Clusters','title');set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color','w');axcopy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clustact = activations; clustwinv= winv; 
% back project only dimensions of interest
load /data/common2/emotion/clusters/GrandICAonScores.mat clustact clustwinv  pcamatall
% dimensions 3 and 5 separate emotions
pcadims = cell2mat(pcadims);
figure;
for nx = 1:length(button)
    onesubj=clustact(:,sum(pcadims(button(1:nx-1)))+1:sum(pcadims(button(1:nx)))); 
    onesubj([2:3,4,5:10],:) = 0;
    %onesubj([3,5,4,6],:) = 0; % looks like nothing at all!!
    backproj = clustwinv*onesubj;

    [weights,sphere,compvars,bias,signs,lrates,activations] = runica(backproj,'pca',1,'extended',1);
    ws = weights*sphere;  winv = pinv(ws);cd (['/data/common2/emotion/',paths{button(nx)}])
    %save BackProjICA.mat winv activations
    emo2 = {' anger',' frustration',' jealousy',' fear' ,' disgust',' grief',' sad',' compassion',' love',' relief',' content',' awe',' happy',' joy',' excited'};
    cols = jet(15);
    row = 2; col = 2;pp = 1;
    subplot(ceil(sqrt(length(button))),round(sqrt(length(button))),nx)
    for e = 1:size(winv,1)
        ph =plot(e,winv(e,1),'.');hold on;
        %ph =plot(winv(e,1),winv(e,2),'.');hold on;
        set(ph,'markersize',20);         set(ph,'color',cols(e,:));
        ph = text(e,winv(e,1),emo2{e});
        %ph = text(winv(e,1),winv(e,2),emo2{e});
        set(ph,'color',cols(e,:)); 
    end;
    mx =  max(winv(:,1));mnx = min(winv(:,1));
    %my =  max(winv(:,2));mny = min(winv(:,2));
    %set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    %for e = 1:size(winv,1)
    %    pl =plot([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)]);
    %    set(pl,'color',cols(e,:))             
    %end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    % set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);

    xlabel('Winv 1'); ylabel('Winv 2'); 
    title(['Subj ',int2str(nx)]);
end;
axcopy;set(gcf,'color','w');
textsc('Back Projection of only dimensions 1 from Grand ICA to data of each subj individually','title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

print  -dpsc2 -Pcoloring 

% from saved data
figure;
for nx = 1:length(button)
    cd (['/data/common2/emotion/',paths{button(nx)}]);
    load BackProjICA.mat winv activations
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
    cols = jet(15);
    row = 2; col = 2;pp = 1;
    subplot(4,5,nx);set(gca,'fontsize',8);
    for e = 1:size(winv,1)
        ph =plot(winv(e,1),winv(e,2),'.');hold on;
        set(ph,'markersize',20);         set(ph,'color',cols(e,:));
        ph = text(winv(e,1),winv(e,2),emo2{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',8);
    end;
    mx =  max(winv(:,1));mnx = min(winv(:,1));
    my =  max(winv(:,2));mny = min(winv(:,2));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    for e = 1:size(winv,1)
        pl =plot([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)]);
        set(pl,'color',cols(e,:));             
    end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    % set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);

    %xlabel('Winv 1'); ylabel('Winv 2'); 
    title(['Subj ',int2str(nx)]);
end;
subplot(4,5,nx+1);set(gca,'fontsize',8);
    for e = 1:size(clustwinv,1)
        ph =plot(clustwinv(e,3),clustwinv(e,5),'.');hold on;
        set(ph,'markersize',20);         set(ph,'color',cols(e,:));
        ph = text(clustwinv(e,3),clustwinv(e,5),emo2{e});
        set(ph,'color',cols(e,:));set(ph,'fontsize',8); 
    end;
    mx =  max(clustwinv(:,3));mnx = min(clustwinv(:,3));
    my =  max(clustwinv(:,5));mny = min(clustwinv(:,5));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    for e = 1:size(clustwinv,1)
        pl =plot([clustwinv(e,3) clustwinv(e,3)],[clustwinv(e,5) clustwinv(e,5)]);
        set(pl,'color',cols(e,:))             
     end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    xlabel('Winv 3'); ylabel('Winv 5'); 
    title('Projected Dimensions');

axcopy;set(gcf,'color','w');
textsc('Back Projection of only dimensions 3 and 5 from Grand ICA to data of each subj individually','title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);


% look at each subj with rotated axes and orientations
nx=8;% subj2
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
 winv(:,2) = winv(:,2)*-1;
 winv(:,1) = winv(:,1)*-1;
 figure;
 for e = 1:size(winv,1)
     ph =plot(winv(e,1),winv(e,2),'.');hold on;
     set(ph,'markersize',20);         set(ph,'color',cols(e,:));
     ph = text(winv(e,1),winv(e,2),emo2{e});
     set(ph,'color',cols(e,:)); set(ph,'fontsize',8);
 end;
 for e = 1:size(winv,1)
     pl =plot([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)]);
     set(pl,'color',cols(e,:));             
 end;

%%%%  correction of winvs to be like projected axes
nx=1;% subj1
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,2) = winv(:,2)*-1;
 activations(2,:) = activations(2,:)*-1;
 save BackProjICA.mat winv activations
nx=2;% subj2
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
%w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
a1 = activations(1,:);a2 = activations(2,:);activations(1,:) = a2;activations(2,:) = a1;
 %winv(:,2) = winv(:,2)*-1; winv(:,1) = winv(:,1)*-1;
 activations(2,:) = activations(2,:)*-1; activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=3;% subj3
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 save BackProjICA.mat winv activations
nx=4;% subj4
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,1) = winv(:,1)*-1;
 activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=5;% subj5
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=6;% subj6
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,1) = winv(:,1)*-1;
 activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=7;% subj7
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
%w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
a1 = activations(1,:);a2 = activations(2,:);activations(1,:) = a2;activations(2,:) = a1;
 save BackProjICA.mat winv activations
nx=8;% subj8
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
save BackProjICA.mat winv activations
nx=9;% subj9
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 save BackProjICA.mat winv activations
nx=10;% subj10
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
save BackProjICA.mat winv activations
nx=11;% subj11
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,2) = winv(:,2)*-1;
activations(2,:) = activations(2,:)*-1;
save BackProjICA.mat winv activations
nx=12;% subj12
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
save BackProjICA.mat winv activations
nx=13;% subj13
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 save BackProjICA.mat winv activations
nx=14;% subj14
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
%w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
a1 = activations(1,:);a2 = activations(2,:);activations(1,:) = a2;activations(2,:) = a1;
 %winv(:,2) = winv(:,2)*-1;
activations(2,:) = activations(2,:)*-1;
save BackProjICA.mat winv activations
nx=15;% subj15
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,2) = winv(:,2)*-1; winv(:,1) = winv(:,1)*-1;
 activations(2,:) = activations(2,:)*-1; activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=16;% subj16
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
a1 = activations(1,:);a2 = activations(2,:);activations(1,:) = a2;activations(2,:) = a1;
 winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=17;% subj17
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 save BackProjICA.mat winv activations
nx=18;% subj18
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
 %winv(:,2) = winv(:,2)*-1;winv(:,1) = winv(:,1)*-1;
 activations(2,:) = activations(2,:)*-1; activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
nx=19;% subj19
cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat 
w1 = winv(:,1); w2 = winv(:,2);winv(:,2) = w1; winv(:,1) = w2;
a1 = activations(1,:);a2 = activations(2,:);activations(1,:) = a2;activations(2,:) = a1;
 winv(:,1) = winv(:,1)*-1;
activations(1,:) = activations(1,:)*-1;
 save BackProjICA.mat winv activations
%%%%%%%%%%%%%%%%%%%%%%%%%%
% from saved data
figure;
for nx = 1:length(button)
    cd (['/data/common2/emotion/',paths{button(nx)}]);
    load BackProjICA.mat winv activations
    emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excited'};
    
    % to plot neg version to interpret neg activations
    %winv(:,1) = winv(:,1)*-1;winv(:,2) = winv(:,2)*-1;
    
    cols = jet(15);
    row = 2; col = 2;pp = 1;
    subplot(4,5,nx);set(gca,'fontsize',8);
    for e = 1:size(winv,1)
        ph =plot(winv(e,1),winv(e,2),'.');hold on;
        set(ph,'markersize',20);         set(ph,'color',cols(e,:));
        ph = text(winv(e,1),winv(e,2),emo2{e});
        set(ph,'color',cols(e,:)); set(ph,'fontsize',8);
    end;
    mx =  max(winv(:,1));mnx = min(winv(:,1));
    my =  max(winv(:,2));mny = min(winv(:,2));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    for e = 1:size(winv,1)
        pl =plot([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)]);
        set(pl,'color',cols(e,:));             
    end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    % set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);

    %xlabel('Winv 1'); ylabel('Winv 2'); 
    title(['Subj ',int2str(nx)]);
end;
subplot(4,5,nx+1);set(gca,'fontsize',8);
    for e = 1:size(winv,1)
        ph =plot(clustwinv(e,3),clustwinv(e,5),'.');hold on;
        set(ph,'markersize',20);         set(ph,'color',cols(e,:));
        ph = text(clustwinv(e,3),clustwinv(e,5),emo2{e});
        set(ph,'color',cols(e,:));set(ph,'fontsize',8); 
    end;
    mx =  max(clustwinv(:,3));mnx = min(clustwinv(:,3));
    my =  max(clustwinv(:,5));mny = min(clustwinv(:,5));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    for e = 1:size(winv,1)
        pl =plot([clustwinv(e,3) clustwinv(e,3)],[clustwinv(e,5) clustwinv(e,5)]);
        set(pl,'color',cols(e,:))             
     end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    xlabel('Winv 3'); ylabel('Winv 5'); 
    title('Projected Dimensions');

axcopy;set(gcf,'color','w');
textsc('Back Projection of only dimensions 3 and 5 from Grand ICA to data of each subj individually','title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
 
% find distances of each emotion (in each subject) from (0,0)
clear dist subjdist
for nx = 1:length(button)
    cd (['/data/common2/emotion/',paths{button(nx)}]);
    load BackProjICA.mat winv activations
    
    % to plot neg version to interpret neg activations
    %winv(:,1) = winv(:,1)*-1;winv(:,2) = winv(:,2)*-1;

    for  em = 1:size(winv,1)
        dist(em) = sqrt(winv(em,1)^2 + winv(em,2)^2);
    end;
    subjdist(:,nx) = dist';
end;
% normalize each subject based on mean
for nx = 1:size(subjdist,2)
    fact(nx) = mean(subjdist(:,nx));
    subjdist(:,nx) = subjdist(:,nx)/fact(nx);
end;
clear epnts subjepnts
for nx = 1:length(button)
    cd (['/data/common2/emotion/',paths{button(nx)}]);
    load BackProjICA.mat winv activations
    
    % to plot neg version to interpret neg activations
    %winv(:,1) = winv(:,1)*-1;winv(:,2) = winv(:,2)*-1;

    for  em = 1:size(winv,1)
        epnts(em,:) = [winv(em,1)/fact(nx),winv(em,2)/fact(nx)];
    end;
    subjepnts(:,:,nx) = epnts;
end;

spiralorder = [5,3,2,1,4,6,12,8,14,9,15,10,13,11,7];
figure;
for em=1:size(subjepnts,1)
    subplot(4,4,em);set(gca,'fontsize',16);
for nx = 1:size(subjepnts,3)
    ph = plot(subjepnts(spiralorder(em),1,nx),subjepnts(spiralorder(em),2,nx),'k.');hold on;
         set(ph,'color',cols(spiralorder(em),:)); set(ph,'markersize',30);   
         set(gca,'ylim',[-2 1.5]); set(gca,'xlim',[-2.3 2.3]);
         plot([0 0],[get(gca,'ylim')],'k');plot([get(gca,'xlim')],[0 0],'k');
         set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
end;
title(emo2{spiralorder(em)});
end;
subplot(4,4,em+1)
for em=1:size(subjepnts,1)
for nx = 1:size(subjepnts,3)
    ph = plot(subjepnts(spiralorder(em),1,nx),subjepnts(spiralorder(em),2,nx),'k.');hold on;
         set(ph,'color',cols(spiralorder(em),:)); set(ph,'markersize',18);           
end;
end;
         plot([0 0],[get(gca,'ylim')],'k');plot([get(gca,'xlim')],[0 0],'k');
         set(gca,'ylim',[-2 1.5]); set(gca,'xlim',[-2.4 2.3]);
title('All Emotions');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
textsc('Back-projection dimensions from all subjects; axes and polarity reversed when necessary to bring all subjects into alignment; Distances normalized within subj by mean distance from (0,0); WINVs ARE REVERSED TO SHOW NEGATIVE WEIGHTINGS','title');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now find factors contributing to the two dimensions
%%   find factor contributions
PCAdim = 2;clear subjfactors
for nx = 1:length(button)
    cd (['/data/common2/emotion/',paths{button(nx)}]);load  BackProjICA.mat     
    normact = activations(PCAdim,:)/std(activations(PCAdim,:));    
    sigfacts = find(normact < -1.5); % chooses factors about z score of 1
    subjfactors{button(nx)} = sigfacts;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now plot dipoles, spectral templates and emotion weightings for only clustered factors
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; clear cols
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
qt = 0;wt = 1; clear dipcols  mnemodiff subjactcell actcell dipsource subjwts subjactcell allnames allidx
for nx = 22:length(gdcomps)
    if ~isempty(subjfactors{nx})
        cpcols{nx} = hsv(length(gdcomps{nx})); clear actcell factwts factname actidx
        EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
        tw = 1;   
        sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
        wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
        ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
        for e = 2:length(numtrials{nx})+1
            emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
        end;
        if wt == 1 
            dipsource{nx,1}(1) = EEG.dipfit.model(1);
        end;
    end; 
    for tp = 1:length(subjfactors{nx})
        for cp = 1:length(nxlists{nx}{subjfactors{nx}(tp)})
            rcp = find(nxlists{nx}{subjfactors{nx}(tp)}(cp) == gdcomps{nx});
            dipsource{nx,tw}(cp) = EEG.dipfit.model(nxlists{nx}{subjfactors{nx}(tp)}(cp));
            dipcols{nx,tw}{cp} = cpcols{nx}(rcp,:);
            actidx{tw}(cp) = rcp;
        end;        
        
        tpwts =  winv(:,subjfactors{nx}(tp))'; 
        if mean(tpwts)<0
            tpwts = tpwts*-1;
            na = activations(subjfactors{nx}(tp),:)*-1;
        else
            na = activations(subjfactors{nx}(tp),:);  
        end;
        actcell{tw} = na;
        factname{tw} = subjfactors{nx}(tp);
        factwts{tw} = tpwts;tw = tw+1; wt = wt+1;
        if wt == 11
            qt = 1;
            break;break
        end;    
    end;
    if ~isempty(subjfactors{nx})
        subjactcell{nx} = actcell;
        subjwts{nx} = factwts;
        allnames{nx} = factname;
        allidx{nx} = actidx;
    end;
    if qt == 1
        break;break;
    end;    
    ALLEEG=[];EEG=[];
end; 
tl = ['Subject Factors contributing to Dimension ',int2str(PCAdim),' of ICA clustering on EmoScores'];
row = 5; col = 6;     
figure; op = 10%28;
for nx = size(dipsource,1):-1:1
    for tw = size(dipsource,2):-1:1
        if ~isempty(dipsource{nx,tw})
        subplot(row,col,op);       
        dipplot(dipsource{nx,tw},'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',dipcols{nx,tw},'projlines','on','projimg','on'); op = op-3; view(60,20); 
        end;
    end;
end;
p=2;
for nx = 1:length(subjactcell)
    actcell = subjactcell{nx};
    if ~isempty( subjactcell{nx})
    for tw = 1:length(actcell)
        if ~isempty(actcell{tw})
            subplot(row,col,p)        
            mxcp = 0;         mncp = 0;  
            for cp = 1:length(allidx{nx}{tw})
                rcp = allidx{nx}{tw}(cp);
                %rcp = find(nxlists{nx}{tw}(cp) == gdcomps{nx});
                ph=plot(freqs,actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp));hold on; 
                set(ph,'color',cpcols{nx}(rcp,:)); set(ph,'linewidth',1.5); 
                mxcp1 = max(actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                mncp1 = min(actcell{tw}(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                if mxcp1 > mxcp
                    mxcp = mxcp1;
                end;
                if mncp1< mncp
                    mncp = mncp1;
                end;        
            end;p = p+1;
            set(gca,'xlim',[0 40]);     set(gca,'ylim',[mncp+mncp*.01 mxcp+mxcp*.01]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on'); set(gca,'fontsize',10);    
            set(gca,'ycolor','c');
            subplot(row,col,p);             
            clear basemat basemat1
            % need to calculate new emomap for each subject
            for e = 2:length(numtrials{nx})+1
                emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
            end;          
            emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
            for e = 1:length(numtrials{nx})-1
                e=emoorder(e);clear newmat
                tempmat = subjwts{nx}{tw}(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
                if e == 1
                    pl=1; for pk = .1:.1:.9
                        basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat1,basemat1-basemat1,'.');hold on;
                    set(ph,'color',[.5 .5 .5]);        set(ph,'markersize',16);
                    pl=1; 
                    for pk = .03:.01:.98
                        basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat,basemat-basemat);hold on;
                    handvec(1,e) = ph;
                    set(ph,'color',[.5 .5 .5]);        set(ph,'linewidth',1.5); 
                    mx = max(basemat-basemat); mn = min(basemat-basemat);
                else
                    pl=1;
                    for pk = .1:.1:.9
                        newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat1,newmat1-basemat1,'.');hold on;
                    set(ph,'color',cols(e,:));        set(ph,'markersize',16);pl=1;
                    for pk = .03:.01:.98
                        newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                    end;
                    ph = plot(basemat,newmat-basemat);hold on;
                    handvec(1,e) = ph;
                    set(ph,'color',cols(e,:));         set(ph,'linewidth',1.5); 
                    mx1 = max(newmat-basemat); mn1 = min(newmat-basemat);
                    if mx1 > mx
                        mx = mx1;
                    end;
                    if mn1 < mn
                        mn = mn1;
                    end;                
                end; 
            end;
            set(gca,'xlim',[basemat(1)-basemat(1)*.02 basemat(end)+basemat(end)*.02]);
            set(gca,'box','off'); set(gca,'fontsize',10); 
            set(gca,'ylim',[mn+mn*.01 mx+mx*.01]);   p = p+2;
            title(['Subj ',int2str(nx),'; Fact ',int2str(allnames{nx}{tw})]);
        end; set(gcf,'color','w');    
    end;
    end;    
end;
ph = textsc(tl,'title'); set(ph,'color','r');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

print  -dpsc2 -Pcoloring 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
