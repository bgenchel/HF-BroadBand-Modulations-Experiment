% Summary scripts for Co-Mod Paper figures
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
gdcomps{36}=[];
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
savedat = 'SpecCoMod'; fullpaths = newpaths;
savedat = 'SpecCoModMuscle'; 
load /data/common1/emotion/AllClustFacs.mat    

strs = {
    'load /data/common1/emotion/DeltaClust.mat',
    'load /data/common1/emotion/ThetaClust.mat',
    'load /data/common1/emotion/AlphaClust.mat',
    'load /data/common1/emotion/BetaClust.mat',
    'load /data/common1/emotion/GammaClust.mat'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IM movie  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp = [2,3,4,5,7,10,12,13,15]; frqlim = [3 125]; nx=2;
im = useims{nx}; % will plot all ims back projected (when 'weights is added)
dim = 1;
weights = fullwts(dim,:);weights = weights/max(weights);
datset = 'sources.set';

for em = 2:length(emos)
  [M] = IMmovie(datset,savedat,fullpaths{nx},im,comp,frqlim,['Cp',int2str(comp),'EmSpc',int2str(dim),emos{em},'.avi'],[em],weights','off',{'r'},{emos{em}});
  close
end;
    
figure;movie(M)

SpecCoModPlot('sources.set',fullpaths{nx},[],[1:15],savedat,[3 128],'n',0,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  FIGURE 2-3 (IM examples)%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frqlim = [3 125]; 
nx = 10;
tmpls = [26,2,7,35,34,1,5,32,12,18,21,39,37,6];
comps = [3,6,19,20,11];
% plot template examples:
SpecCoModPlot('sources.set',fullpaths{nx},comps,tmpls,savedat,frqlim,'n',0,[],[]);
str = ['print /home/julie/Manuscripts/Emotion/figures/nx',int2str(nx),'TemplExamples.eps -depsc -painters -adobe'];eval(str)

% plot mean spectra with back projections
PlotSpecFacEnv('sources.set',savedat,fullpaths{nx},tmpls,comps,[],1,frqlim,1,1,0);
ph=textsc(['Subject ',int2str(nx)],'title');set(ph,'fontsize',16);
str = ['print /home/julie/Manuscripts/Emotion/figures/nx',int2str(nx),'BackProjsPos.eps -depsc -painters'];eval(str)
PlotSpecFacEnv('sources.set',savedat,newpaths{nx},tmpls,comps,[],[],frqlim,1,.99,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  FIGURE 5-7 (clusters)  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colrs = {[0 0 1],[.5 0 .7],[0 1 0],[1 0 0],[0 0 1]};
colrsbk = {[1 .3 1],[1 .6 0],[1 .2 .2],[.5 .5 1],[.2 1 .2]};
for s = 1:length(strs)
    eval(strs{s})
    figure;row=2; col=2;    
    for clust = 1:length(finaltempls)
        sbplot(row,col,clust)
        plot([0 50],[0 0],'k-'); hold on;
        [ph fx lx] = quadplot(freqs,finaltempls{clust}',1,'k');hold on;
        lightcol = colrsbk{s};
        set(ph,'color',lightcol);
        [ph fx lx] = quadplot(freqs,mean(finaltempls{clust},1),2.5,'k');
        set(ph,'color',colrs{s});
        set(gca,'ytick',[-10:5:15]);
        set(gca,'ylim',[min(finaltempls{clust}(:))-.5 max(finaltempls{clust}(:))+.5]);
        set(gca,'ticklength',[.03 .03]);
        if s == length(strs)
            plot([get(gca,'xlim')],[0 0],'k-'); hold on;
        end;
    end;
    str = ['print /home/julie/Manuscripts/Emotion/figures/',strs{s}(28:32),'Templs.eps -depsc -adobe'];eval(str)
end;

for s = 1:length(strs)
    eval(strs{s})
    clear facvec comods justcomps wtsmat1 jcwts denslist
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);        
    figure;pl = 1;row=4; viewnum=[1,2,3];col=length(viewnum);zoom=1.5;
    for clust = 1:length(comods)
        [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},fullpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % last 1 plots solo IMs in black
        pl = pl+length(viewnum);
    end;
    str = ['print /home/julie/Manuscripts/Emotion/figures/',strs{s}(28:32),'DipsLines.tif -dtiff'];eval(str)
end;
%******************************************************
% Plot DENSITY:
densargs = {'mrislices',[63:-12:-25],'mriview','top','geom',[3,3]};%,'cmax',1};% for 12 plots
% **********************
% combine alpha clusters for density plotting
mask = 1; % one will mask by gdcomps
for s =1:length(strs)
    eval(strs{s})
    clear facvec comods justcomps wtsmat1 jcwts denslist
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);
    if s == 3 % for alphas, try putting all together
        denslist{4} = cell(1,length(denslist{1}));        
        for cls = 1:3
            for nx = 1:length(facvec{cls})
                denslist{4}{nx} = [denslist{4}{nx} denslist{cls}{nx}];
                denslist{4}{nx} = [denslist{4}{nx} justcomps{cls}{nx}];
                denslist{4}{nx} = unique(denslist{4}{nx});
            end;    
        end; 
    end;    
    for clust = 1:length(denslist)   
        if mask == 1
            figure;PlotDipoles('sources.set', fullpaths, denslist{clust},[],[],'alldistance','yred','off',densargs,gdcomps); % plots the probability of a dipole        
            ph=textsc([strs{s}(28:32) ,' Cluster ',int2str(clust)],'title');  set(ph,'color','r');        
            str = ['print /home/julie/Manuscripts/Emotion/figures/',strs{s}(28:32),'DensCls',int2str(clust),'MASK.tif -dtiff'];eval(str)
        else
            figure;PlotDipoles('sources.set', fullpaths, denslist{clust},[],[],'alldistance','bred','off',densargs,[]); % plots the probability of a dipole        
            ph=textsc([strs{s}(28:32) ,' Cluster ',int2str(clust)],'title');  set(ph,'color','r');        
            str = ['print /home/julie/Manuscripts/Emotion/figures/',strs{s}(28:32),'DensCls',int2str(clust),'.tif -dtiff'];eval(str)
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure ~10 -- single subject example(s) of IM weights during emos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all pairwise:
cols = jet(15);cols(10,:) = [1 .9 0]; msize = 20;
for nx = 13:35
s = load([fullpaths{nx},savedat,'.mat']);     
sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0); 
ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall
speceig = floatread(s.eigfile,[s.numtrials s.pcs],[],0);
specwts = speceig*winv;  % templates   
winv = specwts;    
figure; row = 3; col=3; pl =  1;
for im = 1:6%size(winv,2)-1
    for imm = im+1:7%size(winv,2)
        if pl > row*col
            textsc(['Subject ',int2str(nx)],'title');
            figure; pl=1;
        end;
        sbplot(row,col,pl);pl = pl+1;  
        for e = 1:length(s.dstrials)
            ph = plot(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im),winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm),'.','markersize',msize - e);
            set(ph,'color',cols(e,:));
            hold on;
        end;
        plot([0 0],[get(gca,'ylim')],'k-');
        plot([get(gca,'xlim')],[0 0],'k-');
        xlabel(['IM ',int2str(im)]);
        ylabel(['IM ',int2str(imm)]);
    end;
end;
textsc(['Subject ',int2str(nx)],'title');
end;
% excellent examples:
nx=2; ims = [4,11,15]; % **
nx=3; ims = [3,4,12];
nx=4; ims = [2];
nx=5; ims = [1,21,31];
nx=6; ims = [1,2,6,19]; % **
nx=7; ims = [2];
nx=8; ims = [];
nx=9; ims = [4,7];  % **%3, ** ,8,10
nx=10; ims = [8,13];
nx=11; ims = [1,4,16]; % **
nx=12; ims = [];
nx=13; ims = [10];
nx=14; ims = [1,2,3,7,8]; % **
nx=15; ims = [1,2]; % **, few trials
nx=16; ims = [2,8,11,13];
nx=17; ims = [4,10,23]; % **
nx=18; ims = [1,18,26]; % **
nx=19; ims = [2,8,13];
nx=20; ims = [1,6,11];
nx=21; ims = [1,4,17]; % ** 1,4
nx=22; ims = [1,4,17];
nx=23; ims = [3,6,9,10]; % **
nx=24; ims = [1,4,5,13]; % **
nx=25; ims = [1,6,10]; % **
nx=26; ims = [3,10]; % **
nx=27; ims = [4,8,23];
nx=28; ims = [2,6,16,18];
nx=29; ims = [5,9,11,14]; % **
nx=30; ims = [4,14];
nx=31; ims = [2];
nx=32; ims = [];
nx=33; ims = [2];
nx=34; ims = [3,4,7,17,22];% **
nx=35; ims = [8,10,35];

figure; msize = 10; row = 3; col=3; pl =  1;
cols = jet(15);cols(10,:) = [1 .9 0]; msize = 20;
s = load([fullpaths{nx},savedat,'.mat']);     
sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0); 
ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall
speceig = floatread(s.eigfile,[s.numtrials s.pcs],[],0);
specwts = speceig*winv;  % templates   
winv = specwts;    
for im = 1:length(ims)-1
    for imm = im+1:length(ims)
        if pl > row*col
            textsc(['Subject ',int2str(nx)],'title');
            figure; pl=1;
        end;msize = 10;
        %figure;
        sbplot(row,col,pl);pl = pl+1;  
        for e = 1:length(s.dstrials)
            ph = plot(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im)),winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm)),'.','markersize',msize);msize = msize-.1;
            %ph = plot(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im))),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm))),'.','markersize',16);
            set(ph,'color',cols(e,:));
            %ph = text(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im))),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm))),emo2{e}); set(ph,'color',cols(e,:));set(ph,'fontsize',16);
            hold on;
        end;
        plot([0 0],[get(gca,'ylim')],'k-');
        plot([get(gca,'xlim')],[0 0],'k-');
        xlabel(['IM ',int2str(ims(im))]);
        ylabel(['IM ',int2str(ims(imm))]);
        textsc(['Subject ',int2str(nx)],'title');
    end;
end;

function PlotScalpMaps(datset,paths,complist,ttl,savettl)

corrcut = .5; useratings = 4;
[subjvalact,valcorr,actcorr,subjidxval,subjidxact] = RegressEmos(emomeans,subjlist,corrcut,useratings);
for nx = 1:35
    if ~isempty(actcorr{nx})
        figure; msize = 5; row = 2; col=2; pl =  1;
        s = load([fullpaths{nx},savedat,'.mat']);     
        sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
        ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall
        ims = valcorr{nx};
        delu = [];
        sd = mean(std(emomeans{nx}'))+1.5*std(std(emomeans{nx}'));
        for im = 1:length(ims)
            if std(emomeans{nx}(ims(im),:)) < sd
                delu = [delu im];
            end;
        end;
        ims(delu) = [];
        for im = 1:length(ims)
            %corrims = [corrims; emomeans{nx}(ims(im),:)];
            if ~isempty(ims)
            sbplot(row,col,pl);pl = pl+1;  
            for e = 1:length(s.dstrials)
                ph = plot(e,winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im)),'.','markersize',msize);
                set(ph,'color',cols(e,:));        hold on;
                ph = plot(e,median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im))),'.','markersize',msize*5);
                set(ph,'color',cols(e,:));
            end;
            ph = plot(emomeans{nx}(ims(im),:),'k-');
            title(['IM ',int2str(ims(im))]);
            plot([get(gca,'xlim')],[0 0],'k-');
            end;
        end;
        textsc(['Subject ',int2str(nx)],'title');
    end;
end;
corrims = []; % for valence:
for xx = 1:size(subjidxval,1)
    if subjidxval(xx,3) < 0
        corrims = [corrims; emomeans{subjidxval(xx,1)}(subjidxval(xx,2),:)*-1];
    else
        corrims = [corrims; emomeans{subjidxval(xx,1)}(subjidxval(xx,2),:)];
    end;    
end;
corrims = []; % for activity:
for xx = 1:size(subjidxact,1)
    if subjidxact(xx,3) < 0
        corrims = [corrims; emomeans{subjidxact(xx,1)}(subjidxact(xx,2),:)*-1];
    else
        corrims = [corrims; emomeans{subjidxact(xx,1)}(subjidxact(xx,2),:)];
    end;    
end;
figure;        
for e = 1:15
    ph = plot(e,corrims(:,e),'.','markersize',10);
    set(ph,'color',cols(e,:));        hold on;
end;
ph = plot(median(corrims,1),'k-');
plot([get(gca,'xlim')],[0 0],'k-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot IM interactions by clusters--------------------:
for nx = 1:35
s = load([fullpaths{nx},savedat,'.mat']);     
sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall

clust1 =3; % plot all sub clusts
clust2 = 5; % vs sub clusts of this clust

eval(strs{clust1})
cidx = finalidx;
eval(strs{clust2})
ccidx = finalidx;
cols = jet(15); msize = 20;
figure; row = 3; col=3; pl =  1;
for c = 1:length(cidx)
    for cc = 1:length(ccidx)
        im = cidx{c}(find(cidx{c}(:,1)==nx),2);
        imm = ccidx{cc}(find(ccidx{cc}(:,1)==nx),2);
        for i = 1:length(im)
            for ii = 1:length(imm)
                if pl > row*col
                    textsc(['Subject ',int2str(nx)],'title');
                    figure; pl=1;
                end;
                sbplot(row,col,pl);pl = pl+1;  
                for e = 1:length(s.dstrials)
                    ph = plot(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im(i)),winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm(ii)),'.','markersize',msize - e);
                    set(ph,'color',cols(e,:));
                    hold on;
                end;
                plot([0 0],[get(gca,'ylim')],'k-');
                plot([get(gca,'xlim')],[0 0],'k-');
                xlabel(['IM ',int2str(im(i))]);
                ylabel(['IM ',int2str(imm(ii))]);
            end;
        end;
    end;
end;
textsc(['Subject ',int2str(nx)],'title');
end;
% for selected ims, show the median weights:

plotsubjs = [2,4,9,11,15,17,18,21,22,23,25,28,29,33,34];%25: 2-7
subjims = {[1,7],[1,2],[1,3],[2,5],[2,4],[1,9],[5,6],[1,3],[1,4],[2,3],[1,7],[1,4],[2,6],[1,3],[1,2]};
figure;row=3; col=3;pl=1;
for nxx = 1:length(plotsubjs)
    nx = plotsubjs(nxx);    ims = subjims{nxx};
    im =ims(1);    imm =ims(2);
    s = load([fullpaths{nx},savedat,'.mat']);     
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall    
    if pl > row*col
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        figure; pl=1;
    end;
    sbplot(row,col,pl); pl=pl+1; set(gca,'fontsize',14);
    for e = 1:length(s.dstrials)
        %ph = plot(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im),winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm),'.','markersize',12);set(ph,'color',cols(e,:));
        
        ph = plot(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im)),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm)),'.','markersize',msize); set(ph,'color',cols(e,:));
        ph = text(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im)),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm)),emo2{e});    set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        hold on;
    end;
    plot([0 0],[get(gca,'ylim')],'k-');plot([get(gca,'xlim')],[0 0],'k-');
    xlabel(['IM ',int2str(im)]);ylabel(['IM ',int2str(imm)]);
    title(['Subj ',int2str(nx)]);
    set(gca,'xlim',[min(winv(:,im))-.04 max(winv(:,im))+.04]);
    set(gca,'ylim',[min(winv(:,imm))-.04 max(winv(:,imm))+.04]);
    %str = ['print /home/julie/Manuscripts/Emotion/figures/Subj',int2str(nx),'ExIMscatter',int2str(im),'-',int2str(imm),'.eps -depsc -painters'];eval(str)
    %str = ['print /home/julie/Manuscripts/Emotion/figures/Subj',int2str(nx),'ExIMmedian',int2str(im),'-',int2str(imm),'.eps -depsc -painters'];eval(str)
end;
     set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    str = ['print /home/julie/Manuscripts/Emotion/figures/SubjExIMscatters1.eps -depsc -painters'];eval(str)

%%%% plot 3 together:------------
nx=31;
s = load([fullpaths{nx},savedat,'.mat']);     
sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
ws = wts*sph;   winv = pinv(ws); clear wts sph ws icamatall

im =1
imm =2
immm =3
figure; % all points, not median
for e = 1:length(s.dstrials)
    pnt1 = winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im);
    pnt2 = winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm);
    pnt3 = winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),immm);
    ph=plot3(pnt1,pnt2,pnt3,'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['IM ',int2str(im)]);ylabel(['IM ',int2str(imm)]);zlabel(['IM ',int2str(immm)]);
str = ['print /home/julie/Manuscripts/Emotion/figures/Subj',int2str(nx),'ExIMscatter',int2str(im),'-',int2str(imm),'-',int2str(immm),'.eps -depsc'];eval(str)

figure; clear pnt1 pnt2 pnt3% just the median of each emo
for e = 1:length(s.dstrials)
    pnt1(1,e) = median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),im));
    pnt2(1,e) = median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),imm));
    pnt3(1,e) = median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),immm));
    ph=plot3(pnt1(1,e),pnt2(1,e),pnt3(1,e),'.');hold on;
    set(ph,'markersize',25);set(ph,'color',cols(e,:));
    ph = text(pnt1(1,e),pnt2(1,e),pnt3(1,e),emo2{e});
    set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
end;
zl = get(gca,'zlim');
for e = 1:length(s.dstrials)
    ph =plot3([pnt1(1,e) pnt1(1,e)],[pnt2(1,e) pnt2(1,e)],[zl(1)  pnt3(1,e)]);
    set(ph,'color',cols(e,:)); set(ph,'linewidth',2)             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['IM ',int2str(im)]);ylabel(['IM ',int2str(imm)]);zlabel(['IM ',int2str(immm)]);
 str = ['print /home/julie/Manuscripts/Emotion/figures/Subj',int2str(nx),'ExIMmedian',int2str(im),'-',int2str(imm),'-',int2str(immm),'.eps -depsc'];eval(str)
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot back-proj to each subject and decompose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common2/emotion/SpectralDecompCombine.mat 
cols = jet(15);cols(10,:) = [1 .9 0];
cls = 6; dim = [1:3];msize = 10;
bp = fullmd{cls}(:,dim)*fullwts{cls}(dim,:);
figure; row = 5; col = 5; pl=1;
for nxx = 1:length(subjlist)
    nx=subjlist(nxx);
    if pl > row*col
        figure; pl=1;
    end;
    os = find(keeptrack{cls}(:,1) == nx);
    dimwts = fullwts{cls}(dim(1),os);
    [val x] = max(abs(dimwts));
    subjim = ktsubj{cls}(os(x),:);
    savespc(nx,1) = dimwts(x);
    dimwts = fullwts{cls}(dim(2),os);
    [val x] = max(abs(dimwts));
    savespc(nx,2) = dimwts(x);
% $$$     sbplot(row,col,pl); pl = pl+1;
% $$$     for e = 1:size(bp,1)
% $$$         ph = plot(e,bp(e,os(x)),'.','markersize',msize); hold on; 
% $$$         set(ph,'color',cols(e,:));
% $$$         ph = text(e,bp(e,os(x)),emo2{e});    
% $$$         set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
% $$$     end;
% $$$     set(gca,'ylim',[-2 2]);
% $$$     title(['S ',int2str(subjim(1)),' IM ',int2str(subjim(2))]);
end;
clear savespc
for cls = 1:5
figure; dim(1) = dimclusts{1}(cls);
dim(2) = dimclusts{2}(cls);pl = (cls-1)*2+1;
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);    
    os = find(keeptrack{cls}(:,1) == nx);
    dimwts = fullwts{cls}(dim(1),os);
    [val x] = max(abs(dimwts));
    savespc(nxx,pl) = max(abs(dimwts));pl=pl+1;
    %iminst = pics{cls}{nx}(find(abs(pics{cls}{nx}(:,2))==keeptrack{cls}(os(x),2)),2);
    %if mean(iminst) < 0
    %    savespc(nxx,pl) = dimwts(x)*-1;pl=pl+1;
    %else
    %    savespc(nxx,pl) = dimwts(x);pl=pl+1;
    %end;
    dimwts = fullwts{cls}(dim(2),os);
    [val xx] = max(abs(dimwts));
    savespc(nxx,pl) = max(abs(dimwts));pl=pl-1;
    %savespc(nxx,pl) = dimwts(xx);pl=pl-1;
    ph = plot(dimwts(x),dimwts(xx),'.','markersize',msize);hold on;
    ph = text(dimwts(x),dimwts(xx),[' ',int2str(nx)]);
    set(ph,'fontsize',18);
    plot([0 0],[get(gca,'ylim')],'-');
    plot([get(gca,'xlim')],[0 0],'-');
end;
title(['Spectral cluster ',int2str(cls),' Dim ',int2str(dim(1)),'-',int2str(dim(2))]);
end;
% md scaling:
dd = pdist(savespc, 'euclidean') ;    
[md,mwts] = cmdscale(dd);
figure; 
for nxx = 1:length(subjlist)
plot3(md(nxx,1),md(nxx,2),md(nxx,3),'.','markersize',msize); hold on;
ph = text(md(nxx,1),md(nxx,2),md(nxx,3),[' ',int2str(subjlist(nxx))]);
set(ph,'fontsize',18);
end;
plot3([0 0],[get(gca,'ylim')],[0 0],'-');
plot3([get(gca,'xlim')],[0 0],[0 0],'-');
plot3([0 0],[0 0],[get(gca,'zlim')],'-');
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['Dim 1']);ylabel(['Dim 2']);zlabel(['Dim 3']);

% re-decompose and plot:--------------------
figure; row = 6; col = 6; pl=1;msize = 15;
flipx = [3,4,9,13,18,19,20,23,24,27,30,34];
flipy = [2,7,9,10,17,20,23,24,25,26,27,30,31,34,35];
for nx = 1:max(keeptrack{cls}(:,1))
    os = find(keeptrack{cls}(:,1) == nx);
    if ~isempty(os)& length(os) > 3
        %onesubj = bp(:,os);
        %[weights,sphere,compvars,bias,signs,lrates,activations] = runica(onesubj,'pca',2,'stop',1e-7,'maxsteps',2000);
        %onewinv = pinv(weights*sphere);
        dimwts = fullwts{cls}(dim(1),os);
        [val x] = max(abs(dimwts));
        onewinv(:,1) = fullmd{cls}(:,dim(1))*fullwts{cls}(dim(1),x);
        dimwts = fullwts{cls}(dim(2),os);
        [val x] = max(abs(dimwts));
        onewinv(:,2) = fullmd{cls}(:,dim(2))*fullwts{cls}(dim(2),x);
        %sbplot(row,col,pl); pl = pl+1;
        if ismember(nx,flipx)
            onewinv(:,1) = onewinv(:,1)*-1;
        end;
        if ismember(nx,flipy)
            onewinv(:,2) = onewinv(:,2)*-1;
        end;
        for e = 1:size(bp,1)
            ph = plot(onewinv(e,1),onewinv(e,2),'.','markersize',msize); hold on; 
            %ph = plot3(onewinv(e,1),onewinv(e,2),onewinv(e,3),'.','markersize',msize); hold on; 
            set(ph,'color',cols(e,:));
            %ph = text(onewinv(e,1),onewinv(e,2),emo2{e});    
            %set(ph,'color',cols(e,:)); set(ph,'fontsize',22); 
        end;
        %zl = get(gca,'zlim');
        %for e = 1:size(onewinv,1)
        %    ph =plot3([onewinv(e,1) onewinv(e,1)],[onewinv(e,2) onewinv(e,2)],[zl(1)  onewinv(e,3)]);
        %    set(ph,'color',cols(e,:)); set(ph,'linewidth',2)             
        %end;
        % xlabel(['Dim 1']);ylabel(['Dim 2']);zlabel(['Dim 3']);
        % set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
        %title(['Subj ',int2str(nx)]);
    end;
end;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get stats on all clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validsbj = [2:21,23:31,33:35];
clear allfacts
for s = 1:length(strs)
    eval(strs{s})
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);        
    % find how many comps/IMs/subjects per cluster
    for cls = 1:length(facvec)
        numcomps = 0; numims = 0; numsubj = 0;
        for nxx = 1:length(validsbj)
            nx = validsbj(nxx);
            if ~isempty(facvec{cls}{nx})
                numcomps = numcomps + length(facvec{cls}{nx});
                numims = numims + length(unique(facvec{cls}{nx}));
                numsubj = numsubj +1;
            end;
        end;
        allfacts{s}(cls,:) = [cls, numsubj, numims, numcomps];
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is the abs and rel RMS the all templates?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

incsubjs = [1:21,23:35];
histrms = zeros(1,0); hirelrms = zeros(1,0);
allhistrms = zeros(1,0); allhirelrms = zeros(1,0);
for nxx = 1:length(incsubjs)
    nx = incsubjs(nxx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = load([fullpaths{nx},savedat,'.mat']);   
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws icamatall

    clear wts sph ws allfacs alltemps
    for tp = 1:size(activations,1)
        clear onebkprj  ws wts data backproj acts plotprj sph allrms x maxval
        for rcp = 1:length(s.complist)
            allrms(rcp) = sqrt(mean(activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp).^2));
        end;
        hirms = allrms(find(allrms > max(allrms)*.5));
        % all templates > 50% of max rms
        allhistrms(end+1:end+length(allrms)) = allrms;
        histrms(end+1:end+length(hirms)) = hirms;
        newrms = allrms/max(allrms);
        newhi = newrms(find(allrms > max(allrms)*.5));
        hirelrms(end+1:end+length(newhi)) = newhi;
        allhirelrms(end+1:end+length(newrms)) = newrms;
    end;
    fprintf('\nsubj %s\n',int2str(nx)); 
end;
figure; 
sbplot(2,2,1); hist(histrms,100); title('RMS: all used templates');
xl = get(gca,'xlim'); set(gca,'xlim',[0 xl(2)]); 
sbplot(2,2,2); hist(hirelrms,100); title('RMS/max(RMS): used templates');
xl = get(gca,'xlim'); set(gca,'xlim',[0 xl(2)]); 
sbplot(2,2,3); hist(allhistrms,100); title('RMS: ALL IC templates');
xl = get(gca,'xlim'); set(gca,'xlim',[0 xl(2)]); 
sbplot(2,2,4); hist(allhirelrms,100); title('RMS/max(RMS):ALL ICs');
xl = get(gca,'xlim'); set(gca,'xlim',[0 xl(2)]); 
set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
figure;
[alldat bins] = hist(allhirelrms,100); 
[dat bins] = hist(hirelrms,bins); 
ph = bar(bins,alldat,'b'); hold on;
ph = bar(bins,dat,'r');
print /home/julie/Manuscripts/CoMod/figures/SuppFigs/TemplRMS.eps -depsc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many button and non-button subjects per cluster?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/AllClustFacsNEW.mat    
incsubjs = [1:21,23:35];
allrep = 0; allno = 0; allfeel = 0;
for nxx = 1:length(incsubjs)
  nx = incsubjs(nxx);
  for im = 1:length(bigwts{nx})
    if nx < 13
      allrep = allrep + length(bigwts{nx}{im});
    elseif nx > 12 & nx < 21
      allno = allno + length(bigwts{nx}{im});
    elseif nx > 21
      allfeel = allfeel + length(bigwts{nx}{im});
    end;
  end;    
end;
button = [1:12]; % repetitive button presses
button = [21:26]; %'only when you feel it' subjects
button = [13:20]; % no button press (apart from the first one)
% input one at a time. Will collect across clusters
clear allbuts allbutsscale
row = 2; col = 3;
for s = 1:length(strs)
    eval(strs{s})
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(gdcomps,fullpaths,finalidx,allbigs,bigwts,orivec);        
    %  Only put this in in the first round to initialize:
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    for cls = 1:length(facvec)
        repbut = 0; nobut = 0; feelbut = 0;
        for nx = 1:length(facvec{cls})
            if ~isempty(facvec{cls}{nx})& nx ~= 22
                if nx < 13
                    repbut = repbut + length(facvec{cls}{nx});
                elseif nx > 12 & nx < 21
                    nobut = nobut + length(facvec{cls}{nx});
                elseif nx > 21
                    feelbut = feelbut + length(facvec{cls}{nx});
                end;
            end;
        end;
        allbuts{s}(cls,:) = [repbut,nobut,feelbut];
        allbutsscale{s}(cls,:) = [repbut/allrep,nobut/allno,feelbut/allfeel]/size(finaltempls{cls},1);
        
        sbplot(row,col,cls)
        ph = bar(1,allbuts{s}(cls,1));hold on;set(ph,'facecolor','b');
        ph = bar(2,allbuts{s}(cls,2));hold on;set(ph,'facecolor','r');
        ph = bar(3,allbuts{s}(cls,3));hold on;set(ph,'facecolor','g');
        ph = bar(4,allbutsscale{s}(cls,1)*100000);set(ph,'facecolor',[.6 .6 .9]);
        ph = bar(5,allbutsscale{s}(cls,2)*100000);set(ph,'facecolor',[.9 .6 .6]);
        ph = bar(6,allbutsscale{s}(cls,3)*100000);set(ph,'facecolor',[.6 .9 .6]);
        title([strs{s}(28:32),' Cluster ',int2str(cls)]);
        set(gca,'xtick',[2,5]);
        set(gca,'xticklabel',{'#cls IMs','#cls IMs/tot#IMs'});
        set(gca,'xlim',[0 7]);
    end;
    sbplot(row,col,cls+1) % legend:
    ph = bar(1,allbuts{s}(cls,1));hold on;set(ph,'facecolor','b');
    ph = bar(2,allbuts{s}(cls,2));hold on;set(ph,'facecolor','r');
    ph = bar(3,allbuts{s}(cls,3));hold on;set(ph,'facecolor','g');
    ph = bar(4,allbutsscale{s}(cls,1)*10);set(ph,'facecolor','b');
    ph = bar(5,allbutsscale{s}(cls,2)*10);set(ph,'facecolor','r');
    ph = bar(6,allbutsscale{s}(cls,3)*10);set(ph,'facecolor','g');
    legend('Repetitive Button (12 subjs)','NO Button (8 subjs)','Feel-it Button (14 subjs)','location','west')
    textsc([strs{s}(28:32),'; Numbers of IM templates in each expt group'],'title');
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    [P(s),ANOVATAB,STATS] = anova1(allbutsscale{s},{'Repetitive','NO button','Feeling it'}) ;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many IC templates per IM?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numics = zeros(1,0);
for nx = 1:length(bigwts)
    for im = 1:length(bigwts{nx})
        numics(1,end+1) = length(bigwts{nx}{im});    
    end;
end;
figure; hist(numics,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%skew,variance and kurtosis of wts histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear vrinc skn kurtos
for nx = 1:35%length(fullpaths)
    s = load([fullpaths{nx},savedat,'.mat']);   
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    ws = wts*sph;      winv = pinv(ws); clear wts sph ws icamatall
    fullwinv = winv; % reassign 'winv' each time
    clear winv
    for im = 1:size(fullwinv,2)
        for t = 1:length(s.dstrials) % break up windows into emos...
            winv = fullwinv(sum(s.dstrials(1:t-1))+1:sum(s.dstrials(1:t)),:);            
            vrinc(nx,im,t) = var(winv(:,im));
            skn(nx,im,t) = skewness(winv(:,im));
            kurtos(nx,im,t) = kurt(winv(:,im));
        end;
    end;
    fprintf('.');
end;
% break up by clusters:
pl=1; clear vars skews kurts
for c = 1:length(strs)
    eval(strs{c})
    for cls = 1:length(finalidx)
        v=zeros(15,0);s=zeros(15,0);k=zeros(15,0);
        for nx = 1:max(finalidx{cls}(:,1))
            ims = unique(abs(finalidx{cls}(find(finalidx{cls}(:,1) == nx),2)));
            vv=[];ss=[];kk=[];
            if ~isempty(ims)
                for e = 1:size(vrinc,3)
                    vv = [vv; vrinc(nx,ims,e)];
                    ss = [ss; skn(nx,ims,e)];
                    kk = [kk; kurtos(nx,ims,e)];
                end;
                v(:,end+1:end+size(vv,2)) = vv;
                s(:,end+1:end+size(ss,2)) = ss;
                k(:,end+1:end+size(kk,2)) = kk;            
            end;
        end;
        vars{pl} = v;
        skews{pl} =s;
        kurts{pl} = k;
        pl = pl+1;
    end;
end;
allmeas = {vars,skews,kurts};  figttl = {'Variance','Skewness','Kurtosis'};
ttl = {'Delta','Theta','AlphaLow','AlphaPeak','AlphaHi','BetaLow','BetaHi','HiLow','Broadband'};
cols = jet(15);cols(10,:) = [1 .9 0];
for m = 1:length(allmeas)
    figure;row = 3; col = 3;
    for cls = 1:length(allmeas{m})    
        sbplot(row,col,cls);
        for e = 1:size(allmeas{m}{cls},1)
            ph = bar(e,mean(allmeas{m}{cls}(e,:),2));hold on;set(ph,'facecolor',cols(e,:)); 
            ph = plot([e e],[mean(allmeas{m}{cls}(e,:),2)-std(allmeas{m}{cls}(e,:)) mean(allmeas{m}{cls}(e,:),2)+std(allmeas{m}{cls}(e,:))],'k-');
        end;
        if cls < 8
            set(gca,'ylim',[-1.5 1.5]);
        else
            set(gca,'ylim',[-5 15]); 
        end;
        title(ttl{cls});set(gca,'xticklabel',[]);
    end;
    textsc([figttl{m},': Mean +- standard dev'],'title');
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Correlate IMs between clusters within subj then avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allvec=cell(1,0);
for s = 1:length(strs)
    eval(strs{s})
    clear facvec comods justcomps wtsmat1 jcwts denslist
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);  
    for cls = 1:length(facvec)
        allvec{end+1} = facvec{cls};
    end;
end;
shuffnum = 200; clear corr bootstats MI MIboot varcorr varstats
for nx = 1:length(fullpaths)
    checkim = [];
    for cls = 1:length(allvec)
        checkim = [checkim allvec{cls}{nx}'];
    end;
    checkim = unique(checkim);
    if ~isempty(checkim)
        [corr{nx}, bootstats{nx},MI{nx},MIboot{nx},varcorr{nx},varstats{nx}] = CorrCoMod(fullpaths{nx},savedat,checkim,shuffnum,'off');
        %save /data/common1/emotion/EmoClustCorrs.mat corr bootstats MI MIboot varcorr varstats allvec
        %save /data/common1/emotion/EmoClustCorrsns.mat corr bootstats MI MIboot varcorr varstats allvec
        %save /data/common1/emotion/EmoClustCorrsbyEMOns.mat corr bootstats MI MIboot varcorr varstats allvec
        %save /data/common1/emotion/EmoClustCorrsbyEMO.mat corr bootstats MI MIboot varcorr varstats allvec
        %load /data/common1/emotion/EmoGammaCorrs.mat corr bootstats MI MIboot varcorr varstats allvec
    end;
    fprintf('\nSubject %s done.',int2str(nx));
end;

for nx = 1:length(bootstats)
    if ~isempty(bootstats{nx})
        for m = 1:size(bootstats{nx},1)
           for mm = 1:size(bootstats{nx},2)
               bootstats{nx}(mm,m,:,:) = bootstats{nx}(m,mm,:,:);
               varstats{nx}(mm,m,:,:) = bootstats{nx}(m,mm,:,:);
               MIboot{nx}(mm,m,:,:) = bootstats{nx}(m,mm,:,:);
           end;
        end;
    end;
end;               

subjlist = [1:length(corr)];
subjlist = [2:21,23:31,33:35];% all but 1,22,32

 [savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(corr,bootstats,subjlist,{'Brain','Muscle','OMT'},emos,allvec,'corr');


 [savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(corr,bootstats,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'corr');

[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(varcorr,varstats,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'varcorr');

[savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(MI,MIboot,subjlist,{'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'},emos,allvec,'mi');

clustlabels = {'Theta','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Gamma'};
%clustlabels = {'Delta','Theta','Alpha','Beta','Gamma'};
% make all clusters and subclusters separate:
pl = 1; clear idxs 
for s = 1:length(strs)
    eval(strs{s});
    for cls = 1:length(finalidx)
        idxs{pl} = finalidx{cls}; pl = pl+1;
    end;
end;


x=5;
for c = 10:10
    clsmat(4,x) = colmeans(c);x=x+1;
end;
% plot sorted correlations
for nx = 1:length(fullpaths)
  for cls1 = 1:length(allvec)-1
    for cls2 = 1:length(allvec)
      keepcorrs = [keepcorrs,corr{nx}(all
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IMs seem not to be correlated, but they do clump, which is a looser association
%%%% I'm going back to emomeans to see if there are any quadrante effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common2/emotion/EmoWeights.mat emomeans emodeciles
clustlabels = {'Theta','Alpha','Beta','Gamma'}; pl = 1;
for s = 1:length(strs) - 1
    eval(strs{s});
    idx1 = []; 
    for cls = 1:length(finalidx)
        idx1 = [idx1; finalidx{cls}];
    end;
    for ss = s+1:length(strs)
        idx2 = [];
        eval(strs{ss})
        for cls = 1:length(finalidx)
            idx2 = [idx2; finalidx{cls}];
        end;
        relmeans1 = []; relmeans2=  [];
        for nx = 1:length(gdcomps)
            ics1 = unique(idx1(find(idx1(:,1) == nx),3)');
            ics2 = unique(idx2(find(idx2(:,1) == nx),3)');
            ics = intersect(ics1,ics2);
            for ic = 1:length(ics)
                ims{1} = idx1(find(idx1(:,1) == nx & idx1(:,3) == ics(ic)),2); % keep orientation and deal in function
                ims{2} = idx2(find(idx2(:,1) == nx & idx2(:,3) == ics(ic)),2); 
                if ~isempty(ims{1}) & ~isempty(ims{2})
                    for im1 = 1:length(ims{1})
                        bigics1 = allbigs{nx}{abs(ims{1}(im1))};
                        for im2 = 1:length(ims{2})
                            bigics2 = allbigs{nx}{abs(ims{2}(im2))};
                            onlyics = intersect(bigics1,bigics2);
                            if ~isempty(onlyics) & ismember(onlyics,ics)                              
                                if ims{1}(im1) < 0
                                    relmeans1 = [relmeans1;emomeans{nx}(abs(ims{1}(im1)),:)*-1];
                                else
                                    relmeans1 = [relmeans1;emomeans{nx}(abs(ims{1}(im1)),:)];
                                end;                                
                                if ims{2}(im2) < 0
                                    relmeans2 = [relmeans2;emomeans{nx}(abs(ims{2}(im2)),:)*-1]; 
                                else
                                    relmeans2 = [relmeans2;emomeans{nx}(abs(ims{2}(im2)),:)]; 
                                end;
                                
                            end;
                        end;
                    end;
                end;
            end;
        end;
        newlabels{pl} = [clustlabels{s},'/',clustlabels{ss}];
        figure; cols = jet(size(relmeans2,2));cols(10,:) = [1 .9 0];
        for e = 1:size(relmeans1,2)
            ph = plot(relmeans1(:,e),relmeans2(:,e),'.');hold on;
            set(ph,'color',cols(e,:));set(ph,'markersize',5); 
            %ph = text(mean(relmeans1(:,e)),mean(relmeans2(:,e)),[' ',emos{e}]);
            %set(ph,'color',cols(e,:));
        end;
        %set(gca,'xlim',[-.08 .08]); set(gca,'ylim',[-.25 .25]);
        plot([0 0],[get(gca,'ylim')],'k-'); plot([get(gca,'xlim')],[0 0],'k-'); 
        title(newlabels{pl}); pl = pl+1;
    end;
end;

        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Density differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = {'Alpha','Beta','Gamma'};

load /data/common4/emotion/CoModAlphaClusts.mat 
comb = [1 2 3]; % combine these clusters for density plot
clear newvec newbigs
for nx = 1:length(facvec{comb(1)})
    newvec{1}{nx} = [facvec{comb(1)}{nx}',facvec{comb(2)}{nx}'];
end;
facvec = newvec;nm = 1;
%%%%%%%%%
load /data/common4/emotion/GammaClusters.mat 
comb = [1 2]; clear newvec  % combine weak and strong co-modulation
newtmpls{2} = [gamatempls{comb(1)};gamatempls{comb(2)}];
for nx = 1:length(facvec{comb(1)})
    newvec{2}{nx} = [facvec{comb(1)}{nx}',facvec{comb(2)}{nx}'];% Co-Mod in slot #2
end;
newvec{1} = facvec{3}; % no co-modulation
newtmpls{1} = gamatempls{3};
facvec = newvec; nm = 3;
%%%%%%%%%
% plot using SnglTrialAnal.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot angle histogram for each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/AllClustFacs.mat 

load /data/common4/emotion/AllCoModAlpha.mat 
load /data/common4/emotion/AllCoModBeta.mat 
load /data/common4/emotion/AllCoModGama.mat 
% need to change 'PlotCoModasDioles()' to advance subplot by 2 instead of 1
    %btstrap = 200; % collect a bunch of angles from randomly connected dipoles within-cluster
    btstrap = gdcomps; % collect random dipoles from all 'gdcomps'
    %btstrap = [];
    [angles,bilats] = PlotCoModAngles(finalidx,gdcomps,fullpaths,bigwts,allbigs,orivec,btstrap);
    str = ['save /data/common4/emotion/alphabilats.mat bilats']; eval(str);
    % clust1: p = .055
    % clust2: p < .01
    % clust3: p < .01
    
     set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
     str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/AlphaCxnDirs.tif -dtiff'];eval(str)
     str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/AlphaCxnDirsWithin.jpg -djpeg'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/AlphaCxnDirs.eps -depsc'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/AlphaCxnDirsWithin.eps -depsc'];eval(str)

    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/BetaCxnDirs.jpg -djpeg'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/BetaCxnDirsWithin.jpg -djpeg'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/BetaCxnDirs.eps -depsc'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/BetaCxnDirsWithin.eps -depsc'];eval(str)
    
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/GammaCxnDirs.jpg -djpeg'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/GammaCxnDirs.eps -depsc'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/GammaCxnDirsWithin.jpg -djpeg'];eval(str)
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/GammaCxnDirsWithin.eps -depsc'];eval(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Density plot of connections (instead)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CoModDensity(finalidx,gdcomps,fullpaths,bigwts,allbigs,orivec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Plot a single trial FFT highly weighted for one fac/comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=2; % cj82:23
comp = 13; % mid occipital:8
templs = [6,14];
frqs = [3 27];
EEG = pop_loadset('awe.set',fullpaths{nx});
s = load([fullpaths{nx},savedat,'.mat']); 
sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
ws = wts*sph;     winv = pinv(ws); clear wts sph ws 
PlotSpecFacEnv('sources.set',savedat,fullpaths{nx},templs,comp,[],[1 1],frqs,1,.99,0);% last 1=mean
    ph=textsc(['Subject ',int2str(nx)],'title');set(ph,'fontsize',16);
    str = ['print /home/julie/Manuscripts/CoMod/nx',int2str(nx),'BackProjsPos',int2str(comp),'.eps -depsc -painters'];eval(str)
    SpecCoModPlot('sources.set',fullpaths{nx},[],templs,savedat,frqs,'n',0,[]);
    str = ['print /data/common4/emotion/Figs/nx',int2str(nx),'BackProjcp',int2str(comp),'.eps -depsc'];eval(str)
    str = ['print /data/common4/emotion/Figs/nx',int2str(nx),'cp',int2str(comp),'Map.jpg -djpeg'];eval(str)


numpnts = 15; clear fac1 fac2 both facwts
[fac1,facwts{1}] = FindWtedFacs(winv,templs(1),templs(2),numpnts);
[fac2,facwts{2}] = FindWtedFacs(winv,templs(2),templs(1),numpnts-2);
[both,facwts{3}] = FindWtedFacs(winv,templs,[],numpnts);
allfacs = {fac1,fac2,both};
shiftmax = 50; xprod = .8;
tmpnts = 1; % gives one second (for validation)
%colrs = {[1 .35 .35],[.5 1 1],[.75 .5 1]};
colrs = {[1 .45 .45],[.5 1 1],[.75 .5 1]};

% now select trials of interest 
%clear seltrials
%%%----------------------------------------------------------------------
% easy way: (one factor at a time) ----------------------
templ = 22; comp = 45;
trials = find(winv(:,templ)>.9)';
%trials = find(winv(:,templ)<-.7)';
[seltrials,setnames] = RetrieveActs(datset,newpaths{nx},comp,trials,savedat);
for tr = 1:size(seltrials,1)
    [pr frqqs] = pwelch(seltrials(tr,:),EEG.srate,EEG.srate/2,EEG.srate*2,EEG.srate);
    pwr(tr,:) = 10*log10(pr');    
end;
figure; trcols = hsv(size(seltrials,1));
    sbplot(2,1,1)
    ph = plot(seltrials','r');hold on;
    plot(mean(seltrials,1),'k','linewidth',2);
    set(gca,'xlim',[1 size(seltrials,2)]);
    title(['Comp ',int2str(comp),', IM ',int2str(templ),'--low trials']);
    sbplot(2,1,2)
    fr = find(frqqs> 3 & frqqs < 80);
    ph = plot(frqqs(fr),pwr(:,fr)','b');hold on;
    ph = plot(frqqs(fr),mean(pwr(:,fr),1),'k','linewidth',2);hold on;
    set(gca,'xlim',[frqqs(fr(1)) frqqs(fr(end))]);
%%%----------------------------------------------------------------------

figure; row = 3; col = 3;pl = 1;
for fc = 2:length(allfacs)
    [seltrials{fc},setnames] = RetrieveActs(datset,fullpaths{nx},comp,allfacs{fc},savedat);
     pwr = zeros(257,0);
    for tr = 1:size(seltrials{fc},1)
        if ~isempty(find(seltrials{fc}(tr,:)))
            [pwr(:,end+1) frqqs] = pwelch(seltrials{fc}(tr,:),EEG.srate,EEG.srate/2,EEG.srate*2,EEG.srate);
        end;
    end;
    pwr(1,:) = []; frqqs(1) = [];    pwr = 10*log10(pwr');
    [pwr] = eegfilt(pwr,EEG.srate,0,40);
    fr = find(frqqs> frqs(1) & frqqs < frqs(2));
    cols = jet(size(seltrials,1));
    sbplot(row,col,pl)
    for tr = 1:size(pwr,1) 
        ph = plot(frqqs(fr),pwr(tr,fr),'k','linewidth',1);hold on; 
        set(ph,'color',colrs{fc});
    end;
    plot(frqqs(fr),mean(pwr(:,fr),1),'k-','linewidth',3);
    set(gca,'xlim',[frqqs(fr(1)) frqqs(fr(end))]);   
    set(gca,'ticklength',[.03 .03]);
    set(gca,'ylim',[-30 2]);
    if fc < 3
    title(['Factor ',int2str(templs(fc))]);
    else
    title('Both');
    end;    
    set(gca,'box','off');
    pl = pl+1;
    sbplot(row,col,[pl pl+1])
    for f = size(seltrials{fc},1):-1:1
        if ~isempty(find(seltrials{fc}(f,:)))
            id = f; break
        end;
    end;               
    [alignmat] = PhaseAlignActs(seltrials{fc}(id,:),seltrials{fc},shiftmax,xprod);
    ph = plot(alignmat','r-','linewidth',.5);hold on;
    set(ph,'color',colrs{fc});
    ph = plot(mean(alignmat,1),'k-','linewidth',2);hold on;
    set(gca,'xlim',[0 size(seltrials{fc},2)]);
    ph = plot([get(gca,'xlim')],[0 0],'k-','linewidth',1.5); 
    set(gca,'ylim',[-5 7]);axis('off')    
    pl = pl+2;
end;
set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
str = ['print /home/julie/Manuscripts/CoMod/nx',int2str(nx),'SelTrialPwrCp',int2str(comp),'Facs',int2str(templs(1)),'-',int2str(templs(2)),'.eps -depsc'];eval(str)

figure; for x=1:size(seltrials{fc},1)
sbplot(4,4,x)
plot(seltrials{fc}(x,:)');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use pre-selected trials to plot scatter plot with.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; colrss = {[1 .2 .2],[.2 1 1],[.75 .5 1]};
ph=plot([0 0],[-4 3],'k-');set(ph,'color',[.5 .5 .5]); hold on;
ph=plot([-3 3],[0 0],'k-');set(ph,'color',[.5 .5 .5]);
plot(winv(:,templs(1)),winv(:,templs(2)),'b.','markersize',10);
set(gca,'box','off');
xlabel(['IM ',templs(1)]); ylabel(['IM ',templs(2)]);
for fc = 1:length(allfacs)
    [indices,x,y] = ScatterWts(winv,templs,length(allfacs{fc}),allfacs{fc},colrss{fc});
end;axis('off')

str = ['print /home/julie/Manuscripts/CoMod/nx',int2str(nx),'ScatterWts',int2str(templs(1)),'-',int2str(templs(2)),'.eps -depsc'];eval(str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scatter plots of 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common2/emotion/DeltaClust.mat
load /data/common2/emotion/ThetaClust.mat 
load /data/common2/emotion/AlphaClust.mat 
load /data/common2/emotion/BetaClust.mat 
load /data/common2/emotion/GammaClust.mat 
SpecCoModPlot('sources.set',fullpaths{nx},[],[1:10],savedat,[3 128],'n',0,[],1);

ims = [5,7];
figure; plot(winv(:,ims(1)),winv(:,ims(2)),'.');
hold on; plot([0 0],[get(gca,'ylim')],'k-');
plot([get(gca,'xlim')],[0 0],'k-');
xlabel(['Dim ',int2str(ims(1))]); 
ylabel(['Dim ',int2str(ims(2))]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Compare waveforms to straight sinusoids %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    tt = [0:tmpnts/(tmpnts*100000):tmpnts];
    if fc == 1
    amp = 4;
    frq = 9.3; 
    shift = 6;
    col = 'm';
    else
    amp = 4;
    frq = 10.15; 
    shift = 5;
    col = 'b';
    end;    
    x = amp * cos(2 * pi * frq * tt + shift);
    tt = tt*256; plot(tt,x,col);


% a straight sinusoid:
 t=[0:.00392:1]*1000;sl=10;
ph =plot(t,seltrials(sl,:),'k-') ;hold on;
    set(ph,'color','g');

plot(t,sin(pi*t+150)*6,'r-');hold on;
set(gca,'xlim',[t(1) t(end)]);
 set(gca,'xlim',[0 650]);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Plot wts/time for each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common2/emotion/DeltaClust.mat
load /data/common2/emotion/ThetaClust.mat 
load /data/common2/emotion/AlphaClust.mat 
load /data/common2/emotion/BetaClust.mat 
load /data/common2/emotion/GammaClust.mat 

[emoorders] = FindEmoOrder(fullpaths,emos);
for clust = 1:length(facvec)
    count = 0;
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust})
            count = count + length(unique(facvec{clust}{nx}));
        end;
    end;
    figure;  labels = 0; % emo labels off
    row = round(sqrt(count)); col = round(sqrt(count)); place = 1;
    %row = 6; col = 5; place = 1;
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust})   
            if row*col < place + length(unique(facvec{clust}{nx}))
                 figure; place=1;
            end;
            [place] = PlotSubjBackProjs(fullpaths{nx},savedat,emos,unique(facvec{clust}{nx}),1,emoorders{nx},row,col,place,labels); % presentation order
            %[place] = PlotSubjBackProjs(fullpaths{nx},savedat,emos,unique(facvec{clust}{nx}),1,emos,row,col,place,labels); % presentation order
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            %str  = ['print /data/common4/emotion/BackProjWts',int2str(nx),'.jpg -djpeg']; eval(str)
            %close
        end
    end;
end;
 figure;[newpl] = PlotSubjBackProjs(fullpaths{nx},savedat,emos,[1:15],1,emoorders{nx},4,4,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Plot wts/time FFTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/AllCoModAlpha.mat 
load /data/common4/emotion/AllCoModBeta.mat 
load /data/common4/emotion/AllCoModGama.mat 
ttls = {'Alpha-Peak','Alpha-Low','Alpha-Hi','Beta-Low','Beta-Hi','','Gamma-Co','Gamma-Part','Gamma-Solo'};
figure; row = 8; col = 15; pl = 1;
for clust = 1:length(facvec)
    PlotClustWtsSpectra(savedat,fullpaths,facvec{clust},1,row,col,pl);
    pl = pl+15;
    %title(ttls{pl}); pl = pl+1;
end;
ph = textsc(['Spectra of Weights over time from Emotion CoMod Decomp'],'title'); set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
xlabel('Log Frequency (Hz)'); ylabel('Log Power (dB)');
print /home/julie/WtsSpectraClustsEmos.jpg -djpeg
print /home/julie/WtsSpectraClusts.jpg -djpeg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Collect stats on ICA decomps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saved in GoodComps.mat are the following variables:
% fms = number of frames used in subject ICA decomposition
% chs = number of channels in corresponding subject
% pc = number of PCA dimensions retained in corresponding subject
% All subjects included except mr72-2
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones

emolengths = zeros(1,0);
for nx = 1:length(fullpaths)
    for em = 1:length(datset)
        EEG = pop_loadset( datset{em},fullpaths{nx},'all');
        emolengths(nx,em) = size(EEG.data,2);
    end;
end;
nchans = zeros(1,0);
for nx = 1:length(fullpaths)
    EEG = pop_loadset( 'sources.set',fullpaths{nx},'all');
    nchans(end+1) = size(EEG.data,1);
end;
min(nchans)
max(nchans)
[mean(nchans),std(nchans)]
2*60*256*2
x=sum(emolengths,2)';
x./(nchans.^2)
min(x./(nchans.^2))
max(x./(nchans.^2))
mean(x./(nchans.^2))
std(x./(nchans.^2))

savedat = 'SpecHPModEmos'; fullpaths = newpaths;
savedat = 'SpecHPModMuscle';  % includes muscle and inf frontal
savedat = 'SpecCoModMoreFreqs'; 
for nx = 1:length(fullpaths)
    s = load([newpaths{nx},savedat,'.mat']); 
    specstats(1,nx) = length(s.rowmeans);
    specpcs(1,nx) = s.pcs;
    [pcared(1,nx) pv]= SpecCoModPvaf(savedat,newpaths{nx},[],[],1,0);
end;
% muscle decomp:
>> mean(pcared) 28.7600

>> max(pcared)   45.8255

>> min(pcared) 18.3262

>> std(pcared) 5.7117
    
for nx = 1:length(fullpaths)
    ngdcomps(1,nx) = length(gdcomps{nx});
end

savename = 'Emo-HP';
for nx = 1:35    
    for bl = 1:nsets(nx)
    EEG = pop_loadset( [savename,'-',int2str(bl),'-',int2str(nchans(nx)),'.set'],oldpaths{nx},'all');
        datlength(nx,bl) = size(EEG.data,2);
    end;
end;
x=sum(datlength,2)';
x./(nchans.^2)
%-------------
% how many comods per IM?
load /data/common1/emotion/AllClustFacs.mat    
load /data/common1/emotion/AllClustFacsWAVE.mat    
load /data/common1/emotion/AllClustFacsMuscle.mat    
allns = [];subjns = [];
for nx = 1:length(allbigs)
    if ~isempty(allbigs{nx})
        clear ncomods
        for im = 1:length(allbigs{nx})
            ncomods(im) = length(allbigs{nx}{im});
        end;
        subjns = [subjns sum(ncomods)];
        allns = [allns ncomods];
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collect emotion orders from all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CoModPwrTrends(savedat,fullpaths{nx},2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Collect emotion orders from all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglab
datset = 'ButtonOnly.set';
clear allord
for nx = 1:length(fullpaths)
    EEG = pop_loadset( datset,fullpaths{nx},'all');
    emord = cell(1,0);
    for ev = 1:length(EEG.event)
        if ~isempty(find(strcmp(EEG.event(ev).type,emos)))
            emord = [emord EEG.event(ev).type];
        end;
    end;
    allord{nx} = emord;
end;
for nx= 1:length(allord)
        fprintf('\n');
    for x = 1:length(allord{nx})
        fprintf('%s\t',allord{nx}{x});
    end;
    if strcmp(allord{nx}{1},'awe')
        if strcmp(allord{nx}{2},'frustration')
            emoidx(nx) = 1;
        elseif strcmp(allord{nx}{2},'fear')
            emoidx(nx) = 2;
        end;
    elseif strcmp(allord{nx}{1},'happy')
        emoidx(nx) = 3;
    elseif strcmp(allord{nx}{1},'compassion')
        emoidx(nx) = 4;
    elseif strcmp(allord{nx}{1},'love')
        emoidx(nx) = 5;
    end;    
end;
for nxx= 1:length(unique(emoidx))
    [x,y] = find(emoidx == nxx);
    for ord = 1:length(y)
        nx = y(ord);
        fprintf('\n');
        for x = 1:length(allord{nx})
            fprintf('%s\t',allord{nx}{x});
        end;
    end;    
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find % variance accounted for by all retained dims of spectra
for nx = 2:length(gdcomps)
    if ~isempty(gdcomps{nx})
        [pcared(1,nx) pv{nx}] = SpecCoModPvaf(savedat,fullpaths{nx},[],[],1,0);    
    end;
    fprintf('\nSubject %s done...\n',int2str(nx));
end;
save /data/common4/emotion/CoModPvafs.mat pcared pv
min(pcared)=  6.18; max(pcared)=23.22;
 clear numims
for nx=1:35 
    for im = 1:length(allbigs{nx})
numims(nx,im) = length(allbigs{nx}{im});     
    
    end;
end;
for nx=2:35
    x=numims(nx,find(numims(nx,:)));
    minnum(1,nx) = min(x);
    maxnum(1,nx) = max(x);
    sumnum(1,nx) = sum(x);
    meannum(1,nx) = mean(x);
    stdnum(1,nx) = std(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
