% specifically for Gamma-focused comod paper:

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
load /data/common1/emotion/AllClustFacsMuscle.mat    
load /data/common1/emotion/GammaGdMsVf.mat % this is the cluster info for joint decomp
% 1:3 = high/low clusters: brain, scalp, ocular-motor
% 4:6 = broadband clusters: brain, scalp, ocular-motor
%savedat = 'SpecHPModEmos'; fullpaths = newpaths;
savedat = 'SpecCoMod'; fullpaths = newpaths;
savedat = 'SpecCoModMuscle';  % includes muscle and inf frontal
savedat = 'SpecCoModMoreFreqs';fullpaths = newpaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  FIGURE 2 (IM examples)%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frqlim = [2 125]; 

nx = 10;tmpls = [26,3,5,7,13,11,12,15,31,20,16,22];
nx = 10;comps = [3,5,6,9,16,19,20,25,22,21,42,92];

nx = 2;tmpls = [3,1,2,16,24,5,20,6,8,9,43,12,39,13];% used in paper
nx = 2;comps = [33,58,70,94,6,5,2,13,10,25];

nx = 9;tmpls = [21,14,15,13,12,19,22,28,11,1,2];% 
nx = 9;comps = [21,14,35,100,5,9,24,28,3,15,11,40,51]; % 
nx = 17;tmpls = [3,7,13,15,18,21,30,34,40];% nx=17
nx = 17;comps = [1,3,8,17,22,15,70,82]; % nx=17
% plot template examples:
SpecCoModPlot('sources.set',newpaths{nx},comps,tmpls,savedat,frqlim,'n',0,[]);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%str = ['print /home/julie/Manuscripts/Gamma/nx',int2str(nx),'TemplExamples.eps -depsc -painters'];eval(str)
%str = ['print /home/julie/Manuscripts/Gamma/nx',int2str(nx),'TemplExamplesLOG.eps -depsc -painters'];eval(str)

% plot mean spectra with back projections
nx=2; comps = [58,5,2];
tmpls = [1,8,20];
for t  =1:length(tmpls)
    PlotSpecFacEnv('sources.set',savedat,newpaths{nx},tmpls(t),comps,[],[],frqlim,1,.99,0);
str = ['print /home/julie/Manuscripts/Gamma/nx',int2str(nx),'BackProjsIM',int2str(tmpls(t)),'.eps -depsc -painters'];eval(str)
end;
ph=textsc(['Subject ',int2str(nx)],'title');set(ph,'fontsize',16);
PlotSpecFacEnv('sources.set',savedat,newpaths{nx},tmpls,comps,[],[],frqlim,1,.99,0);


s = load([newpaths{nx},savedat,'.mat']);     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sph=floatread([newpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([newpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
icamatall = floatread([newpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws icamatall
rcp = 3; tp = 13;
figure;
ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),2); hold on;
%ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); hold on;
set(gca,'xscale','log');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 3D headplots for Fig 3:
plotics = [58,5,2];
EEG = pop_loadset( 'filename', 'sources.set','filepath',fullpaths{2}); 
for ic = 1:length(plotics)
    figure; mypop_headplot(EEG, 0, plotics(ic), int2str(plotics(ic)), [],  'load', '/data/common1/emotion/mi83/sources.spl','electrodes','off');
    str = ['print /home/julie/Manuscripts/Gamma/HeadPlot',int2str(ic),'.tif -dtiff']; eval(str);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral templates:***********
cols = {[.9 .1 .8],[1 .5 0],[.2 1 .2]};
cols = {'r','b','g'};
row=2;%round(sqrt(length(finaltempls))); 
col=2;%ceil(sqrt(length(finaltempls))); 
figure; pl = 1; cls=2;
for cls = 1:length(finaltempls)
    sbplot(row,col,pl); pl = pl+1;
    [realx labelx,han] = quadplot(freqs,finaltempls{cls},2,cols{cls}); hold on;%[0 .75 .75][.2 1 .2]
    [realx labelx,han] = quadplot(freqs,mean(finaltempls{cls},1)',3,'k');%[.16 .5 .3]
    %set(gca,'xlim',[realx(1) realx(end)]);
    set(gca,'ylim',[min(finaltempls{cls}(:)) max(finaltempls{cls}(:))]);
    set(gca,'ticklength',[.03 .03]);
    title([int2str(cls),'-',int2str(size(finaltempls{cls},1))]);      
    %plot([realx(2) realx(2)],[get(gca,'ylim')],'g-');
    %plot([realx(4) realx(4)],[get(gca,'ylim')],'g-');
    plot([get(gca,'xlim')],[0 0],'k-');
    set(gca,'xlim',[1.85 10.75]);
    set(gca,'ylim',[-2.5 10]);
    %set(gca,'xscale','log');
end;
ph=textsc('Template clusters','title');set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/BBclustersCombined.eps -depsc -adobe -painters
%--------------------------------------------------------
% plot the clusters:*************************
clear facvec comods justcomps wtsmat1 jcwts denslist
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);

% plot dipoles
figure;pl = 1;row = 3%length(comods);
 viewnum=[1,2,3];col = 3;%length(viewnum) ;
 zoom= 1.1;
for clust = 2:length(comods)
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
print /home/julie/Manuscripts/Gamma/HighFrqBrainDips.tif -dtiff 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot brain templates and then any OR all associated muscle templs
muscics = cell(length(finaltempls),length(finaltempls));
clustcols = {'r',[.2 .2 1],'g','r',[.2 .2 1],'g','r',[.2 .2 1],'g'};
meancols = {'b','g','r','b','g','r'};
figure; row = 2; col = 2;pl = 1;

clust = 2;    otherclust =1; %  1 for br, 2 for musc

sbplot(row,col,pl);pl = pl+1;
[realx labelx,han] = quadplot(freqs,finaltempls{clust}',2,clustcols{clust}); hold on;%[0 .75 .75][.2 1 .2]
[realx labelx,han] = quadplot(freqs,mean(finaltempls{clust},1)',3,'k');%[.16 .5 .3]
set(gca,'ylim',[min(finaltempls{clust}(:)) max(finaltempls{clust}(:))]);
set(gca,'ticklength',[.03 .03]);
title([int2str(clust),'-',int2str(size(finaltempls{clust},1))]);      
plot([get(gca,'xlim')],[0 0],'k-'); yl = get(gca,'ylim');
set(gca,'xlim',[1.85 10.75]);
set(gca,'ylim',[-6 10]);
%set(gca,'xscale','log');
othertempls = []; keepidx = [];
nxiccell = cell(1,length(gdcomps));

for imm = 1:size(finalidx{clust},1)
    nx = finalidx{clust}(imm,1);
    %currmod = allbigs{nx}{abs(finalidx{clust}(im,2))};
    currmod = [EEG.gdcomps, EEG.ventfrontal, EEG.muscle];

    currmod(find(currmod == finalidx{clust}(im,3))) = [];
    if otherclust == 2 | otherclust == 5| otherclust == 8
        muscmod = currmod(find(ismember(currmod,EEG.muscle)));clust2 = 2;othercol = 'k';% find muscle
    elseif otherclust == 1 | otherclust == 4| otherclust == 7
        muscmod = currmod(find(ismember(currmod,EEG.gdcomps)));clust2 = 1;othercol = 'k';% find gdcomps
    elseif otherclust == 3 | otherclust == 6| otherclust == 9
        muscmod = currmod(find(ismember(currmod,EEG.ventfrontal)));clust2 = 3;othercol = 'k';% find oc motor
    end;
    % find all indices for this subj and IM:
    otheridx = find(kptk(:,1) == finalidx{clust}(im,1) & kptk(:,2) == abs(finalidx{clust}(im,2)));
    % find all ICs for this subj/IM:
    otheric = kptk(find(kptk(:,1) == finalidx{clust}(im,1) & kptk(:,2) == abs(finalidx{clust}(im,2))),3);
    % for the ICs that match the muscle ICs found above, collect the templates:
    othertempls = [othertempls;clustfacs(otheridx(find(ismember(otheric,muscmod))),:)];
    keepidx = [keepidx;kptk(otheridx(find(ismember(otheric,muscmod))),:)];
    nxiccell{nx} = [nxiccell{nx} otheric(find(ismember(otheric,muscmod)))'];
end;
muscics{clust,otherclust} = nxiccell;
sbplot(row,col,pl);pl = pl+1;
[realx labelx,han] = quadplot(origfreqs,othertempls',2,clustcols{clust2}); hold on;% all traces
[realx labelx,han] = quadplot(origfreqs,mean(othertempls,1)',3,othercol);% mean
if size(othertempls,1) == 1
    set(han,'color',clustcols{clust2});
end;
set(gca,'ylim',[min(othertempls(:)) max(othertempls(:))]);
set(gca,'ticklength',[.03 .03]); 
yyl = get(gca,'ylim'); set(gca,'ylim',[min(yyl(1),yl(1)) max(yyl(2),yl(2))]);
title([int2str(clust),'-',int2str(size(othertempls,1))]);
plot([get(gca,'xlim')],[0 0],'k-');    
set(gca,'xlim',[1.85 10.75]);
set(gca,'ylim',[-6 8]);
%set(gca,'xscale','log');
numims(1,clust) = size(othertempls,1);

print /home/julie/Manuscripts/Gamma/BrainCoModMusc.eps -depsc -adobe -painters
print /home/julie/Manuscripts/Gamma/BrainCoModMusc.jpg -djpeg
  
clustnames = {'Brain','Muscle','VentFront'};
figure; row = 3; col = 4;pl = 1;
for clust = 1:3
    keeptempls = cell(1,3);keepidx= cell(1,3);
    clustnxs = unique(finalidx{clust}(:,1));
    for nxx = 1:length(clustnxs)
        nx = clustnxs(nxx);
        EEG = pop_loadset('sources.set' ,newpaths{nx});     
        iclist{1} = EEG.gdcomps; iclist{2} = EEG.muscle; iclist{3} = EEG.ventfrontal; 
        nxims = unique(abs(finalidx{clust}(find(finalidx{clust}(:,1) == nx),2)));
        for imm = 1:length(nxims)
            im = nxims(imm);
            clear sortics
            for typ = 1:3 % brain, muscle, vf
                sortics{typ} = allbigs{nx}{im}(find(ismember(allbigs{nx}{im},iclist{typ})));
                plottempls = [];plotidx = [];
                if ~isempty(sortics{typ})
                    for ic = 1:length(sortics{typ})
                        idx = find(kptk(:,1) == nx & kptk(:,2) == im& kptk(:,3) == sortics{typ}(ic));
                        plottempls = [plottempls; clustfacs(idx,:)];  
                        plotidx = [plotidx; kptk(idx,:)];
                    end;
                end;
                keeptempls{typ} = [keeptempls{typ}; plottempls];
                keepidx{typ} = [keepidx{typ}; plotidx];
            end;
        end;
    end;
    sbplot(row,col,pl);pl = pl+1;
    [han,realx labelx] = quadplot(origfreqs,finaltempls{clust}',2,clustcols{clust}); hold on;% all traces
    [han,realx labelx] = quadplot(origfreqs,mean(finaltempls{clust},1)',3,'k');% mean    
plot([get(gca,'xlim')],[0 0],'k-'); set(gca,'ticklength',[.03 .03]); set(gca,'ylim',[-5 10]);
    title([clustnames{clust},'-',int2str(size(finaltempls{clust},1)),' Templates']);
    for typ = 1:3
        if ~isempty(keeptempls{typ})
            sbplot(row,col,pl);pl = pl+1;
            [han,realx labelx] = quadplot(origfreqs,keeptempls{typ}',2,clustcols{typ}); hold on;% all traces
            [han,realx labelx] = quadplot(origfreqs,mean(keeptempls{typ},1)',3,'k');% mean
plot([get(gca,'xlim')],[0 0],'k-'); set(gca,'ticklength',[.03 .03]);  set(gca,'ylim',[-5 10]);
            title([clustnames{typ},'-',int2str(size(keeptempls{typ},1)),' Templates']);
        else
          pl = pl+1;  
        end;
    end;
end;
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/BrainCoModMusc.eps -depsc -adobe -painters
print /home/julie/Manuscripts/Gamma/BrainCoModMusc.jpg -djpeg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot gamma IM as a power law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=10; im=13; ic = 7; % example 1
nx=17; im = 4; ic = 1; %example 2
nx=17; im = 10; ic = 9; %example 2
nx=6; im = 1; ic = 6; %example 2
nx=6; im = 6; ic = 9; %example 2
nx=30; im = 14; ic = 6; %example 2
nx=29; im = 11; ic = 28; %example 2

s = load([newpaths{nx},savedat,'.mat']);   
sph=floatread([newpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
wts=floatread([newpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0); 
ws = wts*sph; winv = pinv(ws); 
speceig = floatread([newpaths{nx},s.eigfile],[s.numtrials s.pcs],[],0);
specwts = speceig*winv;   
icamatall = floatread([newpaths{nx},savedat,'.fdt'],[s.pcs s.numframes],[],0);    
activations = ws*icamatall;   clear wts sph ws icamatall

figure; ph = logplot(s.freqs,s.meanpwr(ic,:),2,'g');hold on;
ph = logplot(s.freqs,s.meanpwr(ic,:)+activations(im,(ic-1)*length(s.freqs)+1:ic*length(s.freqs)),2,'r');
ph = logplot(s.freqs,s.meanpwr(ic,:)-activations(im,(ic-1)*length(s.freqs)+1:ic*length(s.freqs)),2,'b');
set(gca,'xscale','log'); set(gca,'yscale','log');
title(['S ',int2str(nx),'; IM ',int2str(im)]);

figure; sbplot(2,2,1);
for x = 1:150
    y(x) = 1/x;
end;
plot([1:150],y,'g','linewidth',2); hold on;
set(gca,'xlim',[1 150]);
clear y
for x = 1:150
    y(x) = 5*(1/x);
end;
plot([1:150],y,'r','linewidth',2); hold on;
clear y
for x = 1:150
    y(x) = 1/x^2;
end;
plot([1:150],y,'b','linewidth',2); hold on;
set(gca,'yscale','log');
set(gca,'ylim',[0 12]);
legend('1/f spectrum','5*(1/f)','1/f^2 spectrum');
title(['1/f vs 5/f vs 1/f^2 functions']);
       
sbplot(2,2,2);
for x = 1:150
    y(x) = 1/x;
end;
plot([1:150],y,'g','linewidth',2); hold on;
set(gca,'xlim',[1 150]);
clear y
for x = 1:150
    y(x) = 5*(1/x);
end;
plot([1:150],y,'r','linewidth',2); hold on;
clear y
for x = 1:150
    y(x) = 1/x^2;
end;
plot([1:150],y,'b','linewidth',2); hold on;
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'ylim',[0 12]);
title(['1/f vs 5/f vs 1/f^2 functions; xscale = log']);

sbplot(2,2,3);flims = [4 125];
nx=6; im = 1; ic = 6; %example 1
s = load([newpaths{nx},savedat,'.mat']);   
sph=floatread([newpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
wts=floatread([newpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0); 
ws = wts*sph; winv = pinv(ws); 
speceig = floatread([newpaths{nx},s.eigfile],[s.numtrials s.pcs],[],0);
specwts = speceig*winv;   
icamatall = floatread([newpaths{nx},savedat,'.fdt'],[s.pcs s.numframes],[],0);    
activations = ws*icamatall;   clear wts sph ws icamatall

fr = find(s.freqs > flims(1) & s.freqs < flims(2));
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr),2,'g');hold on;
acts = activations(im,(ic-1)*length(s.freqs)+1:ic*length(s.freqs));
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr)+acts(1,fr),2,'m');
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr)-.2*s.meanpwr(ic,fr),2,'r');
for x = 1:size(s.meanpwr,2)
    newy(x) = s.meanpwr(ic,x) + x^2/10000;
end;
ph = logplot(s.freqs(fr),newy(fr),2,'b');
set(gca,'ylim',[-32 -2]);set(gca,'box','off');
legend('Mean pwr','Mn pwr + IM template','Mean pwr*.2','Mean pwr + x^2 fxn');
title(['S ',int2str(nx),'; IM ',int2str(im)]);

sbplot(2,2,4); flims = [15 125];
fr = find(s.freqs > flims(1) & s.freqs < flims(2));
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr),2,'g');hold on;
acts = activations(im,(ic-1)*length(s.freqs)+1:ic*length(s.freqs));
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr)+acts(1,fr),2,'m');
ph = logplot(s.freqs(fr),s.meanpwr(ic,fr)-.2*s.meanpwr(ic,fr),2,'r');
ph = logplot(s.freqs(fr),newy(fr),2,'b');
set(gca,'xscale','log'); set(gca,'yscale','log');
set(gca,'ylim',[-35 0]);set(gca,'box','off');
title(['log-log zoom of high freqs']);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
print /home/julie/Manuscripts/Gamma/PowerLawFigure.eps -depsc 
print /home/julie/Manuscripts/Gamma/PowerLawFigure.jpg -djpeg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot IM wt examples:------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=2; ims = [1,4];  % **%3, ** ,8,10
nx=9; ims = [4,7];  % **%3, ** ,8,10
nx=17; ims = [10,23]; % **
nx=21; ims = [1,4]; % **
nx=29; ims = [5,9]; % **
nx=34; ims = [7,17];% **

im = 1; imm=2;
cols = jet(15);cols(10,:) = [1 .9 0]; msize = 20;
s = load([newpaths{nx},savedat,'.mat']);
wts = floatread([newpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);
sph = floatread([newpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0);
ws = wts*sph;winv = pinv(ws);
speceig = floatread([newpaths{nx},s.eigfile],[length(s.rowmeans) inf],[],0);
specwts = speceig*winv; winv = specwts;     
%data = floatread([newpaths{nx},savedat,'DAT.fdt'],[s.numrows s.numframes],[],0);

for med = 1:2
    figure;msize = 10;
    for e = 1:length(s.dstrials)
        if med == 1
            ph = plot(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im)),winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm)),'.','markersize',msize);msize = msize-.1;set(ph,'color',cols(e,:));
        else
            ph = plot(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im))),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm))),'.','markersize',32);set(ph,'color',cols(e,:));
            ph = text(median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(im))),median(winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),ims(imm))),emo2{e}); 
            set(ph,'color',cols(e,:));set(ph,'fontsize',16);
        end;        
        hold on;
    end;
    plot([0 0],[get(gca,'ylim')],'k-'); plot([get(gca,'xlim')],[0 0],'k-');
    xlabel(['IM ',int2str(ims(im))]);  ylabel(['IM ',int2str(ims(imm))]);
    textsc(['Subject ',int2str(nx)],'title');
    %str = ['print /home/julie/Manuscripts/Gamma/Subj',int2str(nx),'IMwtEx',int2str(med),'.eps -depsc -painters -adobe']; eval(str);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get stats on all clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validsbj = [2:21,23:31,33:35];
clear allfacts
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
    allfacts(cls,:) = [cls, numsubj, numims, numcomps];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Were any emotions associated with gamma in brain or muscle?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cols = jet(15);cols(10,:) = [1 .9 0]; row = 2; col = 2;
for clust = 2:4%length(finalidx)
    allmeans = [];lastim = [0 0];
    for im = 1:size(finalidx{clust},1)
        currim = [finalidx{clust}(im,1) abs(finalidx{clust}(im,2))];
        if ~isempty(find(currim == lastim)) 
            if length(find(currim == lastim)) < 2 
                if finalidx{clust}(im,2) < 0
                    allmeans = [allmeans; emomeans{finalidx{clust}(im,1)}(abs(finalidx{clust}(im,2)),:)*-1];
                else
                    allmeans = [allmeans; emomeans{finalidx{clust}(im,1)}(finalidx{clust}(im,2),:)];
                end;
            end;
        else
            if finalidx{clust}(im,2) < 0
                allmeans = [allmeans; emomeans{finalidx{clust}(im,1)}(abs(finalidx{clust}(im,2)),:)*-1];
            else
                allmeans = [allmeans; emomeans{finalidx{clust}(im,1)}(finalidx{clust}(im,2),:)];
            end;
        end;        
        lastim = [finalidx{clust}(im,1) abs(finalidx{clust}(im,2))];
    end;
    [P(1,clust),ANOVATAB,STATS] = anova1(allmeans,emos,'on');
    figure;   
    sbplot(row,col,1); 
    for e = 1:size(allmeans,2)
        ph = bar(e,median(allmeans(:,e),1));hold on;set(ph,'facecolor',cols(e,:));
        plot([e e],[median(allmeans(:,e),1)-std(allmeans(:,e),1) median(allmeans(:,e),1)+std(allmeans(:,e),1)],'k-');
    end;set(gca,'xlim',[0 16]);set(gca,'xticklabel',[]);
    title(['Cluster ',int2str(clust),' median']);
    
    sbplot(row,col,2); 
    for e = 1:size(allmeans,2)
        ph = bar(e,mean(allmeans(:,e),1));hold on;set(ph,'facecolor',cols(e,:));
        plot([e e],[mean(allmeans(:,e),1)-std(allmeans(:,e),1) mean(allmeans(:,e),1)+std(allmeans(:,e),1)],'k-');
    end;set(gca,'xlim',[0 16]);set(gca,'xticklabel',[]);
    title(['Cluster ',int2str(clust),' mean']);
    
    sbplot(row,col,3);   
    for e = 1:size(allmeans,2)
        ph = plot([1:size(allmeans,1)],allmeans(:,e),'k-'); hold on;
        set(ph,'color',cols(e,:));
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many IMs significantly differentiate one or more emotions?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some prelimiary plots for reference:
for e = 1:15    
     [freqout] = PlotSpecFacEnv('sources.set',savedat,newpaths{nx},[],gdcomps{nx},e,[],[3 128],1,1,1);
     textsc(['Subj ',int2str(nx),'; All IM back-projections: ',emos{e}],'title');
 end;
subjlist = [2:21,23:31,33:35];% all but 1,22,32
allidx = []; whichcls = []; plotidx = cell(1,5);
for cls = 1:length(finalidx)
    allidx = finalidx{cls};
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        onesubjims = allidx(find(allidx(:,1) == nx),2);
        useims{nx} = unique(abs(onesubjims)');
        [fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,3,nx,useims,'mds',0,['Subj ',int2str(nx)],1);
    end;
end;
%%-----------------------
[emosigs,emoPs] = SigEmoShifts(savedat,newpaths,subjlist,emos);
save /data/common2/emotion/SigEmoShiftsMuscDecomp.mat emosigs emoPs
load /data/common2/emotion/SigEmoShiftsMuscDecomp.mat emosigs emoPs
cols = jet(15);cols(10,:) = [1 .9 0];row = 2; col = 2;
for clust = 1:length(finalidx)
    allsigs = []; lastim = [0 0];
    for im = 1:size(finalidx{clust},1)
        currim = [finalidx{clust}(im,1) abs(finalidx{clust}(im,2))];
        if ~isempty(find(currim == lastim)) 
            if length(find(currim == lastim)) < 2 
                allsigs = [allsigs; emosigs{finalidx{clust}(im,1)}(abs(finalidx{clust}(im,2)),:)];
            end;
        else
            allsigs = [allsigs; emosigs{finalidx{clust}(im,1)}(abs(finalidx{clust}(im,2)),:)];
        end;
        lastim = [finalidx{clust}(im,1) abs(finalidx{clust}(im,2))];
    end;
    [P(1,clust),ANOVATAB,STATS] = anova1(allsigs,emos,'on');
    figure;
    sbplot(row,col,1); 
    for e = 1:size(allsigs,2)
        ph=bar(e,mode(allsigs(:,e),1)); hold on;set(ph,'facecolor',cols(e,:));
    end;set(gca,'xlim',[0 16]);set(gca,'xticklabel',[]);
    title(['Cluster ',int2str(clust),' mode']);
        
    sbplot(row,col,2); 
    for e = 1:size(allsigs,2)
        ph = bar(e,median(allsigs(:,e),1));hold on;set(ph,'facecolor',cols(e,:));
        plot([e e],[median(allsigs(:,e),1)-std(allsigs(:,e),1) median(allsigs(:,e),1)+std(allsigs(:,e),1)],'k-');
    end;set(gca,'xlim',[0 16]);set(gca,'xticklabel',[]);
    title(['Cluster ',int2str(clust),' median']);
    
    sbplot(row,col,3); 
    for e = 1:size(allsigs,2)
        ph = bar(e,mean(allsigs(:,e),1));hold on;set(ph,'facecolor',cols(e,:));
        plot([e e],[mean(allsigs(:,e),1)-std(allsigs(:,e),1) mean(allsigs(:,e),1)+std(allsigs(:,e),1)],'k-');
    end;set(gca,'xlim',[0 16]);set(gca,'xticklabel',[]);
    title(['Cluster ',int2str(clust),' mean']);
    
    sbplot(row,col,4);  
    for e = 1:size(allsigs,2)
    ph = plot([1:size(allsigs,1)],allsigs(:,e),'k-'); hold on;
    set(ph,'color',cols(e,:));
    end;
    textsc('Number of times each emotion was sig diff from any other emotion in a clustered IM','title')
end;
            [P(1,clust),ANOVATAB,STATS] = anova1(allsigs,emos,'on');
            figure;comp = multcompare(STATS,'alpha',.05,'ctype','bonferroni'); 
            comppairs = [];
            for cp = 1:size(comp,1)
                if length(find(comp(cp,[3,5]) == abs(comp(cp,[3,5])))) ~= 1 %(not straddling 0=sig)
                    comppairs = [comppairs;[comp(cp,[1:2])]];
                end;
            end;
