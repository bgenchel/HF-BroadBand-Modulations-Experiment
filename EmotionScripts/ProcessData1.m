% Run 'single-trial' analysis of spectral co-modulation on emotion data
addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
gdcomps{36}=[];
%str = ['load /data/common4/emotion/GoodComps.mat ']; eval(str);
%load /data/common4/emotion/clustfacs.mat 
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
savedat = 'SpecHPModEmos'; fullpaths = newpaths;
%savedat = 'SpecHPModMuscle';  % includes muscle and inf frontal

strs = {
    'load /data/common2/emotion/DeltaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/ThetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/AlphaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/BetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/GammaClust.mat finaltempls finalidx finalmean freqs'};


str = ['load ',fullpaths{nx},savedat,'.mat '];eval(str);    
PlotSpecFacEnv('sources1.set',savedat,fullpaths{nx},[1:pcs],gdcomps{nx},gdcomps{nx},1,[3 45],1,.98,0);

SpecCoModPlot('emo-1-248-NF.set',fullpaths{nx},comps,comps,[1:pcs],savedat,[3 125],'n',0,[]);
PlotSpecFacEnv('sources1.set',savedat,fullpaths{nx},[2,3,5,6,15],gdcomps{nx},gdcomps{nx},1,[3 45],1);
PlotSpecFacEnv('sources1.set',savedat,fullpaths{nx},[2,3,5,6,15],gdcomps{nx},gdcomps{nx},1,[3 45],1);
[pcared pv] = SpecCoModPvaf(savedat,fullpaths{nx},[],1);    
SpecCoModPlot('sources1.set',fullpaths{nx},gdcomps{nx},gdcomps{nx},[1:pcs],savedat,[3 45],'y',0,3);
SpecCoModPlotNew('sources1.set',fullpaths{nx},gdcomps{nx},[1:pcs],savedat,[3 50],'n',0,[]);
SpecCoModPlotNew('Emo-HP-1-200.set',fullpaths{nx},comps,[1:pcs],savedat,[3 125],'n',0,[]);
SpecCoModPlotNew('Emo-HP-1-248.set',fullpaths{nx},comps,[1:pcs],savedat,[3 125],'n',0,[]);
 gdcomps{2} =[    1     2     4     5     6     7     8     9    10    13    14    15    17   20    24]; % no filter
SpecCoModPlotNew('sources.set',fullpaths{nx},gdcomps{nx},gdcomps{nx},savedat,[3 50],'n',0,[],1);

[corr] = PlotWinvCorr(wv1,wv2,locs,locs,30);
 comps

subjlist = [2:21,23:31,33:35];
[emosigs,emoPs] = SigEmoShifts(savedat,newpaths,subjlist,emos);
save /data/common2/emotion/SigEmoShiftsMuscDecomp.mat emosigs emoPs
%%%%%%%%%%%%%%%%%%%%%%%%
nx = 4;
clist{nx} = [2,15];
PlotSpecs('AlldatSpecs.mat',fullpaths{nx}(1:end-5),clist,gdcomps,[1 45]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all templates from each member of a cluster
for clust = 1:length(facvec)
    unqs = 0;
    for nx = 1:length(facvec{clust})
        usedim = [];
        if ~isempty(facvec{clust}{nx})
            for im = 1:length(facvec{clust}{nx})
                if ~ismember(facvec{clust}{nx}(im),usedim)
                    usedim = [usedim facvec{clust}{nx}(im)];
                    unqs = unqs+1;
                end;
            end;        
        end;
    end;
    figure; row = round(sqrt(unqs)); col = ceil(sqrt(unqs)); pl = 1;viewnum = [1,2];zoom=1.1;
    col=9;
    for nx = 1:length(facvec{clust})
        usedim = [];
        if ~isempty(facvec{clust}{nx})
            for im = 1:length(facvec{clust}{nx})
                if ~ismember(facvec{clust}{nx}(im),usedim)
                    usedim = [usedim facvec{clust}{nx}(im)];
                    if ~isempty(comods{clust}{nx}{length(usedim)})
                    sbplot(row,col,pl)     
                    SpecCoModPlotNew('sources1.set',fullpaths{nx},gdcomps{nx},gdcomps{nx},facvec{clust}{nx}(im),savedat,[3 50],'y',0,[]);     
                    ph = text(3,3,['SUBJ ',int2str(nx)]); set(ph,'color','c');set(ph,'fontsize',14);
                    pl = pl+1;
                    clear tmpplot tmpjust tmpwts tmpjwts
                    tmpplot{nx} = {comods{clust}{nx}{length(usedim)}};
                    tmpjust{nx} = justcomps{clust}{nx};
                    tmpwts{nx} = {wtsmat1{clust}{nx}{length(usedim)}};
                    tmpjwts{nx} = jcwts{clust}{nx};
                    
    PlotCoModasDipoles(tmpplot,tmpjust,fullpaths,'sources1.set',row,col,pl,zoom,0,viewnum,tmpwts,tmpjwts,1,[]); pl = pl+2;
                    
                    end;
                end;
            end;        
        end;
    end;axcopy
    textsc(['Cluster ',int2str(clust)],'title');   
     set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
end;
set(gca,'xgrid','off'); set(gca,'xtick',[10:10:40]);set(gca,'xlim',[frs(1) frs(end)]);  
set(gca,'ytick',[0:2:8]);set(gca,'yticklabel',{0,[],4,[],8});
print /data/common4/emotion/Figs/AlphaIMex2.eps -depsc -adobe

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick out only highly weighted epochs to test for cross coherence.
for nx = 31:length(fullpaths)
    str = ['load ',fullpaths{nx},savedat,'.mat '];eval(str);         
    CrsCohSelected(datset,savedat,fullpaths{nx},rmepochs,gdcomps{nx});
    fprintf('\nSubject %s done. \n',int2str(nx));
end;

load /data/common4/emotion/AllCoModAlpha.mat 
load /data/common4/emotion/AllCoModBeta.mat 
load /data/common4/emotion/AllCoModGama.mat 
    for clust = 1:length(facvec)
    for nx = 4:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            str = ['load ',fullpaths{nx},savedat,'.mat '];eval(str);         
            for pos = 1:2
                CrsCohSelected(datset,savedat,fullpaths{nx},rmepochs,gdcomps{nx},pos);
            end;        
        end;      
    end;
    
    end;
%%%% Plot cross coh results
nx=2;
im = 5;  
for cond = 1:3
    PlotCrossCoh('sources1.set',['SelCrsCohIM',int2str(im),'.mat'],fullpaths{nx},gdcomps{nx},1,cond);
end;

% $$$ ic1 = 15; ic2 = 16; pos = 1;
% $$$ s = load([fullpaths{nx},nms{pos}]);
% $$$ figure; imagesc(s.timesout,s.freqsout,s.cohercell{ic1,ic2},[-.2 .2]);
% $$$ set(gca,'yscale','log');
% $$$ set(gca,'ytick',[10:10:s.freqsout(end)]);
% $$$ set(gca,'ydir','norm'); title(nms{pos});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a button press per emotion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
button = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
for nxx = 1:length(button)
    nx = button(nxx);
    EEG = pop_loadset( 'ButtonOnly.set',fullpaths{nx},'all');
    lats = zeros(1,0); evs = cell(1,0);
    for ev = 1:length(EEG.event)
        if ismember(EEG.event(ev).type,emos)
            lats(1,end+1) = EEG.event(ev).latency;
            evs{end+1} = EEG.event(ev).type;
        end;
    end; lats(end+1) = size(EEG.data,2);
    clear pressevents
    for k = 1:length(lats) - 1
        ndat = EEG.data(lats(k):lats(k+1)) - mean(EEG.data(lats(k):lats(k+1)));
        [normdat,outx] = movav(ndat,[],round(length(normdat)*.0082),[],[],[]);
        maxdata = max(normdat);
        thresh = maxdata*.55;  % sets threshold at 10% of max value
        r=1; presscount = 0;
        for g = 1:length(normdat)
            if normdat(1,g)>thresh & normdat(1,g-1)<thresh
                presscount = presscount + 1;
                %EEG.event(end+1).latency = g-50;r = r+1;  % make event 50 frames earlier than threshold
                %EEG.event(end).type = 'press';
                %EEG.event(end).Event_Type = 'Response';
            end;    
        end;
        pressevents(1,k) = presscount;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find % variance accounted for by all retained dims of spectra
for nx = 1:length(gdcomps)
    if ~isempty(gdcomps{nx})
        [pcared(1,nx) pv(nx,:)] = SpecCoModPvaf(savedat,fullpaths{nx},[],1);    
    end;
    fprintf('\nSubject %s done...\n',int2str(nx));
end;
save /data/common4/emotion/CoModPvafs.mat pcared pv
min(pcared)=  6.18; max(pcared)=23.22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a display of all pairwise IM weights distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 2;
clust1 =2; clust2 =5;
eval(strs{clust1})
fi1 = finalidx;
eval(strs{clust2})
fi2 = finalidx;clear finalmean finaltempls finalidx
%fi1 = {subjfacs{9}};
%fi2 = {subjfacs{10}};
for cl1 = 1:length(fi1)
    for cl2 = 1:length(fi2)
        clear quadprobs quadwinv quaddens
        for nx = 1:35
            x=fi1{cl1}(find(fi1{cl1}(:,1) == nx),:);
            y = x(:,2); % all IMS for that subj
            y=unique(y); % count each im only once
            xx=fi2{cl2}(find(fi2{cl2}(:,1) == nx),:);
            yy = xx(:,2); % all IMS for that subj
            yy=unique(yy); % count each im only once
            idxs1 = []; idxs2=[];
            for z=1:length(y)
                for zz = 1:length(yy)
                    idxs1 = [idxs1 y(z)];
                    idxs2 = [idxs2 yy(zz)];
                end;
            end;   
            if ~isempty(idxs1)
                [quadprobs{nx},quadwinv{nx},quaddens{nx}] = WtsDensity(savedat,fullpaths{nx},{idxs1,idxs2},100);
                for pr = 1:length(quadprobs{nx})
                    figure;
                    for q = 1:4
                        sbplot(2,2,q)
                        plot([0 16],[quadwinv{nx}{pr}(1,q) quadwinv{nx}{pr}(1,q)],'k-');hold on; 
                        plot([0 16],[quadprobs{nx}{pr}(1,q) quadprobs{nx}{pr}(1,q)],'m-');hold on; 
                        plot([0 16],[quadprobs{nx}{pr}(2,q) quadprobs{nx}{pr}(2,q)],'m-');hold on; 
                        for e = 1:length(emo2)
                        ph = bar(e,quaddens{nx}{pr}(e,q)); 
                        set(ph,'facecolor',cols(e,:));  
                        ph = text(e,0,emo2{e});
                        set(ph,'rotation',90);
                        end;
                        set(gca,'xlim',[0 16]);
                    end;
                    textsc(['Subject ',int2str(nx)],'title');
                end;
            end;
                %PlotIMwtsInteractions(savedat,fullpaths{nx},{idxs1,idxs2},1,'trial',200); 
            %textsc(['Subject ',int2str(nx),'; T',int2str(cl1),'-G',int2str(cl2)],'title');            
        end;
    end;
end;
% plot emo medians instead

figure; hold on; im = 7;
for em = 1:15
    ph = plot(emomeans{2}(im,em),em,'.'); set(ph,'color',cols(em,:));set(ph,'markersize',20);
          ph = text(emomeans{2}(im,em),em,emo2{em});
      set(ph,'color',cols(em,:)); set(ph,'fontsize',15);
end;
 set(gca,'xlim',[0 18.5])
ph = plot([get(gca,'xlim')],[0 0],'k-');
 set(gca,'fontsize',15);
 ylabel(['Mean IM',int2str(im),' Weight'])
 set(gca,'xticklabel',[]);

find(ismember(allidx(:,3),allidx2(:,3)))
% find indices of all clustered IMs-----------------------
for nx = 1:max(allidx(:,1))
    onesubj = allidx(find(allidx(:,1) == nx),2);
    useims{nx} = unique(onesubj');
end; clear allidx
% useims can be a vector or a cell array of 2 vectors to plot vs each other
PlotIMwtsInteractions(savedat,fullpaths{nx},useims{nx},1);

[distchar,alldatmat,MI,fullxcorr,emxcorr] = DescriptMat(savedat,fullpaths,subjlist,'on');
addpath('/home/jason/mbin/')
[mo,ord] = arr3(MI{2})
%%%% plot the measures by cluster:
for cls = 1:length(finalidx)
allmeas = zeros(6,0);
for nx = 1:35
    emcorrs = zeros(15,size(emxcorr{nx},2),0);
    modes = [];means = [];medians = []; vars = [];skews = [];kurts = [];
    x=finalidx{cls}(find(finalidx{cls}(:,1) == nx),:);
    y = x(:,2); % all IMS for that subj
    y=unique(y); % count each im only once    
    for z = 1:length(y)
        modes = [modes alldatmat{nx}(1,y(z))];
        means = [means alldatmat{nx}(2,y(z))];
        medians = [medians alldatmat{nx}(3,y(z))];
        vars = [vars alldatmat{nx}(4,y(z))];
        skews = [skews alldatmat{nx}(5,y(z))];
        kurts = [kurts alldatmat{nx}(6,y(z))];
        for em = 1:size(emxcorr{nx},3)
            emcorrs(em,:,z) = emxcorr{nx}(y(z),:,em)/emxcorr{nx}(y(z),y(z),em);
        end;
    end;
    en = size(allmeas,2);
    allmeas(1,en+1:en+length(modes)) = modes;
    allmeas(2,en+1:en+length(modes)) = means;
    allmeas(3,en+1:en+length(modes)) =  medians;
    allmeas(4,en+1:en+length(modes)) = vars;
    allmeas(5,en+1:en+length(modes)) = skews;
    allmeas(6,en+1:en+length(modes)) =  kurts;
end;
ttls = {'Mode','Mean','Median','Variance','Skewness','Kurtosis'};
figure; row = 2; col = 3; % no obvious difference in any measure bet clusts
for m = 1:size(allmeas,1)
    sbplot(row,col,m)
    hist(allmeas(m,:),30); hold on;
    ph = plot([0 0],[get(gca,'ylim')],'r-');
    title(ttls{m})
end;
textsc(['Beta Cluster ',int2str(cls)],'title');
end;
%  plot mutual information scores between clusters:-----------------
clust1=2; clust2 = 3; clear allmis
eval(strs{clust1})
fi1 = finalidx;
eval(strs{clust2})
fi2 = finalidx; pl = 1;
figure; pl = 1;
row = round(sqrt(length(fi1)*length(fi2))); 
col = ceil(sqrt(length(fi1)*length(fi2))); % no obvious difference in any measure bet clusts
for cls1=1:length(fi1)
    for cls2 = 1:length(fi2)
        mis = [];
        for nx = 1:35
            x=fi1{cls1}(find(fi1{cls1}(:,1) == nx),:);
            y = x(:,2); % all IMS for that subj
            y=unique(y); % count each im only once
            xx=fi2{cls2}(find(fi2{cls2}(:,1) == nx),:);
            yy = xx(:,2); % all IMS for that subj
            yy=unique(yy); % count each im only once
            subjmis = MI{nx}(y,yy);
            subjmis = reshape(subjmis,[1 size(subjmis,1)*size(subjmis,2)]);
            mis = [mis subjmis];
        end;
        sbplot(row,col,pl); 
        hist(mis,30); hold on;
        ph = plot([0 0],[get(gca,'ylim')],'r-');
        title(['T',int2str(cls1),'-A',int2str(cls2)]);
        allmis{pl} = mis; pl = pl+1;
    end;
end;

        
% weights distribution ideas:
subtract mean and find: correlation
variance correlation: like power coherence
shift in lag between time courses. 
test probability of interaction in partiucular quadrants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot highlighted activations for specified ICs/IMs
nx=2; ics = [9,11,20]; ims = [15,16,19]; stdcut = 2;
PlotCoModActs(datset,fullpaths{nx},ics,ims,savedat,stdcut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = find(freqs > 3);
[peakfreqs peakmags] = PeakFreqs(clustfacs(:,fr),freqs(fr));
[peakfreqs peakmags] = PeakFreqs(mnspecs(:,fr),freqs(fr));
 
    histvals = zeros(0,4);
    histvals1 = zeros(0,1);
    histvals2 = zeros(0,1);
    newkeep = zeros(0,size(kptk,2));
    leftfacs = zeros(0,length(fr));
    leftmns = zeros(0,length(fr));
    noalph = zeros(0,length(fr));
    noharm = zeros(0,length(fr));
    noalphidx = zeros(1,0);
    lenalph = zeros(1,0);
    lenharm = zeros(1,0);
    for sp = 1:size(peakfreqs,1)
        alph = find(peakfreqs(sp,:)  > 8 & peakfreqs(sp,:) < 13);
        if ~isempty(alph)
            harm = find(peakfreqs(sp,:)  > 15 & peakfreqs(sp,:) < 25);
                lenalph(1,end+1) = length(alph);
            if ~isempty(harm)                
                histvals1(end+1) = peakfreqs(sp,harm(1))/peakfreqs(sp,alph(1));
                lenharm(1,end+1) = length(harm);
            else
                noharm(end+1,:) = mnspecs(sp,fr);
            end;
        else
            noalph(end+1,:) = mnspecs(sp,fr);
            noalphidx(end+1) = sp;
        end;
    end;
        if peakfreqs(sp,1) > 8.5 & peakfreqs(sp,1) < 12.5
            if peakfreqs(sp,2) > peakfreqs(sp,1)*2-5 & peakfreqs(sp,2) < peakfreqs(sp,1)*2+5
                histvals(end+1,1) = peakfreqs(sp,1);
                histvals(end,2) = peakfreqs(sp,2)/2;
                histvals(end,3) = peakfreqs(sp,2);
                histvals(end,4) = peakfreqs(sp,1)*2;
                histvals1(end+1) = peakfreqs(sp,2)/peakfreqs(sp,1);
                histvals2(end+1) = peakmags(sp,1) - peakmags(sp,2);
                newkeep(end+1,:) = kptk(sp,:);
            else
                leftfacs(end+1,:) = clustfacs(sp,fr);
                leftmns(end+1,:) = mnspecs(sp,fr);
            end;            
        end;
    end;
    figure; sbplot(2,2,1);plot(freqs(fr),noalph)
    sbplot(2,2,2); plot(freqs(fr),noharm)
    histvals1(find(histvals1==0)) = [];
    histvals2(find(histvals2==0)) = [];
    
    figure; hist(histvals1,50)
    figure; hist(histvals2,50)
    figure; imagesc(histvals)
    plotcomps = cell(1,length(gdcomps));
    for s = 1:length(noalphidx)
        plotcomps{kptk(noalphidx(s),1)}(end+1) = onebig{kptk(noalphidx(s),1)}{kptk(noalphidx(s),2)};
    end;
    for nx = 1:length(plotcomps)
        plotcomps{nx} = unique(plotcomps{nx});
    end;
    
    PlotSomeMaps('sources1.set',fullpaths,plotcomps,'Components with No Alpha Peak');
    densargs = {'mrislices',[63:-8:-25],'mriview','top'};
    figure;PlotDipoles('sources1.set',fullpaths,plotcomps,[],'on','on','alldistance','r','off',densargs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster emotion indices across subjects
subjlist = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
subjlist = [13:20]; % no button press (apart from the first one)
subjlist = [1,2,4:6,8:12,14,17:21,23,25:30,31,33,34,35]; % all 'good' subjects (ones that said they got into it)
subjlist = [1:21,23:34,35];  % all but mr72-2
subjlist = [1,3:9,12,14,16,17,19,21,22,23,24,26,27,29,33]; % females
subjlist = [1,4:6,8,9,12,14,17,19,21,23,26,27,29,33]; % 'good' females
subjlist = [2,10,11,13,15,18,20,25,28,30,31,32,34,35]; % males
subjlist = [2,10,11,18,20,25,28,30,31,34]; % 'good' males
subjlist = [1:35]; clustby = 'e'; % 'e'-emotion, 's'-spectra, 'b'-both
subjlist = [2:21,23:35]; 

% 'e'-emotion, 's'-spectra, 'b'-both, 't' - template
[clustmat,newkeep] = ClustEmoCtx(subjlist,fullpaths,savedat,'t',1);
for p = 1:size(clustmat,1)
    clustmat(p,:) = clustmat(p,:) - mean(clustmat(p,:));
end;

alldist = pdist(clustmat, 'euclidean');
links = linkage(alldist,'ward');%avg>ward same as wtd, but what is that?
figure;[hnd,idx,perm]=  dendrogram(links,30);
figure; clear ctxclusts facclusts
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx))); 
for cls = 1:max(idx)
    onecls = find(idx == cls);
    ctxclusts{cls} = clustmat(find(idx == cls),:);
    facclusts{cls} = newkeep(find(idx == cls),:);
    for nx = 1:length(gdcomps)
        rr = onecls(newkeep(onecls,1) == nx);% find indices for current subj
        facvec{cls}{nx} = newkeep(rr,2)';
    end;
    sbplot(row,col,cls)
    ph = plot(freqs,clustmat(find(idx == cls),:)');hold on;
    set(ph,'color','g');
    ph = plot(freqs,mean(clustmat(find(idx == cls),:),1),'linewidth',2);
    set(ph,'color','b');
    set(gca,'xlim',[3 freqs(end)]); set(gca,'ylim',[-2000 2000]);
    title([int2str(cls),'-',int2str(length(facclusts{cls}))]);
    set(gca,'yticklabel',[]);
end;ctxclusts
axcopy
 figure;ph = plot(freqs,mean(clustmat(find(idx == cls),:),1),'linewidth',2);
 set(gca,'yticklabel',[]);
  axis('off');
 set(ph,'color','m');
print /home/julie/spec7.jpg -djpeg
% plot resulting 'cluster'
%[emomat,kpp] = ClustEmoCtx(subjlist,fullpaths,savedat,'e',1);

figure;row=round(sqrt(length(facclusts))); 
col=ceil(sqrt(length(facclusts)));  cols = jet(15);clear tts
for cl = 1:length(ctxclusts)
    sbplot(row,col,cl)
    stds = std(emomat(find(idx == cl),:));
    ph = plot(mean(emomat(find(idx == cl),:),1),'k-','linewidth',2);hold on;
    set(gca,'xlim',[0 16]);
    plot([get(gca,'xlim')],[0 0],'g-');
    %for em = 1:size(emomat,2)
    %    tts(1,em) = ttest(emomat(find(idx == cl),em),0,.05);
    %end;    
    %ph = plot(emomat(find(idx == cl),:)','k-','linewidth',1.5);hold on;
    for em = 1:size(emomat,2)
        tts = ttest(emomat(find(idx == cl),em),0,.05);
        if tts == 1
            %ph = plot(em,yl - yl*.2,'*');
            %set(ph,'markersize',6); set(ph,'color',[1 .7 0]);        
            ph = plot(em,mean(emomat(find(idx == cl),em),1),'.'); hold on;
            set(ph,'markersize',24); set(ph,'color',cols(em,:));        
        end;    
    end;    
    %for em = 1:length(stds)
    %    ph = plot([em em],[mean(emomat(find(idx == cl),em),1)-stds(em) mean(emomat(find(idx == cl),em),1)+stds(em)],'k-');
    %    set(ph,'color',[1 .7 0]);
    %end; yl = get(gca,'ylim');
    %sigpnts = find(tts);
    %if ~isempty(sigpnts)
    %    for hh = 1:length(sigpnts)
    %        ph = plot([sigpnts(hh) sigpnts(hh)],[yl(2)-yl(2)*.25 yl(2)-yl(2)*.25],'r*');
    %        set(ph,'markersize',8);set(ph,'color',[.8 0 1]);
    %    end;   
    %end;    
    %title([int2str(cl),'-',int2str(length(facclusts{cl}))]);
    title(int2str(cl));
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); 
    mn = min(mean(emomat(find(idx == cl),:))); mx = max(mean(emomat(find(idx == cl),:)));
    set(gca,'ylim',[mn+mn*.1 mx+mx*.1]);
end;axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

% cluster emo vectors within each cluster
for cl = 1:length(ctxclusts)
    newemomat = emomat(find(idx == cl),:);
    alldist = pdist(newemomat, 'correlation');
    links = linkage(alldist,'average');%avg>ward same as wtd, but what is that?
    figure;[hnd,idx2,perm]=  dendrogram(links,2);
    figure;
    row=round(sqrt(max(idx2))); 
    col=ceil(sqrt(max(idx2))); 
    for cls = 1:max(idx2)
        sbplot(row,col,cls)
        plotmean = mean(emomat(find(idx2 == cls),:),1);
        stds = std(emomat(find(idx2 == cls),:));
        for em = 1:size(emomat,2)
            tts(1,em) = ttest(emomat(find(idx2 == cls),em),0,.01);
        end;    
        ph = plot(plotmean,'r-','linewidth',2);hold on;
        set(gca,'xlim',[1 15]);
        plot([get(gca,'xlim')],[0 0],'g-');
        for em = 1:length(stds)
            ph = plot([em em],[plotmean(1,em)-stds(em) plotmean(1,em)+stds(em)],'k-');
            set(ph,'color',[1 .7 0]);
        end;     set(gca,'ylim',[min(plotmean) max(plotmean)]);
        yl = get(gca,'ylim');
        sigpnts = find(tts);
        if ~isempty(sigpnts)
            for hh = 1:length(sigpnts)
                ph = plot([sigpnts(hh) sigpnts(hh)],[yl(2)-yl(2)*.25 yl(2)-yl(2)*.25],'r*');
                set(ph,'markersize',8);set(ph,'color',[1 .5 0]);
            end;   
        end; 
    end;
end;
% Where are these 'factors located'? Need to look at comp(trial) weights
for nxx = 1:length(subjlist)
    nx = subjlist(nxx);
    s = load([fullpaths{nx},savedat,'Stuff.mat']); 
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[25 inf],[],0);    
    ws = wts*sph;    winv = pinv(ws); 
    activations = ws * data;   clear wts sph ws winv meanwts cprms
    for fac = 1:size(activations,1)
        totmeans(nx,fac) = sqrt(mean(activations(fac,:).^2));
        for cp = 1:length(gdcomps{nx})
            cprms(fac,cp) = sqrt(mean(activations(fac,(cp-1)*sum(s.dstrials)+1:cp*sum(s.dstrials)).^2));
            for em = 1:length(s.dstrials)
                meanwts(em,cp,fac) = sqrt(mean(activations(fac,(cp-1)*sum(s.dstrials(1:em))+1:cp*sum(s.dstrials(1:em))).^2));
            end; 
        end;
    end;
    allmeans{nx} = meanwts;
    allcps{nx} = cprms;
    fprintf('.');
end;

emolist = [6,10,11]; % 10,11,13

for clust = 1:length(facvec)
    clear allbigs
    for nx = 1:length(facvec{cls})
        if ~isempty(facvec{cls}{nx})
            colcps = zeros(0,2); colcomps = [];
            for fc = 1:length(facvec{clust}{nx})
            %for emm= 1:length(emolist)
            %    em = emolist(emm);
            %for cp = 1:size(meanwts,2)
            %if meanwts(em,cp,facvec{cls}{nx}) > totmeans(nx,facvec{cls}{nx})
            %    colcps(end+1,:) = [em cp];
            %end;   
            %end;
            %end;
            %allbigs{nx} = gdcomps{nx}(find(max(allcps{nx}(facvec{clust}{nx}(fc),:))));
            colcomps = [colcomps gdcomps{nx}(find(allcps{nx}(facvec{clust}{nx}(fc),:) > totmeans(nx,facvec{clust}{nx}(fc))))];
            %allbigs{nx} = unique(gdcomps{nx}(colcps(find(colcps(:,1)),2)'));
            %impcps{nx} = colcps;
            end;
            allbigs{nx} = unique(colcomps);
        end;
    end;subjidx = zeros(1,0); new=1; clear complist
    for nx = 1:length(allbigs)    
        if ~isempty(allbigs{nx})             
            EEG = pop_loadset('sources.set' ,fullpaths{nx}); 
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            dipsources = EEG.dipfit.model(allbigs{nx}(1));
            tmplist = zeros(1,0);
            for w = 1:length(allbigs{nx})
                dipsources(1,w).posxyz = EEG.dipfit.model(allbigs{nx}(w)).posxyz;
                dipsources(1,w).momxyz = EEG.dipfit.model(allbigs{nx}(w)).momxyz;
                dipsources(1,w).rv = EEG.dipfit.model(allbigs{nx}(w)).rv;
                subjidx(1,end+1) = nx;
                tmplist(1,end+1) = allbigs{nx}(w);
            end;           
            complist{nx}{1} = tmplist;
            if new == 1
                allbesa = dipsources;new = 0;
            else
                allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            end;
            dipsources = []; 
        end; 
    end;
    %PlotCrossLines(complist,fullpaths,'sources.set');
    %ph =textsc(['Dipole Connections for Cluster ',int2str(clust)],'title');
    %set(ph,'color','r');    set(ph,'fontsize',14);
    %str = ['print /data/common4/emotion/Figs/CtxLinesCls',int2str(clust),'.jpg -djpeg'];
    
    
    optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','distance','subjind',subjidx,'methodparam',15);     
    ph =textsc(['Dipole Density for Cluster ',int2str(clust),' Weighted by Subject Proximity'],'title');
    set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/CtxDensCls',int2str(clust),'.jpg -djpeg'];
    eval(str); close; close;
end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DENSITY DIFFERENCE CALCULATION  *****************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find significant differences between modulation factor dipole densities
densargs = {'mrislices',[63:-8:-25],'mriview','top','geom',[4,3],'cmax',.001};% for 12 plots
%densargs = {'mrislices',[55:-9:-9],'mriview','top','geom',[4,3],'cmax',.001}; % for 8 plots
%densargs = {'mrislices',[-33:6:33],'mriview','side','geom',[4,3],'cmax',.001};% for gamma
%densargs = {'mrislices',[35:-10:-35],'mriview','side','geom',[4,3],'cmax',.001};% for gamma
names = {'Alpha','Beta','Gamma'};
clear mnspecs onebig clustenv clustfacs bigwts
%load /data/common4/emotion/CoModAlphaClusts.mat 
load /data/common4/emotion/AllCoModAlpha.mat 
allpairs = {[1 2],[1 3],[2 3]}; nm = 1;

%load /data/common4/emotion/BetaClusters.mat 
load /data/common4/emotion/AllCoModBeta.mat 
allpairs = {[1 2]}; nm = 2;

%load /data/common4/emotion/GammaClusters.mat 
load /data/common4/emotion/AllCoModGama.mat 
allpairs = {[1 2],[1 3],[2 3]};nm = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%******************************************************
for ap = 1:length(allpairs)
    clusts = allpairs{ap};  
    clear complist1 complist2
    for m = 1:length(clusts)
        cls = clusts(m);
        for nx = 1:length(gdcomps)
            pl = 1; usedim = []; jc = zeros(1,0);
            if ~isempty(find(newfacs{cls}(:,1) == nx))
                subtemps = newfacs{cls}(find(newfacs{cls}(:,1) == nx),:);
                for im = 1:size(subtemps,1)                  
                    currim = subtemps(im,2);
                    if ~ismember(currim,usedim)
                        len = length(find(subtemps(:,2) == currim));
                        jc(end+1:end+len) = subtemps(find(subtemps(:,2) == currim),3);
                        usedim = [usedim currim];
                        pl = pl+1;
                    end;
                end;
                if m == 1
                    complist1{nx} = unique(jc);% I think these should be unique...
                else
                    complist2{nx} = unique(jc);
                end;            
            else
                complist1{nx} = [];
                complist2{nx} = [];
            end;
        end;
    end;
    
    shuffnum = 450; % 450 barely works! done for all three clusters
    [minmask maxmask] = SigDensityDiff('sources1.set',fullpaths,complist1,complist2,shuffnum,.01);
    str = ['save /data/common4/emotion/',names{nm},'Clust',int2str(clusts(1)),'-',int2str(clusts(2)),'DensDiff.mat minmask maxmask shuffnum clusts'];eval(str)
%end;
    %str = ['load /data/common4/emotion/',names{nm},'Clust',int2str(clusts(1)),'-',int2str(clusts(2)),'DensDiff.mat '];eval(str)
    %PlotSigDensDiff('sources1.set',fullpaths,complist1,complist2,minmask,maxmask,densargs);
    %ph=textsc([names{nm},' Combine Clusters Diff ',int2str(clusts(1)),'-',int2str(clusts(2))],'title');
    %set(ph,'color','r');
    %tr = ['print /data/common4/emotion/Figs/',names{nm},'DensDiffcombine.jpg -djpeg'];eval(str)
    %str = ['print /data/common4/emotion/Figs/',names{nm},'DensDiff',int2str(clusts(1)),'-',int2str(clusts(2)),'-SIDE.jpg -djpeg'];eval(str)
    clear minmask maxmask
end;
    str = ['print /home/julie/Manuscripts/CoMod/figures/SuppFigs/DensDiffcbar.jpg -djpeg'];eval(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot activation time courses for different factor clusters (aligned to be phase-consistent)
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
load /data/common4/emotion/CoModAlphaClusts.mat 
load /data/common4/emotion/BetaClusters.mat 
load /data/common4/emotion/GammaClusters.mat 

ttls = {'Alpha ','Alpha ','Alpha ','Beta ','Beta ','Gamma ','Gamma ','Gamma '};
cutoff = 1.5; % std from mean of weights distribution
shiftmax = 50; xprod = .3;
for clust = 2:length(facvec)
    howmany = 0;
    for nx = 1:length(facvec{clust})
        howmany = howmany + length(facvec{clust}{nx});
    end;
    figure; row = round(sqrt(howmany)); col = ceil(sqrt(howmany)); pl = 1;
    allalign = zeros(0,256);
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            for fac = 1:length(facvec{clust}{nx})
                factor = facvec{clust}{nx}(fac);
                comps = onebig{nx}{factor};
                [keepepochs,hightrials,highvals] = FindHighWtTrials(datset,fullpaths{nx},savedat,factor,comps,cutoff);
                hightrials = squeeze(hightrials); hightrials = hightrials';
                onemat = zeros(0,length(hightrials));
                for d = 1:length(keepepochs)
                    onemat(end+1:end+size(keepepochs{d},3),:) = squeeze(keepepochs{d})';
                end;         
                [val theone] = max(highvals); % picks act with highest weight from all datasets
                [alignmat] = PhaseAlignActs(hightrials(theone,:),onemat,shiftmax,xprod);
                
                allalign(end+1:end+size(alignmat,1),:) = alignmat;
                
                 %sbplot(6,1,[1 5]); 
                 sbplot(row,col,pl); 
                 imagesc(alignmat(:,1:(size(alignmat,2) - shiftmax)));
                 set(gca,'xlim',[1 199]); pl = pl+1;
                title(['Subj ',int2str(nx),'; Fac ',int2str(facvec{clust}{nx}(fac))]);
                %set(gca,'xticklabel',[]); sbplot(6,1,6);plot(mean(alignmat,1)); set(gca,'xlim',[1 256]);
            end;
        end;
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    textsc(['Alpha Cluster ',int2str(clust)],'title');
    str = ['print /data/common4/emotion/AlphaClust',int2str(clust),'alignmat.eps -depsc'];eval(str)
    str = ['save /data/common4/emotion/AlphaClust',int2str(clust),'alignmat.mat allalign'];eval(str)

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% plot scatter plot and corresponding activations 
nx= 2; 
s = load([fullpaths{nx},savedat,'.mat']); 
sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
ws = wts*sph; winv = pinv(ws); clear wts sph ws

factors = [15,2]; numpnts = 10;comp = 8;
%%%%%%%%%%%%%%

[indices,x,y] = ScatterWts(winv,factors,numpnts);
% !! pause to pick points!
[seltrials,setnames] = RetrieveActs(datset,fullpaths{nx},comp,indices,savedat);
% now plot the selected trials using the x/y coords
% select desired axis
fulax = get(gca,'xlim');set(gca,'xlim',[-max(abs(fulax)) max(abs(fulax))]);
fulax = get(gca,'xlim');fulgn = fulax(2) - fulax(1);
fuly = get(gca,'ylim');set(gca,'ylim',[-max(abs(fuly)) max(abs(fuly))]);
fuly = get(gca,'ylim');fulht = fuly(2) - fuly(1);
corfacx = .025;% because figure axis extends beyond gca
corfacy = .025;% because figure axis extends beyond gca
%zeratx = abs(fulax(1)/fulgn);
%zeraty = abs(fuly(1)/fulht);
plen = .12; %length of activation plots
phgt = .06;
for sl = 1:size(seltrials,1)
    if x(sl) < 0
        lft = (abs(fulax(1)) + x(sl))/fulgn ;
    else
        lft = (abs(fulax(1)) + x(sl))/fulgn - corfacx;
    end;   
    if y(sl) < 0
        bot = (abs(fuly(1)) + y(sl))/fulht  - corfacy;
    else
        bot = (abs(fuly(1)) + y(sl))/fulht - corfacy;
    end;
    if lft > .85
        lft = .845;
    end;    
    rect = [lft,bot,plen,phgt];
    
    axes('position',rect);
    ph = plot(seltrials(sl,:),'k-');
    set(gca,'xlim',[1 size(seltrials,2)]);    set(ph,'linewidth',1);
    %title([setnames{sl},'-',int2str(sl)]);
    title(int2str(sl));
    axis('off')
end;
str = ['print /data/common4/emotion/Figs/Nx2Fc',int2str(factors(1)),'-',int2str(factors(2)),'Cp',int2str(comp),'Acts2.eps -depsc'];eval(str)
%%%%%%
% do pwelch on individual trials
EEG = pop_loadset('awe.set',fullpaths{nx});
clear pwr
for tr = 1:size(seltrials,1)
    [pwr(:,tr) frqs] = pwelch(seltrials(tr,:),EEG.srate,EEG.srate/2,EEG.srate*2,EEG.srate); 
end;
fr = find(frqs> 3&frqs < 25);
figure; cols = jet(size(seltrials,1));
for tr = 1:size(seltrials,1)    
    ph = plot(frqs(fr),pwr(fr,tr),'k');hold on;
    set(ph,'color',cols(tr,:));
end;
hipwr = 10*log10(hipwr);
[hipwr] = eegfilt(hipwr',EEG.srate,0,50);
medpwr = 10*log10(medpwr);
[medpwr] = eegfilt(medpwr',EEG.srate,0,50);
lopwr = 10*log10(lopwr);
[lopwr] = eegfilt(lopwr',EEG.srate,0,50);
s = load([fullpaths{nx},savedat,'.mat']); 
mnpwr = [hipwr(fr);medpwr(fr);lopwr(fr)];mnpwr=mean(mnpwr,1);
figure;
plot(frqs(fr),hipwr(fr),'r');hold on;
plot(frqs(fr),medpwr(fr),'g');
plot(frqs(fr),lopwr(fr),'b');
plot(frqs(fr),s.meanpwr(5,fr),'k--');
set(gca,'xlim',[frqs(fr(1)) frqs(fr(end))]);
str = ['print /data/common4/emotion/Figs/nx23SnglTrFFTs.eps -depsc'];eval(str)

%%%%%%
%%  Plot quadrants as indiv acts and mean (aligned) on a separate figure;
[indices,x,y] = ScatterWts(winv,factors,numpnts,[],[1 0 0]);
% !! pause to pick points!
[seltrials,setnames] = RetrieveActs(datset,fullpaths{nx},comp,indices,savedat);
figure;row = ceil(sqrt(size(seltrials,1)));col = ceil(sqrt(size(seltrials,1)));
cols = jet(numpnts);
for sl = 1:size(seltrials,1)
    sbplot(row,col,sl)
    ph = plot(seltrials(sl,:),'k-');
    set(ph,'color',cols(sl,:));
    set(gca,'xlim',[1 size(seltrials,2)]);    set(ph,'linewidth',2);
    title([setnames{sl},'-',int2str(sl)]);
    %title(int2str(sl));
    axis('off')
end;
sbplot(row,col,sl+1)
ph = plot(mean(seltrials,1),'k-');
set(gca,'xlim',[1 size(seltrials,2)]);    set(ph,'linewidth',1);
title('Straight Mean');     axis('off')

figure;
 t=[0:.00392:1]*1000;sl=10;
ph =plot(t,seltrials(sl,:),'k-') ;hold on;
    set(ph,'color','g');

plot(t,sin(pi*t+150)*6,'r-');hold on;
set(gca,'xlim',[t(1) t(end)]);
 set(gca,'xlim',[0 650]);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% plot scatter plots for all factor combinations
nx= 2; 
s = load([fullpaths{nx},savedat,'Stuff.mat']); 
sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
ws = wts*sph; winv = pinv(ws); clear wts sph ws
figure;row = ceil(sqrt(size(winv,2)));col = ceil(sqrt(size(winv,2)));pl=1;
row = row + row*1.5;col = col+col*1.5;
for fac1 = 1:size(winv,2)-1
    for fac2 = fac1+1:size(winv,2)
    sbplot(row,col,pl) 
    plot(winv(:,fac1),winv(:,fac2),'b.');pl=pl+1;hold on;
    plot([0 0],[get(gca,'ylim')],'r-');
    plot([get(gca,'xlim')],[0 0],'r-');
    title(['Fac ',int2str(fac1),'-Fac ',int2str(fac2)]);
    set(gca,'xlim',[-3 3]);set(gca,'ylim',[-3 3]);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%plot 3D histogram of co-modulation directions
nclust = 3;
for clust = 1:nclust
    clear complist wtsmat frqmat corecomps
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})
            cc = zeros(1,0);
            for fc = 1:length(facvec{clust}{nx})
                tmplist = zeros(1,0);
                for w = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(w),allbigs{nx}{facvec{clust}{nx}(fc)})
                        tmplist(1,end+1) = gdcomps{nx}(w);numcomps=numcomps+1;
                    end;                    
                end;      
                cc(1,end+1) = onebig{nx}{facvec{clust}{nx}(fc)};
                wtsmat{nx}{fc} = bigwts{nx}{facvec{clust}{nx}(fc)};
                %frqmat{nx}{fc} = frqcell{clust}{nx}(fc);
                complist{nx}{fc} = tmplist;
                fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                numfacs = numfacs+1;
            end;
            corecomps{nx} = cc;
            numsubjs(1,end+1) = nx;
        end;
    end; 
    frqmat = [];zoom = 1; wtlims = [1 2];frqlims = [];yon = 1;
    figure; row=2;col=2;pl=1; mxlim = 8.5; gry = .7;
    [angles] = PlotCrossLinesWts(complist,fullpaths,'sources1.set',corecomps,wtsmat,wtlims,frqmat,frqlims,row,col,pl,zoom,0,1); 
    [hh bins] = hist(angles(1,:),60);
    [hhz binz] = hist(angles(2,:),60);
    sbplot(2,2,3);
    ph=plot([-mxlim mxlim],[0 0],'k-');hold on;set(ph,'color',[gry gry gry]);
    ph=plot([0 0],[-mxlim mxlim],'k-');set(ph,'color',[gry gry gry]);
    for ln = 1:length(bins)
        x = cos(bins(ln))*hh(ln);hold on;
        y = sin(bins(ln))*hh(ln);    
        plot([-x x],[-y y],'r-')
    end;
    set(gca,'xlim',[-mxlim mxlim]);set(gca,'ylim',[-mxlim mxlim]);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gca,'box','off');xlabel('X Dir'); ylabel('Y Dir');
    sbplot(2,2,4);
    ph=plot([-mxlim mxlim],[0 0],'k-');hold on;set(ph,'color',[gry gry gry]);
    ph=plot([0 0],[-mxlim mxlim],'k-');set(ph,'color',[gry gry gry]);
    for ln = 1:length(bins)
        yz = cos(binz(ln))*hhz(ln); hold on;   
        z = sin(binz(ln))*hhz(ln);    
        plot([-yz yz],[-z z],'g-')
    end;
    set(gca,'xlim',[-mxlim mxlim]);set(gca,'ylim',[-mxlim mxlim]);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gca,'box','off');xlabel('Y Dir'); ylabel('Z Dir');
    ph=textsc(['Alpha Cluster ',int2str(clust)],'title');
    set(ph,'fontsize',16);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot straight dipoles or density
folder = which('pop_dipfit_settings');folder = folder(1:end-21);
delim  = folder(end);mri = [ folder 'standard_BESA' delim 'avg152t1.mat' ];
mri = load('-mat', mri);mmri = mri.mri;
for clst = 3:2:length(finalclustspecs)
    clear mri ms onebkprj oneorig shuffwtd winv wted data diffboot
    %clust = intclusts(clst); 
    clust = clst;
    new = 1;       denswt = zeros(1,0);    subjidx = zeros(1,0);  
     numsubjs = 0; numfacs = 0; numcomps = 0;plotfacs = zeros(0,99);clear complist
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})             
            for fc = 1:length(facvec{clust}{nx})
                tmplist = zeros(1,0);
                for w = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(w),allbigs{nx}{facvec{clust}{nx}(fc)})
                        tmplist(1,end+1) = gdcomps{nx}(w);numcomps=numcomps+1;
                   end;                    
                end;           
                complist{nx}{fc} = tmplist;
                fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                numfacs = numfacs+1;
            end;
            numsubjs = numsubjs+1;
        end;
    end; 
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})             
            EEG = pop_loadset('sources.set' ,fullpaths{nx}); 
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            for fc = 1:length(facvec{clust}{nx})
                dipsources = EEG.dipfit.model(gdcomps{nx}(1));
                for w = 1:length(allbigs{nx}{facvec{clust}{nx}(fc)})
                    dipsources(1,w).posxyz = EEG.dipfit.model(allbigs{nx}{facvec{clust}{nx}(fc)}(w)).posxyz;
                    dipsources(1,w).momxyz = EEG.dipfit.model(allbigs{nx}{facvec{clust}{nx}(fc)}(w)).momxyz;
                    dipsources(1,w).rv = EEG.dipfit.model(allbigs{nx}{facvec{clust}{nx}(fc)}(w)).rv;
                    subjidx(1,end+1) = nx;
                end;           
                if new == 1
                    allbesa = dipsources;new = 0;
                    nowtall = dipsources;
                else
                    allbesa(end+1:end+size(dipsources,2)) = dipsources; 
                    if fc == 1
                        nowtall(end+1:end+size(dipsources,2)) = dipsources;
                    end;                    
                end;
                dipsources = []; 
            end;
        end; 
    end;
    plotfacs(end+1:end+size(finalclustspecs{clust},1),:) = finalclustspecs{clust};
    figure; 
    sbplot(2,2,1);mydipplot(allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on','spheres','on','dipolelength',0,'dipolesize',20,'color',{'r'});
    view(0,90)
    sbplot(2,2,2);mydipplot(allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on','spheres','on','dipolelength',0,'dipolesize',20,'color',{'r'});
    view(0,0)
    sbplot(2,2,3);mydipplot(allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on','spheres','on','dipolelength',0,'dipolesize',20,'color',{'r'});
    view(90,0)
    subplot(2,2,4)
    plot(freqs,plotfacs);hold on;set(gca,'fontsize',16);
    ph = plot(freqs,mean(plotfacs,1),'k-');set(ph,'linewidth',2.5);
    set(gca,'xgrid','on');    set(gca,'xlim',[freqs(1) freqs(end)]);xlabel('Frequency (Hz)');
    ylabel('Relative Power');    title(['Co-Mod Cluster ',int2str(clust)]);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color','w');
end;

    optdipplot = {nowtall,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15);     
     ph =textsc(['Cls ' ,int2str(clust),'; ',int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components;  UN-weighted'],'title');
set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/Cls',int2str(clust),'Unwted.jpg -djpeg'];
    %eval(str); close; close;
     %[unwt mri] = dipoledensity( optdipplot, 'method','alldistance','methodparam',15);     
    
    optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','distance','methodparam',15,'subjind',subjidx); 
    %figure;dipoledensity( optdipplot, 'method','distance','subjind',subjidx,'methodparam',15,'weight',denswt); 
     ph =textsc(['Cls ' ,int2str(clust),'; ',int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components;  Weighted by subject distance'],'title');set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/Cls',int2str(clust),'SubjDist.jpg -djpeg'];
    %eval(str); close; close;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /data/common4/emotion/clustfacs.mat 
load /data/common4/emotion/PdistCoModsClusts.mat 
figure; 
row=round(sqrt(length(clustspecs))); 
col=ceil(sqrt(length(clustspecs))); 
for cls = 1:length(clustspecs)
    sbplot(row,col,cls)
    ph = plot(freqs,clustspecs{cls}');hold on;
    set(ph,'color','g');
    ph = plot(freqs,mean(clustspecs{cls},1),'linewidth',2);
    set(ph,'color','b');
    set(gca,'xlim',[3 freqs(end)]);
    title([int2str(cls),'-',int2str(size(clustspecs{cls},1))]);
end;
axcopy
ph = textsc(['Clusters of Factor Spectra (largest with *-1) by pdist/linkage/dendrogram'],'title');
set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the envelope summary for spectral factor clusters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cls = 1:2:length(clustenvs)
    figure; row = ceil(sqrt(size(clustenvs{cls},1))); col = round(sqrt(size(clustenvs{cls},1)))+1;
    for fac = 1:size(clustenvs{cls},1)
        sbplot(row,col,fac);
        plot(freqs,clustenvs{cls}(fac,1:99),'r-');hold on;
        plot(freqs,clustenvs{cls}(fac,100:198),'b-');hold on;
        ph = plot(freqs,clustmeans{cls}(fac,:),'k-');set(ph,'linewidth',2.5); 
        set(gca,'xlim',[3 50]);
        title(['Sb ',int2str(facclusts{cls}(fac,1)),' Fac ',int2str(facclusts{cls}(fac,2))]);
    end;
    sbplot(row,col,fac+1);
    plot(freqs,mean(clustenvs{cls}(:,1:99),1),'r-');hold on;
    plot(freqs,mean(clustenvs{cls}(:,100:198),1),'b-');hold on;
    ph = plot(freqs,mean(clustmeans{cls},1),'k-');set(ph,'linewidth',2.5); 
    title('Mean Cluster Envelope/Mean');    
    textsc(['Cluster ',int2str(cls)],'title');set(gca,'xlim',[3 50]);
end;        
%%%%%%%%% Plot individual factors %%%%%%%%%%%%%%%%%%%%%
for cls = 1:2:length(clustenvs)
    figure; row = ceil(sqrt(size(clustenvs{cls},1))); col = round(sqrt(size(clustenvs{cls},1)))+1;
    for fac = 1:size(clustenvs{cls},1)
        sbplot(row,col,fac);
        plot(freqs,clustspecs{cls}(fac,:),'b-');hold on;
        set(gca,'xlim',[3 50]);
        title(['Sb ',int2str(facclusts{cls}(fac,1)),' Fac ',int2str(facclusts{cls}(fac,2))]);
    end;
    sbplot(row,col,fac+1);
    plot(freqs,mean(clustspecs{cls},1),'r-');hold on;
    title('Mean Cluster Envelope/Mean');    set(gca,'xlim',[3 50]);
    textsc(['Cluster ',int2str(cls)],'title');
end;        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard deviation formula%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=zeros(1,0);
for x = 1:length(clustspecs{8}(1,:))
    x=clustspecs{8}(1,x);
    x=x-mean(clustspecs{8}(1,:));
    x=x^2;
    xx(1,end+1) = x;
end;
xx = sum(xx);
xx = xx/(length(clustspecs{8}(1,:)) - 1);
xx = sqrt(xx); xx
% vs RMS:
sqrt(mean(clustspecs{8}(1,:).^2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a subject emotion space by pca on subj winv
emos = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
for nxx = 1:length(incsubjs)
    figure;row=4; col=4;
    nx = incsubjs(nxx);
    str = ['load ',fullpaths{nx},savedat,'Stuff.mat'];eval(str);  
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
    ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws pvf pvfreq
    for fac = 1:size(winv,2)
    sbplot(row,col,fac); cols = jet(15); %%%%%%%%%%%  2D   %%%%%%%%%%%%%%%%%
    for em = 1:length(emos)
        ph = plot(winv(sum(dstrials(1:em-1))+1:sum(dstrials(1:em)),fac),'.'); hold on;
        set(ph,'markersize',5); set(ph,'color',cols(em,:));       
        %ph = text(pc(1,em),pc(2,em),emos{em});
        %set(ph,'color',cols(em,:)); 
    end;  
    end;
    %set(gca,'xgrid','on'); set(gca,'ygrid','on'); 
end;
figure; fac1 = 10; fac2 = 14;
    for em = 1:length(emos)
        ph = plot(winv(sum(dstrials(1:em-1))+1:sum(dstrials(1:em)),fac1),winv(sum(dstrials(1:em-1))+1:sum(dstrials(1:em)),fac2),'.'); hold on;
        set(ph,'markersize',5); set(ph,'color',cols(em,:));       
        %ph = text(pc(1,em),pc(2,em),emos{em});
        %set(ph,'color',cols(em,:)); 
    end;  

figure; cols = jet(15); %%%%%%%%%%%  3D   %%%%%%%%%%%%%%%%%
for em = 1:length(emos)
    ph = plot3(pc(1,em),pc(2,em),pc(3,em),'.'); hold on;
    set(ph,'markersize',24); set(ph,'color',cols(em,:));       
    ph = text(pc(1,em),pc(2,em),pc(3,em),emos{em});
    set(ph,'color',cols(em,:)); 
    set(gca,'xgrid','on'); set(gca,'ygrid','on'); set(gca,'zgrid','on');
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find percent variance accounted for for each back-proj'd fac-comp to original data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/KmeansClustCoMods2.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers

clear pvfreq pvf
for clust = 1:length(facvec)
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            str = ['load ',fullpaths{nx},'SpecCoModStuff.mat numtrials numframes freqs keeptrack rmepochs dstrials pcs savedat comment'];eval(str);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
            data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
            ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws pvf pvfreq
            backprojdat = winv*activations ;
            for fc = 1:length(facvec{clust}{nx}) 
                 x=[1:size(activations,1)];x(facvec{clust}{nx}(1,fc)) = [];
                %acts = activations; acts(facvec{clust}{nx}(1,fc),:) = 0;
                acts = activations; acts(x,:) = 0;
                backproj = winv*acts ;
                for cmp = 1:length(gdcomps{nx})
                    oneorig = backprojdat(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                    %oneorig = reshape(oneorig,1,size(oneorig,1)*size(oneorig,2));
                    onebkprj = backproj(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                    %onebkprj = reshape(onebkprj,1,size(onebkprj,1)*size(onebkprj,2));
                    pvfreq(cmp,:,fc) = 100-((var(onebkprj)./var(oneorig))*100);
                    pvf(cmp,fc) = 100-((sum(var(onebkprj))./sum(var(oneorig)))*100);
                end;                
            end;
            allnxpvfreq{clust}{nx} = pvfreq;
            allnxpvf{clust}{nx} = pvf;
        end;fprintf('.');
    end;
end;
%save /data/common4/emotion/pvafsALL.mat allnxpvf allnxpvf
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% BETTER!:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use ratio of back-proj var and total var instead OR RMS back projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALLEEG=[];EEG=[];
for clust = 1:length(facvec)
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            clear onebkprj winv ws wts data backproj activations acts plotprj sph oneorig eigvec sv backprojdat
            str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
            data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
            ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws pvf pvfreq
            backprojdat = winv*activations ; clear cmpscale rmsscale
            for fc = 1:length(facvec{clust}{nx}) 
                 x=[1:size(activations,1)];x(facvec{clust}{nx}(1,fc)) = [];
                acts = activations; acts(x,:) = 0;
                backproj = winv*acts ;
                for cmp = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(cmp),allbigs{nx}{facvec{clust}{nx}(fc)})
                        oneorig = backprojdat(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                        onebkprj = backproj(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                        cmpscale(cmp,fc) = var(onebkprj)/var(oneorig);
                        oneact = activations(facvec{clust}{nx}(1,fc),length(freqs)*(cmp-1)+1:length(freqs)*cmp);
                        oneact = sqrt(mean(oneact.^2));
                        rmswinv = sqrt(mean(winv(:,facvec{clust}{nx}(1,fc)).^2));
                        rmsscale(cmp,fc) = oneact*rmswinv;
                    else
                        cmpscale(cmp,fc) = 0;
                        rmsscale(cmp,fc) = 0;
                    end;                    
                end;                
            end;
            allrmss{clust}{nx} = rmsscale;
            allcmpscales{clust}{nx} = cmpscale;
        end;fprintf('.');
    end;
end;
save /data/common4/emotion/ProjWtCell.mat allcmpscales allrmss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  using the Perc pwr from each factor/component, weight a density plot for each cluster
load /data/common4/emotion/ProjWtCell.mat allcmpscales allrmss
intclusts = [1:13,15:19];%  [2,4,5,6,8,9,11,13,15,16,19];
folder = which('pop_dipfit_settings');folder = folder(1:end-21);
delim  = folder(end);mri = [ folder 'standard_BESA' delim 'avg152t1.mat' ];
mri = load('-mat', mri);mmri = mri.mri;
for clst = 1:length(intclusts)
    clear mri ms onebkprj oneorig shuffwtd winv wted data diffboot
    clust = intclusts(clst);  new = 1;       denswt = zeros(1,0);    subjidx = zeros(1,0);  
    clust = intclusts(clst);  numsubjs = 0; numfacs = 0; numcomps = 0;clear complist
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})             
            for fc = 1:length(facvec{clust}{nx})
                tmplist = zeros(1,0);
                for w = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(w),allbigs{nx}{facvec{clust}{nx}(fc)})
                        tmplist(1,end+1) = gdcomps{nx}(w);numcomps=numcomps+1;
                   end;                    
                end;           
                complist{nx}{fc} = tmplist;
                fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                numfacs = numfacs+1;
            end;
            numsubjs = numsubjs+1;
        end;
    end; 
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})             
            EEG = pop_loadset('sources.set' ,fullpaths{nx}); 
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            for fc = 1:1%length(facvec{clust}{nx})
                dipsources = EEG.dipfit.model(gdcomps{nx}(1));
                for w = 1:length(gdcomps{nx})
                    dipsources(1,w).posxyz = EEG.dipfit.model(gdcomps{nx}(w)).posxyz;
                    dipsources(1,w).momxyz = EEG.dipfit.model(gdcomps{nx}(w)).momxyz;
                    dipsources(1,w).rv = EEG.dipfit.model(gdcomps{nx}(w)).rv;
                    %denswt(1,end+1) = allcmpscales{clust}{nx}(w,fc);                    
                    denswt(1,end+1) = allrmss{clust}{nx}(w,fc);
                    subjidx(1,end+1) = nx;
                end;           
                if new == 1
                    allbesa = dipsources;new = 0;
                    nowtall = dipsources;
                else
                    allbesa(end+1:end+size(dipsources,2)) = dipsources; 
                    if fc == 1
                        nowtall(end+1:end+size(dipsources,2)) = dipsources;
                    end;                    
                end;
                dipsources = []; 
            end;
        end; 
    end;
    optdipplot = {nowtall,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15);     
     ph =textsc(['Cls ' ,int2str(clust),'; ',int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components;  UN-weighted'],'title');
set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/PdistUNwtCls',int2str(clust),'.jpg -djpeg'];
    eval(str); close; close;
     %[unwt mri] = dipoledensity( optdipplot, 'method','alldistance','methodparam',15);     
    
    optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',denswt); 
    %figure;dipoledensity( optdipplot, 'method','distance','subjind',subjidx,'methodparam',15,'weight',denswt); 
     ph =textsc(['Cls ' ,int2str(clust),'; ',int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components;  Weighted by RMS'],'title');set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/PdistWtdCls',int2str(clust),'.jpg -djpeg'];
    eval(str); close; close;
end;

    
    for shuff = 1:500
        denswt1 = shuffle(denswt);
        [wted mri] =dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',denswt1); 
        shuffwtd(:,:,:,shuff) = wted;
    end;
    for x = 1:size(shuffwtd,1)
        for y = 1:size(shuffwtd,2)
            for z = 1:size(shuffwtd,3)
                flims = squeeze(shuffwtd(x,y,z,:));
                valims = sort(flims); 
                diffboot(x,y,z,1) = valims(end-5);
                %diffboot(x,y,z,2) = valims(2);   % do only one side (upper)          
            end;
        end;
        fprintf('.');
    end;
    ms = wted;
    ms(find(ms<diffboot(:,:,:,1))) = 0;
    mri3dplot(ms,mmri);set(gcf,'color','w');
    ph = textsc(['Bootstrap masked weighted density; Cluster ',int2str(clust)],'title');set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/SigFacDensity',int2str(clust),'.jpg -djpeg'];
    eval(str);  close;    
end;

    
    fprintf(['Cluster ',int2str(clust),' done...']);
    %unwt = unwt/max(max(max(unwt)));
    %wted = wted/max(max(max(wted)));
    %%  opt subtraction of densities
    for zz = 1:size(wted,3)
        diffdens(:,:,zz) = unwt(:,:,zz) - wted(:,:,zz);
        for bs = 1:100
            shuffun = unwt(:,:,zz);  
            shuffun = shuffle(shuffun,1);shuffun = shuffle(shuffun,2);
            shuffwt = unwt(:,:,zz);  
            shuffwt = shuffle(shuffwt,1);shuffwt = shuffle(shuffwt,2);
            boots(:,:,bs) = shuffun - shuffwt;
        end;clear diffboot
        for x = 1:size(boots,1)
            for y = 1:size(boots,2)
                flims = squeeze(boots(x,y,:));
                valims = sort(flims); 
                diffboot(x,y,2) = valims(end-1);
                diffboot(x,y,1) = valims(2);
                
                %diffboot(x,y,2) = max(flims);
                %diffboot(x,y,1) = min(flims);
            end;
        end;
        ms = diffdens(:,:,zz);
        ms(find(ms>diffboot(:,:,1)&ms<diffboot(:,:,2))) = 0;
        maskdiff(:,:,zz) = ms;
        fprintf('.');
    end;
    figure;
    for z = 1:91
        subplot(10,10,z)
        imagesc(maskdiff(:,:,z)',[-.7 .7]);
        set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    end;
    %mymriplot(maskdiff,mmri);set(gcf,'color','w');
    mymriplot(diffdens,mmri);set(gcf,'color','w');
    ph = textsc(['unwted minus wted density; Cluster ',int2str(clust)],'title');set(ph,'color','r');    set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/EmoFacClustDensDiff',int2str(clust),'.jpg -djpeg'];
    eval(str);  close;    
end;
,'plotargs',{'cmax',.6,'cmap',jet}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Great, now make 'PlotCrossLines.m' plots to show specifically which components are co-modulated (densities do not reveal this because the densities could be from isolated components)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First find highest weighted comps for each factor of interest (for each cluster)
load /data/common4/emotion/clustfacs.mat  
load /data/common4/emotion/PdistCoModsClusts.mat 

combineclust = [1,3,27,30,33,35,43,57,81]; % alpha no shift
combineclust = [7,9,55,71,73,79,95]; % alpha shift
combineclust = [11,15,19,23,85]; % gamma
combineclust = [5,13,45,51,87]; % 20 Hz
combineclust = [21,25,75,91]; % 30 Hz

intclusts = [69:2:100];
numsubjs = zeros(1,0); numfacs = 0; numcomps = 0;plotfacs = zeros(0,99);clear complist
for clst = 1:length(combineclust)
    %clust = combineclust(clst);  
    clust = clst;  
numsubjs = zeros(1,0); numfacs = 0; numcomps = 0;plotfacs = zeros(0,99);clear complist
    for nx = 1:length(gdcomps)    
        if ~isempty(facvec{clust}{nx})             
            for fc = 1:length(facvec{clust}{nx})
                tmplist = zeros(1,0);
                for w = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(w),allbigs{nx}{facvec{clust}{nx}(fc)})
                        tmplist(1,end+1) = gdcomps{nx}(w);numcomps=numcomps+1;
                   end;                    
                end;           
                complist{nx}{fc} = tmplist;
                fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                numfacs = numfacs+1;
            end;
            numsubjs(1,end+1) = nx;
        end;
    end; 
    plotfacs(end+1:end+size(clustspecs{clust},1),:) = clustspecs{clust};
%end;
numsubjs = unique(numsubjs); numsubjs = length(numsubjs);
    PlotCrossLines(complist,fullpaths,'sources.set');
    subplot(2,2,4)
    %plot(freqs,clustenvs{clust}(:,1:99),'r-');hold on; % for envelope
    %plot(freqs,clustenvs{clust}(:,100:198),'b-');hold on;
    %ph = plot(freqs,mean(clustenvs{clust}(:,1:99),1),'k-');set(ph,'linewidth',2.5);
    %ph = plot(freqs,mean(clustenvs{clust}(:,100:198),1),'k-');set(ph,'linewidth',2.5);
    %plot(freqs,clustspecs{clust});hold on;set(gca,'fontsize',16);
    %ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',2.5);
    plot(freqs,plotfacs);hold on;set(gca,'fontsize',16);
    ph = plot(freqs,mean(plotfacs,1),'k-');set(ph,'linewidth',2.5);
    %set(gca,'ylim',[-3 8]); 
    set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);xlabel('Frequency (Hz)');
    ylabel('Relative Power');
    title(['Co-Mod Cluster ',int2str(clust)]);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
     ph =textsc([int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components'],'title');
     set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/AlphaClustsLines.jpg -djpeg'];eval(str)
    close
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find percent trials in each emotion (to find shifts)

[allposperc] = PercDistShift(incsubjs,savedat,fullpaths,emos,0);
allpercs = zeros(0,length(emos));  kt = zeros(0,2);
for nx = 1:length(allposperc)
    if ~isempty(allposperc{nx})
        allpercs(end+1:end+15,:) = allposperc{nx};
        kt(end+1:end+15,:) = [repmat(nx,[15 1]) [1:15]'];        
    end;
end;
alldist = pdist(allpercs, 'correlation'); % euc better than seuc
links = linkage(alldist,'complete');%ward best
figure;[hnd,idx,perm]=  dendrogram(links,40);
figure; clear ctxclusts facclusts
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx))); 
for cls = 1:max(idx)
    onecls = find(idx == cls);
    emoclusts{cls} = allpercs(onecls,:);
    facclusts{cls} = kt(onecls,:);  
    for nx = 1:length(gdcomps)
        facvec{cls}{nx} = kt(onecls(kt(onecls,1) == nx),2)';
    end;
    sbplot(row,col,cls)
    ph = plot([1:15],emoclusts{cls}');hold on;
    set(ph,'color','g');
    ph = plot([1:15],mean(emoclusts{cls},1),'linewidth',2);
    set(ph,'color','r');    set(gca,'xlim',[1 15]);
    title([int2str(cls),'-',int2str(size(emoclusts{cls},1))]);
    plot([get(gca,'xlim')],[0 0],'k-.','linewidth',2);
end;emoclusts
axcopy
% plot sig diffs of mean
figure; clear ctxclusts facclusts
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx)));  cols = jet(15);
for cls = 1:max(idx)
    sbplot(row,col,cls)
    ph = plot([1:15],mean(emoclusts{cls},1),'linewidth',2); hold on;
    set(ph,'color','k');    set(gca,'xlim',[1 15]);
    yl = max(abs(mean(emoclusts{cls},1)))+max(abs(mean(emoclusts{cls},1)))*.1;
    set(gca,'ylim',[-yl yl]);
    for em = 1:size(emoclusts{cls},2)
        tts = ttest(emoclusts{cls}(:,em),0,.05);
        if tts == 1
            %ph = plot(em,yl - yl*.2,'*');
            %set(ph,'markersize',6); set(ph,'color',[1 .7 0]);        
            ph = plot(em,mean(emoclusts{cls}(:,em),1),'.');
            set(ph,'markersize',24); set(ph,'color',cols(em,:));        
        end;    
    end;    
    title([int2str(cls),'-',int2str(size(emoclusts{cls},1))]);
    plot([get(gca,'xlim')],[0 0],'k-.','linewidth',2);
end;axcopy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternatively, plot the emotion vectors for spectral clusters
[allposperc, allnegperc] = PercDistShift(incsubjs,savedat,fullpaths,emos,1);
allpercs = zeros(0,length(emos));  kt = zeros(0,2);
for nx = 1:length(allposperc)
    if ~isempty(allposperc{nx})
        for p = 1:size(allposperc{nx},1)
        allpercs(end+1,:) = allposperc{nx}(p,:);
        kt(end+1,:) = [nx p];        
        allpercs(end+1,:) = allposperc{nx}(p,:)*-1;
        %allpercs(end+1,:) = allnegperc{nx}(p,:);
        kt(end+1,:) = [nx p];   
        end;
    end;
end;
figure;  
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx))); 
for cls = 1:max(idx)
    onecls = find(idx == cls);
    assocemos{cls} = allpercs(onecls,:);
    sbplot(row,col,cls)
    ph = plot([1:15],assocemos{cls}');hold on;
    set(ph,'color','g');
    ph = plot([1:15],mean(assocemos{cls},1),'linewidth',2);
    set(ph,'color','r');    set(gca,'xlim',[0 16]);
    title([int2str(cls),'-',int2str(size(assocemos{cls},1))]);
    plot([get(gca,'xlim')],[0 0],'k-.','linewidth',2);
end;
axcopy
% plot sig diffs of mean
figure; clear ctxclusts facclusts
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx)));  cols = jet(15);
for cls = 1:max(idx)
    sbplot(row,col,cls)
    ph = plot([1:15],mean(assocemos{cls},1),'linewidth',2); hold on;
    set(ph,'color','k');    set(gca,'xlim',[0 16]);
    yl = max(abs(mean(assocemos{cls},1)))+max(abs(mean(assocemos{cls},1)))*.1;
    set(gca,'ylim',[-yl yl]);
    for em = 1:size(assocemos{cls},2)
        tts = ttest(assocemos{cls}(:,em),0,.01);
        if tts == 1
            %ph = plot(em,yl - yl*.2,'*');
            %set(ph,'markersize',6); set(ph,'color',[1 .7 0]);        
            ph = plot(em,mean(assocemos{cls}(:,em),1),'.');
            set(ph,'markersize',24); set(ph,'color',cols(em,:));        
        end;    
    end;    
    title([int2str(cls),'-',int2str(size(assocemos{cls},1))]);
    ph = plot([get(gca,'xlim')],[0 0],'k-.','linewidth',2);
    set(ph,'color',[.5 .5 .5]);
end;axcopy

for cls = 1:max(idx)
[P,atab,stats]  = anova1(assocemos{cls},emos,'off');
figure; mcomp = multcompare(stats, 'alpha',.05);
textsc(['Cluster ',int2str(cls),'; P == ',num2str(P)],'title');
end;

, 'ctype','bonferroni'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot of mean spectrum plus back-projs of factors (envelope)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load /data/common4/emotion/KmeansClustCoMods.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers

for clust = 1:length(facvec)
    p=1;
    for nx = 1:length(gdcomps{nx})
        if ~isempty(facvec{clust}{nx})
            str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
            data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
            ws = wts*sph;    activations = ws*data;    winv = pinv(ws); 
            cols = jet(length(gdcomps{nx}));
            for fac = 1:length(facvec{clust}{nx})                
                figure; row = ceil(sqrt(length(gdcomps{nx}))); col = ceil(sqrt(length(gdcomps{nx}))); pl = 1;
                sbplot(row,col,pl)
                ph = plot(freqs,clustspecs{clust}(p,:),'k-'); p = p+1;
                set(ph,'linewidth',2);set(ph,'color',[.9 .5 0]);
                set(gca,'xlim',[3 50]);
                col = ceil(sqrt(length(gdcomps{nx})));
                x=[1:size(activations,1)];x(facvec{clust}{nx}(1,fac)) = [];                
                acts = activations; acts(x,:) = 0;
                backproj = winv*acts ; clear plotprj
                for cmp = 1:length(gdcomps{nx})
                    sbplot(row,col,pl)
                    onebkprj = backproj(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                    for trls = 1:size(onebkprj,1)
                        plotprj(trls,:) = onebkprj(trls,:) + meanpwr(cmp,:);
                    end;
                    ph = plot(freqs,plotprj,'k-'); hold on;
                    set(ph,'color',cols(cmp,:));
                    ph = plot(freqs,meanpwr(cmp,:),'k-'); 
                    set(ph,'linewidth',2);set(gca,'xlim',[3 50]);
                    fprintf('.');pl = pl+1; set(gca,'xticklabel',[]);
                    title(['Sb ',int2str(nx),';Fc ',int2str(facvec{clust}{nx}),';Cp ',int2str(gdcomps{nx}(cmp))]);
                end;
                str = ['print /data/common4/emotion/Figs/Cls',int2str(clust),'Sb',int2str(nx),'Fc',int2str(facvec{clust}{nx}(fac)),'BackProj.jpg -djpeg'];eval(str)
                close
            end;
        end;
    end;   
end;
% easier
intclusts = [1:30];
for clst = 1:length(intclusts)
    clust = intclusts(clst);  
    for nx = 1:length(gdcomps)
        if ~isempty(facvec{clust}{nx})
            str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
            for fac = 1:length(facvec{clust}{nx})               
                
                SpecCoModPlot('sources.set',fullpaths{nx},gdcomps{nx},savedat,[0 50],pcs,numtrials,numframes,freqs,'y',1,facvec{clust}{nx}(fac));
                str = ['print /data/common4/emotion/Cls',int2str(clust),'Sb',int2str(nx),'Fc',int2str(facvec{clust}{nx}(fac)),'BackProj.jpg -djpeg'];eval(str);  close                
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot cluster summary figure with highest var comp from each factor
for clust = 1:length(facvec)
    howmany = 0;
    for nx = 1:length(facvec{clust})
        if ~isempty(facvec{clust}{nx})
            howmany = howmany + length(facvec{clust}{nx});
        end;
    end;
    figure; pl = 1;row = round(sqrt(howmany)); col = ceil(sqrt(howmany));
    cols = jet(howmany);
    for nx = 1:length(gdcomps)
        if ~isempty(facvec{clust}{nx})
            clear onebkprj winv ws wts data backproj activations acts plotprj sph
            str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
            data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
            ws = wts*sph;    activations = ws*data;    winv = pinv(ws); 
            for fac = 1:length(facvec{clust}{nx})
                x=[1:size(activations,1)];x(facvec{clust}{nx}(1,fac)) = [];                
                acts = activations; acts(x,:) = 0;
                backproj = winv*acts ; clear plotprj
                cmp = find(cell2mat(onebig{nx}(facvec{clust}{nx}(fac))) == gdcomps{nx});
                onebkprj = backproj(:,length(freqs)*(cmp-1)+1:length(freqs)*cmp); 
                for trls = 1:size(onebkprj,1)
                    plotprj(trls,:) = onebkprj(trls,:) + meanpwr(cmp,:);
                end;
                subplot(row,col,pl)
                ph = plot(freqs,plotprj,'k-'); hold on; 
                set(ph,'color',cols(pl,:));pl = pl+1;
                ph = plot(freqs,meanpwr(cmp,:),'k-'); 
                set(ph,'linewidth',2);set(gca,'xlim',[3 50]);
                fprintf('.'); set(gca,'xticklabel',[]);
                title(['Sb ',int2str(nx),';Fc ',int2str(facvec{clust}{nx}(fac)),';Cp ',int2str(gdcomps{nx}(cmp))]);
            end;
        end;
    end;
    str = ['print /data/common4/emotion/Cls',int2str(clust),'BackProjs.jpg -djpeg'];eval(str);  close                
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find percent of positive trials for each subject and cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 0;subjlist = [1:35];
[allposperc] = PercDistShift(subjlist,savedat,fullpaths,emos,ploton);
load /data/common4/emotion/AllPercs.mat allposperc

clear subjfacs
for nx = 1:length(allposperc)
    if ~isempty(allposperc{nx})
        for em = 1:size(allposperc{nx},1)
            subjfacs{em}{nx} = find(allposperc{nx}(em,:) > .2); 
        end;        
    end;
end;
cols = jet(15);
for em = 1:size(allposperc{nx},1)
    figure; pl = 1; row= 4; col = 4;    
    for nx = 1:length(subjfacs{em})
        if ~isempty(subjfacs{em}{nx})            
            sbplot(row,col,pl)
            for ff = 1:length(subjfacs{em}{nx})
                for e = 1:length(emos)
                    ph=plot(e,allposperc{nx}(e,subjfacs{em}{nx}(ff)),'k.'); hold on; 
                    set(ph,'markersize',15);
                    set(ph,'color',cols(e,:));
                    ph = text(e,.1,emos{e});hold on;
                    set(ph,'rotation',90);set(ph,'color',cols(e,:));
                end;
                plot(allposperc{nx}(:,subjfacs{em}{nx}(ff)),'k-'); hold on; 
                title(['Sb ',int2str(nx),' Fac ',int2str(subjfacs{em}{nx}(ff))]);
            end;pl = pl+1;
        end;
    end;
    textsc(['EMOTION: ',emos{em}],'title');
end;

                
pcmat = zeros(15,0);
for nx = 1:length(allposperc)
    if ~isempty(allposperc{nx})
        pcmat(:,end+1:end+size(allposperc{nx},2)) = allposperc{nx};
    end;
end;

pcdims = 10;% run ica on emoscores
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(pcmat,'pca',pcdims,'extended',1,'stop',1e-7);
 ws = weights*sphere; winv = pinv(ws);

% plot each dimension separately:
cols = jet(15);cols(10,:) = [.9 .9 0];
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure;  
for wv = 1:size(winv,2)
    subplot(round(sqrt(size(winv,2))),ceil(sqrt(size(winv,2))),wv)
    for e = 1:size(winv,1)
    ph=plot(e,winv(e,wv),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    ph = text(e,winv(e,wv),emo2{e});
    set(ph,'color',cols(e,:)); 
    end;
end;axcopy
textsc(['All comps; all subjs; PCA to ',int2str(size(winv,2))],'title');
%%%%%%%%%%%
%%  Plot 3 Dims vs each other:
cols = jet(15);
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure; % just 3  dims vs each other
c1 = 5; c2 = 2; c3 = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Run a decomposition of winvs with 15 indicators for which emotion it is. 
load /data/common4/emotion/EmoClustALLGood.mat pcamatall button allpcs
load /data/common4/emotion/PCASpacegdsubjs.mat ws winv activations button allpcs

emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
str = ['load /data/common4/emotion/GoodComps.mat gdcomps numsets gdchan paths fullpaths']; eval(str);
alldat = zeros(15,0);emoind = zeros(15,0); clear mnemodiff forstats
for nxs = 1:length(incsubjs)
    nx = incsubjs(nxs);
    str = ['load ',fullpaths{nx},'SpecCoModStuff.mat numtrials numframes freqs keeptrack rmepochs dstrials pcs savedat comment'];eval(str);  
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
    ws = wts*sph;    activations = ws*data;    winv = pinv(ws); clear wts sph ws allfacs alltemps winv2
    
    % find RMS of each factor activation and multiply with winv
    %for f = 1:size(winv,2)
    %    rms = sqrt(mean(activations(f,:).^2)); 
    %    winv(:,f) = winv(:,f)*rms;
    %end;
    tpind = ones(15,size(winv,1)); tpind = tpind*-1;
    for e = 1:length(emos)  % start with 2 for2 straight nums (not diffs)
        tpind(e,sum(dstrials(1:e-1))+1:sum(dstrials(1:e))) = 1; %tempmat = sort(tempmat);
        %tpind(e,sum(dstrials(1:e-1))+1:sum(dstrials(1:e))) = 300/dstrials(e); %tempmat = sort(tempmat);
        tpind(e,:) = tpind(e,:) - mean(tpind(e,:));
    end;     
    %emoind(:,end+1:end+size(winv,1)) = tpind;%instead of diff
    alldat = zeros(15,0);
    alldat(:,end+1:end+size(winv,1)) = winv'; 
    alldat(end+1:end+size(tpind,1),:) = tpind;
    fprintf('\n One More SUBJECT Done: %i',nx);clear sph  wts ws activations icamatall
    [weights,sphere,compvars,bias,signs,lrates,activations] = runica(alldat,'stop',1e-7,'maxsteps',1000,'pca',15);
    winv2 = pinv(weights*sphere);
    figure; cols=jet(15);
    for tmp = 1:size(winv2,2)
        sbplot(round(sqrt(size(winv2,2))),ceil(sqrt(size(winv2,2))),tmp)
        ph = plot([1:15],winv2(16:end,tmp),'k.-','linewidth',2,'markersize',10);
        xl = get(gca,'ylim'); set(gca,'xlim',[0 16]);hold on;
        set(gca,'ylim',xl);
        for em = 1:length(emos)
            ph = text(em,xl(1),emos{em});
            set(ph,'rotation',90);set(ph,'color',cols(em,:));
        end; 
        plot([get(gca,'xlim')],[0 0],'k:');
    end;
figure; imagesc(winv2)
 figure; imagesc(alldat)

[pc,eigvec,sv] = runpca(totdat,15);
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(pc,'extended',1,'stop',1e-7,'maxsteps',2000);
figure; imagesc(winv*eigvec);
save /data/common4/emotion/contextdecompstuff.mat weights sphere activations pc eigvec totdat
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find spectrum of trial wts for all clustered factors
% 2 Hz sampling rate for 'trials' (50% overlap of 1 sec epochs)
load /data/common4/emotion/KmeansClustCoMods.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers
intclusts = [2,4,5,6,8,9,11,13,15,16,19];
figure; row = 3; col =4; pl = 0;
for clst =1:length(intclusts)
    clust = intclusts(clst);
    if ~isempty(clustspecs{clust})
    pl = pl+1;
    subplot(row,col,pl)
    %plot(freqs(cfr),clustspecs{clust}(:,cfr));hold on;
    %ph = plot(freqs(cfr),mean(clustspecs{clust}(:,cfr),1),'k-');set(ph,'linewidth',1.5);
    plot(freqs,clustspecs{clust});hold on;
    ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',1.5);
    set(gca,'ylim',[-5 15]); set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
    title(['Cluster ',int2str(clst)]);
    end;
end;
axcopy
textsc(['Emotion Spectral Co-Modulation Factor Clusters; Kmeans to ',int2str(size(allsums,2))],'title');
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
% dstrials give LAST trial for each emotion.

for clust = 1:length(clustspecs)
    for nx = 1:length(gdcomps)
        if ~isempty(facvec{clust}{nx})
            str = ['load ',fullpaths{nx},'SpecCoModStuff.mat numtrials numframes freqs keeptrack rmepochs dstrials pcs savedat comment'];eval(str);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
            ws = wts*sph;    winv = pinv(ws); clear wts sph ws allfacs alltemps dspwr 
            for fc = 1:length(facvec{clust}{nx})
                %facwts = winv(:,facvec{clust}{nx}(fc)); % for all trials
                for ds = 1:length(dstrials) 
                    if dstrials(ds) > 70
                        %[pwr(:,fc),frq] = pwelch(facwts,32,16,128,2) ; 
                        %pwr(:,fc) = 10*log10(pwr(:,fc));
                        facwts = winv(sum(dstrials(1:ds-1))+1:sum(dstrials(1:ds)),facvec{clust}{nx}(fc));%by emo
                        [pwr,frq] = pwelch(facwts,32,16,128,2) ;%by emo
                        pwr = 10*log10(pwr);
                        dspwr(fc,:,ds) = pwr';%by emo
                    end;
                end;                
            end;
            %facmod{clust}{nx} = pwr'; % for all trials
            facmod{clust}{nx} = dspwr;%by emo
        end;
        fprintf('.');
    end;
end;
%save /data/common4/emotion/WtsSpectraAllEmos.mat facmod frq
%save /data/common4/emotion/WtsSpectraByEmo.mat facmod frq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/WtsSpectraAllEmos.mat facmod frq meanall
% plot wts spectra for all trials together
intclusts = [2,4,5,8,9,11,13,15,18];
figure; row = round(sqrt(length(intclusts))); col = ceil(sqrt(length(intclusts))); fr = find(frq > .05);
for clst = 1:length(intclusts)
    clust = intclusts(clst);        
    if ~isempty(facmod{clust})
    colspecs = [];
    sbplot(row,col,clst)
    for nx = 1:length(facmod{clust})
        if ~isempty(facmod{clust}{nx})
            for fc = 1:size(facmod{clust}{nx},1)
            ph= plot(frq(fr),facmod{clust}{nx}(fc,fr)-meanall(1,fr),'r');hold on;
            colspecs(end+1,:) = facmod{clust}{nx}(fc,fr)-meanall(1,fr);
            end;
        end;
    end;
    ph = plot(frq(fr),mean(colspecs,1),'k-'); set(ph,'linewidth',2);
                %set(gca,'ylim',[0 .6]);
                set(gca,'xlim',[frq(fr(1)) frq(fr(end))]);
    title(['Cluster ',int2str(clust)]);                    
    end;
end;axcopy
    ph=textsc('Spectral decomposition of wts from Co-Mod decomp; in dB; Mean of all spectra subtracted from all','title');set(ph,'fontsize',14);axcopy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/WtsSpectraByEmo.mat facmod frq meanall
% plot wts spectra for all emos separately (subplotted by subj)
intclusts = [2,4,5,8,9,11,13,15,18];
cols = jet(15); fr = find(frq > .05);
for clst = 1:length(facmod)
    clust = intclusts(clst);    leng = 0;
    for m = 1:length(facmod{clust})
        if ~ isempty(facmod{clust}{m})
            len = size(facmod{clust}{m},1);
            leng = leng+len;
        end;
    end;    
    figure; row = round(sqrt(leng)); col = ceil(sqrt(leng)); pl=1;
    for nx = 1:length(facmod{clust})
        if ~isempty(facmod{clust}{nx})
            for fc = 1:size(facmod{clust}{nx},1)
                sbplot(row,col,pl);  colspecs = [];
                for em = 1:size(facmod{clust}{nx},3)
                    %ph= plot(frq(fr),facmod{clust}{nx}(fc,fr,em)-meanall(1,fr),'r');hold on;
                    ph= plot(frq(fr),facmod{clust}{nx}(fc,fr,em),'r');hold on;
                    set(ph,'color',cols(em,:));
                    %colspecs(end+1,:) = facmod{clust}{nx}(fc,fr,em)-meanall(1,fr);
                    colspecs(end+1,:) = facmod{clust}{nx}(fc,fr,em);
                end;
                ph = plot(frq(fr),mean(colspecs,1),'k-'); set(ph,'linewidth',1.5);pl = pl+1;
                %set(gca,'ylim',[0 .6]);
                set(gca,'xlim',[frq(fr(1)) frq(fr(end))]);
                title(['Subj ',int2str(nx),'-Fac ',int2str(facvec{clust}{nx}(fc))]);                    
            end;
        end;
    end;
    ph=textsc(['Spectral Template Cluster ',int2str(clust)],'title');set(ph,'fontsize',14);axcopy
end;
% plot wts spectra for all emos separately(subplotted by emo)
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
intclusts = [2,4,5,8,9,11,13,15,18]; colspecs = cell(1,15);
cols = jet(15);cols(10,:) = [.9 .9 0]; fr = find(frq > .05);
for clst = 1:length(intclusts)
    clust = intclusts(clst);   
    %figure;row = 4; col = 4; pl=1;
    for nx = 1:length(facmod{clust})
        if ~isempty(facmod{clust}{nx})
            for fc = 1:size(facmod{clust}{nx},1)
                for em = 1:size(facmod{clust}{nx},3)
                    %sbplot(row,col,em);  
                    %ph= plot(frq(fr),facmod{clust}{nx}(fc,fr,em),'r');hold on;
                    %set(ph,'color',cols(em,:));
                    colspecs{em} = [colspecs{em}; facmod{clust}{nx}(fc,fr,em)-meanall(1,fr)];
                    %colspecs{em} = [colspecs{em}; facmod{clust}{nx}(fc,fr,em)];
                end;
            end;
        end;
    end;
    for emo = 1:15
        sbplot(row,col,emo);  
        ph= plot(frq(fr),colspecs{emo},'r');hold on; set(ph,'color',cols(emo,:));
        ph = plot(frq(fr),mean(colspecs{emo},1),'k-'); set(ph,'linewidth',2);pl = pl+1;
        set(gca,'ylim',[-10 10]);
        set(gca,'xlim',[frq(fr(1)) frq(fr(end))]);
        title(emos{emo});     
    end;    
    ph=textsc(['Spectral Template Cluster ',int2str(clust)],'title');set(ph,'fontsize',14);axcopy
    fprintf('.'); clf
    clustemospecs{clust} = colspecs;
end;

figure; row = round(sqrt(length(intclusts)));col = ceil(sqrt(length(intclusts))); pl = 1;
for clust = 1:length(clustemospecs)
    if ~isempty(clustemospecs{clust})
        sbplot(row,col,pl);  
        for emo = 1:15
            ph= plot(frq(fr),mean(clustemospecs{clust}{emo},1),'r');hold on; set(ph,'color',cols(emo,:));
            %ph = plot(frq(fr),mean(clustemospecs{clust}{emo},1),'k-'); set(ph,'linewidth',2);pl = pl+1;
            %set(gca,'ylim',[-10 10]);
            set(gca,'xlim',[frq(fr(1)) frq(fr(end))]);
            
        end; 
        title(['Clust ',int2str(clust)]); pl = pl+1;
    end;
end;
    ph=textsc(['Mean Wts Spectra from each emotion (color-coded) for each spectral cluster '],'title');set(ph,'fontsize',14);axcopy

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/EmoValence.mat emoval emos comment  


% for each freq, correlate behavioral rating with mean across factors
%%  NO CORRELATION IN ANY CLUSTER AT ANY FREQ IS SIG BY PERMUTATION TEST
figure; row = 4; col = 5; cols = jet(25);pl = 1;
emocorr = emoval/std(emoval);
for clust = 1:length(facmod)
    if ~isempty(facmod{clust})
        colspecs = [];
        for nx = 1:length(facmod{clust})
            if ~isempty(facmod{clust}{nx})
                for fc = 1:size(facmod{clust}{nx},1)
                    for em = 1:size(facmod{clust}{nx},3)
                        colspecs(end+1,:,em) = facmod{clust}{nx}(fc,:,em);
                    end;
                end;
            end;
        end;
        for fr = 1:length(frq)
            forcorr = squeeze(mean(colspecs(:,fr,:),1))';
            forcorr = forcorr/std(forcorr);            
            [corr(1,fr),indx,indy,corrs] = matcorr(emocorr,forcorr);
            for shf = 1:100
            [fakecorr(shf,fr),indx,indy,corrs] = matcorr(emocorr,shuffle(forcorr));
            end;
        end;
        mask(1,2) = max(max(fakecorr));
        mask(1,1) = min(min(fakecorr));
        subplot(row,col,pl)
        ph=plot(frq,corr,'k-'); set(ph,'color',cols(clust,:));pl = pl+1;hold on;
        set(ph,'linewidth',2);
        ph=title(['Cluster ',int2str(clust)]); set(ph,'fontsize',14);
        %set(gca,'ylim',[0 1]);
        plot([get(gca,'xlim')],[mask(1,1) mask(1,1)],'r:');
        plot([get(gca,'xlim')],[mask(1,2) mask(1,2)],'r:');
    end;
end;
ph=textsc(['Correlation of behavioral emotion scores with MEAN factor wts spectra at each freq; dashed lines are outer limits of permutation test (100 times) '],'title');set(ph,'fontsize',14);axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot all factors in a cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intclusts = [2,4,5,8,9,11,13,15,18]; % interesting factors
eeglab 
% change to include a super-impose option
emo = 11;
PlotSpecCoModAcrSubj(fullpaths,savedat,gdcomps,subjfacs{emo},freqs,12,['Spectral Decomp cluster of emotion data; EMOTION ',emos{emo}]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcnum = 16;
%[pc,eigvec,sv] = runpca(clustfacs',pcnum);
%x=alldip/mean(std(alldip));  

%icadat = [pc*10;x'];
%icadat = pc;

[weights,sphere,compvars,bias,signs,lrates,activations] = runica(clustfacs','stop',1e-7,'maxsteps',5000,'pca',16);
ws = weights*sphere;  winv = pinv(ws);
spectemps = winv;

%spectemps = eigvec*winv(1:pcnum,:);
figure; row = round(sqrt(pcnum)); col = ceil(sqrt(pcnum))+1;
for clust = 1:size(spectemps,2)
    subplot(row,col,clust)
    plot(freqs,spectemps(:,clust),'k-');
    set(gca,'xlim',[freqs(1) freqs(end)]);set(gca,'xgrid','on');
end;axcopy

% which factors are highly weighted?
clear facvechi facveclo th tl
for clust = 1:size(activations,1)
    hifacs = find(activations(clust,:)>mean(activations(clust,:))+1.5*std(activations(clust,:)));
    lofacs = find(activations(clust,:)<mean(activations(clust,:))-1.5*std(activations(clust,:)));
    th{clust} = clustfacs(hifacs,:);
    tl{clust} = clustfacs(lofacs,:); 
    for nx = 1:length(gdcomps)
        facvechi{clust}{nx} = kptk(hifacs(kptk(hifacs,1) == nx),2)';
        facveclo{clust}{nx} = kptk(lofacs(kptk(lofacs,1) == nx),2)';
    end;
end;
% easy plot
figure; clust = 8;  pl = 1;
row = round(sqrt(size(th{clust},1)));
col = round(sqrt(size(th{clust},1)))+1;
for rr = 1:size(th{clust},1)
    sbplot(row,col,pl)
    ph =plot(freqs,th{clust}(rr,:),'k-');pl = pl+1;hold on;
    set(gca,'fontsize',7);set(gca,'box','off');
    set(gca,'xlim',[freqs(1) freqs(end)]); set(gca,'ylim',[-2 10]);   
    plot([10 10],[get(gca,'ylim')],'r-');
    set(gca,'xtick',[5:5:freqs(end)]);  set(gca,'xticklabel',{[] 10 [] [] [] 30 [] [] [] []}); 
    set(gca,'xgrid','on');

end;
figure; clust = 8;  pl = 1;
row = round(sqrt(size(tl{clust},1)));
col = round(sqrt(size(tl{clust},1)))+1;
for rr = 1:size(tl{clust},1)
    sbplot(row,col,pl)
    ph =plot(freqs,tl{clust}(rr,:),'k-');pl = pl+1;hold on;
    set(gca,'fontsize',7);set(gca,'box','off');
    set(gca,'xlim',[freqs(1) freqs(end)]); set(gca,'ylim',[-2 10]);   
    plot([10 10],[get(gca,'ylim')],'r-');
    set(gca,'xtick',[5:5:freqs(end)]);  set(gca,'xticklabel',{[] 10 [] [] [] 30 [] [] [] []}); 
    set(gca,'xgrid','on');
end;


    
clust = 1;
PlotSpecCoModAcrSubj(fullpaths,gdcomps,facvec{clust},freqs,12,['Spectral Decomp cluster of emotion data; Cluster ',int2str(clust)]);


%%%%%% 
% check that all dipoles are under 15%
clear allrvs rvs
for nx = 1:length(gdcomps)    
    EEG = pop_loadset('sources.set', fullpaths{nx});    
    rvs = {EEG.dipfit.model(gdcomps{nx}).rv};
    rvs = cell2mat(rvs);
    if find(rvs > .15)
        nocomps{nx} = gdcomps{nx}(find(rvs > .15));
        newcomps{nx} = gdcomps{nx};
        newcomps{nx}(find(rvs > .15)) = [];
    end;
    allrvs{nx} = rvs;
    ALLEEG=[];EEG=[];
end;

nocomps =
 
  Columns 1 through 8
 
    [34]    [17]     []    [10,14]     []    [23]     []     []
 
  Columns 9 through 14
 
    [1x3 double]    [1x2 double]    [1x4 double]     []    [15]     []
 
  Columns 15 through 22
 
    [1x2 double]    [25]    [1x3 double]     []    [23]    [19]    [14]    [29]
 
  Columns 23 through 30
 
    [1x3 double]     []     []    [19]    [1x6 double]     []    [34]     []
 
  Columns 31 through 35
 
     []     []    [23]     []    [1x3 double]
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tried a sliding correlation matrix for clustering... not as good as kmeans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try doing a sliding correlation analysis to cluster comps
load /data/common4/emotion/clustfacs.mat clustfacs alldip kptk freqs
clear corr
for cp = 1:size(clustfacs,1)-1
    for ocp = cp+1:size(clustfacs,1)
        %for slide = 9:11
    [corr(cp,ocp),indx,indy,corrs] = matcorr(clustfacs(cp,14:28),clustfacs(ocp,14:28));   
        %end;
    end;
    fprintf('.');
end;
figure; imagesc(corr(:,:,3)); colorbar;
%find high corr
clear hicorr
for cp = 1:size(clustfacs,1)-1
    for ocp = cp+1:size(clustfacs,1)
        %if ~isempty(find(corr(cp,ocp,:) > .8))
            hicorr(cp,ocp) = max(abs(corr(cp,ocp,1))); 
        %end;
    end;
    fprintf('.');
end;
%figure; imagesc(hicorr); colorbar;

clsted = ones(1,size(hicorr,2));p=1;clear corrclust
%for corcut = .97:-.01:.92
corcut = .99;
for cp = 1:size(hicorr,1)
    if clsted(cp) > 0
        ocomps = [];oocomps = [];
        if ~isempty(find(hicorr(cp,:) > corcut))|~isempty(find(hicorr(:,cp) > corcut))
            oths = find(hicorr(cp,:) > corcut);
            for xp = 1:length(oths)
                ocomps = [ocomps find(hicorr(:,oths(xp))> corcut)'];
                if ~isempty(find(hicorr(oths(xp),:) > corcut))|~isempty(find(hicorr(oths(xp),:) > corcut))
                    ooths = find(hicorr(oths(xp),:) > corcut);
                    for xxp = 1:length(ooths)
                        oocomps = [oocomps find(hicorr(:,ooths(xxp))> corcut)'];
                    end;
                end;
            end;
            rwhicorr = find(abs(hicorr(cp,:)) > corcut);
            colhicorr = find(abs(hicorr(:,cp)) > corcut)';
            %rwhicorr(find(clsted))= [];
            %colhicorr(find(clsted))= [];
            %oocomps(find(clsted))= [];
            %ocomps(find(clsted))= [];
            corrclust{p} = unique([cp,rwhicorr,colhicorr,unique(ocomps),unique(oocomps)]); 
            p = p+1;
            %hicorr(cp,:) = 0;
            %hicorr(unique(ocomps),:) = 0;
            %hicorr(:,unique(ocomps)) = 0;
            %hicorr(unique(oocomps),:) = 0;
            %hicorr(:,unique(oocomps)) = 0;
            clsted(cp) = 0;
            clsted(unique(ocomps)) = 0;
            clsted(unique(oocomps)) = 0;
            clsted(rwhicorr) = 0;
            clsted(colhicorr) = 0;
        end;        
    end;
end;
length( find(clsted))
length(corrclust)

check=0;corcut = .98;thresh = .99;
newnotclsted = length( find(clsted));
while thresh > .8
for cp = 1:size(hicorr,1)
    if  clsted(cp) > 0
       newclust = 0;
        [val mostcorr] = max(hicorr(cp,:));
        if val > thresh
             check = check+1;
             for x = 1:length(corrclust)
                if ismember(mostcorr,corrclust{x})
                    corrclust{x} = [corrclust{x},cp];
                    %hicorr(cp,:) = 0;hicorr(:,cp) = 0;
                    clsted(cp) = 0; newclust= 0;
                end;
            end;
            if newclust == 1
                ocomps = [];oocomps = [];
                if ~isempty(find(hicorr(cp,:) > corcut))|~isempty(find(hicorr(:,cp) > corcut))
                    oths = find(hicorr(cp,:) > corcut);
                    for xp = 1:length(oths)
                        ocomps = [ocomps find(hicorr(:,oths(xp))> corcut)'];
                        if ~isempty(find(hicorr(oths(xp),:) > corcut))|~isempty(find(hicorr(oths(xp),:) > corcut))
                            ooths = find(hicorr(oths(xp),:) > corcut);
                            for xxp = 1:length(ooths)
                                oocomps = [oocomps find(hicorr(:,ooths(xxp))> corcut)'];
                            end;
                        end;
                    end;
                    rwhicorr = find(hicorr(cp,:) > corcut);
                    colhicorr = find(hicorr(:,cp) > corcut)';
                    corrclust{p} = unique([cp,rwhicorr,colhicorr,unique(ocomps),unique(oocomps)]); 
                    p = p+1;
                    %hicorr(cp,:) = 0;
                    %hicorr(unique(ocomps),:) = 0;
                    %hicorr(:,unique(ocomps)) = 0;
                    %hicorr(unique(oocomps),:) = 0;
                    %hicorr(:,unique(oocomps)) = 0;
                    clsted(unique(ocomps)) = 0;
                    clsted(unique(oocomps)) = 0;
                    clsted(find(hicorr(cp,:) > corcut)) = 0;
                    clsted(find(hicorr(:,cp) > corcut)) = 0; 
                end;   
            end;
        end;        
    end;
end;
notclsted = length( find(clsted));
if notclsted == newnotclsted
    thresh = thresh - .05;
end;
newnotclsted = notclsted
end;
length( find(clsted))
length(corrclust)


clear clustspecs
for cl = 1:length(corrclust)
    clustspecs{cl} = clustfacs(corrclust{cl},:);
end;
figure; row = 6; col =7; pl = 0;
for clust =1:length(clustspecs)
    pl = pl+1;
    subplot(row,col,pl)
    %plot(freqs(cfr),clustspecs{clust}(:,cfr));hold on;
    %ph = plot(freqs(cfr),mean(clustspecs{clust}(:,cfr),1),'k-');set(ph,'linewidth',1.5);
    plot(freqs,clustspecs{clust});hold on;
    ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',1.5);
    set(gca,'ylim',[-5 15]); set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
end;
subplot(row,col,pl+1)
plot(freqs,clustfacs(find(clsted),:));hold on;
ph = plot(freqs,mean(clustfacs(find(clsted),:),1),'k-');set(ph,'linewidth',1.5);
set(gca,'ylim',[-5 15]); set(gca,'xgrid','on');
set(gca,'xlim',[freqs(1) freqs(end)]);
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Try clustering on back-proj envelope%%%%%%%%%%%%%%%%%%%%%%%
pcnum = 5;
[weights,sphere,compvars,bias,signs,lrates,activations]  = runica(clustfacs,'extended',1,'pca',pcnum,'stop',1e-7,'maxsteps',1000);
winv = pinv(weights*sphere);

 reps=5;optk = 16;
[kout, C,sumd, allsums] = kmeans(winv,optk,'replicates',reps); 
clear clustidxs clustspecs facvec outliers
kkout = kout;
for cl = 1:optk
    oneclust = find(kout == cl);
    fout = allsums(oneclust)/std(allsums(oneclust));
    outliers{cl} = find(fout > 9)';
    kkout(oneclust(outliers{cl})) = 0;
    oneclust = find(kkout == cl);
    %clustspecs{cl} = clustenv(oneclust,:);
    %clustspecs{cl} = normfacs(oneclust,:);
    clustspecs{cl} = clustfacs(oneclust,:);
    for nx = 1:length(gdcomps)
        facvec{cl}{nx} = kptk(oneclust(kptk(oneclust,1) == nx),2)';
    end;
end;
outliers
figure; row = 4; col =5; pl = 0;
for clust =1:length(clustspecs)
    if ~isempty(clustspecs{clust})
    pl = pl+1;
    subplot(row,col,pl)
    %plot(freqs(cfr),clustspecs{clust}(:,cfr));hold on;
    %ph = plot(freqs(cfr),mean(clustspecs{clust}(:,cfr),1),'k-');set(ph,'linewidth',1.5);

    plot(freqs,clustspecs{clust}(:,1:99),'r-');hold on;
    %plot(freqs,clustspecs{clust}(:,100:198),'b-');hold on;
    ph = plot(freqs,mean(clustspecs{clust}(:,1:99),1),'k-');set(ph,'linewidth',2.5);
    %ph = plot(freqs,mean(clustspecs{clust}(:,100:198),1),'k-');set(ph,'linewidth',2.5);
    %set(gca,'ylim',[-5 10]); 
    set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
    title(['Cluster ',int2str(clust)]);
    end;
end;
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 


save /data/common4/emotion/KmeansClustCoMods3.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers % 3 is the clustering of the envelopes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cluster factor spectra by kmeans   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; row = 21; col = 21;pl = 0;
for fc = 1:size(clustfacs,1)
    pl = pl+1;
    sbplot(row,col,pl)
    plot(freqs,clustfacs(fc,:));hold on;
    set(gca,'ylim',[-5 15]); set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
    plot([10 10],[get(gca,'ylim')],'r-');hold on;
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    fprintf('.');
end;



 cfr = find(freqs > 4 & freqs < 40);
for x=1:size(clustfacs,1)
    %normfacs(x,:) = clustfacs(x,:)/std(clustfacs(x,cfr));
    normfacs(x,:) = clustfacs(x,:)/sqrt(mean(clustfacs(x,cfr).^2));    
end;

pcnum = 25; cfr = find(freqs > 3 & freqs < 50);
[pc,eigvec,sv] = runpca(normfacs(:,:)',pcnum);

 %[optk,centr,clst,Cg] = Kmeangrp(pc',35,3,1);

 reps=5;optk = 20;
[kout, C,sumd, allsums] = kmeans(pc',optk,'replicates',reps); 
clear clustidxs clustspecs facvec outliers
kkout = kout;
for cl = 1:optk
    oneclust = find(kout == cl);
    fout = allsums(oneclust)/std(allsums(oneclust));
    outliers{cl} = find(fout > 5)';
    kkout(oneclust(outliers{cl})) = 0;
    oneclust = find(kkout == cl);
    clustspecs{cl} = clustfacs(oneclust,:);
    %clustspecs{cl} = normfacs(oneclust,:);
    for nx = 1:length(gdcomps)
        facvec{cl}{nx} = kptk(oneclust(kptk(oneclust,1) == nx),2)';
    end;
end;
outliers
figure; row = 4; col =5; pl = 0;
for clust =1:length(clustspecs)
    if ~isempty(clustspecs{clust})
    pl = pl+1;
    subplot(row,col,pl)
    %plot(freqs(cfr),clustspecs{clust}(:,cfr));hold on;
    %ph = plot(freqs(cfr),mean(clustspecs{clust}(:,cfr),1),'k-');set(ph,'linewidth',1.5);
    plot(freqs,clustspecs{clust});hold on;
    ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',1.5);
    set(gca,'ylim',[-5 10]); set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);
    title(['Cluster ',int2str(clust)]);
    end;
end;
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
%%  DONE!!!  :
% save /data/common4/emotion/KmeansClustCoMods.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers
% save /data/common4/emotion/KmeansClustCoMods2.mat clustfacs clustspecs kptk facvec kout allsums freqs outliers
