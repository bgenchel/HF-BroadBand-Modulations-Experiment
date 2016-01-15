%%  For use in clustering spectral factor templates across subjects

eeglab
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
str = ['load /data/common4/emotion/GoodComps.mat gdcomps numsets gdchan paths fullpaths']; eval(str);
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; % for all new ones
frqlim = [0 50];overlap = 2; % 50% overlap(1 sec epochs)
savedat = 'SpecCoMod';pc=15;
nx=2;str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cluster spectral co-mod templates across subjects
savedat = 'SpecCoMod';
clustfacs = []; clustenv = []; mnspecs = []; kptk = [];  pl = 1; 
incsubjs = [2:21,23:35]; 
for nxx = 1:length(gdcomps)
    clear data activations winv  
    nx = incsubjs(nxx);
    str = ['load ',fullpaths{nx},savedat,'Stuff.mat '];eval(str);  
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.sph'],[numtrials numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(pcs),'.wts'],[pcs numtrials],[],0); 
    data = floatread([fullpaths{nx},savedat,'.fdt'],[numtrials numframes],[],0);    
    ws = wts*sph;    activations = ws*data;    clear wts sph ws  
    for tp = 1:size(activations,1)
        clear allrms alltemps
        for rcp = 1:length(gdcomps{nx})
            alltemps(rcp,:) = activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp); 
            allrms(rcp) = sqrt(mean(activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp).^2));
        end;
        clustfacs = [clustfacs;alltemps];
        kptk = [kptk;[repmat(nx,[size(alltemps,1) 1]) repmat(tp,[size(alltemps,1) 1]) gdcomps{nx}']];
        mnspecs = [mnspecs;meanpwr]; 
        [mxv mxidx] = max(allrms);
        mxspec = alltemps(mxidx,:);
        mxval = mxspec * mxspec';
        for cp = 1:size(alltemps,1)
            scaleby(1,cp) = mxspec * alltemps(cp,:)';
        end;
        
    end;
    fprintf('.');
end;
load /data/common4/emotion/ClustAllFacs.mat clustfacs mnspecs kptk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Try clustering factor spectra with pdist %%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/ClustAllFacs.mat 
%%%% the following seems to be an optimal combination

pl = 1; clear newkeep normfacs
for f = 1:size(clustfacs,1)
    %normfacs(pl,:) = clustfacs(f,:)/sqrt(mean(clustfacs(f,:).^2)) ;
    normfacs(pl,:) = clustfacs(f,:);
    newkeep(pl,:) = kptk(f,:);pl = pl+1;
    %normfacs(pl,:) = (clustfacs(f,:)/sqrt(mean(clustfacs(f,:).^2)))*-1;
    normfacs(pl,:) = clustfacs(f,:)*-1;
    newkeep(pl,:) = kptk(f,:);pl = pl+1; % double so both orientations represented
end;

pcnum =10; % reducing dimensions good. 10 maybe?
[weights,sphere,compvars,bias,signs,lrates,activations]  = runica(normfacs,'extended',1,'pca',pcnum,'stop',1e-7,'maxsteps',1000);
winv = pinv(weights*sphere);

alldist = pdist(winv, 'euclidean'); % euc better than seuc
links = linkage(alldist,'ward');%ward best
figure;[hnd,idx,perm]=  dendrogram(links,50);
figure; clear ctxclusts facclusts facvec
row=round(sqrt(max(idx))); 
col=ceil(sqrt(max(idx))); 
for cls = 1:max(idx)
    onecls = find(idx == cls);
    ctxclusts{cls} = normfacs(onecls,:);
    facclusts{cls} = newkeep(onecls,:);  
    clustspecs{cls} = normfacs(onecls,:);
    rr = onecls(newkeep(onecls,1) == nx);% find indices for current subj
    facvec{cls} = newkeep(rr,2)';
    sbplot(row,col,cls)
    ph = plot(freqs,clustspecs{cls}');hold on;
    set(ph,'color','g');
    ph = plot(freqs,mean(clustspecs{cls},1),'linewidth',2);
    set(ph,'color','b');
    set(gca,'xlim',[3 freqs(end)]);
    title([int2str(cls),'-',int2str(length(facclusts{cls}))]);
end;ctxclusts
axcopy

save /data/common4/emotion/PdistallFacs.mat clustfacs clustspecs kptk facvec freqs idx newkeep

load /data/common4/emotion/ClustAllFacs.mat 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  using the Perc pwr from each factor/component, weight a density plot for each cluster

intclusts = [1:2:50];% 
nx=2;
EEG = pop_loadset('sources.set' ,fullpaths{nx}); 
for clst = 1:length(intclusts)
    clust = intclusts(clst);  
    new = 1;       denswt = zeros(1,0);    subjidx = zeros(1,0);  
     numsubjs = 0; numfacs = 0; numcomps = 0;clear complist
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
        for w = 1:size(facclusts{clust},1)
            dipsources(1,w).posxyz = EEG.dipfit.model(facclusts{clust}(w,3)).posxyz;
            dipsources(1,w).momxyz = EEG.dipfit.model(facclusts{clust}(w,3)).momxyz;
            dipsources(1,w).rv = EEG.dipfit.model(facclusts{clust}(w,3)).rv;
            subjidx(1,end+1) = nx;
        end;           
        if new == 1
            allbesa = dipsources;new = 0;
        else
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
        end;
        dipsources = []; 
    end;
    optdipplot = {allbesa,'gui','off','image','mri','coordformat','spherical','normlen','on'};
    figure;dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',denswt); 
     
    ph =textsc(['Subject ',int2str(nx),'; Cluster ',int2str(clust)],'title');set(ph,'color','r');    
    set(ph,'fontsize',14);
    %str = ['print /data/common4/emotion/Figs/Sb',int2str(nx),'FacClust',int2str(clust),'.jpg -djpeg'];
    %eval(str); close; close;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intclusts = [1:30];
for clst = 1:length(intclusts)
    clust = intclusts(clst);  numfacs = 0; numcomps = 0;clear complist
            for fc = 1:length(facclusts{clust})
                fac = facclusts{clust}(fc,2);
                tmplist = zeros(1,0);
                for w = 1:length(gdcomps{nx})
                    if ismember(gdcomps{nx}(w),facclusts{clust}(:,3))
                        tmplist(1,end+1) = gdcomps{nx}(w);numcomps=numcomps+1;
                   end;                    
                end;           
                complist{nx}{fc} = tmplist;
                fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                numfacs = numfacs+1;
            end;
    PlotCrossLines(complist,fullpaths,'sources.set');
    subplot(2,2,4)
    %plot(freqs,clustspecs{clust}(:,1:99),'r-');hold on; % for envelope
    %plot(freqs,clustspecs{clust}(:,100:198),'b-');hold on;
    %ph = plot(freqs,mean(clustspecs{clust}(:,1:99),1),'k-');set(ph,'linewidth',2.5);
    %ph = plot(freqs,mean(clustspecs{clust}(:,100:198),1),'k-');set(ph,'linewidth',2.5);
    plot(freqs,clustspecs{clust});hold on;set(gca,'fontsize',16);
    ph = plot(freqs,mean(clustspecs{clust},1),'k-');set(ph,'linewidth',2.5);
    %set(gca,'ylim',[-3 8]); 
    set(gca,'xgrid','on');
    set(gca,'xlim',[freqs(1) freqs(end)]);xlabel('Frequency (Hz)');
    ylabel('Relative Power');
    title(['Co-Mod Cluster ',int2str(clust)]);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
     ph =textsc([int2str(numsubjs),'/33 Subjs; ',int2str(numfacs),' Factors; ',int2str(numcomps),' Components'],'title');
     set(ph,'fontsize',14);
    str = ['print /data/common4/emotion/Figs/CoModClsLines',int2str(clust),'.jpg -djpeg'];eval(str)
    close
end;
