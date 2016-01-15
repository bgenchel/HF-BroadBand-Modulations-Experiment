% Epoch datasets on button press in subjects who pressed more than once and look for pre/post changes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/home/julie/MatlabScripts/emotion')
DataInfo    % this matlab files loads all subject info needed
fullpaths=newpaths;

repbutton = [2:12];  ttl = 'Repetitive button presses';% really subj 1 too
feelbutton = [21,23:26];  ttl = 'only when feeling it';
nobutton = [13:20,27:31,33:35]; ttl = ' no button press';% (apart from the first one)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPOLE density diffs between button groups?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buttgrps = {repbutton,feelbutton,nobutton};
buttls = {'repetitive','feelingit','nobutton'};

% All ICs from each subject in each group
clear clustcps
for clust = 1:3
    button = buttgrps{clust};
    clustcps{clust} = cell(1,35);
    for nxx = 1:length(button)
        nx = button(nxx);
        clustcps{clust}{nx} = gdcomps{nx};
    end;
end
vw = {'top','side','rear'};
dipviews = {[0 0 1],[1 0 0],[0 -1 0]};
slices = {[64:-8:-24],[-55:10:55],[65:-14:-89]};
geo = {'geom',[3,4]};

[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);

for c1=1:2
    for c2 = c1+1:length(clustcps)
        [minmask maxmask] = SigDensityDiff('sources.set','sources.set',fullpaths,fullpaths,clustcps{c1},clustcps{c2},[],[],200,.01,'dips');    
        
        for v = 1:2 
            densargs = ['mrislices',slices{v},'mriview',vw{v}, geo];                    
            [densdiff] = PlotSigDensDiff('sources.set','sources.set',fullpaths,fullpaths,clustcps{c1},clustcps{c2},[],[],minmask,maxmask,densargs,'dips');
            
            ph = textsc([buttls{c1},'-',buttls{c2},' dipole density difference'],'title');set(ph,'color','r');
            eval(['print ',fullpaths{nx}(1:end-5),buttls{c1},'-',buttls{c2},'DensDiff','_',vw{v},'.jpg -djpeg']);
        end;
    end;
end;

% ICs with gamma IMs only:--------------
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(ffs,allbigs,bigwts,orivec);
buttls = {'repetitive','feelingit','nobutton'};

for c1=1:2
    for c2 = c1+1:length(denslist)
        [minmask maxmask] = SigDensityDiff('sources.set','sources.set',fullpaths,fullpaths,denslist{c1},denslist{c2},[],[],200,.01,'dips');    
        
        for v = 1:2 
            densargs = ['mrislices',slices{v},'mriview',vw{v}, geo];                    
            [densdiff] = PlotSigDensDiff('sources.set','sources.set',fullpaths,fullpaths,denslist{c1},denslist{c2},[],[],minmask,maxmask,densargs,'dips');
            
            ph = textsc([buttls{c1},'-',buttls{c2},' dipole density diff ALPHA'],'title');set(ph,'color','r');
            eval(['print ',fullpaths{nx}(1:end-5),buttls{c1},'-',buttls{c2},'AlphaIMDensDiff','_',vw{v},'.jpg -djpeg']);
        end;
    end;
end;
densdiff = abs(densdiff);
save /home/julie/tmpdens.mat densdiff

list_brain_areas('/home/julie/tmpdens.mat',0,['/home/julie/tmpRegions.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IM differences around button press?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
button = repbutton;
button = feelbutton;

for nxx = 1:length(button)
    nx = button(nxx);
    for e = 1:length(datset)
        EEG = pop_loadset( 'filename', datset{e}, 'filepath', fullpaths{nx});
        if length(find(ismember({EEG.event.type},'press'))) > 4
            figure; pl = 1; row = round(sqrt(length(gdcomps{nx}))); col = ceil(sqrt(length(gdcomps{nx})));
            figure; pl = 1; row = length(gdcomps{nx}); col = 1;
            
            %EEG = pop_epoch( EEG, {  'press'  }, [-6  6], 'newname', 'Button Epochs', 'epochinfo', 'yes');
            for ic = 1:length(gdcomps{nx})
                alphapwr = abs(hilbert(eegfilt(EEG.icaact(gdcomps{nx}(ic),:),EEG.srate,8,13))).^2;
                alphapwr = 10*log10(alphapwr);
                 [alphapwr,outx] = movav(alphapwr,[1:length(alphapwr)],EEG.srate);
                 evtms = cell2mat({EEG.event(find(ismember({EEG.event.type},'press'))).latency});
                sbplot(row,col,pl); pl = pl+1;
                plot([1:length(alphapwr)],alphapwr);hold on;
                set(gca,'xlim',[1 length(alphapwr)]); set(gca,'ylim',[min(alphapwr) max(alphapwr)]);
                yl = get(gca,'ylim'); 
                for ev = 1:length(evtms)
                    ph = plot([evtms(ev) evtms(ev)],[yl(1) yl(2)],'r-','linewidth',2);
                end;
                title(['IC ',int2str(gdcomps{nx}(ic))]);
% $$$                 sbplot(row,col,pl); pl = pl+1;
% $$$                 erpimage( EEG.icaact(gdcomps{nx}(ic), :), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), ['Act, IC ',int2str(gdcomps{nx}(ic))], 2, 1 , 'yerplabel','', 'erp', 'on','cbar', 'off', 'topo', { mean(EEG.icawinv(:,gdcomps{nx}(ic)),2) EEG.chanlocs EEG.chaninfo },'limits' , [-500 1000 NaN NaN NaN NaN NaN NaN NaN]);
% $$$                 sbplot(row,col,pl); pl = pl+1;
% $$$                 erpimage( EEG.icaact(gdcomps{nx}(ic), :), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), ['Pwr, IC ',int2str(gdcomps{nx}(ic))], 2, 1 , 'yerplabel','', 'erp', 'on', 'cbar', 'off', 'plotamps', 'on', 'coher',[8 12 0.01] ,'topo', { mean(EEG.icawinv(:,gdcomps{nx}(ic)),2) EEG.chanlocs EEG.chaninfo },'limits' , [-1000 1500 NaN NaN NaN NaN NaN NaN NaN] );
            end;
            textsc(['Subject ',int2str(nx),', ',emos{e}],'title');
        end;
    end;
end;

% Look at IMs instead:
savedat = 'SpecCoModMoreFreqs'; 
load /data/common1/emotion/AllClustFacs.mat    

for nxx = 1:length(button)
    nx = button(nxx);
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;    winv = pinv(ws); 
    clear wts sph ws 
    speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  % templates   
    winv = specwts;  clear delpoints emeans  
    allevs = [];
    for e = 1:length(datset)
        EEG = pop_loadset( 'filename', datset{e}, 'filepath', fullpaths{nx});
        epsize = .5; % one second
        for evtm = epsize/2*EEG.srate:round(EEG.srate/s.overlap):size(EEG.data,2)-epsize/2*EEG.srate
            % create events to make overlapping 1 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';% for string event codes  
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-epsize epsize]);% was-1 1 recently
        EEG = pop_select(EEG,'notrial',s.rmepochs{e});
        presstrial = [];
        for ep = 1:length(EEG.epoch)
            if ismember('press',EEG.epoch(ep).eventtype)
                presstrial = [presstrial ep];
            end;
        end;
        allevs = [allevs presstrial+(s.keeptrack(e,1)-1)];
    end;
    edat = epoch(winv',allevs,[-8 8]);
    tms = [-8:1:8];
    figure; pl = 1; row = 6; col = 6;    
    for ic = 1:length(gdcomps{nx})
        for im = 1:length(allbigs{nx})
            if ismember(gdcomps{nx}(ic),allbigs{nx}{im})
                if pl > row*col
                    textsc(['Subject ',int2str(nx)],'title');
                    figure; pl = 1;
                end;
                sbplot(row,col,pl); pl = pl+1;
                imagesc(tms,[1:size(edat,3)],squeeze(edat(im,:,:))');hold on;
                plot([0 0],[get(gca,'ylim')],'k-');
                title(['IC ',int2str(gdcomps{nx}(ic)),', IM ',int2str(im)]);
            end;
        end;
    end;
    textsc(['Subject ',int2str(nx)],'title');
end;

savedat = 'SpecCoModWavePrePost';
datset = {'anger.set','frustration.set','jealousy.set','fear.set' ,'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set','prebase.set','postbase.set'}; % for all new ones

button = [2:18,20,21,23:31,33:35];
figure; row=6; col=6;
clear higrats
for nxx = 1:length(button)
    nx = button(nxx);
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;    winv = pinv(ws); 
    clear wts sph ws 
    speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  % templates   
    winv = specwts;  clear delpoints emeans  
    emowts = median(winv(1:s.keeptrack(15,2),:));
    bslnwts = median(winv(s.keeptrack(16,1):end,:));
    allrats = []; clear icimidx icimgrats
    for ic = 1:length(gdcomps{nx})
        imgrat = [];imidx = [];
        for im = 1:length(allbigs{nx})
            if ismember(gdcomps{nx}(ic),allbigs{nx}{im})
                imgrat = [imgrat bslnwts(1,im)-emowts(1,im)];
                imidx = [imidx im];
                allrats = [allrats bslnwts(1,im)-emowts(1,im)];
            end;
        end;
        icimgrats{ic} = imgrat;
        icimidx{ic} = imidx;
    end;
    for ic = 1:length(gdcomps{nx})
        highs = find(abs(icimgrats{ic}) > median(allrats)+std(allrats)*2);
        imdir = single(icimgrats{ic}(highs)./abs(icimgrats{ic}(highs))); % orientation of IM
        higrats{nx}{ic} = icimidx{ic}(highs).*imdir;
    end;
    %sbplot(row,col,nx);    plot(allrats)
end;
clustcps{1} = cell(1,35);
addtmpl = [];addidx = []; addmeans = [];
for nx = 1:length(higrats)
    for ic = 1:length(higrats{nx})
        if ~isempty(higrats{nx}{ic})
            clustcps{1}{nx} = [clustcps{1}{nx} gdcomps{nx}(ic)];
            for im = 1:length(higrats{nx}{ic})
                currim = find(kptk(:,1) == nx & abs(kptk(:,2)) == abs(higrats{nx}{ic}(im)) & kptk(:,3) == gdcomps{nx}(ic));
                imdir = kptk(currim,2)/abs(kptk(currim,2))*(higrats{nx}{ic}(im)/abs(higrats{nx}{ic}(im))); % orientation of IM
                addtmpl = [addtmpl; clustfacs(currim,:)*imdir];
                addmeans = [addmeans; mnspecs(currim,:)];
                addidx = [addidx; abs(kptk(currim,:))];
            end;
        end;
    end;
end
figure; PlotDipoleClusters('sources.set', fullpaths,gdcomps,clustcps,1,2,2,1,'Hi diff bet rest & imagery',[1,2,3,4],[1 0 0],'off');
  

figure; quadplot(s.freqs,addtmpl);
figure; quadplot(s.freqs,mean(addtmpl,1),2);


[deltaclust, thetaclust,alphaclust,betaclust,gamaclust,freqs] = SortModTempls(addtmpl,addidx,addmeans,s.freqs,[3 128],[],'quad');

for x = 1:size(alphaclust{1},1)
    [val idx] = max(abs(alphaclust{1}(x,find(s.freqs>7&s.freqs<14))));
    if alphaclust{1}(x,idx+length(find(s.freqs<=7)))< 0 
        alphaclust{3}(x,2) = alphaclust{3}(x,2)*-1;
    end;
end;
for x = 1:size(betaclust{1},1)
    [val idx] = max(abs(betaclust{1}(x,find(s.freqs>13&s.freqs<35))));
    if betaclust{1}(x,idx+length(find(s.freqs<=13)))< 0 
        betaclust{3}(x,2) = betaclust{3}(x,2)*-1;
    end;
end;
for x = 1:size(gamaclust{1},1)
    [val idx] = max(abs(gamaclust{1}(x,find(s.freqs>63))));
    if gamaclust{1}(x,idx+length(find(s.freqs<=63)))< 0 
        gamaclust{3}(x,2) = gamaclust{3}(x,2)*-1;
    end;
end;
clear finalidx finaltmpls
finalidx{1} = alphaclust{3}(find(alphaclust{3}(:,2) > 0),:);
finaltmpls{1} = alphaclust{1}(find(alphaclust{3}(:,2) > 0),:);

finalidx{2} = alphaclust{3}(find(alphaclust{3}(:,2) < 0),:);
finaltmpls{2} = alphaclust{1}(find(alphaclust{3}(:,2) < 0),:);

finalidx{3} = betaclust{3}(find(betaclust{3}(:,2) > 0),:);
finaltmpls{3} = betaclust{1}(find(betaclust{3}(:,2) > 0),:);

finalidx{4} = betaclust{3}(find(betaclust{3}(:,2) < 0),:);
finaltmpls{4} = betaclust{1}(find(betaclust{3}(:,2) < 0),:);

finalidx{5} = gamaclust{3}(find(gamaclust{3}(:,2) > 0),:);
finaltmpls{5} = gamaclust{1}(find(gamaclust{3}(:,2) > 0),:);

finalidx{6} = gamaclust{3}(find(gamaclust{3}(:,2) < 0),:);
finaltmpls{6} = gamaclust{1}(find(gamaclust{3}(:,2) < 0),:);
            
[facvec,comods,wtsmat1,justcomps,jcwts,denslist] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);

row = length(comods);
viewnum=[1,2,3];col = length(viewnum)  + 1;
zoom= 1.6;
figure;pl = 1;
for clust = 1:length(comods)
    sbplot(row,col,pl); pl = pl+1;
    ph = quadplot(s.freqs,finaltmpls{clust});
    [angles] = PlotCoModasDipoles(comods{clust},justcomps{clust},newpaths,'sources.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,[]); % next to last 1 plots solo IMs in black
    pl = pl+length(viewnum);
end;
