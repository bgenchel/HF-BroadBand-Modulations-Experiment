% this script is to find peak freqs in spectra, normalize and plot on a besa head

% find power spectrum kb77 and all rest voltage corrected for nonWM
%eeglab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coefs = cwt(EEG.icaact(1,:),[2:.5:20],'morl','plot');
%%%%  DONE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear allspec
for nx = 34:length(paths)
    if ~isempty(gdcomps{nx})
        EEG = pop_loadset( 'prebase.set',['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);
    cd (['/data/common2/emotion/',paths{nx},'/ersps/']);
    sph=floatread(['/data/common2/emotion/',paths{nx},sphs{nx}],sphsize{nx}); 
    wts=floatread(['/data/common2/emotion/',paths{nx},wtss{nx}],wtssize{nx}); 
         EEG.icaweights=wts;EEG.icasphere=sph;EEG.icawinv=[];
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
       
        % special for continuous dataset%%%%%%%%%%%%%%%
    for evtm = 256:512:size(EEG.data,2)-256  % go up by 3 sec to create 3 sec epochs
    EEG.event(end+1) =  EEG.event(1);% appends events to the end
        EEG.event(end).latency = evtm;
        EEG.event(end).type = 'fake';        
    end;
    EEG = pop_epoch( EEG,{'fake'} , [-1 1]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
    EEG = pop_rmbase( EEG,[-1000 1000]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_rejkurt(EEG,0,gdcomps{nx} ,4,4,0,1);        
    EEG = pop_jointprob(EEG,0,gdcomps{nx} ,4,4,0,1);
    EEG = pop_rejkurt(EEG,0,gdcomps{nx},5,5,0,1);        
    EEG = pop_jointprob(EEG,0,gdcomps{nx} ,5,5,0,1);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tm = find(EEG.times > -700 & EEG.times < 700);  % input time to be analyzed!
        fixcompspec = zeros(length(gdcomps{nx}),257);
        for comp = 1:length(gdcomps{nx})%size(EEG.data,1)
            [pwr freqs] = pwelch(EEG.icaact(gdcomps{nx}(comp),:),256,128,512,EEG.srate);
            %[pwr freqs] = pwelch(reshape(EEG.icaact(gdcomps{nx}(comp),tm,:),1,length(tm)*size(EEG.icaact,3)),256,128,512,EEG.srate);
            fixcompspec (comp,:) = pwr';
        end;        
        allspec{nx} = fixcompspec;
        ALLEEG=[];EEG=[];
    end;
end;
comment = 'uses pwelch method to estimate power spectral density for EEG.icaact of all good comps for each subject. comps saved in compact matrix so row does not equal component, but gdcomps{nx}(row)==component; pwelch params = data,256,128,512,EEG.srate;prebase.set';
save /data/common2/emotion/clusters/allPreBasespecs.mat freqs allspec comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************************************************************************
load /data/common1/stern/eeg/Sternall/allchanspecs.mat freqs allspec comment
chanspec = allspec;
load /data/common1/stern/eeg/Sternall/fixcompspc.mat 
fixspec = allspec;
load /data/common1/stern/eeg/Sternall/alltaskspec.mat 
taskspec = allspec;
load /data/common1/stern/eeg/Sternall/mem0-1000spec.mat 
memspec = allspec;
load /data/common1/stern/eeg/Sternall/ig0-1000spec.mat 
igspec = allspec;
load /data/common1/stern/eeg/Sternall/probe0-1000spec.mat 
probespec = allspec;
conds = {chanspec, fixspec, taskspec,memspec, igspec, probespec};
fr = find(freqs > 1  & freqs < 15);
clear allsubjpeaks allconds
peaksurround = 3; %larger number gets only broad peaks, tall, sharp peaks can be missed if too large
for cond = 1:length(conds)
    for nx = 1:length(allspec)
        if ~isempty(gdcomps{nx})
            compspec = conds{cond}{nx};clear allmax
            if cond == 1
                X = mean(conds{cond}{nx}(:,fr),1);pl = 1;  mxma = [];
                for n = peaksurround+1:length(X)-peaksurround
                    if X(n) > X(n-[1:peaksurround])& X(n) > X(n+[1:peaksurround])%
                        mxma(1,pl) = n; pl = pl+1;                
                    end;
                end;
                allmax = freqs((fr(1)-1) + mxma)';
            else        
                for a = 1:length(gdcomps{nx})
                    X = compspec(a,fr);pl = 1;  mxma = [];
                    for n = peaksurround+1:length(X)-peaksurround
                        if X(n) > X(n-[1:peaksurround])& X(n) > X(n+[1:peaksurround])%
                            mxma(1,pl) = n; pl = pl+1;                
                        end;
                    end;
                    allmax{a} = freqs((fr(1)-1) + mxma)';
                end;
            end;
            allsubjpeaks{nx} = allmax;
        end;
    end;
    allconds{cond} = allsubjpeaks;
end;
%%%%%%%%%%%%%%
%  compare peaks between conditions
%  1=meanchan; 2=fix,3=taskbase, 4=mem,5=ig,6=probe;
frange = [8 13];
cond1 = 3;      cond2 = 4;
pl = 1; plotcomps = cell(1,length(gdcomps)); 
for nx = 1:length(allconds{cond1})
    speaks1 = allconds{cond1}{nx};clear allmax
    speaks2 = allconds{cond2}{nx};clear allmax
    for a = 1:length(gdcomps{nx})
        cpeaks1 = speaks1{a};
        cpeaks2 = speaks2{a};already = 0;
        for c = 1:length(cpeaks1)
            if cpeaks1(c) > frange(1) & cpeaks1(c) <frange(2)
                if ~isempty(find(cpeaks2 >frange(1) & cpeaks2<frange(2)))
                    idx = find(cpeaks2 >frange(1) & cpeaks2<frange(2));
                    if abs(cpeaks1(c)-cpeaks2(idx)) > freqs(6)-freqs(4)
                        if already < 1
                            plotcomps{nx}(end+1) = gdcomps{nx}(a);
                            already = 1;
                        end;
                    end;
                else
                    if already < 1
                        plotcomps{nx}(end+1) = gdcomps{nx}(a);
                        already = 1;
                    end;
                end;                    
            end;
        end;
    end;
end;

%%%%%%%%%%%
% Plot all the comps with diff peaks super-imposed
cols = {'k','c','m','b','g','r'};
figure; row = 4; col = 4; pl=1;
for nx = 1:length(allspec)
    if ~isempty(plotcomps{nx})
        for a = 1:length(plotcomps{nx})
            if pl > row*col
                axcopy; pl=1; figure;
                subplot(row, col, pl);
            else                
                subplot(row, col, pl);
            end;
            for cond = 1:length(conds)
                if cond == 1
                    X = mean(conds{cond}{nx}(:,fr),1);
                    X = X/std(X);
                else
                    compspec = conds{cond}{nx};
                    X = compspec(find(gdcomps{nx} == plotcomps{nx}(a)),fr);
                    X = X/std(X);
                end;
                ph = plot(freqs(fr),X,'b');  hold on;
                set(ph,'color',cols{cond});set(ph,'linewidth',1.5);
            end;pl = pl+1;
            set(gca,'xlim',[3 freqs(fr(end))]); set(gca,'xgrid','on'); 
            title(['Sb: ',int2str(nx),'; cp: ',int2str(plotcomps{nx}(a))]);
            %set(gca,'yticklabel',[]);
        end;
    end;
end;axcopy;
%%%%%%%%%%%%%%
% Plot mean chan spec with all alpha comps of each subj
cond = 3; % task baseline
figure;pl = 1;
for nx = 1:length(paths)
    if ~isempty(gdcomps{nx})
    cols = jet(length(gdcomps{nx}));
    subplot(5,5,pl)
    for comp = 1:length(gdcomps{nx})
        compspec = conds{cond}{nx};
        X = compspec(comp,fr);
        X = X/std(X);
        ph = plot(freqs(fr),X,'b');  hold on; 
        set(ph,'color',cols(comp,:));set(ph,'linewidth',1.5);
    end;pl = pl+1;
    X = mean(conds{cond}{nx}(:,fr),1);
    X = X/std(X);
    ph = plot(freqs(fr),X,'k');  hold on;set(ph,'linewidth',2);
    set(gca,'xlim',[3 freqs(fr(end))]); set(gca,'xgrid','on'); 
    title(['Sb: ',int2str(nx),'; cp: ',int2str(gdcomps{nx}(comp))]);
    end;
end;axcopy
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
 subplot(2,3,5)
 cond=4;
         compspec = conds{cond}{nx};
         X = compspec(comp,fr);
         X = X/std(X);
         ph = plot(freqs(fr),X,'b');  hold on;
         set(ph,'color',cols(comp,:));set(ph,'linewidth',1.5);
     set(gca,'xlim',[3 freqs(fr(end))]); set(gca,'xgrid','on'); 
     title(['Sb: ',int2str(nx),'; cp: ',int2str(gdcomps{nx}(comp))]);


%%%%%%%%%%%%%%
% collect data of interest into a single matrix for histogram plotting

cond =2;  % 1=meanchan; 2=fix,3=taskbase, 4=mem,5=ig,6=probe;
allfreqs = zeros(1,0);pl=1;
for nx = 1:length(conds{cond})
    allmax = allconds{cond}{nx};
    for a = 1:length(allmax)
        allfreqs(1,end+1:end+length(allmax{a})) = allmax{a};
        for b = 1:length(allmax{a})
        keeptrack(pl,:) = [nx,gdcomps{nx}(a)];pl = pl+1;
        end;
    end;
end;
figure; hist(allfreqs,30);
% currently 38 diff freqs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot dipole density for a given freq bin
findfreqs = sort(unique(allfreqs));
for ff = 1:length(findfreqs)
    cps = find(allfreqs == findfreqs(ff));
    relcomps = keeptrack(cps,:); 
    cpoi = cell(1,length(gdcomps));   
    for w = 1:length(cps)
    cpoi{relcomps(w,1)}(end+1) = relcomps(w,2);
    end;
    clustcps{ff} = cpoi;
    realfreq(ff) = freqs(find(freqs == findfreqs(ff)));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at dipole locations;
% for looking at dipoles later:
eeglab
paths = {'/ap82/Sternberg/','/cj82/Sternberg/','/ds76/Sternberg/','/ec81/Sternberg/','/jo74/Sternberg/','/ke70/Sternberg/','/km81/Sternberg/','/mk79/Sternberg/','/nf68/Sternberg/','/tp62/Sternberg/','/ds80/Sternberg/','/kb77/Sternberg/','/cz84/Sternberg/','/gm84/Sternberg/','/ts79/Sternberg/','/ny84/Sternberg/','/ft84/Sternberg/','/gv84/Sternberg/','/ka83/Sternberg/','/cy82/Sternberg/','/jb84/Sternberg/','/rd81/Sternberg/','/km81/Sternberg2/','/jo74/Shortst/','/bt78/','/as78/'};

for nx = 1:length(paths)  
    EEG = pop_loadset('sources.set' ,['/data/common1/stern/eeg/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/scott/matlab');
for clust = 1:26%length(clustcps)
    numsubj = 0;
    for x = 1:length(gdcomps)
        if ~isempty(clustcps{clust}{x})
            numsubj = numsubj+1;
        end;
    end;    
   allbesa =  EEG.dipfit.model(1); colset = zeros(1,length(clustcps{clust}));
    for ss = 1:length(clustcps{clust})
        if ~isempty(clustcps{clust}{ss})        
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            dipsources = EEG.dipfit.model(clustcps{clust}{ss}(1));
            for w = 1:length(clustcps{clust}{ss})
                dipsources(1,w) = EEG.dipfit.model(clustcps{clust}{ss}(w));
            end;           
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            colset(ss) = length(dipsources);
            dipsources = [];
        end; 
    end;
    allbesa(1) = [];
    %% Assigns each subject a different color; and makes subjidx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcols = jet(length(clustcps{clust})); 
    plotcol = cell(0); pl = 0;clear subjidxmat
    for ss = 1:length(clustcps{clust})
        if colset(ss) ~= 0
            for nn = 1:colset(ss)
                plotcol{pl+nn} = dcols(ss,:);
                subjidxmat(pl+nn) = ss;
            end;
            pl = pl+nn;
        end;
    end;
    optdipplot = {allbesa,'dipolelength',0,'gui','off','dipolesize',18,'image','mri','spheres','on','color',plotcol};
    figure; mydipoleentropy( optdipplot,  subjidxmat,'distance',10);  cbar;
    %figure; dipoledensity( optdipplot, 'subjind',subjidxmat, 'method','distance','methodparam',10, 'mri_view' ,'top', 'gridsize',16); 
 set(gcf,'color','w');    
ph = textsc(['Sternberg subject dipole density (10); ',num2str(realfreq(clust)),' Hz Only;  ',int2str(numsubj),' out of 23 subj'],'title');   
set(ph,'fontsize',16);
savename = ['print -djpeg /home/julie/Scott/PeakFreqs/SubjDipDens-',num2str(realfreq(clust)),'Hz.jpg'];
    eval(savename)
    close
end;  % for density plots

 figure; dipplot(allbesa, 'dipolesize',20,'image','mri','spheres','on','gui','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all dipoles for each subject and plot subjects individually with color-coded dipoles
lofreq = 4;hifreq = 13;
cols = jet(length(find(freqs(fr)> lofreq & freqs(fr)< hifreq))); 
dipsize = [8:1:length(find(freqs(fr)> lofreq & freqs(fr)< hifreq))+7];
fs = freqs(fr);fs = fs(find(fs > lofreq & fs < hifreq));
figure;sp=1;
for nx = 1:12%length(paths)
    if ~isempty(gdcomps{nx})
    allbesa =  EEG.dipfit.model(1); pl = 1;clear freqcol freqsize
        EEG = eeg_retrieve(ALLEEG, nx); CURRENTSET = nx;
        sd = find(keeptrack(:,1) == nx);
        frqs = allfreqs(sd);
        comps = keeptrack(sd,2);
        dipsources = EEG.dipfit.model(comps(1));
        for w = 1:length(comps)
            if frqs(w) > lofreq & frqs(w) < hifreq
            dipsources(1,pl) = EEG.dipfit.model(comps(w));
            freqcol{pl} = cols(find(freqs(fr) == frqs(w))-length(find(freqs(fr)<lofreq)),:);
            freqsize(pl) = dipsize(find(freqs(fr) == frqs(w))-length(find(freqs(fr)<lofreq))); pl = pl+1;
            end;
        end;           
        allbesa(end+1:end+size(dipsources,2)) = dipsources; 
        dipsources = [];       allbesa(1) = []; 
        subplot(3,4,sp)
        %dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',freqcol,'projlines','on','projimg','on'); view(60,30); camzoom(1);sp = sp+1;clear allbesa
        dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',freqcol); view(0,90); camzoom(1.1);sp = sp+1;clear allbesa
ph=text(-75,100,5,['Subject: ',int2str(nx)]);set(ph,'color','r');
    end;
end;
set(gcf,'color','w');  
ph=colorbar;  set(ph,'ytick',[0:.058:1]); set(ph,'yticklabel',[fs(1:end)]);
set(ph,'fontsize',7);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%  Plot scalp maps for each freq bin
 pfs = find(realfreq > 8 & realfreq < 13);
 for ff = 1:length(pfs)
     howmany = 0;figure; pl=1
     for xx = 1:length(clustcps)
         howmany = howmany+length(clustcps{pfs(ff)}{xx});
     end; row = ceil(sqrt(howmany));col = ceil(sqrt(howmany));    
     for nx = 1:length(gdcomps)
         if ~isempty(clustcps{pfs(ff)})
             EEG = eeg_retrieve(ALLEEG, nx); CURRENTSET = nx;
             for cp = 1:length(clustcps{pfs(ff)}{nx})
                 subplot(row,col,pl)
                 topoplot(EEG.icawinv(:,clustcps{pfs(ff)}{nx}(cp)),EEG.chanlocs,'electrodes','off');
                 pl = pl+1;title(['S ',int2str(nx),' cp ',int2str(clustcps{pfs(ff)}{nx}(cp))]);
             end;
         end;
     end;
     textsc(['All scalp maps with peak frequency: ',num2str(realfreq(pfs(ff))),' Hz'],'title');
 end;