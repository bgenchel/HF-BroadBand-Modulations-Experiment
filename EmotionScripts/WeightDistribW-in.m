% explores the within subject weights for each trial (and trial distibutions for each emotion) to see if anything is more telling of similar emotions other than the straight mean.

load /data/common2/emotion/clusters/subjdims.mat % subjdims numtrials pcadims wtsfile sphfile comment
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones

zscorecut = 2; %zscore cutoff for positive and neg weighted trials
figure;percdiff = 5; % (currently absolute percent) percent diff between pos and neg percents to be counted as a sig factor for that emotion.
clear sadhappfacs happyorsadall
for nx = 1:length(gdcomps)
    clear sph wts icamatall winv wtdis ws
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20        
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    activations = (wts*sph)*icamatall;    winv = pinv(wts*sph);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end; 
    for fac = 1:size(winv,2)
        wtdis = winv(:,fac)/std(winv(:,fac));
        postrials = find(wtdis > zscorecut);    
        negtrials = find(wtdis < -zscorecut); clear tottrials etrials percNEGtrials percPOStrials
        for e = 1:length(emomap)-1
            tottrials = length(emomap(e):emomap(e+1)-1);
            etrials = length(find(negtrials >= emomap(e) & negtrials < emomap(e+1)));
            percNEGtrials(1,e) = etrials/tottrials*100;clear etrials
            etrials = length(find(postrials >= emomap(e) & postrials < emomap(e+1)));
            percPOStrials(1,e) = etrials/tottrials*100;
        end;
        allnegperc{fac} = percNEGtrials;
        allposperc{fac} = percPOStrials;
    end;

    cols = jet(length(emomap)-1);clf
    %figure; 
    row = round(sqrt(size(winv,2))); col = ceil(sqrt(size(winv,2))); 
    for fac = 1:size(winv,2)
        subplot(row,col,fac)
        for e = 1:length(emomap)-1
            ph = bar(e,allnegperc{fac}(1,e)); hold on;
            set(ph,'facecolor',cols(e,:));
            ph = text(e,0,emos{e}); set(ph,'rotation',90);
        end;    
        yl = get(gca,'ylim'); if yl(2) < 6
            set(gca,'ylim',[0 6]); end;
            set(gca,'xlim',[0 16]);set(gca,'xtick',[1:15]);set(gca,'xticklabel',[]);
            title(['Factor ',int2str(fac)]);
    end;
    textsc(['Subj ',int2str(nx),' Percent of total trials with NEGATIVE weights below zscore of -',int2str(zscorecut)],'title');
    
    cols = jet(length(emomap)-1);clf
    %figure; 
    row = round(sqrt(size(winv,2))); col = ceil(sqrt(size(winv,2))); 
    for fac = 1:size(winv,2)
        subplot(row,col,fac)
        for e = 1:length(emomap)-1
            ph = bar(e,allposperc{fac}(1,e)); hold on;
            set(ph,'facecolor',cols(e,:));
            ph = text(e,0,emos{e}); set(ph,'rotation',90);
        end;    
        yl = get(gca,'ylim'); if yl(2) < 6
            set(gca,'ylim',[0 6]); end;
            set(gca,'xlim',[0 16]);set(gca,'xtick',[1:15]);set(gca,'xticklabel',[]);
            title(['Factor ',int2str(fac)]);
    end;
    textsc(['Subj ',int2str(nx),' Percent of total trials with POSITIVE weights above norm score of ',int2str(zscorecut)],'title');

    clear efacs
    for e = 1:length(emomap)-1
        majfacs = []; m=1;
        for fac = 1:size(winv,2)
            if  allposperc{fac}(1,e) > percdiff  %- allnegperc{fac}(1,e)
                majfacs(1,m) = fac; m = m+1;
            end;
            if  allnegperc{fac}(1,e) > percdiff  % - allposperc{fac}(1,e)
                majfacs(1,m) = -fac; m = m+1;
            end;
        end;
        efacs{e} =  majfacs;
    end;
    subjefacs{nx} = efacs;
    
    %%% find specific diff between happy and sad (sad = 7, happy = 13)
        efacs = []; happyorsad = []; m=1; pd = 10; %pd = percent diff cutoff
        for fac = 1:size(winv,2)
            if  (allposperc{fac}(1,7)  - allposperc{fac}(1,13)) > pd
                happyorsad(1,m) = fac;% sad
                efacs(1,m) = fac; m = m+1;
            elseif (allposperc{fac}(1,13)  - allposperc{fac}(1,7)) > pd
                happyorsad(2,m) = fac; % happy
                efacs(1,m) = fac; m = m+1;                
            end;
            if  (allnegperc{fac}(1,7) - allnegperc{fac}(1,13))> pd
                happyorsad(1,m) = -fac;  %sad
                efacs(1,m) = -fac; m = m+1;
            elseif (allnegperc{fac}(1,13) - allnegperc{fac}(1,7))> pd
                happyorsad(2,m) = -fac; %happy
                efacs(1,m) = -fac; m = m+1;                
            end;
        end;
        sadhappfacs{nx} =  efacs;
        happyorsadall{nx} = happyorsad;
        allnegperc{fac}(1,13)
        allnegperc{fac}(1,7)
        allposperc{fac}(1,13)
        allposperc{fac}(1,7)
        
    fprintf('\nSubject %s Done...',int2str(nx));    
end;
clear allnegperc percNEGtrials allposperc percPOStrials majfacs efacs tottrials  etrials sph wts winv icamatall activations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find differences in weights between two emotions
em1 = 7; % sad
em2 = 13; % happy
alpha = .0000000001;

clear sadhappfacs happyorsadall multfact subjefacs
for nx = 1:length(gdcomps)
    fprintf('\nWorking on subject %s ... ',int2str(nx));
    clear sph wts icamatall winv wtdis ws activations 
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20        
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    activations = (wts*sph)*icamatall;    winv = pinv(wts*sph);    emomap = ones(1,1);
    clear icamatall sph wts     
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end; clear mean1 mean2 p
    for fac = 1:size(winv,2)
        wtdis = winv(:,fac)';
        tottrials(1,1) = length(emomap(em1):emomap(em1+1));
        tottrials(1,2) = length(emomap(em2):emomap(em2+1));
        mean1(1,fac) = mean(wtdis(emomap(em1):emomap(em1+1)));
        mean2(1,fac) = mean(wtdis(emomap(em2):emomap(em2+1)));
        [hh p(1,fac)] = ttest2(wtdis(emomap(em1):emomap(em1+1)),wtdis(emomap(em2):emomap(em2+1)),.01,'both');
        %[d(1,fac),P(1,fac),stats] =  manova1([wtdis(emomap(em1):emomap(em1+1)),wtdis(emomap(em2):emomap(em2+1))]',[ones(1,tottrials(1,1)),ones(1,tottrials(1,2))*2]',.01);
        submean(1,fac) = mean1(1,fac) - mean2(1,fac);         
    end;
    multfact{nx}(1,:) = mean1;
    multfact{nx}(2,:) = mean2;          
    subjefacs{nx} =  find(p<alpha);    
end;
%%%%%%%%%%%%
% Plot dipoles and spectra for each emotion given by percent highly weighted trials above
% rearrange factors into emotion clusters
load /data/common2/emotion/nxlists.mat nxlists
             
freqs = [1:.5:50];fr = find(freqs > 3 & freqs<36); % freqs to plot   
figure;pl = 1; row = 6; col = 6;
for nx = 1:length(gdcomps)
    ALLEEG=[];EEG=[];
    clear sph wts icamatall winv wtdis ws activations
    addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    activations = (wts*sph)*icamatall; 
    if ~isempty(subjefacs{nx})
        for f = 1:length(subjefacs{nx}) 
            fac = subjefacs{nx}(f);
            dipsources = [];dipsources = EEG.dipfit.model(nxlists{nx}{fac}(1));
            ccs1 = winter(length(nxlists{nx}{fac}));   clear cols 
            ccs2 = autumn(length(nxlists{nx}{fac}));   clear cols 
            if pl == row*col+1
                set(gcf,'color','w');
                textsc(['Factors with significantly different weights for Sad and Happy'],'title');
                set(gcf,'PaperOrientation','landscape');
                set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                figure; pl = 1;
            end;                
            for c = 1:length(nxlists{nx}{fac})
                dipsources(1,c) = EEG.dipfit.model(nxlists{nx}{fac}(c));
                cols1{c} = ccs1(c,:);
                cols2{c} = ccs2(c,:);
            end;           
            subplot(row,col,pl)
            dipplot(dipsources,'image','mri','gui','off','normlen','on','dipolesize',30,'dipolelength',0 ,'color',cols2,'projlines','on','projimg','on','spheres','on');pl = pl+1;  view(60,30);
            subplot(row,col,pl)
            for c = 1:length(nxlists{nx}{fac})
                cp = find(nxlists{nx}{fac}(c) == gdcomps{nx});
                oneact = activations(fac,length(freqs)*(cp-1)+1:length(freqs)*cp);
                ph = plot(freqs(fr),oneact(fr)*multfact{nx}(1,f),'k-'); hold on;
                set(ph,'color',cols1{c}); set(ph,'linewidth',2);
                ph = plot(freqs(fr),oneact(fr)*multfact{nx}(2,f),'k-'); hold on;
                set(ph,'color',cols2{c}); set(ph,'linewidth',2);
            end; set(gca,'xlim',[3 35]);
            %ph = title(emos{e}); set(ph,'color','r');
            ph = title(['S ',int2str(nx),' Fac ',int2str(fac)]); set(ph,'color','r');
            set(gca,'xgrid','on'); pl = pl+1;                
        end;            
    end; 
end;
set(gcf,'color','w');
textsc(['Factors with significantly different weights for Sad and Happy'],'title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find highest variance components for each spectral factor:

for nx = 1:length(gdcomps)
    clear sph wts icamatall winv wtdis ws
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20        
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    activations = (wts*sph)*icamatall;    clear newact
    for tp = 1:size(activations,1)
        for cp = 1:length(gdcomps{nx})
            newact(tp,cp,:) = activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp);
        end;      
    end;
    for tp = 1:size(activations,1)
        for ff = 1:length(freqs)
            tmpact = newact(tp,:,ff);
            tvar(tp,ff) = var(tmpact);
        end;
    end;
    % to Normalize by lowest variance freqencies for each template
    for tp = 1:size(activations,1)
        mvar = min(tvar(tp,:));
        lowvar{tp} = find(tvar(tp,:) >= mvar & tvar(tp,:) <= mvar+mvar*.4); % find freqs with low var
        tmpcps = newact(tp,:,lowvar{tp}); tmpcps = squeeze(tmpcps); % makes a cp X freqs
        newact(tp,:,:) = newact(tp,:,:)/mean(std(tmpcps'));   % mean of std over all comps        
    end;
    
    % find most variant components by zscore cutoff
    clear tplists
    for tp = 1:size(activations,1)
        cplist = zeros(0);forstd = newact(tp,:,:); forstd = squeeze(forstd);% makes a cp X freqs
        cut = mean(std(forstd'))+std(std(forstd'));  % includes comps from mean to .3 *std
        for cp = 1:length(gdcomps{nx})
            if std(forstd(cp,:)) > cut
                cplist(end+1) = gdcomps{nx}(cp);
            end;
        end;
        tplists{tp} = cplist;
    end;
    nxlists{nx} = tplists;
    
    fprintf('\n SUBJECT %i Done',nx);
end;
clear newact sbsph sbwts forsts sph wts icamatall activations tmpact tvar mvar lowvar tmpcps cplist tplists
save /data/common2/emotion/nxlists.mat nxlists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dipoles and spectra for each emotion given by percent highly weighted trials above
% rearrange factors into emotion clusters
load /data/common2/emotion/nxlists.mat nxlists
                
freqs = [1:.5:50];fr = find(freqs > 3 & freqs<36); % freqs to plot   
figure;pl = 1; 
for nx = 1:length(gdcomps)
    ALLEEG=[];EEG=[];
    clear sph wts icamatall winv wtdis ws activations
    addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
    if length(gdcomps{1}) < 20
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    end;
    activations = (wts*sph)*icamatall; 
    e=1;
    %for e = 1:length(emos)
        if ~isempty(subjefacs{nx}{e})
            for f = 1:length(subjefacs{nx}{e}) 
                fac = subjefacs{nx}{e}(f);
                if fac < 0
                    mult = -1; fac = abs(fac);
                else
                    mult = 1;
                end;               
                dipsources = [];dipsources = EEG.dipfit.model(nxlists{nx}{fac}(1));
                ccs = jet(length(nxlists{nx}{fac}));   clear cols 
                if pl == row*col+1
                    set(gcf,'color','w');
                    textsc([emos{e},': Highly weighted factors for given emotions'],'title');
                    set(gcf,'PaperOrientation','landscape');
                    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                    figure; pl = 1;
                end;                
                for c = 1:length(nxlists{nx}{fac})
                    dipsources(1,c) = EEG.dipfit.model(nxlists{nx}{fac}(c));
                    cols{c} = ccs(c,:);
                end;           
                subplot(row,col,pl)
                dipplot(dipsources,'image','mri','gui','off','normlen','on','dipolesize',30,'dipolelength',0 ,'color',cols,'projlines','on','projimg','on','spheres','on');pl = pl+1;  view(60,30);
                subplot(row,col,pl)
                for c = 1:length(nxlists{nx}{fac})
                    cp = find(nxlists{nx}{fac}(c) == gdcomps{nx});
                    oneact = activations(fac,length(freqs)*(cp-1)+1:length(freqs)*cp);
                    ph = plot(freqs(fr),oneact(fr)*mult,'k-'); hold on;
                    set(ph,'color',cols{c}); set(ph,'linewidth',2);
                end; set(gca,'xlim',[3 35]);
                %ph = title(emos{e}); set(ph,'color','r');
                ph = title(['Subj ',int2str(nx)]); set(ph,'color','r');
                set(gca,'xgrid','on'); pl = pl+1;                
            end;            
        end; 
    end;
    set(gcf,'color','w');
    %textsc(['Subject ',int2str(nx),' Highly weighted factors for given emotions'],'title');
    textsc([emos{e},' Highly weighted factors for all subjs'],'title');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dipoles and spectra for factors specified by sad-happy differential
load /data/common2/emotion/nxlists.mat nxlists
figure;row = 6; col = 6;pl=1;fr = find(freqs>3 & freqs<35);
for nx = 10:length(gdcomps)
    if ~isempty(sadhappfacs{nx})        
        ALLEEG=[];EEG=[];
        clear sph wts icamatall winv wtdis ws activations
        addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
        EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
        rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
        sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
        wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[pcadims{nx} subjdims{nx}(1)]); 
        if length(gdcomps{1}) < 20
            icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
        else
            icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
        end;
        activations = (wts*sph)*icamatall; 
        for f = 1:length(sadhappfacs{nx}) 
            fac = sadhappfacs{nx}(f);
            if fac < 0
                mult = -1; fac = abs(fac);
            else
                mult = 1;
            end;              
            dipsources = [];dipsources = EEG.dipfit.model(nxlists{nx}{fac}(1));
            ccs = jet(length(nxlists{nx}{fac}));   clear cols 
            if pl == row*col+1
                set(gcf,'color','w');
                textsc(['Highly weighted factors for sad vs happy'],'title');
                set(gcf,'PaperOrientation','landscape');
                set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                figure; pl = 1;
            end;                
            for c = 1:length(nxlists{nx}{fac})
                dipsources(1,c) = EEG.dipfit.model(nxlists{nx}{fac}(c));
                cols{c} = ccs(c,:);
            end;           
            subplot(row,col,pl)
            dipplot(dipsources,'image','mri','gui','off','normlen','on','dipolesize',30,'dipolelength',0 ,'color',cols,'projlines','on','projimg','on','spheres','on');pl = pl+1;  view(60,30);
            subplot(row,col,pl)
            for c = 1:length(nxlists{nx}{fac})
                cp = find(nxlists{nx}{fac}(c) == gdcomps{nx});
                oneact = activations(fac,length(freqs)*(cp-1)+1:length(freqs)*cp);
                ph = plot(freqs(fr),oneact(fr)*mult,'k-'); hold on;
                set(ph,'color',cols{c}); set(ph,'linewidth',2);
            end; set(gca,'xlim',[3 35]);
            %ph = title(emos{e}); set(ph,'color','r');
            if happyorsadall{nx}(1,f) ~= 0
                ph = title(['Subj ',int2str(nx),' sad']); 
            elseif happyorsadall{nx}(2,f) ~= 0
                ph = title(['Subj ',int2str(nx),' happy']); 
            end;  set(ph,'color','r');          
            set(gca,'xgrid','on'); pl = pl+1;                
        end;            
    end; 
end;
set(gcf,'color','w');
textsc([' Highly weighted factors for sad vs happy'],'title');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect high variance spectra for each relevant factor as well as 3D loc for clustering.
button = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
button = [13:21,23:34]; % no button press (apart from the first one) and no mr72-2
button = [1,2,4:6,8:12,14,17:21,23,25:30,31,33,34]; % all 'good' subjects (ones that said they got into it)
button = [1:21,23:34];  % not mr72-2
button = [1,3:9,12,14,16,17,19,21,22,23,24,26,27]; % females
button = [1,4:6,8,9,12,14,17,19,21,23,26,27]; % 'good' females
button = [2,10,11,13,15,18,20,25]; % males
button = [2,10,11,18,20,25,34]; % 'good' males

eeglab
freqs = [1:.5:50];fr = find(freqs > 3 & freqs<36); % freqs to plot   
alldips = zeros(0,6); allspec = zeros(0,length(fr)); pp = 1;
for nx = 1:length(button)
    ALLEEG=[];EEG=[];
    clear sph wts icamatall winv wtdis ws activations
    addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{button(nx)}],[pcadims{button(nx)} subjdims{button(nx)}(1)]); 
    if length(gdcomps{1}) < 20
        icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    else
        icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    end;
    activations = (wts*sph)*icamatall; 
    for e = 1:length(emos)
        if ~isempty(subjefacs{button(nx)}{e})
            for f = 1:length(subjefacs{button(nx)}{e}) 
                fac = subjefacs{button(nx)}{e}(f);
                if fac < 0
                    mult = -1; fac = abs(fac);
                else
                    mult = 1;
                end;    clear dipsources           
                if ~isempty(nxlists{button(nx)}{fac})
                    for c = 1:length(nxlists{button(nx)}{fac})
                        tp = EEG.dipfit.model(nxlists{button(nx)}{fac}(c)).posxyz;
                        dipsources(c,1:3) = tp(1,:);
                        dipsources(c,4:6) = tp(2,:);
                        clear tp
                    end;  
                    alldips(end+1:end+c,:) = dipsources;
                    clear allact
                end;
                if ~isempty(nxlists{button(nx)}{fac})
                    for c = 1:length(nxlists{button(nx)}{fac})
                        cp = find(nxlists{button(nx)}{fac}(c) == gdcomps{button(nx)});
                        allact(c,:) = activations(fac,length(freqs)*(cp-1)+1:length(freqs)*cp);   
                        keeptrack(pp,:) = [button(nx),e,f,gdcomps{button(nx)}(cp)]; pp = pp+1;
                    end; 
                    allspec(end+1:end+c,:) = allact(:,fr);
                end;
            end;            
        end; 
    end;
    fprintf('\n SUBJECT %i Done.\n',button(nx));
end;
clear activations icamatall wts sph dipsources allact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 addpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
ALLEEG=[];EEG=[];
for nx = 1:length(paths)  
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
    rmpath('/data/common/matlab/eeglab/plugins/dipfit2.0/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decompose only one emotion at a time:
figure;row = 8; col =8;  qq = 1;
nclusts = 4;
for e = 1:15
    if qq > 64 
        set(gcf,'color','w');figure;qq = 1;
    end;    
    oneemo = find(keeptrack(:,2) == e);
    % decompose dipoles with PCA
    [pc,eigvec,sv] = runpca(alldips(oneemo,:));
    clustacts{1} = eigvec'; clustwinv{1} = pc'; clustsv{1} = sv(find(sv));

    % decompose spectra with PCA
    [pc,eigvec,sv] = runpca(allspec(oneemo,:));
    clustacts{2} = eigvec'; clustwinv{2} = pc'; clustsv{2} = sv(find(sv));

    % determine # dims to use by using all dims accounting for 5% or more of data
    for q = 1:length(clustsv)
        if ~isempty(clustsv{q})
            maxdims(1,q) =  max(find(clustsv{q}/max(clustsv{q})>.05));
        end;
    end;

    forclust = zeros(size(clustacts{1},2),0);
    for ei = 1:length(clustacts)
        if ~isempty(clustacts{ei})
            forclust(:,end+1:end+maxdims(ei)) = clustacts{ei}(1:maxdims(ei),:)'/mean(std(clustacts{ei}(1:maxdims(ei),:),2));
        end;
    end;
 [kout, C,sumd, allsums] = kmeans(allspec(oneemo,:),nclusts,'replicates',5 );
 %[weights,sphere,compvars,bias,signs,lrates,activations] = runica(x,'stop',1e-7,'maxsteps',2000,'extended',1,'pca',3);
 
 %[kout, C,sumd, allsums] = kmeans(forclust,nclusts,'replicates',5 );
  %[kout, C,sumd, allsums] = kmeans(clustacts{2}(1:65,:)',nclusts,'replicates',5 );
   clear clustcps
    for clust = 1:size(C,1)
        cp = find(kout == clust); 
        cp(find(abs(zscore(allsums(cp,clust)))>20)) = [];
        % deleting > ? stds cleans up signif 
        cpidx{clust} = cp;
        kc = keeptrack(oneemo,:);
        relcomps = kc(cp,:); 
        cpoi = cell(1,length(gdcomps));   
        for w = 1:length(cp)
            cpoi{relcomps(w,1)}(end+1) = relcomps(w,4);
        end;
        clustcps{clust} = cpoi;
    end;
    emoclusts{e} = clustcps; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot resulting cluster dipoles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cols = jet(15);
    for clust = 1:4%nclusts       
        numsubj = 0;numdips = 0;
        for x = 1:length(gdcomps)
            if ~isempty(clustcps{clust}{x})
                numsubj = numsubj+1;
                numdips = numdips+length(clustcps{clust}{x});
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
        subplot(row,col,qq);
        dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(e,:)},'projlines','on','projimg','on');qq = qq+1;
        %if orient == 1
        %ph=text(60,-100,145,['Clust: ',int2str(clust)]);
        %else
        %    ph=text(60,-100,145,['Clust: -',int2str(clust)]);    
        %end;
        %set(ph,'color','y'); set(ph,'fontsize',14);
        view(60,20);
        %textsc(['Spectral Co-Modulation Clusters Dipoles: ',emos{e}],'title');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot resulting cluster spectra
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cp = find(kout == clust); 
        cp(find(abs(zscore(allsums(cp,clust)))>20)) = [];
        % deleting > ? stds cleans up signif 
        cpidx{clust} = cp;
        kc = keeptrack(oneemo,:);
        %relcomps = kc(cp,:); 
        specidx = oneemo(cp);
        plotspec = allspec(specidx,:);
        subplot(row,col,qq)
        ph = plot(freqs(fr),plotspec,'k-'); hold on;
        set(ph,'color',[.5 .5 .5]);set(ph,'linewidth',.75);    
        
        ph = plot(freqs(fr),mean(plotspec,1)); qq = qq+1;
        set(ph,'color',cols(e,:)); set(ph,'linewidth',2.5);
        set(gca,'xlim',[3 35]);set(gca,'xgrid','on');
        if clust == 1
            ph = title(emos{e});set(ph,'color','r');set(ph,'fontsize',16);
        end;
        textsc(['Spectral Co-Modulation Clusters: ',emos{e}],'title');
    end;
        fprintf('\n EMOTION %s Done.\n',emos{e});
end;
set(gcf,'color','w');
