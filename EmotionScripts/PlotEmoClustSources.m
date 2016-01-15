% plots sources weighted highly in ICA clustering as a page of heads, one for each emotion

eeglab
paths = {'/data/common2/emotion/tl81/','/data/common2/emotion/mi83/','/data/common2/emotion/ms82/','/data/common2/emotion/js75/','/data/common2/emotion/kw78/','/data/common2/emotion/tv81/','/data/common2/emotion/sr81/','/data/common2/emotion/an70/','/data/common2/emotion/sg75/'};
for nx = 1:length(paths)
EEG = pop_loadset( 'sources.set', paths{nx});
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;


%allemocell = { gdcompsawe,gdcompsfrustration, gdcompsjoy, gdcompsanger,gdcompssad,gdcompshappy,gdcompsfear,gdcompslove,gdcompsjealousy, gdcompscompassion,gdcompscontent, gdcompsgrief,gdcompsrelief,gdcompsdisgust,gdcompsexcite};

allsubjcell = {subj1gdcomps,subj2gdcomps,subj3gdcomps,subj4gdcomps,subj5gdcomps,subj6gdcomps,subj7gdcomps,subj8gdcomps,subj9gdcomps};
emos = {'prebase','awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
%emoset = [2,4,10,8,14,12,6,9,7,1,11,5,13,3,15];
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
figure; sb=1;
for emo = 1:length(emos)
    allsources =  EEG.dipfit.model(1);
    for nx = 1:length(allsubjcell)
        % find sources for gdcompsa set
        if ~isempty(allsubjcell{nx}{emoset(emo)})        
            EEG = eeg_retrieve(ALLEEG, nx); CURRENTSET = nx;
            dipsources= EEG.dipfit.model(1);
            for w = 1:length(allsubjcell{nx}{emoset(emo)})
                dipsources(1,w)= EEG.dipfit.model(allsubjcell{nx}{emoset(emo)}(w));
            end;           
            allsources(end+1:end+size(dipsources,2)) = dipsources;
            dipsources=[];    
        end;
    end;
    allsources(1) = [];
    if ~isempty(allsources)        
        subplot(4,5,sb);
        dipplot(allsources,'dipolelength',0,'gui','off','dipolesize',30,'color',{'r'},'image','mri','projimg','on','spheres','on');  %  ,'projlines','on','projcol',{'y'}
        view(50,20);
        camzoom(.92);
        ph=text(-80,0,100,emos(emoset(emo)));
        set(ph,'color','y');
        set(ph,'fontsize',15);
    end;sb = sb+1; 
end;
ph = textsc('Highly weighted components in Below-Baseline Alpha/Beta Power;','title');
set(ph,'fontsize',15);
set(ph,'color','r');
set(gcf,'color','w');


emos = {'prebase','awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
allsubjcell = {subj1gdcomps,subj2gdcomps,subj3gdcomps,subj4gdcomps,subj5gdcomps,subj6gdcomps,subj7gdcomps,subj8gdcomps,subj9gdcomps};
emoset = [1,3,5,11,9,15,13,7,10,8,2,12,6,14,4,16,17];
figure; sb=1; nx = 8;
EEG = eeg_retrieve(ALLEEG, nx); CURRENTSET = nx;
for emo = 1:length(emos)
    allsources =  EEG.dipfit.model(1);
    if ~isempty(allsubjcell{nx}{emoset(emo)})        
            dipsources= EEG.dipfit.model(1);
            for w = 1:length(allsubjcell{nx}{emoset(emo)})
                dipsources(1,w)= EEG.dipfit.model(allsubjcell{nx}{emoset(emo)}(w));
            end;           
            allsources(end+1:end+size(dipsources,2)) = dipsources;
            dipsources=[];    
    end;
    allsources(1) = [];
    if ~isempty(allsources)        
        subplot(4,5,sb);
        dipplot(allsources,'dipolelength',0,'gui','off','dipolesize',30,'color',{'r'},'image','mri','projimg','on','spheres','on');  %  ,'projlines','on','projcol',{'y'}
        view(50,20);
        camzoom(.92);
        ph=text(-80,0,100,emos(emoset(emo)));
        set(ph,'color','y');
        set(ph,'fontsize',15);
    end;sb = sb+1; 
end;
ph = textsc(['Highly weighted comps in Subject ',int2str(nx),', Above-Baseline Alpha/Beta'],'title');
set(ph,'fontsize',15);
set(ph,'color','r');
set(gcf,'color','w');




allsubjcomps = {subj1gdcomps};
emos = {'awe','frustration','joy','anger','happy','sad','love','fear' ,'compassion','jealousy','content','grief','relief','disgust','excite'};
emoset = [2,4,10,8,14,12,6,9,7,1,11,5,13,3,15];
figure; sb=1;
for emo = 1:length(emos)
    allsources =  EEG.dipfit.model(1);
    for nx = 1:length(allsubjcomps)
        % find sources for gdcompsa set
        if ~isempty(allsubjcomps{nx}{emoset(emo)})        
            EEG = eeg_retrieve(ALLEEG, nx); CURRENTSET = nx;
            dipsources= EEG.dipfit.model(1);
            for w = 1:length(allsubjcomps{nx}{emoset(emo)})
                dipsources(1,w)= EEG.dipfit.model(allsubjcomps{nx}{emoset(emo)}(w));
            end;           
            allsources(end+1:end+size(dipsources,2)) = dipsources;
            dipsources=[];    
        end;
    end;
    allsources(1) = [];
    if ~isempty(allsources)        
        subplot(4,4,sb);
        dipplot(allsources,'dipolelength',0,'gui','off','dipolesize',25,'color',{'r'},'image','mri','projimg','on','projlines','on','projcol',{'y'},'spheres','on');
        view(50,20);
        camzoom(.92);
        ph = title(emos(emoset(emo)));
        set(ph,'color','y');
        set(ph,'fontsize',15);
    end;sb = sb+1;
end;
ph = textsc('Highly weighted components in Cluster 1; Positive Weights;','title');
set(ph,'fontsize',15);
set(ph,'color','r');
set(gcf,'color','w');


 figure; nx=2;
for cp = 1:length(gdcomps{nx})
figure; nx=2;
dipplot(EEG.dipfit.model(gdcomps{nx}(cp)),'dipolelength',0,'gui','off','dipolesize',25,'color',{'r'},'image','mri','projimg','on','projlines','on','projcol',{'y'},'spheres','on');
view(50,30); camzoom(.8)
        ph=text(-80,0,100,int2str(gdcomps{nx}(cp)));
        set(ph,'color','y');
        set(ph,'fontsize',15);
set(gcf,'color','w');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

end;

for cp = 1:length(gdcomps{nx})
figure; nx=2;
    topoplot(EEG.icawinv(:,gdcomps{nx}(cp)),EEG.chanlocs,  'electrodes', 'off');
end;
