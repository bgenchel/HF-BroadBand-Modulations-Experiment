% plots three views of dipole clusters from ClustbyScalp.m
eeglab
paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/'};
ALLEEG=[];
for nx = 1:length(paths)
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
eeglab redraw;

nclust = 15;   col = 6;row = 4; cols = jet((nclust+1)/2);
figure; pl = 1; cc=1;
for clust = 1:nclust
    if clust == 9
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        set(gcf,'color','w');
        figure;pl = 1;cc=1;
    end;    
if clust == 1 
subj1 = [3  13  18  27]; 
subj2 = [21  29]; 
subj3 = []; 
subj4 = []; 
subj5 = [18  34]; 
subj6 = [1   6  17]; 
subj7 = [1  6  8]; 
subj8 = [1]; 
subj9 = [18  29  36]; 
subj10 = [2  16]; 
subj11 = [12  14  21]; 
subj12 = [8  25]; 
subj13 = [13  23]; 
subj14 = []; 
subj15 = [13  20]; 
subj16 = []; 
subj17 = []; 
subj18 = [11]; 
subj19 = [13  24  26]; 
subj20 = [14]; 
subj21 = [32]; 
end;
  %  K means dipole location Cluster: 2 ; Number of clusters requested: 15
if clust == 2 
subj1 = [19  34]; 
subj2 = [9]; 
subj3 = [12]; 
subj4 = [25]; 
subj5 = [11]; 
subj6 = [20]; 
subj7 = [7]; 
subj8 = []; 
subj9 = []; 
subj10 = []; 
subj11 = []; 
subj12 = [12]; 
subj13 = []; 
subj14 = [17]; 
subj15 = [11]; 
subj16 = [16]; 
subj17 = []; 
subj18 = [17]; 
subj19 = [18]; 
subj20 = [13]; 
subj21 = [20]; 
end;
  %  K means dipole location Cluster: 3 ; Number of clusters requested: 15
if clust == 3 
subj1 = [16  28]; 
subj2 = [15  16]; 
subj3 = [13  17]; 
subj4 = [22]; 
subj5 = [12  25]; 
subj6 = [10]; 
subj7 = []; 
subj8 = []; 
subj9 = []; 
subj10 = [9  26]; 
subj11 = [9]; 
subj12 = [20  29]; 
subj13 = []; 
subj14 = [13  24]; 
subj15 = []; 
subj16 = [8  15]; 
subj17 = [14]; 
subj18 = []; 
subj19 = [15  17]; 
subj20 = [16]; 
subj21 = []; 
end;
  %  K means dipole location Cluster: 4 ; Number of clusters requested: 15
if clust == 4 
subj1 = []; 
subj2 = [2  17]; 
subj3 = [15]; 
subj4 = [8]; 
subj5 = [4  28]; 
subj6 = []; 
subj7 = []; 
subj8 = []; 
subj9 = [11  19  24]; 
subj10 = [27]; 
subj11 = []; 
subj12 = [15]; 
subj13 = [3  19  20  25  30]; 
subj14 = []; 
subj15 = [1]; 
subj16 = []; 
subj17 = [2  7  8]; 
subj18 = [2  9]; 
subj19 = [23]; 
subj20 = [3  19]; 
subj21 = [2   9  17]; 
end;
  %  K means dipole location Cluster: 5 ; Number of clusters requested: 15
if clust == 5 
subj1 = []; 
subj2 = [7   8  24]; 
subj3 = []; 
subj4 = []; 
subj5 = []; 
subj6 = [3]; 
subj7 = []; 
subj8 = []; 
subj9 = [10  42]; 
subj10 = []; 
subj11 = [4   5   7  16]; 
subj12 = []; 
subj13 = []; 
subj14 = []; 
subj15 = [14  21]; 
subj16 = [5]; 
subj17 = []; 
subj18 = [3  7]; 
subj19 = [5]; 
subj20 = []; 
subj21 = [3   4  10  14  18  23]; 
end;
  %  K means dipole location Cluster: 6 ; Number of clusters requested: 15
if clust == 6 
subj1 = [5  10  25]; 
subj2 = [4]; 
subj3 = [7 10]; 
subj4 = [1   4   5  10  14]; 
subj5 = [5   6  22]; 
subj6 = [23]; 
subj7 = [10]; 
subj8 = [8]; 
subj9 = [3   6   7  14]; 
subj10 = [3   8  13]; 
subj11 = [1   2   8  13]; 
subj12 = [17  21]; 
subj13 = [5  8 10]; 
subj14 = [10]; 
subj15 = [2   6  15  16  23  37]; 
subj16 = [4]; 
subj17 = [10]; 
subj18 = [5  8]; 
subj19 = [10  14  21  22]; 
subj20 = [7]; 
subj21 = []; 
end;
  %  K means dipole location Cluster: 7 ; Number of clusters requested: 15
if clust == 7 
subj1 = []; 
subj2 = [5]; 
subj3 = []; 
subj4 = []; 
subj5 = []; 
subj6 = []; 
subj7 = []; 
subj8 = []; 
subj9 = [37]; 
subj10 = [5]; 
subj11 = [15]; 
subj12 = []; 
subj13 = [1  2]; 
subj14 = []; 
subj15 = [4]; 
subj16 = [9]; 
subj17 = [1]; 
subj18 = []; 
subj19 = []; 
subj20 = [6]; 
subj21 = [1  11  21  22  24]; 
end;
  %  K means dipole location Cluster: 8 ; Number of clusters requested: 15
if clust == 8 
subj1 = [11  21  30]; 
subj2 = [6  19]; 
subj3 = []; 
subj4 = [17]; 
subj5 = [24]; 
subj6 = [16]; 
subj7 = []; 
subj8 = []; 
subj9 = []; 
subj10 = [7  33]; 
subj11 = [34]; 
subj12 = []; 
subj13 = [6  15]; 
subj14 = [4  11  14]; 
subj15 = [38]; 
subj16 = []; 
subj17 = [12]; 
subj18 = [10]; 
subj19 = [11]; 
subj20 = []; 
subj21 = [15]; 
end;
  %  K means dipole location Cluster: 9 ; Number of clusters requested: 15
if clust == 9 
subj1 = [7   8   9  22]; 
subj2 = [12]; 
subj3 = []; 
subj4 = [2]; 
subj5 = [10]; 
subj6 = [15]; 
subj7 = []; 
subj8 = []; 
subj9 = [5]; 
subj10 = []; 
subj11 = [22]; 
subj12 = []; 
subj13 = []; 
subj14 = [16  20  21]; 
subj15 = []; 
subj16 = [3  18]; 
subj17 = [16]; 
subj18 = []; 
subj19 = []; 
subj20 = [12]; 
subj21 = [19]; 
end;
  %  K means dipole location Cluster: 10 ; Number of clusters requested: 15
if clust == 10 
subj1 = []; 
subj2 = []; 
subj3 = []; 
subj4 = []; 
subj5 = []; 
subj6 = [8]; 
subj7 = [22]; 
subj8 = [11]; 
subj9 = []; 
subj10 = []; 
subj11 = []; 
subj12 = []; 
subj13 = [27]; 
subj14 = []; 
subj15 = []; 
subj16 = []; 
subj17 = []; 
subj18 = []; 
subj19 = []; 
subj20 = []; 
subj21 = [6  12]; 
end;
  %  K means dipole location Cluster: 11 ; Number of clusters requested: 15
if clust == 11 
subj1 = [1  15  17  20]; 
subj2 = []; 
subj3 = [22]; 
subj4 = [19]; 
subj5 = [1  14  23]; 
subj6 = [11]; 
subj7 = []; 
subj8 = [12]; 
subj9 = [8]; 
subj10 = [6]; 
subj11 = [23  32]; 
subj12 = []; 
subj13 = []; 
subj14 = [12  19]; 
subj15 = [40]; 
subj16 = [2 10]; 
subj17 = [11]; 
subj18 = [18  19]; 
subj19 = [4  27  30]; 
subj20 = [21]; 
subj21 = []; 
end;
  %  K means dipole location Cluster: 12 ; Number of clusters requested: 15
if clust == 12 
subj1 = [14]; 
subj2 = [22]; 
subj3 = [5]; 
subj4 = []; 
subj5 = []; 
subj6 = [2]; 
subj7 = [15]; 
subj8 = []; 
subj9 = [13  17]; 
subj10 = [14]; 
subj11 = [17  18]; 
subj12 = [2  26]; 
subj13 = [9  14]; 
subj14 = [1  3]; 
subj15 = [7  18  25]; 
subj16 = []; 
subj17 = [3  19]; 
subj18 = []; 
subj19 = [7  19]; 
subj20 = [1 10]; 
subj21 = [5  16]; 
end;
  %  K means dipole location Cluster: 13 ; Number of clusters requested: 15
if clust == 13 
subj1 = [26  29]; 
subj2 = [26]; 
subj3 = [18]; 
subj4 = [9  12]; 
subj5 = [16]; 
subj6 = []; 
subj7 = []; 
subj8 = []; 
subj9 = [9  16  20]; 
subj10 = [15]; 
subj11 = [29]; 
subj12 = [14]; 
subj13 = [4  17]; 
subj14 = [2  9]; 
subj15 = []; 
subj16 = [11]; 
subj17 = []; 
subj18 = [16]; 
subj19 = [8  33]; 
subj20 = [15]; 
subj21 = []; 
end;
  %  K means dipole location Cluster: 14 ; Number of clusters requested: 15
if clust == 14 
subj1 = [12]; 
subj2 = [13  18]; 
subj3 = [14  19]; 
subj4 = [13  16]; 
subj5 = [9  19]; 
subj6 = [14]; 
subj7 = []; 
subj8 = [9]; 
subj9 = [22]; 
subj10 = []; 
subj11 = [10]; 
subj12 = []; 
subj13 = []; 
subj14 = [18  23  28]; 
subj15 = []; 
subj16 = [7  25]; 
subj17 = [13]; 
subj18 = [20]; 
subj19 = [12  20  28]; 
subj20 = [18]; 
subj21 = []; 
end;
  %  K means dipole location Cluster: 15 ; Number of clusters requested: 15
if clust == 15 
subj1 = [4   6  36]; 
subj2 = [1]; 
subj3 = [6]; 
subj4 = [3]; 
subj5 = [2   3  15  17]; 
subj6 = [5  9]; 
subj7 = [2  5  9]; 
subj8 = [2  6]; 
subj9 = [1   2  23  30  32]; 
subj10 = [12]; 
subj11 = [6  11]; 
subj12 = [5  7]; 
subj13 = []; 
subj14 = [6  8]; 
subj15 = [3  5]; 
subj16 = [17]; 
subj17 = [6]; 
subj18 = [12  14]; 
subj19 = [6]; 
subj20 = [11]; 
subj21 = [8]; 
end;
gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21};

    allbesa =  EEG.dipfit.model(1);
    for ss = 1:length(gdcomps)
        if ~isempty(gdcomps{ss})        
            EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
            dipsources= EEG.dipfit.model(gdcomps{ss}(1));
            for w = 1:length(gdcomps{ss})
                dipsources(1,w)= EEG.dipfit.model(gdcomps{ss}(w));
            end;           
            allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            dipsources=[];
        end; 
    end;
    allbesa(1) = [];
    subplot(row,col,pl)
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(cc,:)});pl = pl+1;
    ph =title(['Cluster ',int2str(clust)]);
    set(ph,'color','w');
    set(ph,'fontsize',14);
    subplot(row,col,pl)
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(cc,:)});pl = pl+1; view(90,0)
    subplot(row,col,pl)
    dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',25,'dipolelength',0,'spheres','on','color',{cols(cc,:)});pl = pl+1; cc = cc+1;view(0,0)   
end;
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color','w');
