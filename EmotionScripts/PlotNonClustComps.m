% finds non clustered comps from all subjs and plots
% call up sources for future plotting
%  eeglab
paths = {'/data/common1/stern/eeg/ap82/Sternberg/','/data/common1/stern/eeg/cj82/Sternberg/','/data/common1/stern/eeg/ds76/Sternberg/','/data/common1/stern/eeg/ec81/Sternberg/','/data/common1/stern/eeg/jo74/Sternberg/','/data/common1/stern/eeg/ke70/Sternberg/','/data/common1/stern/eeg/km81/Sternberg/','/data/common1/stern/eeg/mk79/Sternberg/','/data/common1/stern/eeg/nf68/Sternberg/','/data/common1/stern/eeg/tp62/Sternberg/','/data/common1/stern/eeg/ds80/Sternberg/','/data/common1/stern/eeg/kb77/Sternberg/','/data/common1/stern/eeg/cz84/Sternberg/','/data/common1/stern/eeg/gm84/Sternberg/','/data/common1/stern/eeg/ts79/Sternberg/','/data/common1/stern/eeg/ny84/Sternberg/','/data/common1/stern/eeg/ft84/Sternberg/','/data/common1/stern/eeg/gv84/Sternberg/','/data/common1/stern/eeg/ka83/Sternberg/','/data/common1/stern/eeg/cy82/Sternberg/','/data/common1/stern/eeg/jb84/Sternberg/','/data/common1/stern/eeg/rd81/Sternberg/','/data/common1/stern/eeg/km81/Sternberg2/','/data/common1/stern/eeg/jo74/Shortst/','/data/common1/stern/eeg/bt78/','/data/common1/stern/eeg/as78/'};
for nx = 1:length(paths)
        EEG = pop_loadset( 'sources.set', paths{nx});
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end;
EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
eeglab redraw;
%
%%  All "good" compe
ap82 = [3:10,12,14:20,23:25,27:29,34,35];% ap82
cj82 = [4,6,7,9:14 16:18,20 22];%cj82
ds76 = [6,7,8,10:14,16:18,20 23 26 29 30 31];%ds76
ec81 = [3,5:10,12,13,16,17,20,22:28,30,32:33,37:39];  %ec81
jo74 = [2 3 4 6 7 8 9 11 12 14 16,17,19:21,23 25 26 28 35 37];%jo74  +31 is good
ke70 = [6 8,9,12 14 19 26];%ke70 % 42 wasn't in original analysis
km81=  [2,4:16,18:20,22,23,26,30,31,33,38:40]; % km81 not noise comps
mk79 = [5:9 11,13:19,22,23 28 30,31];%mk79 
nf68 = [2,5:12,15,18:24,28,35];  %nf68
tp62 = [3,4,8:10,15,16];  %tp62
ds80 = [2:5,7:16,18,23,25]; %ds80
kb77 = [2:4,6:12,14,16:20,23,31,33,42,44]; %kb77
subj13 = [5,7,9,17,18,32];%cz84
subj14 = [3,5:12,14,16,17,18];%gm84
subj15 = [5,6,10,12,13,16:19,21,23,40,47];%ts79  24  outside of head
subj16 = [4:10,12:17,19:24,32,36]; %ny84: 25 is outside head, but good activity
subj17 = [8:10,12,14:23,25:30,33,35,40,41,42,57];%ft84
subj18 = [11,13,14,15,22,24,25,30];%gv84
subj19 = [6,7,9:11,13:18,20,22:25,27,28,31];%ka83
subj20 = [4,5,7:12,14:15,18:24,26:28,32,36,37,43,51];%cy82
subj21 = [5:7,9:12,14:20,23,24,26,29,30,33,35]; %jb84
subj22 = [4,11,12,14,16,19,20,21,26,32];%rd81
subj23 = [3,5:7,9,10,12:21,23,24,26:28,32,34,44];%km81_2
subj24 = [1,3:12,18:21,26,28,30,31,32];%jo74_2
subj25 = [4,6,7,10,11,12,13,18,19,40,52];%bt78
subj26 = [4,5,7,9,10,12:17,19,21]; % as78
gdcomps = {ap82,cj82,ds76,ec81,jo74,ke70,km81,mk79,nf68,tp62,ds80,kb77,subj13,subj14,subj15,subj16,subj17,subj18,subj19,subj20,subj21,subj22,subj23,subj24 subj25 subj26};

% list of each cluster with a,b,c,d, etc
ap82a= [3   6   9 8  15  24  25  35]; 
cj82a= []; 
ds76a= [11  12 18  20   26 30];  %      % OCCIPITAL
ec81a= [13  16  27]; 
jo74a= [4   7  11  14  21  37];  % 8 and 17 were way more like parietal, so removed
ke70a= []; 
km81a= [7   8  11  13]; 
mk79a= [19]; 
nf68a= [8  9  12]; 
tp62a= [8]; 
ds80a= [15 23]; 
kb77a= [ 3 7  11  14  31]; 
subj13a= [];    % cz84
subj14a= [];    % gm84
subj15a= [];    % ts79
subj16a= [12  21];    % ny84
subj17a= [8  10  12  22  25];    % ft84
subj18a= [14];    % gv84
subj19a= [];    % ka83
subj20a= [23 28];    % cy82
subj21a= [];    % jb84
subj22a= [26];    % rd81
subj23a= [];    % km81_2 
subj24a= [];%[8  20];    % jo74_2 
subj25a= [11  19];    % bt78 
subj26a= [12  13  15];    % as78   % 13 works for probe, not mem
gdcompsa = {ap82a,cj82a,ds76a,ec81a,jo74a,ke70a,km81a,mk79a,nf68a,tp62a,ds80a,kb77a,subj13a,subj14a,subj15a,subj16a,subj17a,subj18a,subj19a,subj20a,subj21a,subj22a,subj23a,subj24a subj25a subj26a};

ap82b= [5];       %  out of 18  components are bilateral
cj82b= [9]; 
ds76b= [7]; 
ec81b= [28];                  %   **P1 Occipital Cluster**
jo74b= [16]; 
ke70b= [6]; 
km81b= [30]; 
mk79b= [6]; 
nf68b= []; 
tp62b= [4]; 
ds80b= [2 11]; 
kb77b= []; 
subj13b= [9];    % cz84
subj14b= [14];    % gm84
subj15b= [10];    % ts79
subj16b= [4];    % ny84
subj17b= [];    % ft84
subj18b= [];    % gv84
subj19b= [];    % ka83
subj20b= [9  19];    % cy82
subj21b= [];    % jb84
subj22b= [20];    % rd81
subj23b= [];%[24];    % km81_2 
subj24b= [];%[7  12];    % jo74_2 
subj25b= [];    % bt78 
subj26b= [4];    % as78 
gdcompsb = {ap82b,cj82b,ds76b,ec81b,jo74b,ke70b,km81b,mk79b,nf68b,tp62b,ds80b,kb77b,subj13b,subj14b,subj15b,subj16b,subj17b,subj18b,subj19b,subj20b,subj21b,subj22b,subj23b,subj24b subj25b subj26b};



ap82c= [6  17  20]; 
cj82c= [4   7  11  12  16];    % sub 10 Hz cluster PARIETAL
ds76c= []; 
ec81c= []; 
jo74c= [8  17]; 
ke70c= [9]; 
km81c= [2   5   9  10  14  16  23  38  40]; 
mk79c= []; 
nf68c= [6  7  10  15 19  20]; 
tp62c= []; 
ds80c= []; 
kb77c= [2 8  10 16 19]; 
subj13c= [];    % cz84
subj14c= [3   9  18];    % gm84
subj15c= [];    % ts79
subj16c= [5 7  14  20  22];    % ny84
subj17c= [16  23  28  29 30];    % ft84
subj18c= [];    % gv84
subj19c= [9  18  24  28];    % ka83
subj20c= [7   8  15  21];    % cy82
subj21c= [5   6   7  9  10  14  15  16  17  29  30  35];    % jb84
subj22c= [];    % rd81
subj23c= [];%[5  13  15  16  19  26];    % km81_2 
subj24c= [];%[9  21];    % jo74_2 
subj25c= [];    % bt78 
subj26c= [13];    % as78   % 13 works for mem, not probe
gdcompsc = {ap82c,cj82c,ds76c,ec81c,jo74c,ke70c,km81c,mk79c,nf68c,tp62c,ds80c,kb77c,subj13c,subj14c,subj15c,subj16c,subj17c,subj18c,subj19c,subj20c,subj21c,subj22c,subj23c,subj24c subj25c subj26c};

ap82d= [29]; 
cj82d= [13  18]; 
ds76d= [10  13]; 
ec81d= [37]; 
jo74d= [12  23  25];                        % motor cluster
ke70d= [14]; 
km81d= [20]; 
mk79d= [19]; 
nf68d= [10]; 
tp62d= [15  16]; 
ds80d= []; 
kb77d= []; 
subj13d= [];    % cz84
subj14d= [16];    % gm84
subj15d= [17];    % ts79
subj16d= [23];    % ny84
subj17d= [21];    % ft84
subj18d= [22];    % gv84
subj19d= [];    % ka83
subj20d= [26];    % cy82
subj21d= [24];    % jb84
subj22d= [];    % rd81
subj23d= [];%[12];    % km81_2 
subj24d= [];%[18  26];    % jo74_2 
subj25d= [10];    % bt78 
subj26d= [14];    % as78 
gdcompsd = {ap82d,cj82d,ds76d,ec81d,jo74d,ke70d,km81d,mk79d,nf68d,tp62d,ds80d,kb77d,subj13d,subj14d,subj15d,subj16d,subj17d,subj18d,subj19d,subj20d,subj21d,subj22d,subj23d,subj24d subj25d subj26d};

ap82e= [10  18  12];     % updated march 10th
cj82e= [14  10]; 
ds76e= []; 
ec81e= []; 
jo74e= [3  26]; 
ke70e= [19]; 
km81e= [26 4]; 
mk79e= []; 
nf68e= [24]; 
tp62e= []; 
ds80e= [4 10  5]; 
kb77e= [12 20]; 
subj13e= [];    % cz84
subj14e= [6 7  5];    % gm84
subj15e= [];    % ts79
subj16e= [10 8];    % ny84
subj17e= [17];    % ft84
subj18e= [];    % gv84
subj19e= [];    % ka83
subj20e= [18  12];    % cy82
subj21e= [];    % jb84
subj22e= [16];    % rd81
subj23e= [];    % km81_2 
subj24e= [];%[5];    % jo74_2 
subj25e= [6];    % bt78 
subj26e= [9];    % as78 
gdcompse = {ap82e,cj82e,ds76e,ec81e,jo74e,ke70e,km81e,mk79e,nf68e,tp62e,ds80e,kb77e,subj13e,subj14e,subj15e,subj16e,subj17e,subj18e,subj19e,subj20e,subj21e,subj22e,subj23e,subj24e subj25e subj26e};

% All components anterior to central sulcus
ap82f= [10,12,18];  % all are rv <10% (except jo35->close) and none outside brain
cj82f= [10,14];       
ds76f= [8,16,31];
ec81f= [3,7,9,12,20];
jo74f= [2,3,26,28,35];
ke70f= [19];%42 (not computed)
km81f= [2,4,26,31];
mk79f= [5,8,14,18] ;
nf68f= [18 20 21 24];
tp62f= [3];
ds80f= [3:5,10,14,18,25];
kb77f= [4,12,20,33];
subj13f = [5,17,18,32];%cz84
subj14f = [5:8];%gm84
subj15f = [5,13,19,21,47];%ts79
subj16f = [6,8,10,13,19]; %ny84: 25 is outside head, but good activity
subj17f = [9,17,18,20,26,33,41,42,57];%ft84
subj18f = [13,25,30];%gv84
subj19f = [6,7,11,13];%ka83
subj20f = [4,5,12,18,43];%cy82
subj21f = [11,15,24,26,30,33]; %jb84
subj22f = [4,11,16];%rd81
subj23f = [6,7,14,23];%km81_2
subj24f = [3,5,19,28];%jo74_2
subj25f = [6];%bt78
subj26f = [5,9,10,16,17]; % as78
gdcompsf = {ap82f,cj82f,ds76f,ec81f,jo74f,ke70f,km81f,mk79f,nf68f,tp62f,ds80f,kb77f,subj13f,subj14f,subj15f,subj16f,subj17f,subj18f,subj19f,subj20f,subj21f,subj22f,subj23f,subj24f subj25f subj26f};

% only choose ONE of the following (both alter gdcomps %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds comps NOT clustered
for nx =1: length(gdcomps)
    gdcomps{nx}(ismember(gdcomps{nx},gdcompsa{nx}))=[];    
    gdcomps{nx}(ismember(gdcomps{nx},gdcompsb{nx}))=[];
    gdcomps{nx}(ismember(gdcomps{nx},gdcompsc{nx}))=[];
    gdcomps{nx}(ismember(gdcomps{nx},gdcompsd{nx}))=[];
    gdcomps{nx}(ismember(gdcomps{nx},gdcompse{nx}))=[];    
    gdcomps{nx}(ismember(gdcomps{nx},gdcompsf{nx}))=[];    
end;
% Finds comps  clustered
for nx =1: length(gdcompsa)
union1 = union(gdcompsa{nx},gdcompsb{nx});
union1 = union(union1,gdcompsc{nx});
union1 = union(union1,gdcompsd{nx});
union1 = union(union1,gdcompse{nx});
gdcomps{nx} = union1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% following is comps from ICA on Alpha Ignore comps %%%%%%%%%%%%%%%%%%
%EEG.sources.component  % to find --->
%%%% enter subject sources desired------------------------------------->
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a= [];b=[];c=[];d=[];e=[];f=[];g=[];h=[];k=[];l=[];m=[];n=[];o=[];p=[];q=[];r=[];s=[];t=[];u=[];v=[];ww=[];x=[];yy=[];z=[];aa=[] ;bb=[];   
    for ss = 1:length(gdcomps)
        EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
        cmpindx=[];
        if ~isempty(gdcomps(ss))
            if ~isempty(EEG.dipfit)
                if ss==25
                    aa= gdcomps{ss};                 
                end;
                if ss==26
                    bb= gdcomps{ss};                 
                end;
            else        
                y=[EEG.sources.component];
                for num=1:length(gdcomps{ss})  %finds indexs of interest
                    pp = find(y==gdcomps{ss}(num));
                    cmpindx(1,end+1:end+length(pp)) = pp;
                end;
                if ss==1
                    a= cmpindx;                 
                end;
                if ss==2
                    b= cmpindx;                 
                end;
                if ss==3
                    c= cmpindx;                 
                end;
                if ss==4
                    d= cmpindx;                 
                end;
                if ss==5
                    e= cmpindx;                 
                end;
                if ss==6
                    f= cmpindx;                 
                end;
                if ss==7
                    g= cmpindx;                 
                end;
                if ss==8
                    h= cmpindx;                 
                end;
                if ss==9
                    k= cmpindx;                 
                end;
                if ss==10
                    l= cmpindx;                 
                end;
                if ss==11
                    m= cmpindx;                 
                end;
                if ss==12
                    n= cmpindx;                 
                end;
                if ss==13
                    o= cmpindx;                 
                end;
                if ss==14
                    p= cmpindx;                 
                end;
                if ss==15
                    q= cmpindx;               
                end;
                if ss==16
                    r= cmpindx; 
                end;
                if ss==17
                    s= cmpindx;                 
                end;
                if ss==18
                    t= cmpindx;                 
                end;
                if ss==19
                    u= cmpindx;                 
                end;
                if ss==20
                    v= cmpindx;                 
                end;
                if ss==21
                    ww= cmpindx;                
                end;
                if ss==22
                    x= cmpindx;                 
                end;
                if ss==23
                    yy= cmpindx;                
                end;
                if ss==24
                    z= cmpindx;                 
                end;
            end;            
        end;
    end;

    %%%%%%%% LIST OF COMP INDEXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear tmparray1 tmparray2 tmparray3 tmparray4 tmparray5 tmparray6 tmparray7 tmparray8 tmparray9 tmparray10 tmparray11 tmparray12 tmparray13 tmparray14 tmparray15 tmparray16 tmparray17 tmparray18 tmparray19 tmparray20 tmparray21 tmparray22 tmparray23 tmparray24 tmparray25 tmparray26
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(a)
        tmparray1(1,w)=EEG.sources(1,a(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 2); CURRENTSET = 2;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(b)
        tmparray2(1,w)=EEG.sources(1,b(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 3); CURRENTSET = 3;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(c)
        tmparray3(1,w)=EEG.sources(1,c(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 4); CURRENTSET = 4;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(d)
        tmparray4(1,w)=EEG.sources(1,d(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 5); CURRENTSET = 5;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(e)
        tmparray5(1,w)=EEG.sources(1,e(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 6); CURRENTSET = 6;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(f)
        tmparray6(1,w)=EEG.sources(1,f(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 7); CURRENTSET = 7;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(g)
        tmparray7(1,w)=EEG.sources(1,g(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 8); CURRENTSET = 8;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(h)
        tmparray8(1,w)=EEG.sources(1,h(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 9); CURRENTSET = 9;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(k)
        tmparray9(1,w)=EEG.sources(1,k(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 10); CURRENTSET = 10;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(l)
        tmparray10(1,w)=EEG.sources(1,l(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 11); CURRENTSET = 11;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(m)
        tmparray11(1,w)=EEG.sources(1,m(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 12); CURRENTSET = 12;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(n)
        tmparray12(1,w)=EEG.sources(1,n(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 13); CURRENTSET = 13;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(o)
        tmparray13(1,w)=EEG.sources(1,o(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 14); CURRENTSET = 14;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(p)
        tmparray14(1,w)=EEG.sources(1,p(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 15); CURRENTSET = 15;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(q)
        tmparray15(1,w)=EEG.sources(1,q(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 16); CURRENTSET = 16;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(r)
        tmparray16(1,w)=EEG.sources(1,r(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 17); CURRENTSET = 17;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(s)
        tmparray17(1,w)=EEG.sources(1,s(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 18); CURRENTSET = 18;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(t)
        tmparray18(1,w)=EEG.sources(1,t(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 19); CURRENTSET = 19;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(u)
        tmparray19(1,w)=EEG.sources(1,u(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 20); CURRENTSET = 20;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(v)
        tmparray20(1,w)=EEG.sources(1,v(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 21); CURRENTSET = 21;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(ww)
        tmparray21(1,w)=EEG.sources(1,ww(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 22); CURRENTSET = 22;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(x)
        tmparray22(1,w)=EEG.sources(1,x(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 23); CURRENTSET = 23;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(yy)
        tmparray23(1,w)=EEG.sources(1,yy(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 24); CURRENTSET = 24;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(z)
        tmparray24(1,w)=EEG.sources(1,z(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 25); CURRENTSET = 25;
        tmparray25=EEG.dipfit.model(aa); 
    EEG = eeg_retrieve(ALLEEG, 26); CURRENTSET = 26;
        tmparray26= EEG.dipfit.model(bb);

    %%%%%%%%%%% Put all the tmparrays together %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% must comment out assignments that have empty comp indexes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear allbesa allbesa2
    if ~isempty(a)
        allbesa(1,1:length(a)) = tmparray1;
    end;
    if ~isempty(b)
        allbesa(1,length(a)+1:length(a)+length(b)) = tmparray2;
    end;
    if ~isempty(c)
        allbesa(1,length(a)+length(b)+1:length(a)+length(b)+length(c)) = tmparray3;
    end;
    if ~isempty(d)
        allbesa(1,length(a)+length(b)+length(c)+1:length(a)+length(b)+length(c)+length(d)) = tmparray4;%ec
    end;
    if ~isempty(e)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+1:length(a)+length(b)+length(c)+length(d)+length(e)) = tmparray5; %jo
    end;
    if ~isempty(f)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)) = tmparray6; %ke
    end;
    if ~isempty(g)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)) = tmparray7;%km
    end;
    if ~isempty(h)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)) = tmparray8;%mk
    end;
    if ~isempty(k)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)) = tmparray9;%nf
    end;
    if ~isempty(l)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)) = tmparray10;%tp
    end;
    if ~isempty(m)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)) = tmparray11;%ds8
    end;
    if ~isempty(n)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)) = tmparray12;%kb
    end;
    if ~isempty(o)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)) = tmparray13;
    end;
    if ~isempty(p)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)) = tmparray14;
    end;
    if ~isempty(q)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)) = tmparray15;
    end;
    if ~isempty(r)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)) = tmparray16;
    end;
    if ~isempty(s)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)) = tmparray17;
    end;
    if ~isempty(t)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)) = tmparray18;
    end;
    if ~isempty(u)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)) = tmparray19;
    end;
    if ~isempty(v)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)) = tmparray20;
    end;
    if ~isempty(ww)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)) = tmparray21;%kb
    end;
    if ~isempty(x)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)) = tmparray22;
    end;
    if ~isempty(yy)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)) = tmparray23;
    end;
    if ~isempty(z)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)+length(z)) = tmparray24;
    end;
    if ~isempty(aa)
        allbesa2(1,1:length(aa)) = tmparray25;
    end;
    if ~isempty(bb)
        allbesa2(1,length(aa)+1:length(aa)+length(bb)) = tmparray26;
    end;
    nonclust = allbesa;
    nonclust2 = allbesa2;
figure;  color = {[1 .88 0]};  pl = 4; % 1,4
subplot(2,3,pl)
dipplot(nonclust,'dipolesize',22,'dipolelength',0,'color',color,'image','mri','gui','off','view',[1 0 0]);hold on;
dipplot(nonclust2,'dipolesize',22,'dipolelength',0,'color',color,'image','mri','gui','off','view',[1 0 0]);camzoom(.5);
subplot(2,3,pl+1)
dipplot(nonclust,'dipolesize',22,'dipolelength',0,'color',color,'gui','off','image','mri');hold on;
dipplot(nonclust2,'dipolesize',22,'dipolelength',0,'color',color,'image','mri','gui','off');camzoom(.5);
subplot(2,3,pl+2)
dipplot(nonclust,'dipolesize',22,'dipolelength',0,'color',color,'image','mri','gui','off','view',[0 -1 0]);hold on;
dipplot(nonclust2,'dipolesize',22,'dipolelength',0,'color',color,'image','mri','gui','off','view', [0 -1 0]);camzoom(.5);

% get separate allbesas for clusters
nclust = 5;   clear avgmem avgig avgdiff avgmat oneplot oneplot2
lim =2;  % color limits on ERSPs
for clust = 1:nclust
        if clust==1     %-1 +4
        ap82= [5]; 
        cj82= [9]; 
        ds76= [7]; 
        ec81= [28];                  %   **P1 Occipital Cluster**
        jo74= [16]; 
        ke70= [6]; 
        km81= [30]; 
        mk79= [6]; 
        nf68= []; 
        tp62= [4]; 
        ds80= [2 11]; 
        kb77= []; 
        subj13= [9];    % cz84
        subj14= [14];    % gm84
        subj15= [10];    % ts79
        subj16= [4];    % ny84
        subj17= [];    % ft84
        subj18= [];    % gv84
        subj19= [];    % ka83
        subj20= [9  19];    % cy82
        subj21= [];    % jb84
        subj22= [20];    % rd81
        subj23= [];%[24];    % km81_2 
        subj24= [];%[7  12];    % jo74_2 
        subj25= [];    % bt78 
        subj26= [4];    % as78 
    end; 
    if clust==2     % -4 and -5
        ap82= [3   6   9 8  15  24  25  35]; 
        cj82= []; 
        ds76= [11  12 18  20   26 30];  %      % OCCIPITAL
        ec81= [13  16  27]; 
        jo74= [4   7 11  14  21  37]; 
        ke70= []; 
        km81= [7   8  11  13]; 
        mk79= [19]; 
        nf68= [8  9  12]; 
        tp62= [8]; 
        ds80= [15 23]; 
        kb77= [ 3 7  11  14  31]; 
        subj13= [];    % cz84
        subj14= [];    % gm84
        subj15= [];    % ts79
        subj16= [12  21];    % ny84
        subj17= [8  10  12  22  25];    % ft84
        subj18= [14];    % gv84
        subj19= [];    % ka83
        subj20= [23 28];    % cy82
        subj21= [];    % jb84
        subj22= [26];    % rd81
        subj23= [];    % km81_2 
        subj24= [];%[8  20];    % jo74_2 
        subj25= [11  19];    % bt78 
        subj26= [12  15];    % as78   % 13 works for probe, not mem
    end;
    if clust == 3   
        ap82= [6  17  20]; 
        cj82= [4   7  11  12  16];    % sub 10 Hz cluster PARIETAL
        ds76= []; 
        ec81= []; 
        jo74= [8  17]; 
        ke70= [9]; 
        km81= [2   5   9  10  14  16  23  38  40]; 
        mk79= []; 
        nf68= [6  7  10  15 19  20]; 
        tp62= []; 
        ds80= []; 
        kb77= [2 8  10 16 19]; 
        subj13= [];    % cz84
        subj14= [3   9  18];    % gm84
        subj15= [];    % ts79
        subj16= [5 7  14  20  22];    % ny84
        subj17= [16  23  28  29 30];    % ft84
        subj18= [];    % gv84
        subj19= [9  18  24  28];    % ka83
        subj20= [7   8  15  21];    % cy82
        subj21= [5   6   7  9  10  14  15  16  17  29  30  35];    % jb84
        subj22= [];    % rd81
        subj23= [];%[5  13  15  16  19  26];    % km81_2 
        subj24= [];%[9  21];    % jo74_2 
        subj25= [];    % bt78 
        subj26= [13];    % as78   % 13 works for mem, not probe
    end;
    if clust == 4  
        ap82= [29]; 
        cj82= [13  18]; 
        ds76= [10  13]; 
        ec81= [37]; 
        jo74= [12  23  25];                        % motor cluster
        ke70= [14]; 
        km81= [20]; 
        mk79= [19]; 
        nf68= [10]; 
        tp62= [15  16]; 
        ds80= []; 
        kb77= []; 
        subj13= [];    % cz84
        subj14= [16];    % gm84
        subj15= [17];    % ts79
        subj16= [23];    % ny84
        subj17= [21];    % ft84
        subj18= [22];    % gv84
        subj19= [];    % ka83
        subj20= [26];    % cy82
        subj21= [24];    % jb84
        subj22= [];    % rd81
        subj23= [];%[12];    % km81_2 
        subj24= [];%[18  26];    % jo74_2 
        subj25= [10];    % bt78 
        subj26= [14];    % as78 
    end;
    if clust==5
        ap82= [10  18  12]; 
        cj82= [14  10]; 
        ds76= [];                          % Frontal
        ec81= []; 
        jo74= [3  26]; 
        ke70= [19]; 
        km81= [26]; 
        mk79= []; 
        nf68= [18  24]; 
        tp62= []; 
        ds80= [4 10  5  3]; 
        kb77= [12]; 
        subj13= [];    % cz84
        subj14= [7  5];    % gm84
        subj15= [19];    % ts79
        subj16= [10 8];    % ny84
        subj17= [];    % ft84
        subj18= [];    % gv84
        subj19= [11];    % ka83
        subj20= [18  12   5];    % cy82
        subj21= [];    % jb84
        subj22= [16];    % rd81
        subj23= [];    % km81_2 
        subj24= [];%[5];    % jo74_2 
        subj25= [6];    % bt78 
        subj26= [9 16];    % as78 
    end;
    allcomps = cat(2,ap82,cj82);
    allcomps = cat(2,allcomps,ds76);
    allcomps = cat(2,allcomps,ec81);
    allcomps = cat(2,allcomps,jo74);
    allcomps = cat(2,allcomps,ke70);
    allcomps = cat(2,allcomps,km81);
    allcomps = cat(2,allcomps,mk79);
    allcomps = cat(2,allcomps,nf68);
    allcomps = cat(2,allcomps,tp62);
    allcomps = cat(2,allcomps,ds80);
    allcomps = cat(2,allcomps,kb77);
    allcomps = cat(2,allcomps,subj13);
    allcomps = cat(2,allcomps,subj14);
    allcomps = cat(2,allcomps,subj15);
    allcomps = cat(2,allcomps,subj16);
    allcomps = cat(2,allcomps,subj17);
    allcomps = cat(2,allcomps,subj18);
    allcomps = cat(2,allcomps,subj19);
    allcomps = cat(2,allcomps,subj20);
    allcomps = cat(2,allcomps,subj21);
    allcomps = cat(2,allcomps,subj22);
    allcomps = cat(2,allcomps,subj23);
    allcomps = cat(2,allcomps,subj24);
    allcomps = cat(2,allcomps,subj25);
    allcomps = cat(2,allcomps,subj26);
    gdcomps = {ap82,cj82,ds76,ec81,jo74,ke70,km81,mk79,nf68,tp62,ds80,kb77,subj13,subj14,subj15,subj16,subj17,subj18,subj19,subj20,subj21,subj22,subj23,subj24,subj25,subj26};
    
    a= [];b=[];c=[];d=[];e=[];f=[];g=[];h=[];k=[];l=[];m=[];n=[];o=[];p=[];q=[];r=[];s=[];t=[];u=[];v=[];ww=[];x=[];yy=[];z=[];  aa=[];bb=[];
    for ss = 1:length(gdcomps)
        EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
        cmpindx=[];
        if ~isempty(gdcomps(ss))
            if ~isempty(EEG.dipfit)
                if ss==25
                    aa= gdcomps{ss};                 
                end;
                if ss==26
                    bb= gdcomps{ss};                 
                end;
            else        
                y=[EEG.sources.component];
                for num=1:length(gdcomps{ss})  %finds indexs of interest
                    pp = find(y==gdcomps{ss}(num));
                    cmpindx(1,end+1:end+length(pp)) = pp;
                end;
                if ss==1
                    a= cmpindx; 
                end;
                if ss==2
                    b= cmpindx; 
                end;
                if ss==3
                    c= cmpindx; 
                end;
                if ss==4
                    d= cmpindx; 
                end;
                if ss==5
                    e= cmpindx; 
                end;
                if ss==6
                    f= cmpindx; 
                end;
                if ss==7
                    g= cmpindx; 
                end;
                if ss==8
                    h= cmpindx; 
                end;
                if ss==9
                    k= cmpindx; 
                end;
                if ss==10
                    l= cmpindx; 
                end;
                if ss==11
                    m= cmpindx; 
                end;
                if ss==12
                    n= cmpindx; 
                end;
                if ss==13
                    o= cmpindx; 
                end;
                if ss==14
                    p= cmpindx; 
                end;
                if ss==15
                    q= cmpindx; 
                end;
                if ss==16
                    r= cmpindx; 
                end;
                if ss==17
                    s= cmpindx; 
                end;
                if ss==18
                    t= cmpindx; 
                end;
                if ss==19
                    u= cmpindx; 
                end;
                if ss==20
                    v= cmpindx; 
                end;
                if ss==21
                    ww= cmpindx; 
                end;
                if ss==22
                    x= cmpindx; 
                end;
                if ss==23
                    yy= cmpindx; 
                end;
                if ss==24
                    z= cmpindx; 
                end;
            end;            
        end;
    end;
    %%%%%%%% LIST OF COMP INDEXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear tmparray1 tmparray2 tmparray3 tmparray4 tmparray5 tmparray6 tmparray7 tmparray8 tmparray9 tmparray10 tmparray11 tmparray12 tmparray13 tmparray14 tmparray15 tmparray16 tmparray17 tmparray18 tmparray19 tmparray20 tmparray21 tmparray22 tmparray23 tmparray24 tmparray25 tmparray26
    EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(a)
        tmparray1(1,w)=EEG.sources(1,a(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 2); CURRENTSET = 2;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(b)
        tmparray2(1,w)=EEG.sources(1,b(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 3); CURRENTSET = 3;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(c)
        tmparray3(1,w)=EEG.sources(1,c(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 4); CURRENTSET = 4;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(d)
        tmparray4(1,w)=EEG.sources(1,d(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 5); CURRENTSET = 5;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(e)
        tmparray5(1,w)=EEG.sources(1,e(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 6); CURRENTSET = 6;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(f)
        tmparray6(1,w)=EEG.sources(1,f(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 7); CURRENTSET = 7;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(g)
        tmparray7(1,w)=EEG.sources(1,g(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 8); CURRENTSET = 8;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(h)
        tmparray8(1,w)=EEG.sources(1,h(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 9); CURRENTSET = 9;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(k)
        tmparray9(1,w)=EEG.sources(1,k(w));
    end;
    EEG = eeg_retrieve(ALLEEG, 10); CURRENTSET = 10;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(l)
        tmparray10(1,w)=EEG.sources(1,l(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 11); CURRENTSET = 11;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(m)
        tmparray11(1,w)=EEG.sources(1,m(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 12); CURRENTSET = 12;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(n)
        tmparray12(1,w)=EEG.sources(1,n(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 13); CURRENTSET = 13;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(o)
        tmparray13(1,w)=EEG.sources(1,o(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 14); CURRENTSET = 14;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(p)
        tmparray14(1,w)=EEG.sources(1,p(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 15); CURRENTSET = 15;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(q)
        tmparray15(1,w)=EEG.sources(1,q(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 16); CURRENTSET = 16;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(r)
        tmparray16(1,w)=EEG.sources(1,r(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 17); CURRENTSET = 17;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(s)
        tmparray17(1,w)=EEG.sources(1,s(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 18); CURRENTSET = 18;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(t)
        tmparray18(1,w)=EEG.sources(1,t(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 19); CURRENTSET = 19;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(u)
        tmparray19(1,w)=EEG.sources(1,u(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 20); CURRENTSET = 20;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(v)
        tmparray20(1,w)=EEG.sources(1,v(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 21); CURRENTSET = 21;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(ww)
        tmparray21(1,w)=EEG.sources(1,ww(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 22); CURRENTSET = 22;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(x)
        tmparray22(1,w)=EEG.sources(1,x(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 23); CURRENTSET = 23;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(yy)
        tmparray23(1,w)=EEG.sources(1,yy(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 24); CURRENTSET = 24;
    %EEG.sources = rmfield(EEG.sources, 'besaextori');
    for w=1:length(z)
        tmparray24(1,w)=EEG.sources(1,z(w)); 
    end;
    EEG = eeg_retrieve(ALLEEG, 25); CURRENTSET = 25;
    tmparray25=EEG.dipfit.model(aa); 
    EEG = eeg_retrieve(ALLEEG, 26); CURRENTSET = 26;
    tmparray26= EEG.dipfit.model(bb);

    %%%%%%%%%%% Put all the tmparrays together %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% must comment out assignments that have empty comp indexes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear allbesa allbesa2
    if ~isempty(a)
        allbesa(1,1:length(a)) = tmparray1;
    end;
    if ~isempty(b)
        allbesa(1,length(a)+1:length(a)+length(b)) = tmparray2;
    end;
    if ~isempty(c)
        allbesa(1,length(a)+length(b)+1:length(a)+length(b)+length(c)) = tmparray3;
    end;
    if ~isempty(d)
        allbesa(1,length(a)+length(b)+length(c)+1:length(a)+length(b)+length(c)+length(d)) = tmparray4;%ec
    end;
    if ~isempty(e)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+1:length(a)+length(b)+length(c)+length(d)+length(e)) = tmparray5; %jo
    end;
    if ~isempty(f)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)) = tmparray6; %ke
    end;
    if ~isempty(g)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)) = tmparray7;%km
    end;
    if ~isempty(h)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)) = tmparray8;%mk
    end;
    if ~isempty(k)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)) = tmparray9;%nf
    end;
    if ~isempty(l)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)) = tmparray10;%tp
    end;
    if ~isempty(m)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)) = tmparray11;%ds8
    end;
    if ~isempty(n)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)) = tmparray12;%kb
    end;
    if ~isempty(o)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)) = tmparray13;
    end;
    if ~isempty(p)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)) = tmparray14;
    end;
    if ~isempty(q)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)) = tmparray15;
    end;
    if ~isempty(r)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)) = tmparray16;
    end;
    if ~isempty(s)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)) = tmparray17;
    end;
    if ~isempty(t)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)) = tmparray18;
    end;
    if ~isempty(u)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)) = tmparray19;
    end;
    if ~isempty(v)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)) = tmparray20;
    end;
    if ~isempty(ww)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)) = tmparray21;%kb
    end;
    if ~isempty(x)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)) = tmparray22;
    end;
    if ~isempty(yy)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)) = tmparray23;
    end;
    if ~isempty(z)
        allbesa(1,length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)+1:length(a)+length(b)+length(c)+length(d)+length(e)+length(f)+length(g)+length(h)+length(k)+length(l)+length(m)+length(n)+length(o)+length(p)+length(q)+length(r)+length(s)+length(t)+length(u)+length(v)+length(ww)+length(x)+length(yy)+length(z)) = tmparray24;
    end;
    if ~isempty(aa)
        allbesa2(1,1:length(aa)) = tmparray25;
    end;
    if ~isempty(bb)
        allbesa2(1,length(aa)+1:length(aa)+length(bb)) = tmparray26;
    end;
    oneplot{1,clust}  = allbesa;
    if exist('allbesa2')
        oneplot2{1,clust} = allbesa2;
    end;
end;    
    figure; pl = 1;        % Plot Ignore-Mem = Diff Plots  %%%%%%%%%%%%%
color = {{[0 .5 1]}  {'r'} {'g'} {[1 .56 0]} {[1 0 .8]} };% 1=light blue; orange,pink
for p=1:5
    subplot(2,3,pl)
    [sources X Y Z XE YE ZE] = dipplot(oneplot{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[1 0 0],'color',color{p},'image','mri');hold on;
%    if ~isempty(oneplot2{p})
    [sources X Y Z XE YE ZE] = dipplot(oneplot2{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[1 0 0],'color',color{p},'image','mri');camzoom(.6);
%    end;
    subplot(2,3,pl+1)
    [sources X Y Z XE YE ZE] = dipplot(oneplot{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[0 0 1],'color',color{p},'image','mri');hold on;
%    if ~isempty(oneplot2{p})
    [sources X Y Z XE YE ZE] = dipplot(oneplot2{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[0 0 1],'color',color{p},'image','mri');camzoom(.6);
%    end;
    subplot(2,3,pl+2)
    [sources X Y Z XE YE ZE] = dipplot(oneplot{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[0 -1 0],'color',color{p},'image','mri');hold on;
%    if ~isempty(oneplot2{p})
    [sources X Y Z XE YE ZE] = dipplot(oneplot2{p},'dipolesize',22,'dipolelength',0,'gui','off','view',[0 -1 0],'color',color{p},'image','mri');camzoom(.6);
%    end;    
end;
