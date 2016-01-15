% clusters the time-irrelevant cross coherence during different emotions. Takes the mean coh over 'time'

subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];%  jo74
subj2 = [1,3:23,25:30,35];  % tl81  23 and 35 are pulse
subj3 = [1,2,4:9,12:24,26,30,31,33,34,46,51];   % mi83   14, 23 are pulse
subj4 = [10,12,13,15,17,18,21,22,23,28,33,36];   % ms82
gdcomps = {subj1, subj2, subj3, subj4};

paths = {'/data/common1/emotion/jo74/','/data/common1/emotion/tl81/','/data/common1/emotion/mi83/','/data/common1/emotion/ms82/'};
datpaths = {'/data/common1/emotion/jo74/crosscomps/','/data/common1/emotion/tl81/crosscomps/','/data/common1/emotion/mi83/crosscomps/','/data/common1/emotion/ms82/crosscomps/'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
emlen = length(emos);
datsets = {'prebase.mat','awe.mat','frustration.mat','joy.mat','anger.mat','happy.mat','sad.mat','love.mat','fear.mat' ,'compassion.mat','jealousy.mat','content.mat','grief.mat','relief.mat','disgust.mat','excite.mat','postbase.mat'};

nx = 2;cd (datpaths{nx})
load  prebase.mat
fr = find (freqsout<50);
tm = find (timesout>-1500 & timesout<1500);
allcohs = zeros(length(fr),0);
for nx = 2:length(gdcomps)
    cd (datpaths{nx})
    for e = 1:length(emos)
        load(datsets{e})
        for w=1:length(gdcomps{nx})-1
            a = gdcomps{nx}(w);
            b = gdcomps{nx}(w+1:end);
            twocoh = onecoh{a};
            bootmat1 = oneboot{a};   %chooses a 3D array from cell array
            minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
            minboot = repmat (minboot,[1,size(timesout,2),1]); %makes (63,200,71) matrix of minboot
            twocoh(find(twocoh <= minboot)) = 0;%finds and zeros cohmats less than minmask
            twocoh = twocoh(fr,tm,b);
            twocoh = mean(twocoh,2);   % means over times; makes freqs X 1 X comps
            twocoh = squeeze(twocoh);  % makes a freqs X comps
            for mn = 1:size(twocoh,2)
                twocoh(:,mn) = twocoh(:,mn)-mean(twocoh(:,mn));
            end;            
            allcohs(:,end+1:end+size(twocoh,2)) = twocoh;
        end;        
    fprintf('\nOne more Emotion done: %s',emos{e});
    end;
    fprintf('\nOne more subj done: %s',int2str(nx));
    clear cohercell crsbootcell oneemo oneboot
end;


 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(allcohs,'pca',10,'stop',1e-7,'maxsteps',3000); 
 % puts in freqs X emo*comp pairs
 
numcomp =size(weights,1);
winv = pinv(weights*sphere);
%%%%%%%%%%%%%%%%%%%% Plot distributions of all ICA comps %%%%%%%%%%%%%%%%%
figure;pl=1;
for comp = 1:numcomp
    subplot(ceil(sqrt(numcomp)),ceil(sqrt(numcomp)),pl)
    plot(freqsout(fr),winv(:,comp));pl = pl+1;
    set(gca,'xtick',[5:5:50]);
    set(gca,'xgrid','on');
    set(gca,'xlim',[3 50]);
end;



% Find > zero activations for clusters 1:12
poscut = 5;
for n = 1:size(activations,1)
    figure;plot(activations(n,:)); title(int2str(n));
    tmpact = activations(n,:)/std(activations(n,:));
    hiact{n} = find(tmpact > poscut);
end;
% find sum of positive activations for all components
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};
subjname = {'subj1','subj2','subj3','subj4'};
printname = {'jo74','tl81','mi83','ms82'};
fid = fopen('/data/common1/emotion/CrossCohICAClusts','a');%
for clust = 1:size(activations,1)
    numtot = 0;
        fprintf(fid,'\n\nICA CLUSTER: %s POSITIVE   0-50Hz; mean over time input to ICA; PCA to 10; Activation cut off is std(activations(clust,:))',int2str(clust));
        fprintf(fid,'\nActivation Cutoff: %g ',poscut);
        %NEGATIVE POSITIVE
    sc = hiact{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
        if iseven(length(gdcomps{nx}))
            ntrials = length(gdcomps{nx})*floor((length(gdcomps{nx})-1)/2);
        else
            ntrials = length(gdcomps{nx})*(length(gdcomps{nx})-1)/2;
        end;       
        for em = 1:length(emos)
            ft = sc(find(sc>=numtot+1 & sc<=numtot+ntrials));
            gdcompsa =[]; gdcompsb=[];
            for rr = 1:length(ft)
                prevcp = 0;                
                for cmp = 1:length(gdcomps{nx})-1
                    if ft(rr)>= numtot + prevcp+1 & ft(rr) <= numtot + prevcp +length(gdcomps{nx})-cmp
                        gdcompsa(rr) = gdcomps{nx}(cmp);
                        gdcompsb(rr) = gdcomps{nx}(ft(rr)-(numtot+prevcp)+cmp);
                    end;                    
                    prevcp = prevcp + length(gdcomps{nx})-cmp;
                end;
            end;
            numtot = ntrials+numtot; 
            gdemocompsa{em} = gdcompsa; 
            gdemocompsb{em} = gdcompsb;
        end;
        gdsubjcompsa{nx} = gdemocompsa;
        gdsubjcompsb{nx} = gdemocompsb;
        end;
    end;
    for pemo = 1:length(emos)
        fprintf(fid,'\n%s  Emotion: %s ','%    ',emos{pemo});
        for psub = 1:length(gdcomps)
            fprintf(fid,'\n%sa= [%s]; %s   %s',subjname{psub},int2str(gdsubjcompsa{psub}{pemo}),' % ',printname{psub});
            fprintf(fid,'\n%sb= [%s]; %s   %s',subjname{psub},int2str(gdsubjcompsb{psub}{pemo}),' % ',printname{psub});
        end;
            fprintf(fid,'\ngdcomps%sa = {subj1a,subj2a,subj3a,subj4a};',emos{pemo});       
            fprintf(fid,'\ngdcomps%sb = {subj1b,subj2b,subj3b,subj4b};',emos{pemo});       
    end;    
end;
fclose(fid);

for x=1:17
gdemocomps{x}{1}
end;
% make component lists for each subj with sum above certain cutoff
printname = {'ap82','cj82','ds76','ec81','jo74','ke70','km81','mk79','nf68','tp62','ds80','kb77','subj13','subj14','subj15','subj16','subj17','subj18','subj19','subj20','subj21','subj22','subj23','subj24','subj25','subj26'};
fid = fopen('/data/common1/stern/eeg/Sternall/SnglTrialClusters-MLonly','a');%
multfact = 1.5;
for clust = 1:toclust
        fprintf(fid,'\nICA CLUSTER: %s POSITIVE   0-1400ms; 0-20Hz MEMORIZE single trial input;ALL 25 SUBJ (not gv). (April 30); activations > 0; PCA to 12; score based on SUM of > zero weightings (activations); MIDLINE FRONTAL COMPS ONLY',int2str(clust));
        fprintf(fid,'\nCutoff SCORE: %g ',thresh);
        %NEGATIVE POSITIVE
    sumcell =zeros(1,0);
        for nx = 1:length(gdcomps)
            for em = 1:length(emos)
                sumcell(1,end+1:end+size(allclustsums{clust}{nx}{em},2)) = allclustsums{clust}{nx}{em};
            end;
        end;
        thresh =  mean(sumcell)+multfact*(std(sumcell));   
        for nx = 1:length(gdcomps)
            if ~isempty(gdcomps{nx})
                sums = allclustsums{clust}{nx};
                abth = find(sums > thresh);
                subcomps = gdcomps{nx}(abth);
                clustcomps{nx} = subcomps;
                allvals{nx} = sums(abth);
                allsubjperc{clust} = allvals;
                allclustcomps{clust} = clustcomps;    
                fprintf(fid,'\n%s= [%s]; %s Sums:  %s',printname{nx},int2str(clustcomps{nx}),' % ',int2str(allvals{nx}));
            end;
        end;
    end;
fclose(fid);

poscut = 1;
negcut = -1;
for comp_number = 1:gdcomp;
    z= winv(:,comp_number)/std(winv(:,comp_number));
    [val idx] = sort (z(1,:));                     
    allz = z(find(z>poscut));
    if pos == 1
        sigindx = find(val> poscut);  % POS weights .05
    else
        sigindx = find(val<= negcut);  % NEG weights
    end;
    sigcomps =  idx(sigindx); 
    contrib_comps=zeros(1,length(sigcomps));
    for y = 1:length(sigcomps)
        contrib_comps(1,y) = sigcomps(length(sigcomps)+1 -y);
    end;
    fid = fopen('/data/common1/emotion/CrossCohClusts','a');%AlphaPCA
    fprintf(fid,'\nCross Coh ICA Cluster: %s All crosses: 0-50 Hz; mean over random time epochs ',int2str(comp_number));
    if pos == 1
        fprintf(fid,'\n norm score > %s',num2str(poscut));
    else
        fprintf(fid,'\n norm score < %s',num2str(negcut));
    end;
    
    s1idx = contrib_comps<=s1num;
    s2idx= contrib_comps>=s1num+1  &contrib_comps<=s1num+s2num;
    s3idx= contrib_comps>=s1num+s2num+1  &contrib_comps<=s1num+s2num+s3num;
    s4idx= contrib_comps>= s1num+s2num+s3num+1 &contrib_comps<=s1num+s2num+s3num+s4num;
    s5idx= contrib_comps>=s1num+s2num+s3num+s4num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num;
    %s6idx= contrib_comps>= s1num+s2num+s3num+s4num+s5num+1 &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num;
    %s7idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num;
    %s8idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+s7num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num;
    %s9idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num;
    %s10idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num;
    %s11idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num+s11num;
    %s12idx= contrib_comps>=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num+s11num+1  &contrib_comps<=s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num+s11num+s12num;
    sigcomps=[]; sigindx=[]; s1comps=[]; s2comps=[]; s3comps=[]; s4comps=[]; s5comps=[]; s6comps=[]; s7comps=[]; s8comps=[]; s9comps=[]; s10comps=[];s11comps=[];s12comps=[];
    s1comps = contrib_comps(s1idx);
    for u = 1:length(contrib_comps(s2idx))
        y=contrib_comps(s2idx);
        y=y(u);
        s2comps(1,u) = y - s1num;
    end;
    for u = 1:length(contrib_comps(s3idx))
        y=contrib_comps(s3idx);
        y=y(u);
        s3comps(1,u) = y - (s1num+s2num);
    end;
    for u = 1:length(contrib_comps(s4idx))
        y=contrib_comps(s4idx);
        y=y(u);
        s4comps(1,u) = y - (s1num+s2num+s3num);
    end;
    for u = 1:length(contrib_comps(s5idx))
    y=contrib_comps(s5idx);
    y=y(u);
    s5comps(1,u) = y - (s1num+s2num+s3num+s4num);
    end;
    %for u = 1:length(contrib_comps(s6idx))
    %y=contrib_comps(s6idx);
    %y=y(u);
    %s6comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num);
    %end;
    %for u = 1:length(contrib_comps(s7idx))
    %y=contrib_comps(s7idx);
    %y=y(u);
    %s7comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num+s6num);
    %end;
    %for u = 1:length(contrib_comps(s8idx))
    %y=contrib_comps(s8idx);
    %y=y(u);
    %s8comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num+s6num+s7num);
    %end;
    %for u = 1:length(contrib_comps(s9idx))
    %y=contrib_comps(s9idx);
    %y=y(u);
    %s9comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num);
    %end;
    %for u = 1:length(contrib_comps(s10idx))
    %y=contrib_comps(s10idx);
    %y=y(u);
    %s10comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num);
    %end;
    %for u = 1:length(contrib_comps(s11idx))
    %y=contrib_comps(s11idx);
    %y=y(u);
    %s11comps(1,u) = y - (s1num+s2num+s3num+s4num+s5num+s6num+s7num+s8num+s9num+s10num);
    %end;
    %for u = 1:length(contrib_comps(s12idx))
    %y=contrib_comps(s12idx);
    %y=y(u);
    %s12comps(1,u) = y - (s1num+s2num+dsnum+s4num+s5num+s6num+s7num+s8num+s9num+s10num+s11num);
    %end;
    %%%%%%%%%%%%%% To find which crosses these numbers refer to %%%%%%%%%%%%%%%%%%%%%%%%%%
    % load new subject assignments-crs from one comp generated from script above
    %% gives a comp1 X comp2 X subject 3D matrix of results------------------------------
    subj = {s1comps s2comps s3comps s4comps};% s5comps s6comps s7comps s8comps s9comps s10comps s11comps s12comps
    orig = {subj1 subj2 subj3 subj4};
    clear crosses
    for sc = 1:length(subj)
        nums = subj{sc};
        org = orig{sc};
        m = length(org);
        if ~isempty(nums)
            for n=1:length(nums)
                if nums(n) <= m-1
                    crosses{sc}(n,1) = org(1);
                    crosses{sc}(n,2) = org(nums(n)+1);
                end;
                if nums(n) >= m &  nums(n) <= (m-1)+(m-2)
                    crosses{sc}(n,1) = org(2);
                    crosses{sc}(n,2) = org((nums(n)-(m-1))+2);
                end;
                if nums(n) >= (m-1)+(m-2)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)
                    crosses{sc}(n,1) = org(3);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)))+3);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)
                    crosses{sc}(n,1) = org(4);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)))+4);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)
                    crosses{sc}(n,1) = org(5);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)))+5);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)
                    crosses{sc}(n,1) = org(6);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)))+6);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)
                    crosses{sc}(n,1) = org(7);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)))+7);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)
                    crosses{sc}(n,1) = org(8);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)))+8);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)
                    crosses{sc}(n,1) = org(9);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)))+9);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)
                    crosses{sc}(n,1) = org(10);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)))+10);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)
                    crosses{sc}(n,1) = org(11);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)))+11);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)
                    crosses{sc}(n,1) = org(12);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)))+12);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13);
                    crosses{sc}(n,1) = org(13);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)))+13);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)
                    crosses{sc}(n,1) = org(14);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)))+14);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)
                    crosses{sc}(n,1) = org(15);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)))+15);
                end;
                if nums(n) >=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)
                    crosses{sc}(n,1) = org(16);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)))+16);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)
                    crosses{sc}(n,1) = org(17);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)))+17);
                end;
                if nums(n) >=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)
                    crosses{sc}(n,1) = org(18);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)))+18);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)
                    crosses{sc}(n,1) = org(19);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)))+19);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+1 &  nums(n) <= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)
                    crosses{sc}(n,1) = org(20);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)))+20);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)
                    crosses{sc}(n,1) = org(21);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)))+21);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)
                    crosses{sc}(n,1) = org(22);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)))+22);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)
                    crosses{sc}(n,1) = org(23);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)))+23);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)
                    crosses{sc}(n,1) = org(24);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)))+24);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)+1 &  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)+(m-25)
                    crosses{sc}(n,1) = org(25);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)))+25);
                end;
                if nums(n) >= (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)+(m-25)+1&  nums(n) <=  (m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)+(m-25)+(m-26)
                    crosses{sc}(n,1) = org(26);
                    crosses{sc}(n,2) = org((nums(n)-((m-1)+(m-2)+(m-3)+(m-4)+(m-5)+(m-6)+(m-7)+(m-8)+(m-9)+(m-10)+(m-11)+(m-12)+(m-13)+(m-14)+(m-15)+(m-16)+(m-17)+(m-18)+(m-19)+(m-20)+(m-21)+(m-22)+(m-23)+(m-24)+(m-25)))+26);
                end;
            end; % for n (length(nums))
        end; % for if
    end; % for sc
    fprintf(fid,'\nScores: %s ',num2str(allz));
    for nx = 1:length(gdcomps)
        if isempty(crosses{nx})
        fprintf(fid,'\n%sa= []; ',pname{nx});
        fprintf(fid,'\n%sb= []; ',pname{nx});
           else 
        fprintf(fid,'\n%sa= [%s]; ',pname{nx},int2str(crosses{nx}(:,1)'));
        fprintf(fid,'\n%sb= [%s]; ',pname{nx},int2str(crosses{nx}(:,2)'));
        end;
    end;
    fclose(fid);
end;

