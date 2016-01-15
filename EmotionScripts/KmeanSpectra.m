% clusters smoothed, newtimef spectral output using k means across emotions
%***  Need to run ClustSnglTrPower.m first to get ContDataERSPs.mat
% early button press sessions (not instructed to only press on surge)
 

subjspecs = {'Subj1trials.fdt','Subj2trials.fdt','Subj3trials.fdt','Subj4trials.fdt','Subj5trials.fdt','Subj6trials.fdt','Subj7trials.fdt','Subj8trials.fdt','Subj9trials.fdt','Subj10trials.fdt','Subj11trials.fdt','Subj12trials.fdt','Subj13trials.fdt','Subj14trials.fdt','Subj15trials.fdt','Subj16trials.fdt','Subj17trials.fdt','Subj18trials.fdt','Subj19trials.fdt','Subj20trials.fdt','Subj21trials.fdt','Subj22trials.fdt','Subj23trials.fdt'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};

%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqs = [1:.5:50];  reps = 5; % number of kmeans restarts

%for nx = 1:length(gdcomps)
nx=1;clear allsums kout
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    for cmp = 1:length(gdcomps{nx})
        kmat = zeros(0,length(freqs)); 
        for k = 1:length(Alllongersps)
            oneemo = Alllongersps{k}; % freqsXtimesXcomponent
            oneemo = mean(oneemo,2);oneemo = squeeze(oneemo);outdata = (oneemo(:,gdcomps{nx}(cmp)));
            %[outdata,outx] = movav(oneemo(:,:,gdcomps{nx}(cmp)),[1:size(oneemo,2)],20);
            kmat(end+1:end+size(outdata,2),:) = outdata';
            %subjtrials(k) = size(outdata,2);
        end;
        [optk(cmp),centr,clst,Cg] =Kmeangrp(kmat,5,reps,0); % find optimal cluster number, don't plot
        [kout(:,cmp), C,sumd, allsums{cmp}] = kmeans(kmat,optk(cmp),'replicates',reps);
        fprintf('\n One More Component Done: %i of %i',cmp,length(gdcomps{nx}));
        kmatall(:,:,cmp) = kmat;
    end;
    %allkouts{nx} = kout;
    %optks{nx} = optk;
    %numtrials{nx} = subjtrials;  clear Alllongersps  totrials icamat  onecmp oneemo kout centr clst Cg optk
    fprintf('\n One More SUBJECT Done: %i\n',nx);
%end;
comment = 'Ran Kmeangrp.m which finds optimal number of kmeans clusters to use. Then ran kmeans with that value. Shold include:kout optk subjtrials comment ';
save KmeansClusterInfo.mat kout optk subjtrials comment kmatall

%%%%%%%%
% find emotion assignments of clusters  USING AVGERAGED TRIALS
cd (['/data/common2/emotion',paths{nx},'/ersps/']);  
cols(1,:) = [1,0,0];cols(2,:) = [0,0,1];cols(3,:) = [0,1,0];cols(4,:) = [.5,.5,.5];cols(5,:) = [.5,.5,.5];
clustcomps = find(optk<5);  thresh = 800;
row = ceil(sqrt(length(clustcomps))); col = floor(sqrt(length(clustcomps)));
figure;pl = 1; thresh = 40; 
for cp = 1:length(clustcomps)
    cmp = clustcomps(cp);
    subplot(row,col,pl);
    
    for e = 1:size(kout,1)
        ph = plot(freqs,kmatall(e,:,cmp),'k-'); hold on;
        if allsums{cmp}(e,kout(e,cmp)) < thresh
            set(ph,'color',cols(kout(e,cmp),:));
        else        
            set(ph,'color',cols(4,:));
        end;     
    end;pl = pl+1; 
    title(int2str(gdcomps{nx}(cmp))); set(gca,'box','off');set(gca,'fontsize',6);
    yl = get(gca,'ylim'); 
    abthresh = find(allsums{cmp}(:,1)>thresh);
    clustno = find(kout(:,cmp)==1);cut = intersect(abthresh,clustno);
    for q = 1:length(cut)
        clustno(clustno == cut(q)) = [];
    end;    
    emo1 = emos(clustno);
    ph=text(1,yl(1)+5,emo1');set(ph,'fontsize',6);
    abthresh = find(allsums{cmp}(:,2)>thresh);
    clustno = find(kout(:,cmp)==2);cut = intersect(abthresh,clustno);
    for q = 1:length(cut)
        clustno(clustno == cut(q)) = [];
    end;    
    emo2 = emos(clustno);
    ph=text(22,yl(1)+5,emo2');set(ph,'fontsize',6);
    if optk(cmp) > 2
    abthresh = find(allsums{cmp}(:,3)>thresh);
    clustno = find(kout(:,cmp)==3);cut = intersect(abthresh,clustno);
    for q = 1:length(cut)
        clustno(clustno == cut(q)) = [];
    end;    
    emo3 = emos(clustno);
    ph=text(42,yl(1)+5,emo3');set(ph,'fontsize',6);    
    end;        
end;
 


%%%%%%%%

%for nx = 1:length(gdcomps)
freqs = [1:.5:50];  reps = 2; % number of kmeans restarts
nx = 5 ;    clear allsums kout
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    for cmp = 1:length(gdcomps{nx})
        kmat = zeros(0,length(freqs)); 
        for k = 1:length(Alllongersps)
            oneemo = Alllongersps{k}; % freqsXtimesXcomponent
            [outdata,outx] = movav(oneemo(:,:,gdcomps{nx}(cmp)),[1:size(oneemo,2)],20);
            trialmat = zeros(0,size(outdata,1));
            for trl = 1:5:size(outdata,2)-11
                onetrl = mean(outdata(:,trl:trl+10),2);
                trialmat(end+1,:) = onetrl';
            end;
            kmat(end+1:end+size(trialmat,1),:) = trialmat;
            subjtrials(k) = size(trialmat,1);
        end;
        [optk(cmp),centr,clst,Cg] =Kmeangrp(kmat,5,reps,0); % find optimal cluster number, don't plot
        [kout(:,cmp), C,sumd, allsums{cmp}] = kmeans(kmat,optk(cmp),'replicates',reps);
        fprintf('\n One More Component Done: %i of %i',cmp,length(gdcomps{nx}));
        kmatall(:,:,cmp) = kmat;
    end;
    %allkouts{nx} = kout;
    %optks{nx} = optk;
    %numtrials{nx} = subjtrials;  clear Alllongersps  totrials icamat  onecmp oneemo kout centr clst Cg optk
    fprintf('\n One More SUBJECT Done: %i\n',nx);
%end;
comment = 'Ran Kmeangrp.m which finds optimal number of kmeans clusters to use. Then ran kmeans with that value. Shold include:kout optk subjtrials comment ';
save KmeansClusterInfo.mat kout optk subjtrials comment kmatall

% find emotion assignments of clusters  USING SINGLE TRIALS
 cd (['/data/common2/emotion',paths{nx},'/ersps/']);  load KmeansClusterInfo.mat
 cols(1,:) = [1,0,0];cols(2,:) = [0,0,1];cols(3,:) = [0,1,0];cols(4,:) = [.5,.5,.5];
emomap = ones(1,1);
for e = 2:length(subjtrials)+1
        emomap(1,e) = emomap(e-1) + subjtrials(e-1); % marks where each emotion STARTS
end;
clustcomps = find(optk<5);  thresh = 1000; % upper distance limit of plotted spectra
for cp = 1:length(clustcomps)
figure;pl = 1;  cmp = clustcomps(cp);
for e = 1:length(subjtrials)
    subplot(6,6,pl);clear clustwts
    clustvals = kout(emomap(e):emomap(e+1)-1,cmp);
    onlywts = allsums{cmp}(emomap(e):emomap(e+1)-1,:);
    for k = 1:optk(cmp)
        clustwts(k) = length(find(clustvals == k));
    end;            
    hist(clustvals,optk(cmp)); set(gca,'xlim',[0 optk(cmp)+1]); pl = pl+1;
    title(emos{e});
    %[s r] = max(clustwts); % r gives position(cluster)
    realspec = kmatall(emomap(e):emomap(e+1)-1,:,cmp);
    subplot(6,6,pl)
    for w = 1:optk(clustcomps(cmp))
        %plotspec = realspec(find(clustvals==w),:);
        wtdist = find(onlywts(find(clustvals==w),w)<thresh);
        plotspec = realspec(wtdist,:);
        if ~isempty(plotspec)
            ph = plot(freqs,mean(plotspec),'r-');hold on;
            set(ph,'color',cols(w,:));
        end;
    end;
pl = pl+1; set(gca,'ylim',[-10 5]);   
    set(gca,'xlim',[1 50]); set(gca,'xgrid','on');
end;
end;

