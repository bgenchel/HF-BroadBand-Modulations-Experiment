% clusters coherence from crossf across different components using ICA
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
%subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];

%  Pull out features of crosses
%load  /data/common1/emotion/jo74/RndEpochsXcoh.mat
%load /home/julie/EmoTF/jo74CrsEmos.mat 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
% freqsout(1:20)
%    3.0000    3.7500    4.5000    5.2500    6.0000    6.7500    7.5000
%    8.2500    9.0000    9.7500   10.5000   11.2500   12.0000   12.7500
%   13.5000   14.2500   15.0000   15.7500   16.5000   17.2500

frqrange = find(freqsout < 30);
allcrs  =zeros(length(subj1),length(frqrange)*length(emos)); m=1;
for cp1 = 1:length(subj1)-1
    for cp2 = 2:length(subj1)
        tmpcrs  =zeros(1,0);
        for e = 1:length(emos)
            cohercell = emocoh{e};
            onecoh = cohercell{1,subj1(cp1)};
            crsbootcell = emoboot{e};
            bootmat1 = crsbootcell {1,subj1(cp1)};   %chooses a 3D array from cell array
            minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
            minmask = repmat (minboot,[1,size(timesout,2),1]); 
            onecoh(find(onecoh <= minmask)) = 0;%finds and zeros cohmats less than minmask
            onecoh = onecoh(frqrange,:,subj1(cp2));
            oneval = mean(onecoh,2);
            tmpcrs(1,end+1:end+size(oneval,1)) = oneval';
        end;
        allcrs(m,:) = tmpcrs; m=m+1;
    end;
fprintf('\n cp1 number done: %i',cp1);
end;

%  225 'chan' X  360 pnts (number of crosses X (freq*emo) matrix)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Run PCA
[pc,eigvec,sv] = runpca(allcrs');% 
%%% pc = #comps X #comps : gives you weights
%%% eigvec= totalength X #comps: gives you the ersp comp maps once re-transformed
%[pc,score,latent,tsquare]=princomp(pcamatall); %to see contribs to data
%figure; pareto(latent);

pcaresults = zeros (length(emos)-4,length(fr),size(weights,1));
for n = 1:size(pc,2)
tmpcomp = eigvec (:,n);  %makes a totalength X 1
tmpcomp = reshape(tmpcomp,length(fr),length(emos)-4);  %Good
tmpcomp = tmpcomp';
pcaresults(:,:,n) = tmpcomp;
end; 
% NOW PLOT THE RESULTS %%
figure;
for k=1:9
subplot(3,3,k);
imagesc(freqsout(frqrange),[1:size(pcaresults,1)],pcaresults(:,:,k));
end;
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run ICA

[weights,sphere,compvars,bias,signs,lrates,activations] = runica (allcrs,'pca',10);

% 328 steps for all ?, fr 3-30
icadata = zeros (length(emos),length(frqrange),size(weights,1));
for n = 1:size(weights,1)
    tmpcomp = activations (n,:);
    tmpcomp = reshape(tmpcomp,length(frqrange),length(emos));  %Good
    tmpcomp = tmpcomp';
    icadata(:,:,n) = tmpcomp;
end; 

%%%%  Image the results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image as color
figure;pl=1;
for comp = 1:size(weights,1)
    subplot(3,4,pl)
    imagesc([1:size(icadata,1)],freqsout(frqrange),icadata(:,:,comp)',[-10 10]);hold on;
    title(int2str(comp)); 
    set(gca,'ticklength',[.03 .03]);
    pl=pl+1;
end;

figure;pl=1;
for comp = 1:gdcomp
    subplot(4,4,pl)
    for n = 1:length(subj1)
        plot(freqs(fr),icadata(n,:,comp));hold on;
    end;
    set(gca,'xlim',[3 30]);
    title(int2str(comp)); 
    pl=pl+1;
end;


%%%%%%%%%%%%%%%%%%%
winv = pinv(weights*sphere);
figure;
for n = 1:size(weights,1)
    subplot(4,4,n)
    g=bar(winv(n,:));hold on;
    %set(gca,'ylim',[-.02 .02]);
    set(gca,'xlim',[0 17]);
    set(gca,'xtick',[1:16]);
    set(gca,'xgrid','on');
    set(gca,'fontsize',7);
        %set(g,'color',cb);
end;

%%%%%%%%%%%%%%%%%%%
%%%%% plot winvs against each other:
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
p1 = icadata(:,10:12,1);
p1=mean(p1,2);
p2 = icadata(:,10:12,9);
p2=mean(p2,2);
p3 = icadata(:,11:13,6);
p3=mean(p3,2);
p4 = icadata(:,15:17,6);
p4=mean(p4,2);
p5 = icadata(:,10:12,6);
p5=mean(p5,2);
p6 = icadata(:,11:13,6);
p6=mean(p6,2);
p7 = icadata(:,24:36,6);
p7=mean(p7,2);
p8 = icadata(:,10:12,6);
p8=mean(p8,2);
p9 = icadata(:,10:12,6);
p9=mean(p9,2);
p10 = icadata(:,10:12,6);
p10=mean(p10,2);

figure;
plot3(p3,p5,p7,'.');hold on;
for k=1:14
text(p3(k),p5(k),p7(k),emos(k));
end;
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');



%  Try clustering within a single emo to get a footprint of cross coh

%  Pull out features of crosses
%load  /data/common1/emotion/jo74/RndEpochsXcoh.mat
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
% freqsout(1:20)
%    3.0000    3.7500    4.5000    5.2500    6.0000    6.7500    7.5000
%    8.2500    9.0000    9.7500   10.5000   11.2500   12.0000   12.7500
%   13.5000   14.2500   15.0000   15.7500   16.5000   17.2500

frqrange = find(freqsout < 30);
for e = 1:length(emos)
    m=1;  clear tmpcrs 
    for cp1 = 1:length(subj1)-1
        for cp2 = 2:length(subj1)
            cohercell = emocoh{e};
            onecoh = cohercell{1,subj1(cp1)};
            crsbootcell = emoboot{e};
            bootmat1 = crsbootcell {1,subj1(cp1)};   %chooses a 3D array from cell array
            minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
            minmask = repmat (minboot,[1,size(timesout,2),1]); 
            %onecoh(find(onecoh <= minmask)) = 0;%finds and zeros cohmats less than minmask
            onecoh = onecoh(frqrange,:,subj1(cp2));
            onecoh = reshape(onecoh,1,size(onecoh,1)*size(onecoh,2));
            tmpcrs(m,:) = onecoh; m=m+1;
        end;
        allcrs{e} = tmpcrs; 
    end;
    fprintf('\n emo number done: %i',e);
end;
% Run ICA

[weights,sphere,compvars,bias,signs,lrates,activations] = runica (allcrs{e},'pca',10);

% 439 steps for awe unmasked, fr 3-30
tm = length(timesout);
icadata = zeros (tm,length(frqrange),size(weights,1));
for n = 1:size(weights,1)
    tmpcomp = activations (n,:);
    tmpcomp = reshape(tmpcomp,length(frqrange),tm);  %Good
    tmpcomp = tmpcomp';
    icadata(:,:,n) = tmpcomp;
end; 

%%%%  Image the results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image as color
figure;pl=1;
for comp = 1:size(weights,1)
    subplot(3,4,pl)
    imagesc(timesout,freqsout(frqrange),icadata(:,:,comp)',[-10 10]);hold on;
    title(int2str(comp)); 
    set(gca,'ticklength',[.03 .03]);
    pl=pl+1;
end;
