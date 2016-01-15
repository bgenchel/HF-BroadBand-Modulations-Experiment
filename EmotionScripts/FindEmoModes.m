% within subj finds clusters of spectra for each emotion and categorizes differently (ie anger1,anger2...)
% Input header from PlotEmoClusters.m

figure; pp = 1; 
for nx = 2:2%length(gdcomps)
    mnforica = zeros(15,0);
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    figure;
    for e = 1:length(numtrials{nx})-1
        clear trwts
        for trl = 1:numtrials{nx}(e)
            tpwts =  winv(emomap(e)+(trl-1),:); 
            [x trwts(trl)] = max(abs(tpwts));
        end;
        subplot(4,4,e)
        hist(trwts);set(gca,'ylim',[0 100]);
        emofacts{e} = trwts;
    end; 
    fprintf('\n One More SUBJECT Done: %i',nx);
    [weights,sphere,compvars,bias,signs,lrates,activations] = runica(mnforica,'pca',3);
    ws = weights*sphere; winv = pinv(ws);



figure; pl = 1;
reps= 2;
for nx = 2:2%length(gdcomps)
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    icamatall = zeros(0,length(freqs)*length(gdcomps{nx}));
    for k = 1:length(Alllongersps)
        oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
        for cmp = 1:length(gdcomps{nx})
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            trialmat = zeros(0,size(oneemo,1));
            for trl = 1:2:size(onecmp,2)-5 % was 2; cut down to every 3 for subj 9 and 23 b/c too much data
                onetrl = mean(onecmp(:,trl:trl+5),2);
                trialmat(end+1,:) = onetrl';
            end;
            [optk,centr,clst,Cg] =Kmeangrp(trialmat,5,reps,0); % find optimal cluster number, don't plot
            [kout, C,sumd, allsums{cmp}] = kmeans(trialmat,optk,'replicates',reps);
            for pt = 1:optk
                ptmat(pt,:) = mean(trialmat(kout==pt,:));
            %subplot(4,4,pl)            
            %plot(freqs,mean(trialmat(kout==pt,:)));pl = pl+1;
            %title([int2str(k),' ',int2str(cmp),' ',int2str(pt)]);
            end;
                cpmat(end+1:end+size(ptmat,1),:) = ptmat;               
        end;
        kmat(end+1:end+size(cpmat,1),:) = cpmat;
    end;
    %floatwrite(icamatall, ['/data/common2/emotion/clusters/',subjspecs{nx}]); 
    subjdims{nx} = size(icamatall);
    numtrials{nx} = subjtrials;    
    fprintf('\n One More SUBJECT Done: %i\n',nx);

    [kout(:,cmp), C,sumd, allsums{cmp}] = kmeans(kmat,optk(cmp),'replicates',reps);
    
end;


figure;
subplot(2,2,1)
plot(freqs,trialmat(kout==1,:),'k');hold on;
subplot(2,2,2)
plot(freqs, trialmat(kout==2,:),'k');

figure;
subplot(2,2,1)
plot(freqs,x(1:100,:),'k');hold on;
subplot(2,2,2)
plot(freqs, y(1:100,:),'k');

