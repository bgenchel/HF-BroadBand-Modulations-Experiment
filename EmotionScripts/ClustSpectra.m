% collects spectra from all comps/emos, subtracts baseline and runs PCA or ICA to cluster
eeglab

gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21};

datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-245.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph245pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','sph245pc100.sph','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[245 245],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 ],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250]};

subjspecs = {'Subj1Specs.fdt','Subj2Specs.fdt','Subj3Specs.fdt','Subj4Specs.fdt','Subj5Specs.fdt','Subj6Specs.fdt','Subj7Specs.fdt','Subj8Specs.fdt','Subj9Specs.fdt','Subj10Specs.fdt','Subj11Specs.fdt','Subj12Specs.fdt','Subj13Specs.fdt','Subj14Specs.fdt','Subj15Specs.fdt','Subj16Specs.fdt','Subj17Specs.fdt','Subj18Specs.fdt','Subj19Specs.fdt','Subj20Specs.fdt','Subj21Specs.fdt'};
emoset = {'prebase.set','awe.set','frustration.set','joy.set','anger.set','happy.set','sad.set','love.set','fear.set' ,'compassion.set','jealousy.set','content.set','grief.set','relief.set','disgust.set','excite.set','postbase.set'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for nx = 1:length(gdcomps)
    cd (['/data/common2/emotion',paths{nx},'/ersps/']);load ContDataERSPs.mat
    icamatall = zeros(length(freqs),0);
    for k = 1:length(Alllongersps)
        oneemo = Alllongersps{k}; % freqsXtimesXcomponent   
            trialmat = zeros(0,size(oneemo,1));
        for cmp = 1:length(gdcomps{nx})
            onecmp = oneemo(:,:,gdcomps{nx}(cmp));
            for trl = 1:4:size(oneemo,2)-3
                onetrl = mean(onecmp(:,trl:trl+3),2);
                trialmat(end+1,:) = onetrl';
            end;
        fprintf('\n One More Component Done: %i of %i',cmp,length(gdcomps{nx}));
        end;
        icamatall(:,end+1:end+size(trialmat,1)) = trialmat';
        fprintf('\n One More Emotion Done: %i\n',k);
    end;
end;

        [weights,sphere,compvars,bias,signs,lrates,activations] = runica(icamatall,'pca',25,'maxsteps',2000,'stop',1e-7);
ws = weights*sphere;
winv = pinv(ws);

figure;row = ceil(sqrt(size(winv,2))); col = ceil(sqrt(size(winv,2))); pl = 1;
for tp = 1:size(winv,2)
        subplot(row,col,pl)
        plot(freqs,winv(:,tp)); pl = pl+1;
        set(gca,'xlim',[0 50]);
        set(gca,'ylim',[-5 5]);set(gca,'xgrid','on');       
    end;
end;



% to delete uncorrelated trials for ICA clustering
load /data/common1/emotion/Corridx.mat  
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
    load /data/common1/emotion/jo74/emXsbcptr.mat
alltrials = zeros(length(subj1)*length(freqs),14);
for em = 1:14
    load /data/common1/emotion/jo74/emXsbcptr.mat
    oneemo = alltremo{em};  % pulls out alltrials (TF X comps*trials)
    clear alltremo
    numtrial = size(oneemo,2)/length(subj1);
    allt = zeros(0,1);
    for cp = 1:length(subj1)
        onecomp = oneemo(:,(cp-1)*numtrial+1:cp*numtrial);
        onecomp = onecomp';
        onecomp = reshape(onecomp,size(onecomp,1),length(freqs),length(tm));
        onecomp = shiftdim(onecomp,1);  % makes a freqs X times X trials
      % to keep correlated trials
        onecomp = onecomp(:,:,allcorridx{em}{cp});
        onecomp = mean(onecomp,2);
        onecomp = mean(onecomp,3);
        onecomp = squeeze(onecomp);  % freqs X 1
        allt(end+1:end+size(onecomp,1),1) = onecomp;
        clear onecomp
    end;
    alltrials(:,em) = allt;
    fprintf('\nDone Emotion: %s',emos{em});    
end;

% run PCA on alltrials

[pc,eigvec,sv] = runpca(alltrials');

% pc is emo X comps*spectra
% eigvec is emo(pcs) X emo(inputs)

figure;
ph =plot3(eigvec(1,:),eigvec(2,:),eigvec(3,:),'k.');hold on;
text(eigvec(1,:),eigvec(2,:),eigvec(3,:),emos);

set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');


% run ICA on alltrials

[weights,sphere,compvars,bias,signs,lrates,activations] =runica(alltrials,'pca',14);

winv = weights* sphere;
