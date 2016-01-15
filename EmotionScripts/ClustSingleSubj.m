% early button press sessions (not instructed to only press on surge)
subj1 = [1,3:22,25:30,34,36];  % tl81  23,35 =pulse; 
subj2 = [1,2,4:9,12,13,15:19,21,22,24,26,29];% mi83  14,23 = pulse; 31,34,46,51:muscle
subj3 = [5,6,7,10,12,13,14,15,17,18,19,22];   % ms82  11,37 = pulse; 32 = muscle?
subj4 = [1:5,8:10,12:14,16,17,19,22,25];  % js75  7 is pulse
subj5 = [1:6,9:12,14:19,22:25,28,34]; % kw78
%  instructed to press only at peak of emotion
subj6 = [1:3,5,6,8:11,14:17,20,23];% jo82 
subj7 = [1,2,5:10,15,22];% kl80 
subj8 = [1,2,6,8,9,11,12];% ar81 
subj9 = [1,2,3,5:11,13,14,16:20,22:24,29,30,32,36,37,42]; % eb79  But no button presses actually recorded
subj10 = [2,3,5:8,10:12,16:18,21,23,24,27,28]; % dg75  
subj11 = [1,2,4:18,21:23,29,32,34];  % an82  29 maybe
subj12 = [2,5,7,8,12,14,15,17,20,21,25,26,29];  % jw84
% NO button press sessions
subj13 = [1:6,8,9,10,13:15,17,19,20,23,25,27,30];  % tv81  11 = pulse?? 
subj14 = [1:4,6,8:14,16:21,23,24,28];  % sr81   22,27 = pulse  19 maybe
subj15 = [1:7,11,13:16,18,20,21,23,25,37,38,40];  % an70  9,12 = pulse
subj16 = [2:5,7:11,15:18,25];  % sg75  14 = pulse
subj17 = [1:3,6:8,10:14,16,19];% mr72
subj18 = [2,3,5,7:12,14,16:20];% dk74
subj19 = [4:8,10:15,17:24,26:28,30,33];% dn86  
subj20 = [1,3,6,7,10:16,18,19,21];% mr71
subj21 = [1:6,8:12,14:24,32];  % md85
%  instructed to press only at peak of emotion
subj22 = [1:3,6,7,9,10,11,13:15,17:25,28,29];  % mr72, second session
subj23 = [1:3,7:9,11:15,18:20,22,25]; %cj82   19 is 24%, 20 is 19%,22 outside brain
subj24 = [2:4,6:13,15,17,21,24];  % kc66
subj25 = [9,11,13,14,16,18,19,21,23,24];  % ts79
subj26 = [1,4,6,7,8,11,12,14:20];  % es76
gdcomps = {subj1, subj2, subj3, subj4 subj5, subj6, subj7, subj8,subj9,subj10, subj11, subj12, subj13 subj14, subj15, subj16, subj17,subj18,subj19,subj20, subj21, subj22, subj23, subj24, subj25, subj26};
paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/','/mr72-2/','/cj82/','/kc66/','/ts79/','/es76/'};

datsets = {'emo-1-241.set','emo-1-248.set','emo-1-238.set','emo-1-253.set','emo-1-250.set','emo-1-244.set','emo-1-248.set','emo-1-231.set','emo-1-250.set','emo-1-243.set','emo-1-251.set','emo-1-251.set','emo-1-180.set','emo-1-249.set','emo-1-250.set','emo-1-246.set','emo-1-237.set','emo-1-250.set','emo-1-250.set','emo-1-251.set','emo-1-250.set','emo-1-246.set','emo-1-248.set','emo-1-240.set','emo-1-254.set','emo-1-246.set'};

paths = {'/tl81/','/mi83/','/ms82/','/js75/','/kw78/','/jo82/','/kl80/','/ar81/','/eb79/','/dg75/','/an82/','/jw84/','/tv81/','/sr81/','/an70/','/sg75/','/mr72/','/dk74/','/dn86/','/mr71/','/md85/','/mr72-2/','/cj82/','/kc66/','/ts79/','/es76/'};


sphs = {'sph241.sph','sph248.sph','sph238-110.sph','sph253pc100.sph','sph250pc110.sph','sph244pc100.sph','sph248pc100.sph','sph231pc100.sph','sph250pc100.sph','sph243pc100.sph','sph251pc100.sph','sph251pc100.sph','sph180-90.sph','sph249pc100.sph','sph250pc100.sph','sph246pc100.sph','sph237pc100.sph','sph250pc100.sph','sph250pc100.sph','sph251pc100.sph','sph250pc100.sph','sph246pc100.sph','sph248pc100.sph','sph240pc100.sph','sph254pc100.sph','sph246pc100.sph'};
wtss = {'wts241.wts','wts248.wts','wts238-110.wts','wts253pc100.wts','wts250pc110.wts','wts244pc100.wts','wts248pc100.wts','wts231pc100.wts','wts250pc100.wts','wts243pc100.wts','wts251pc100.wts','wts251pc100.wts','wts180-90.wts','wts249pc100.wts','wts250pc100.wts','wts246pc100.wts','wts237pc100.wts','wts250pc100.wts','wts250pc100.wts','wts251pc100.wts','wts250pc100.wts','wts246pc100.wts','wts248pc100.wts','wts240pc100.wts','wts254pc100.wts','wts246pc100.wts'};

%             tl81       mi83      ms82     js75     kw78     jo82      kl80      ar81        eb79     dg75    an82 jw84   tv81      sr81     an70      sg75    mr72     dk74       dn86      mr71     md85     
sphsize = {[241 241],[248 248],[238 238],[253 253],[250 250],[244 244],[248 248],[231 231],[250 250],[243 243],[251 251],[251 251],[180 180],[249 249],[250 250],[246 246],[237 237],[250 250],[250 250],[251 251],[250 250],[246 246],[248 248],[240 240],[254 254],[246 246]};
wtssize = {[160 241],[160 248],[110 238],[100 253],[110 250],[100 244],[100 248],[100 231],[100 250],[100 243],[100 251],[100 251], [90 180],[100 249],[100 250],[100 246],[100 237],[100 250],[100 250],[100 251],[100 250],[100 246],[100 248],[100 240],[100 254],[100 246]};
 
subjspecs = {'Subj1trials.fdt','Subj2trials.fdt','Subj3trials.fdt','Subj4trials.fdt','Subj5trials.fdt','Subj6trials.fdt','Subj7trials.fdt','Subj8trials.fdt','Subj9trials.fdt','Subj10trials.fdt','Subj11trials.fdt','Subj12trials.fdt','Subj13trials.fdt','Subj14trials.fdt','Subj15trials.fdt','Subj16trials.fdt','Subj17trials.fdt','Subj18trials.fdt','Subj19trials.fdt','Subj20trials.fdt','Subj21trials.fdt','Subj22trials.fdt','Subj23trials.fdt','Subj24trials.fdt','Subj25trials.fdt','Subj26trials.fdt'};
emos = {'prebase','awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite','postbase'};

sphfile = {'Subj1pc10.sph','Subj2pc10.sph','Subj3pc10.sph','Subj4pc10.sph','Subj5pc10.sph','Subj6pc10.sph','Subj7pc10.sph','Subj8pc10.sph','Subj9pc10.sph','Subj10pc10.sph','Subj11pc10.sph','Subj12pc10.sph','Subj13pc10.sph','Subj14pc10.sph','Subj15pc10.sph','Subj16pc10.sph','Subj17pc10.sph','Subj18pc10.sph','Subj19pc10.sph','Subj20pc10.sph','Subj21pc10.sph','Subj22pc10.sph','Subj23pc10.sph','Subj24pc10.sph','Subj25pc10.sph','Subj26pc10.sph'};
wtsfile = {'Subj1pc10.wts','Subj2pc10.wts','Subj3pc10.wts','Subj4pc10.wts','Subj5pc10.wts','Subj6pc10.wts','Subj7pc10.wts','Subj8pc10.wts','Subj9pc10.wts','Subj10pc10.wts','Subj11pc10.wts','Subj12pc10.wts','Subj13pc10.wts','Subj14pc10.wts','Subj15pc10.wts','Subj16pc10.wts','Subj17pc10.wts','Subj18pc10.wts','Subj19pc10.wts','Subj20pc10.wts','Subj21pc10.wts','Subj22pc10.wts','Subj23pc10.wts','Subj24pc10.wts','Subj25pc10.wts','Subj26pc10.wts'};
load /data/common2/emotion/clusters/subjdims.mat % subjdims  numtrials comment
freqs = [1:.5:50]; 
%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);
clear newact
    for tp = 1:size(winv,2)
        for cp = 1:length(gdcomps{nx})
            newact(tp,cp,:) = activations(tp,length(freqs)*(cp-1)+1:length(freqs)*cp);
        end;      
    end;
    for tp = 1:size(winv,2)
        for ff = 1:length(freqs)
            tmpact = newact(tp,:,ff);
            tvar(tp,ff) = var(tmpact);
        end;
    end;
    % to Normalize by lowest variance freqencies for each template
    for tp = 1:size(winv,2)
        mvar = min(tvar(tp,:));
        lowvar{tp} = find(tvar(tp,:) >= mvar & tvar(tp,:) <= mvar+mvar*.4); % find freqs with low var
        tmpcps = newact(tp,:,lowvar{tp}); tmpcps = squeeze(tmpcps); % makes a cp X freqs
        newact(tp,:,:) = newact(tp,:,:)/mean(std(tmpcps'));   % mean of std over all comps        
    end;
    
    % find most variant components by zscore cutoff
    clear tplists
    for tp = 1:size(winv,2)
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
    
    fprintf('\n One More SUBJECT Done: %i',nx);
end;
eeglab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% now plot dipoles, spectral templates and emotion weightings
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; clear cols
%cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
cols = jet(15);
row = 10; col = 5; 
for nx = 18:length(gdcomps)
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    figure; op = 46;cpcols = hsv(length(gdcomps{nx}));
    for tp = size(winv,2):-1:1
        subplot(row,col,op);        clear dipcols
        for cp = 1:length(nxlists{nx}{tp})
            rcp = find(nxlists{nx}{tp}(cp) == gdcomps{nx});
            dipcols{cp} = cpcols(rcp,:);
        end;        
        dipplot(EEG.dipfit.model(nxlists{nx}{tp}),'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',dipcols,'projlines','on','projimg','on'); op = op-5; view(60,20);   camzoom(1.7)
    end;p=2;
    for tp = 1:size(winv,2)
        coi = nxlists{nx}{tp};
        tpwts =  winv(:,tp)'; 
        for pn = 1:2
            if pn == 1
                tpwts = tpwts(find(tpwts>0));
                newact = activations(tp,:);
            else
                tpwts = tpwts(find(tpwts<0));
                newact = activations(tp,:)*-1;  
            end;
            subplot(row,col,p)        
            mxcp = max(newact);         mncp = min(newact);  
            for cp = 1:length(nxlists{nx}{tp})
                rcp = find(nxlists{nx}{tp}(cp) == gdcomps{nx});
                ph=plot(freqs,newact(length(freqs)*(rcp-1)+1:length(freqs)*rcp));hold on; 
                set(ph,'color',cpcols(rcp,:)); set(ph,'linewidth',1.5); 
                mxcp1 = max(newact(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                mncp1 = min(newact(length(freqs)*(rcp-1)+1:length(freqs)*rcp)); 
                if mxcp1 > mxcp
                    mxcp = mxcp1;
                end;
                if mncp1< mncp
                    mncp = mncp1;
                end;    
            end;p = p+1;
            set(gca,'xlim',[0 40]);     set(gca,'ylim',[mncp+mncp*.01 mxcp+mxcp*.01]);   
            set(gca,'xtick',[5:5:40]);  set(gca,'xticklabel',{[] 10 [] 20 [] 30 [] 40}); 
            set(gca,'xgrid','on'); set(gca,'fontsize',10);    
            set(gca,'ycolor','c');
            subplot(row,col,p); 
            clear basemat basemat1
            emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
            pl = 1; clear tpemoscores
            % Specs are collected twice (once is *-1)
            for e = 2:length(numtrials{nx})-1  % start with 2 for straight nums (not diffs)
                e=emoorder(e);
                tpwts =  winv(:,tp)';                
                tempmat = tpwts(emomap(e):emomap(e+1)-1); 
                if pn == 1
                    tpemoscores(e-1,pl) = sum(tempmat(find(tempmat>0)))/length(tempmat);
                else                    
                tpemoscores(e-1,pl) = abs(sum(tempmat(find(tempmat<0))))/length(tempmat);
                end;
            end; ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(row,col,p)        
            tpemoscores(find(newsubj{nx}==0),pl) = 0; % zero out 'bad' emotions
            for e = 1:size(tpemoscores,1)
                if tpemoscores(e,:) ~= 0
                    ph =bar(e,tpemoscores(e,:));hold on;
                    set(ph,'facecolor',cols(e,:));
                    ph=text(e,0,emo2{e}); set(ph,'fontsize',7);
                    set(ph,'rotation',90);
                    set(gca,'ylim',[0 max(tpemoscores)+max(tpemoscores)*.01]);
                end;
            end;set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]); set(gca,'box','off');
          ph=title(['Factor ',int2str(tp)]); set(ph,'fontsize',7); p = p+1;
        end; p = p+1;% to pn
    end; set(gcf,'color','w');    
    ph=textsc(['Subject ',int2str(nx),'      All ICA Factors with Weights Split between Positive and Negative'],'title');set(ph,'fontsize',14);set(ph,'color','r');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    print  -dpsc2 -Pcoloring 
    close
    ALLEEG=[];EEG=[];
end;

% NO PLOT FOR STATISTICS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    for tp = 1:size(winv,2)
        coi = nxlists{nx}{tp};
        tpwts =  winv(:,tp)'; 
        for pn = 1:2
            emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
            tpemoscores = zeros(0,1); grpvar = zeros(0,1);
            for e = 2:length(numtrials{nx})-1  % start with 2 for straight nums (not diffs)
                e=emoorder(e);
                tpwts =  winv(:,tp)';                
                tempmat = tpwts(emomap(e):emomap(e+1)-1); 
                if pn == 1
                    tpemoscores1(end+1:end+length(find(tempmat>0)),1) = tempmat(find(tempmat>0))';
                else                    
                    tpemoscores2(end+1:end+length(find(tempmat<0)),1) = tempmat(find(tempmat<0))';
                end;
            end; ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %            tpemoscores(find(newsubj{nx}==0),pl) = 0; % zero out 'bad' emotions
        end;
            if size(tpemoscores,1)> 3
                [F,pr,df,ss,ms,varcomp,varprop] = anova(tpemoscores,grpvar,10);
                sigfacts(tp,pn) = pr;            
            end;
    end;
    sigfacts
    sigsubjs{nx} = sigfacts;
end;

