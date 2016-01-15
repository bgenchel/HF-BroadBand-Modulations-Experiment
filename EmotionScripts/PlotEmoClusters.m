% finds comps with max variance for each factor, then plots those dipoles with spectral templates

% input subject lists from /data/common2/emotion/Complist
eeglab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sphfile = {'Subj1pc10.sph','Subj2pc10.sph','Subj3pc10.sph','Subj4pc10.sph','Subj5pc10.sph','Subj6pc10.sph','Subj7pc10.sph','Subj8pc10.sph','Subj9pc10.sph','Subj10pc10.sph','Subj11pc10.sph','Subj12pc10.sph','Subj13pc10.sph','Subj14pc10.sph','Subj15pc10.sph','Subj16pc10.sph','Subj17pc10.sph','Subj18pc10.sph','Subj19pc10.sph','Subj20pc10.sph','Subj21pc10.sph','Subj22pc10.sph','Subj23pc10.sph','Subj24pc10.sph','Subj25pc10.sph','Subj26pc10.sph'};
%wtsfile = {'Subj1pc10.wts','Subj2pc10.wts','Subj3pc10.wts','Subj4pc10.wts','Subj5pc10.wts','Subj6pc10.wts','Subj7pc10.wts','Subj8pc10.wts','Subj9pc10.wts','Subj10pc10.wts','Subj11pc10.wts','Subj12pc10.wts','Subj13pc10.wts','Subj14pc10.wts','Subj15pc10.wts','Subj16pc10.wts','Subj17pc10.wts','Subj18pc10.wts','Subj19pc10.wts','Subj20pc10.wts','Subj21pc10.wts','Subj22pc10.wts','Subj23pc10.wts','Subj24pc10.wts','Subj25pc10.wts','Subj26pc10.wts'};
sphfile = {'Subj1pc15.sph','Subj2pc15.sph','Subj3pc15.sph','Subj4pc15.sph','Subj5pc15.sph','Subj6pc15.sph','Subj7pc15.sph','Subj8pc15.sph','Subj9pc15.sph','Subj10pc15.sph','Subj11pc15.sph','Subj12pc15.sph','Subj13pc15.sph','Subj14pc15.sph','Subj15pc15.sph','Subj16pc15.sph','Subj17pc15.sph','Subj18pc15.sph','Subj19pc15.sph','Subj20pc15.sph','Subj21pc15.sph','Subj22pc15.sph','Subj23pc15.sph','Subj24pc15.sph','Subj25pc15.sph','Subj26pc15.sph','Subj27pc15.sph'};
wtsfile = {'Subj1pc15.wts','Subj2pc15.wts','Subj3pc15.wts','Subj4pc15.wts','Subj5pc15.wts','Subj6pc15.wts','Subj7pc15.wts','Subj8pc15.wts','Subj9pc15.wts','Subj10pc15.wts','Subj11pc15.wts','Subj12pc15.wts','Subj13pc15.wts','Subj14pc15.wts','Subj15pc15.wts','Subj16pc15.wts','Subj17pc15.wts','Subj18pc15.wts','Subj19pc15.wts','Subj20pc15.wts','Subj21pc15.wts','Subj22pc15.wts','Subj23pc15.wts','Subj24pc15.wts','Subj25pc15.wts','Subj26pc15.wts','Subj27pc15.wts'};
load /data/common2/emotion/clusters/subjdims.mat % subjdims  numtrials comment
freqs = [1:.5:50]; 

%%%%%%%%% END VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
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

% now plot dipoles, spectral templates and emotion weightings
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; clear cols
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
row = 5; col = 6;  clear mnemodiff
for nx = 13:21%length(gdcomps)
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    figure; op = 28;cpcols = hsv(length(gdcomps{nx}));
    for tp = size(winv,2):-1:1
        subplot(row,col,op);        clear dipcols
        for cp = 1:length(nxlists{nx}{tp})
            rcp = find(nxlists{nx}{tp}(cp) == gdcomps{nx});
            dipcols{cp} = cpcols(rcp,:);
        end;        
        dipplot(EEG.dipfit.model(nxlists{nx}{tp}),'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',dipcols,'projlines','on','projimg','on'); op = op-3; view(60,20);  
    end;p=2;
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; 
        if mean(tpwts)<0
            tpwts = tpwts*-1;
            newact = activations(tp,:)*-1;
        else
            newact = activations(tp,:);  
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
        for e = 1:length(numtrials{nx})-1
            e=emoorder(e);clear newmat
            tempmat = tpwts(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
            if e == 1
                pl=1; for pk = .1:.1:.9
                    basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat1,basemat1-basemat1,'.');hold on;
                set(ph,'color',[.5 .5 .5]);        set(ph,'markersize',16);
                pl=1; 
                for pk = .03:.01:.98
                    basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat,basemat-basemat);hold on;
                handvec(1,e) = ph;
                set(ph,'color',[.5 .5 .5]);        set(ph,'linewidth',1.5); 
                mx = max(basemat-basemat); mn = min(basemat-basemat);
            else
                pl=1;
                for pk = .1:.1:.9
                    newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                ph = plot(basemat1,newmat1-basemat1,'.');hold on;
                set(ph,'color',cols(e,:));        set(ph,'markersize',16);pl=1;
                for pk = .03:.01:.98
                    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
               ph = plot(basemat,newmat-basemat);hold on;
                handvec(1,e) = ph;
                set(ph,'color',cols(e,:));         set(ph,'linewidth',1.5); 
                mx1 = max(newmat-basemat); mn1 = min(newmat-basemat);
                if mx1 > mx
                    mx = mx1;
                end;
                if mn1 < mn
                    mn = mn1;
                end;                
            end; 
            set(gca,'xlim',[basemat(1)-basemat(1)*.02 basemat(end)+basemat(end)*.02]);
            set(gca,'box','off'); set(gca,'fontsize',10); title(['Factor ',int2str(tp)]);
        end;
        set(gca,'ylim',[mn+mn*.01 mx+mx*.01]);   p = p+2;
    end; set(gcf,'color','w');    
    ph=textsc(['Subject ',int2str(nx),'; All Comp Templates and Percentile Comparisons with Prebaseline Weights for all Factors'],'title');set(ph,'fontsize',12);set(ph,'color','r');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    print  -dpsc2 -Pcoloring 
    close
    ALLEEG=[];EEG=[];
end;


% simple version, no plot
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
clear mnemodiff forstats
for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[15 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; 
        if mean(tpwts)<0
            tpwts = tpwts*-1;
            newact = activations(tp,:)*-1;
        else
            newact = activations(tp,:);  
        end;
        for cp = 1:length(nxlists{nx}{tp})
            rcp = find(nxlists{nx}{tp}(cp) == gdcomps{nx});
        end;
        clear basemat basemat1
        for e = 1:length(numtrials{nx})-1
            e=emoorder(e);clear newmat
            tempmat = tpwts(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
            if e == 1
                pl=1;                 
                for pk = .1:.1:.9
                    basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;pl=1;
                for pk = .03:.01:.98
                    basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
            else
                pl=1;
                for pk = .1:.1:.9
                    newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;pl=1;
                for pk = .03:.01:.98
                    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                forstats(:,e,tp) = newmat;
                mnemodiff(tp,e-1,nx) = mean(newmat-basemat);%%****  keeps mean scores
            end; 
        end;
    end; 
    subjstats{nx} = forstats;
    fprintf('\n One More SUBJECT Done: %i',nx);
end;
% compare emoscores across emos and do stats
alpha = .0000001; clear alldiffs allmasks mask diffmat
for tp = 1:10
    [P,ANOVATAB,STATS] = anova1(forstats(:,1:15,tp),nopem,'off');
    %figure;
    [compare,means,H] = multcompare(STATS,alpha,'off','bonferroni');
    for q = 1:size(compare,1)
        diffmat(compare(q,1),compare(q,2)) = compare(q,4);
        if  compare(q,5)>0 & compare(q,3)<0
            mask(compare(q,1),compare(q,2)) = 1;
        elseif compare(q,5)<0 & compare(q,3)>0
            mask(compare(q,1),compare(q,2)) = 1;    
        else        
            mask(compare(q,1),compare(q,2)) = 0;    
        end;    
    end;
    alldiffs{tp} = diffmat;
    allmasks{tp} = mask;
    sigtp(tp) = P;
end;
figure;
for tp = 1:length(alldiffs)
    subplot(3,4,tp)
    plotdiffs = alldiffs{tp};
    plotdiffs(find(allmasks{tp}))=0;
    imagesc(abs(plotdiffs),[-2 2]);
    for w = 1:size(alldiffs{tp},1)
        ph =text(1,w,emos{w+1});set(ph,'fontsize',7);
    end;  
    set(gca,'xtick',[2:2:14]); set(gca,'xgrid','on');title(['Factor ',int2str(tp)]);
end;axcopy
textsc(['Mulitple Comparisons between emotions; Only SIGNIFICANT voxels (p < ',num2str(alpha),')'],'title');
subplot(3,4,tp+1);    imagesc(alldiffs{tp},[-1 1]);colorbar
% for each template: sig diff emotion pairs
%In this example the confidence interval does not contain 0.0, so the difference is significant at the 0.05 level. If the confidence interval did contain 0.0, the difference would not be significant at the 0.05 level.


% for each subject, plot bar graph instead of lines to show emo dependency
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
cols(1,1:3) = [.5 .5 .5];cols(2:length(emos)-1,:) = jet(length(emos)-2);cols(end+1,:) = [0 0 0];
row = 5; col = 6;  clear mnemodiff
for nx = 1:6%length(gdcomps)
    EEG = pop_loadset( 'sources.set', ['/data/common2/emotion/',paths{nx}]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); 
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;
    figure; op = 28;cpcols = hsv(length(gdcomps{nx}));clear mnemodiff
    for tp = size(winv,2):-1:1
        subplot(row,col,op);        clear dipcols
        for cp = 1:length(nxlists{nx}{tp})
            rcp = find(nxlists{nx}{tp}(cp) == gdcomps{nx});
            dipcols{cp} = cpcols(rcp,:);
        end;        
        dipplot(EEG.dipfit.model(nxlists{nx}{tp}),'image','mri','gui','off','normlen','on','dipolesize',35,'dipolelength',0,'spheres','on','color',dipcols,'projlines','on','projimg','on'); op = op-3; view(60,20);  
    end;p=2;
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; 
        if mean(tpwts)<0
            tpwts = tpwts*-1;
            newact = activations(tp,:)*-1;
        else
            newact = activations(tp,:);  
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
        clear basemat basemat1
        emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; clear cols
        for em = 1:length(numtrials{nx})-1
            e=emoorder(em);clear newmat
            tempmat = tpwts(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
            if e == 1
                pl=1; for pk = .1:.1:.9
                    basemat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                pl=1; 
                for pk = .03:.01:.98
                    basemat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
            else
                pl=1;
                for pk = .1:.1:.9
                    newmat1(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                pl=1;
                for pk = .03:.01:.98
                    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
                end;
                mnemodiff(tp,e) = mean(newmat-basemat);%%****  keeps mean scores
            end; 
        end;
        subplot(row,col,p); 
emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'};
%emoorder = [4,2,10,8,14,12,6,9,7,13,11,1,5,3,15];  % no pre or post baseline
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; clear cols
cols = jet(size(mnemodiff,2));mxbar = max(mnemodiff(tp,:));mnbar = min(mnemodiff(tp,:));
        for e = 2:size(mnemodiff,2)
            ph =bar(e-1,mnemodiff(tp,e));hold on;
            set(ph,'facecolor',cols(e,:));
            ph=text(e-1,0,emos{emoorder(e)}); set(ph,'fontsize',7);
            set(ph,'rotation',90);
        end;set(gca,'xlim',[0 16]); set(gca,'ylim',[mnbar+mnbar*.01 mxbar+mxbar*.01]);
        set(gca,'xticklabel',[]); set(gca,'box','off');p = p+2;     
    end; set(gcf,'color','w');    
    ph=textsc(['Subject ',int2str(nx),'; All Comp Templates and Percentile Comparisons with Prebaseline Weights for all Factors'],'title');set(ph,'fontsize',12);set(ph,'color','r');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    print  -dpsc2 -Pcoloring 
    close
    ALLEEG=[];EEG=[];
end;
