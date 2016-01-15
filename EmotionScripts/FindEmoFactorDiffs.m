% to search subject ICA clusters for side of weights that best separates emotions
 
%Input heading from PlotEmoClusters.m
% input subject ratings from DescriptStats.m, use last one that zeros out bad emos
emoorder = [1,5 3 11 9 15  13 7  10  8 14 12  2  6 4 16,17]; 
realorder = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'};
emo2 = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excited'};

for nx = 1:length(gdcomps)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{nx}],[subjdims{nx}(1) subjdims{nx}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{nx}],[10 subjdims{nx}(1)]); 
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{nx}],[subjdims{nx}(1) subjdims{nx}(2)]);
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{nx})+1
        emomap(1,e) = emomap(e-1) + numtrials{nx}(e-1); % marks where each emotion STARTS
    end;pl = 1;
    tw = 1;     clear mnforica 
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; 
        clear tpforica
        for e = 2:length(numtrials{nx})-1 % start at 2 to skip prebase
            e=emoorder(e);pp=1;clear newmat
            tempmat = tpwts(emomap(e):emomap(e+1)-1); 
            tempmat = zscore(tempmat)';
            hivals(e-1) = mean(tempmat(find(tempmat>.35)))*(length(tempmat(find(tempmat>.35))))/length(tempmat);
           lovals(e-1) = mean(tempmat(find(tempmat<-.35)))*(length(tempmat(find(tempmat>.35))))/length(tempmat);
        end;
        lodist = mean(pdist(lovals'));
        hidist = mean(pdist(hivals'));
        if lodist - hidist > .01
            tpside(tp) = -1;
            mnforica(:,tp) = lovals';
        elseif hidist - lodist > .01
            tpside(tp) = 1;
            mnforica(:,tp) = hivals';
        end;     
        %mnforica(:,tw) = hivals';
        %mnforica(:,tw+1) = hivals';tw = tw+2;
    end; 
    % remove the mean from each factor
    for q=  1:size(mnforica,2)
        mnforica(:,q) = mnforica(:,q)/mean(mnforica(:,q));
    end;
   [weights,sphere,compvars,bias,signs,lrates,activations] = runica(mnforica,'pca',3,'extended',1,'stop',1e-8);
    ws = weights*sphere; winv = pinv(ws);
% winv = activations';
figure;
for e = 1:size(winv,1)
        ph =plot3(winv(e,1),winv(e,2),winv(e,3),'.');hold on;
        set(ph,'markersize',20);
        set(ph,'color',cols(e,:));
        ph = text(winv(e,1),winv(e,2),winv(e,3),enames{e});
        set(ph,'color',cols(e,:)); 
    end;pp = pp+1;
    mx =  max(winv(:,1));mnx = min(winv(:,1));
    my =  max(winv(:,2));mny = min(winv(:,2));
    mz =  max(winv(:,3));mnz = min(winv(:,3));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
    zl = get(gca,'zlim');
    for e = 1:size(winv,1)
        pl =plot3([winv(e,1) winv(e,1)],[winv(e,2) winv(e,2)],[zl(1)  winv(e,3)]);
        set(pl,'color',cols(e,:))             
    end;
    set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
    % set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);

    xlabel(['Winv ',int2str(1)]); ylabel(['Winv ',int2str(2)]); zlabel(['Winv ',int2str(3)]); 
   ph = title(['Subject ',int2str(nx)]); set(ph,'fontsize',12);
set(gcf,'color','w');

figure;  % for pca to 2
for e = 1:size(winv,1)
        ph =plot(winv(e,1),winv(e,2),'.');hold on;
        set(ph,'markersize',20);
        set(ph,'color',cols(e,:));
        ph = text(winv(e,1),winv(e,2),enames{e});
        set(ph,'color',cols(e,:)); 
    end;pp = pp+1;
    mx =  max(winv(:,1));mnx = min(winv(:,1));
    my =  max(winv(:,2));mny = min(winv(:,2));
    set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);
    set(gca,'xgrid','on');  set(gca,'ygrid','on');
    xlabel(['Winv ',int2str(1)]); ylabel(['Winv ',int2str(2)]); 
   ph = title(['Subject ',int2str(nx)]); set(ph,'fontsize',12);
set(gcf,'color','w');
 
   [weights,sphere,compvars,bias,signs,lrates,activations] = runica(mnforica,'pca',5,'extended',1,'stop',1e-8);
    ws = weights*sphere; winv = pinv(ws);
figure;  row = 5; col = 4;pp = 1;  % for multiple subjects
for cb = 1:size(winv,2)-2
    for bc = cb+1:size(winv,2)-1
        for cbc = bc+1:size(winv,2)
            subplot(row,col,pp)
            for e = 1:size(winv,1)
                ph =plot3(winv(e,cb),winv(e,bc),winv(e,cbc),'.');hold on;
                set(ph,'markersize',20);
                set(ph,'color',cols(e,:));
                ph = text(winv(e,cb),winv(e,bc),winv(e,cbc),enames{e});
                set(ph,'color',cols(e,:)); 
            end;
            mx =  max(winv(:,cb));mnx = min(winv(:,cb));
            my =  max(winv(:,bc));mny = min(winv(:,bc));
            mz =  max(winv(:,cbc));mnz = min(winv(:,cbc));
            set(gca,'xlim',[mnx mx]); set(gca,'ylim',[mny my]);set(gca,'zlim',[mnz mz]);
            zl = get(gca,'zlim');
            for e = 1:size(winv,1)
                pl =plot3([winv(e,cb) winv(e,cb)],[winv(e,bc) winv(e,bc)],[zl(1)  winv(e,cbc)]);
                set(pl,'color',cols(e,:))             
            end;
            set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
             set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'zticklabel',[]);
            
            xlabel(['Winv ',int2str(cb)]); ylabel(['Winv ',int2str(bc)]); zlabel(['Winv ',int2str(cbc)]); 
            pp = pp+1;
        end;
    end;
end;axcopy
set(gcf,'color','w');
