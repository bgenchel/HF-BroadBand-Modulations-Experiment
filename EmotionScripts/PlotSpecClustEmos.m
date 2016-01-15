% plots the mean emotion weights for each of the ~8 spectral clusters
%
%
%

function PlotSpecClustEmos(savedat,fullpaths,emomeans)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Plot the clusters as means of spectral clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % emomeans -- {nx}(IM,emo)
    emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
    msize = 8;
    cols = jet(15);cols(10,:) = [.9 .9 0];
    q = load('/data/common4/emotion/clustfacs.mat'); 
    clustcols(1,:) = [1 0 0];clustcols(2,:) = [1 .3 .3];clustcols(3,:) = [1 .6 .6];
    clustcols(4,:) = [0 0 1];clustcols(5,:) = [.3 .3 1];
    clustcols(6,:) = [0 1 0];clustcols(7,:) = [.3 1 .3];clustcols(8,:) = [.6 1 .6];
    
    row = 4; col = 6;
    figure; pl = 1;      clust=0;
    for ds = 1:3
        if ds == 1
            %r = load('/data/common4/emotion/CoModAlphaClusts.mat');
            r = load('/data/common4/emotion/AllCoModAlpha.mat');
            %specs = r.alphatempls;
        elseif ds == 2
            %r = load('/data/common4/emotion/BetaClusters.mat');
            r = load('/data/common4/emotion/AllCoModBeta.mat');
            %specs = r.betatempls;
        elseif ds == 3
            %r = load('/data/common4/emotion/GammaClusters.mat');
            r = load('/data/common4/emotion/AllCoModGama.mat');
            %specs = r.gamatempls;
        end;            
        specs = r.newtempls;
        for cls = 1:length(r.facvec)
            clust = clust+1;
            submeans = zeros(0,15);
            for nx = 1:length(r.facvec{cls})
                if ~isempty(r.facvec{cls}{nx})
                    submeans(end+1:end+length(r.facvec{cls}{nx}),:) = emomeans{nx}(r.facvec{cls}{nx},:);
                    %onesubcls(end+1:end+size(clustspecs{cls}{nx},1),:) = clustspecs{cls}{nx}; % templates
                    %submeans(end+1:end+length(clustlist{cls}{nx}),:) = emomeans{nx}(clustlist{cls}{nx},:);
                end;
            end;
            
            sbplot(row,col,[pl pl+1])
            ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
            ph = plot([1:size(emomeans{nx},2)],mean(submeans,1),'-','linewidth',2);hold on;
            set(ph,'color',[.5 .5 .5]);
            ph = plot([1:size(emomeans{nx},2)],mean(submeans,1),'.','markersize',msize);hold on;
            set(ph,'color',[.5 .5 .5]);
            for em = 1:size(submeans,2)
                ph = plot(em,submeans(:,em),'.','markersize',msize);hold on;
                set(ph,'color',cols(em,:));
                ph = plot(em,mean(submeans(:,em),1),'.','markersize',msize*3);
                set(ph,'color',cols(em,:));
            end;        
            set(gca,'ylim',[-.5 .5]); set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
            title(['Spectral Cluster ',int2str(clust),' -- Emotion Weightings (Original Winv)']);
            pl = pl+2;            
            sbplot(row,col,pl)
            ph = plot([r.frs(1) r.frs(end)],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
            ph = plot(r.frs,specs{cls},'k-');     hold on;      
            set(ph,'color',clustcols(clust,:));
            set(gca,'xgrid','on');
            set(gca,'xlim',[r.frs(1) r.frs(end)]);pl = pl+1; 
            title([int2str(size(specs{cls},1)),' members']);
            %str = ['print /data/common4/emotion/KmeansClusts',int2str(cls),'.eps -depsc']; eval(str);
        end;
    end;     
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    axcopy
     
    textsc('Emotion weightings for 8 spectral clusters (no further clustering)','title');
