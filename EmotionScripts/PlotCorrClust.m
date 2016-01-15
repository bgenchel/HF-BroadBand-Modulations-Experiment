% Takes output of CorrEmoMeans() and plots as members of co-mod clusters
% creates a cell array 'clustlist' which is a:
% {EmoDimClusts}{SpecClusts}{nx}(IM idxs) array
%
%
%
% polarvec -- only necessary for ica clusters; give a 1 or -1 depending on the IM weight in the dimension







function PlotCorrClust(gdcomps,fullpaths,clustmeans,clusttrack,emomeans,polarvec)
    
    if ~exist('polarvec')
        polarvec = ones(1,1000);
    end;
    
    clear clustlist clustspecs clustemos
    r = load('/data/common4/emotion/CoModAlphaClusts.mat'); % to get r.frs
    nxs = cell(1,length(r.facvec{1}));
    for clust = 1:length(clusttrack)
        for cls = 1:8 % 8 spectral clusters
            clustlist{clust}{cls} = nxs;
            clustspecs{clust}{cls} = nxs;
            clustdims{clust}{cls} = nxs;
            clustemos{clust}{cls} = nxs;
            fullvec{clust}{cls} = nxs;
        end;
    end;
    for clust = 1:length(clusttrack)
        for nx = 1:length(r.facvec{1})
            cls = 1; 
            % find members of alpha clusters:    
            r = load('/data/common4/emotion/CoModAlphaClusts.mat');
            nowdims = clusttrack{clust}(find(clusttrack{clust}(:,1) == nx),2);
            for fc = 1:length(r.facvec) % go through all spec clusters
                for im = 1:length(nowdims)
                    nowfac =nowdims(im) ;
                    if ~isempty(find(ismember(nowfac,r.facvec{fc}{nx})))
                        clustlist{clust}{cls}{nx}(end+1) = nowfac; 
                        clustemos{clust}{cls}{nx}(end+1,:) = clustmeans{clust}(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac),:);
                        fullvec{clust}{cls}{nx}(end+1) = polarvec(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac));
                    end;
                end;
                if ~isempty(clustlist{clust}{cls}{nx})
                    for memb = 1:length(clustlist{clust}{cls}{nx})
                        curridx = find(r.alphafacs{fc}(:,1) == nx & r.alphafacs{fc}(:,2) == abs(clustlist{clust}{cls}{nx}(memb)));
                        clustspecs{clust}{cls}{nx}(end+1,:) = r.alphatempls{fc}(curridx,:);
                    end;
                end;
                cls = cls+1;
            end;
            % find members of beta clusters:    
            r = load('/data/common4/emotion/BetaClusters.mat');
            nowdims = clusttrack{clust}(find(clusttrack{clust}(:,1) == nx),2);
            for fc = 1:length(r.facvec) % go through all spec clusters
                for im = 1:length(nowdims)
                    nowfac =nowdims(im) ;
                    if ~isempty(find(ismember(nowfac,r.facvec{fc}{nx})))
                        clustlist{clust}{cls}{nx}(end+1) = nowfac; 
                        clustemos{clust}{cls}{nx}(end+1,:) = clustmeans{clust}(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac),:);
                        fullvec{clust}{cls}{nx}(end+1) = polarvec(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac));
                    end;
                end;
                if ~isempty(clustlist{clust}{cls}{nx})
                    for memb = 1:length(clustlist{clust}{cls}{nx})
                        curridx = find(r.betafacs{fc}(:,1) == nx & r.betafacs{fc}(:,2) == abs(clustlist{clust}{cls}{nx}(memb)));
                        clustspecs{clust}{cls}{nx}(end+1,:) = r.betatempls{fc}(curridx,:);
                    end;
                end;
                cls = cls+1;
            end;
            
            % find members of gamma clusters:    
            r = load('/data/common4/emotion/GammaClusters.mat');
            nowdims = clusttrack{clust}(find(clusttrack{clust}(:,1) == nx),2);
            for fc = 1:length(r.facvec) % go through all spec clusters
                for im = 1:length(nowdims)
                    nowfac =nowdims(im) ;
                    if ~isempty(find(ismember(nowfac,r.facvec{fc}{nx})))
                        clustlist{clust}{cls}{nx}(end+1) = nowfac; 
                        clustemos{clust}{cls}{nx}(end+1,:) = clustmeans{clust}(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac),:);
                        fullvec{clust}{cls}{nx}(end+1) = polarvec(find(clusttrack{clust}(:,1) == nx&clusttrack{clust}(:,2) == nowfac));
                    end;
                end;
                if ~isempty(clustlist{clust}{cls}{nx})
                    for memb = 1:length(clustlist{clust}{cls}{nx})
                        curridx = find(r.gamafacs{fc}(:,1) == nx & r.gamafacs{fc}(:,2) == abs(clustlist{clust}{cls}{nx}(memb)));
                        clustspecs{clust}{cls}{nx}(end+1,:) = r.gamatempls{fc}(curridx,:);
                    end;
                end;
                cls = cls+1;
            end;
        end;
        fprintf('Cluster %s done.\n',int2str(clust));
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Plot the clusters as means of spectral clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % emomeans -- {nx}(IM,emo)
    emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
    msize = 16; lw = 2.5; fs = 16;
    cols = jet(15);cols(10,:) = [.9 .9 0];
    q = load('/data/common4/emotion/clustfacs.mat'); 
    row = 5; col = 4;
    clustcols(1,:) = [.8 0 0];clustcols(2,:) = [1 0 0];clustcols(3,:) = [1 .4 .4];
    clustcols(4,:) = [0 0 1];clustcols(5,:) = [.4 .4 1];
    clustcols(6,:) = [0 .8 0];clustcols(7,:) = [0 1 0];clustcols(8,:) = [.4 1 .4];

    for cls = 1:length(clustlist)
        figure; pl = 1;      pg = 1;  row = 5; col = 4; zoom = 1.3; 
        sbplot(row,col,[pl pl+2])
        ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
        ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'-','linewidth',2);hold on;
        set(ph,'color',[.5 .5 .5]);
        ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'.','markersize',msize*2);hold on;
        set(ph,'color',[.5 .5 .5]);
        for em = 1:size(clustmeans{cls},2)
            ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
            set(ph,'color',cols(em,:));
            %ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2);
            %set(ph,'color',cols(em,:));
            ph = text(em,-1.5,emos{em});
            set(ph,'rotation',90);
            set(ph,'color',cols(em,:));             
        end;        
        set(gca,'ylim',[-1.6 1.5]); set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
        ph = title([int2str(cls),'Emotion Weights -- ',int2str(size(clustmeans{cls},1)),' IMs across 35 subjects']);
                set(ph,'fontsize',fs);
        pl = pl+4;
        
        for sc = 1:length(clustlist{cls})% sepctral clusters
            if pl > (row * col) - col
                %axcopy;
                %ph = textsc(['Emotion Weights -- ',int2str(size(clustmeans{cls},1)),' IMs across 35 subjects'],'title');
                %textsc(['Cluster ',int2str(cls),' Emotion Space; ',int2str(size(clustmeans{cls},1)),' Total IMs'],'title');
                %set(ph,'fontsize',fs);
               set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                str = ['print /data/common4/emotion/ICAClust',int2str(cls),'pg',int2str(pg),'.eps -depsc']; eval(str);
                 pg = pg+1;
                figure; pl = 1;     row = 6;zoom = 1.5;
            end;            
               
            onesubcls1 = zeros(0,length(r.frs));
            onesubcls2 = zeros(0,length(r.frs));
            submeans1 = zeros(0,15);
            submeans2 = zeros(0,15);
            numsubjs1 = zeros(1,0); 
            numsubjs2 = zeros(1,0);
            clear complist1 wtsmat1 frqmat1 corecomps1 complist2 wtsmat2 frqmat2 corecomps2
            for nx = 1:length(clustlist{cls}{sc})
                if ~isempty(clustlist{cls}{sc}{nx})
                    for fc = 1:length(clustlist{cls}{sc}{nx})
                        if fullvec{cls}{sc}{nx}(fc) > 0
                            onesubcls1(end+1,:) = clustspecs{cls}{sc}{nx}(fc,:);
                            submeans1(end+1,:) = clustemos{cls}{sc}{nx}(fc,:);
                            cc1 = zeros(1,0);
                            tmplist = zeros(1,0);
                            for w = 1:length(gdcomps{nx})
                                if ismember(gdcomps{nx}(w),q.allbigs{nx}{clustlist{cls}{sc}{nx}(fc)})
                                    tmplist(1,end+1) = gdcomps{nx}(w);
                                end;                    
                            end;      
                            cc1(1,end+1) = q.onebig{nx}{clustlist{cls}{sc}{nx}(fc)};
                            wtsmat1{nx}{fc} = q.bigwts{nx}{clustlist{cls}{sc}{nx}(fc)};
                            frqmat1{nx}{fc} = 10;
                            %frqmat1{nx}{fc} = frqcell{clust}{nx}(fc);
                            complist1{nx}{fc} = tmplist;
                            fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                            numsubjs1(1,end+1) = nx;
                            corecomps1{nx} = cc1;
                        else
                            onesubcls2(end+1,:) = clustspecs{cls}{sc}{nx}(fc,:)*-1;
                            submeans2(end+1,:) = clustemos{cls}{sc}{nx}(fc,:);
                            cc2 = zeros(1,0);
                            tmplist = zeros(1,0);
                            for w = 1:length(gdcomps{nx})
                                if ismember(gdcomps{nx}(w),q.allbigs{nx}{clustlist{cls}{sc}{nx}(fc)})
                                    tmplist(1,end+1) = gdcomps{nx}(w);
                                end;                    
                            end;      
                            cc2(1,end+1) = q.onebig{nx}{clustlist{cls}{sc}{nx}(fc)};
                            wtsmat2{nx}{fc} = q.bigwts{nx}{clustlist{cls}{sc}{nx}(fc)};
                            frqmat2{nx}{fc} = 10;
                            %frqmat2{nx}{fc} = frqcell{clust}{nx}(fc);
                            complist2{nx}{fc} = tmplist;
                            fprintf('\n%s  %s\n',int2str(nx),int2str(tmplist));
                            numsubjs2(1,end+1) = nx;
                            corecomps2{nx} = cc2;
                        end;
                        % separate pos and neg orientation IMs and get dipole info
                    end;
                end; 
            end;
            if ~isempty(onesubcls1)
                numsubjs1 = unique(numsubjs1); numsubjs1 = length(numsubjs1);
                frqmat1 = []; wtlims = [1 2];frqlims = [];yon = 1;
                viewnum = 2;
                [allspots] = PlotCrossLinesWts(complist1,fullpaths,'sources1.set',corecomps1,wtsmat1,wtlims,frqmat1,frqlims,row,col,pl,zoom,0,[],viewnum);
                pl = pl+2;               
                sbplot(row,col,pl)
                ph = plot([1:size(emomeans{nx},2)],submeans1,'-');    set(ph,'color',[.4 .4 .4]);hold on;
                for em = 1:size(submeans1,2)
                    ph = plot(em,submeans1(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                end;        
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.6 .6 .6]);hold on;        
                ph = plot([1:size(emomeans{nx},2)],mean(submeans1,1),'-','linewidth',2);hold on;
                set(ph,'color',[.25 .25 .25]);
                ph = plot([1:size(emomeans{nx},2)],mean(submeans1,1),'.','markersize',msize*1.25);hold on;
                set(ph,'color',[.25 .25 .25]);
                set(gca,'ylim',[-1.5 1.5]); set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                ph = title(['Spec Cls ',int2str(sc)]);
                set(ph,'fontsize',fs-2);
                pl = pl+1;            
                sbplot(row,col,pl)
                ph = plot([r.frs(1) r.frs(end)],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
                for memb = 1:size(onesubcls1,1)
                    ph = plot(r.frs,onesubcls1(memb,:),'k-','linewidth',lw);hold on;      
                    set(ph,'color',clustcols(sc,:));
                end;
                set(gca,'xgrid','on');
                set(gca,'xlim',[r.frs(1) r.frs(end)]);set(gca,'xticklabel',[]);pl = pl+1; 
                 ph = title([int2str(numsubjs1),' sbjs, ',int2str(size(submeans1,1)),' IMs']);               
                set(ph,'fontsize',fs-2);
            end;            
            if ~isempty(onesubcls2)
                numsubjs2 = unique(numsubjs2); numsubjs2 = length(numsubjs2);
                frqmat2 = []; wtlims = [1 2];frqlims = [];yon = 1;
                viewnum = 2;
                [allspots] = PlotCrossLinesWts(complist2,fullpaths,'sources1.set',corecomps2,wtsmat2,wtlims,frqmat2,frqlims,row,col,pl,zoom,0,[],viewnum);
                pl = pl+2;               
                sbplot(row,col,pl)
                ph = plot([1:size(emomeans{nx},2)],submeans2,'-');    set(ph,'color',[.4 .4 .4]);hold on;
                for em = 1:size(submeans2,2)
                    ph = plot(em,submeans2(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                end;        
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.6 .6 .6]);hold on;        
                ph = plot([1:size(emomeans{nx},2)],mean(submeans2,1),'-','linewidth',2);hold on;
                set(ph,'color',[.25 .25 .25]);
                ph = plot([1:size(emomeans{nx},2)],mean(submeans2,1),'.','markersize',msize*1.25);hold on;
                set(ph,'color',[.25 .25 .25]);
                set(gca,'ylim',[-1.5 1.5]); set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                ph = title(['Spec Cls ',int2str(sc)]);
                set(ph,'fontsize',fs-2);
                pl = pl+1;            
                sbplot(row,col,pl)
                ph = plot([r.frs(1) r.frs(end)],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
                for memb = 1:size(onesubcls2,1)
                    ph = plot(r.frs,onesubcls2(memb,:),'k-','linewidth',lw);     hold on;      
                    set(ph,'color',clustcols(sc,:));
                end;
                set(gca,'xgrid','on');
                set(gca,'xlim',[r.frs(1) r.frs(end)]);set(gca,'xticklabel',[]);pl = pl+1; 
                ph = title([int2str(numsubjs2),' sbjs, ',int2str(size(submeans2,1)),' IMs']);               
                set(ph,'fontsize',fs-2);
            end;            
        end;
        axcopy;
        %ph=textsc(['Emotion Weighting -- ',int2str(size(clustmeans{cls},1)),' IMs across subjects'],'title');
        %textsc(['Cluster ',int2str(cls),' Emotion Space; ',int2str(size(clustmeans{cls},1)),' Total IMs'],'title');
        %set(ph,'fontsize',fs);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        str = ['print /data/common4/emotion/ICAClust',int2str(cls),'pg',int2str(pg),'.eps -depsc']; eval(str);
        
    end;
    

