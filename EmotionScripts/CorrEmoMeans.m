% takes the mean weights for emotions from all subjects and clusters by correlation/distance?
%
% [clustmeans, clusttrack, collmeans] = CorrEmoMeans(savedat,whichfacs,gdcomps,emos,emomeans,fullpaths,corrcut,nclusts,plotall,clustmeth)
%
% INPUTS:
% 
% emomeans -- {subj}(IM,em)
% clustmeth -- 'k' for kmeans, 'd' for pdist/linkage, 'i' for ICA clustering, 'p' for PCA clustering


function [clustmeans, clusttrack, collmeans, keeptrack, polarvec] = CorrEmoMeans(savedat,gdcomps,emos,emomeans,fullpaths,corrcut,nclusts,plotall,clustmeth)
    
    cols = jet(15);cols(10,:) = [.9 .9 0];
    
    clear corr indx indy
    for nx1 = 1:length(emomeans) - 1
      if ~isempty(emomeans{nx1})
        for nx2 = nx1 + 1:length(emomeans)
          if ~isempty(emomeans{nx2})
            [corr(nx1,nx2,:),indx(nx1,nx2,:),indy(nx1,nx2,:),corrs] = matcorr(emomeans{nx1},emomeans{nx2});
          end;
        end;
      end;
    end;
    
    % collect highly correlated dims from emomeans
    imidx = cell(35,35);
    dimidx = cell(1,length(emomeans));
    for nx = 1:size(corr,1)
        keepcuts1 = zeros(1,0);
        for nxx = nx+1:length(emomeans)
        keepcuts2 = zeros(1,0);
            if ~isempty(find(abs(corr(nx,nxx,:)) > corrcut))
                addons = find(abs(corr(nx,nxx,:)) > corrcut);
                keepcuts1(1,end+1:end+length(addons)) = squeeze(indx(nx,nxx,find(abs(corr(nx,nxx,:)) > corrcut)));
                keepcuts2(1,end+1:end+length(addons)) = squeeze(indy(nx,nxx,find(abs(corr(nx,nxx,:)) > corrcut)));
            end;
            %dimidx{nxx}(1,end+1:end+length(unique(keepcuts2))) = unique(keepcuts2);
            dimidx{nxx}(1,end+1:end+length(keepcuts2)) = keepcuts2;
        end;
        dimidx{nx}(1,end+1:end+length(keepcuts1)) = keepcuts1;
        %dimidx{nx}(1,end+1:end+length(unique(keepcuts1))) = unique(keepcuts1);
    end;
    
    minim = 4;%8
    collmeans = zeros(0,15); keeptrack = zeros(0,2);
    for nx = 1:length(dimidx)
      if ~isempty(dimidx{nx})
        %dimidx{nx} = unique(dimidx{nx});
        for dm = 1:length(unique(dimidx{nx}))
          if length(find(dimidx{nx} == dm)) > minim
            keeptrack(end+1,:) = [nx, dm];
            collmeans(end+1,:) = emomeans{nx}(dm,:);            
          end;
        end;
      end;
    end;            
    
    if clustmeth == 'k'
        % Now cluster collmeans by kmeans with a single reseed (stable after that):
        stds= 0;
        [idx, C, SUMD, dists]  = kmeans(collmeans,nclusts,'distance','correlation','replicates',25);
        clear seedmat
        for cls = 1:max(idx)
            onecls = find(idx == cls);
            % do std cutoff for all clusts
            currdists = dists(onecls,cls);
            cutval = mean(currdists) + std(currdists)*stds;
            onecls(find(currdists > cutval)) = [];
            seedmat(cls,:) = mean(collmeans(onecls,:),1);
        end;
        [idx, C, SUMD, dists]  = kmeans(collmeans,nclusts,'start',seedmat);

        stds= 1.5; % new, less strict standard dev cutoff
        figure;  msize = 10;
        row=round(sqrt(max(idx))); 
        col=ceil(sqrt(max(idx)));  clear keepcut keepouts clustmeans clusttrack
        for cls = 1:max(idx)
            onecls = find(idx == cls);
            % do std cutoff for all clusts
            currdists = dists(onecls,cls);
            cutval = mean(currdists) + std(currdists)*stds;
            onecls(find(currdists > cutval)) = [];
            clustmeans{cls} = collmeans(onecls,:); 
            clusttrack{cls} = keeptrack(onecls,:);
            if ~isempty(clustmeans{cls})            
                sbplot(row,col,cls)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;            
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',cols(em,:));
                end;        
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(cls),'; ',int2str(size(clustmeans{cls},1)),' members']);
            end;
        end;
        ph=textsc(['Emotion Dimension KMEANS Clusters; Corr Cut: ',num2str(corrcut)],'title');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        
    elseif clustmeth == 'i'
        [weights,sphere,compvars,bias,signs,lrates,activations] = runica(collmeans','pca',3,'stop',1e-7);
        ws = weights*sphere; winv = pinv(ws);
        figure; msize = 10;polarvec = zeros(1,0);
        row=round(sqrt(size(winv,2))); 
        col=ceil(sqrt(size(winv,2)));  clear clustmeans clusttrack
        for cls = 1:size(winv,2)
            onecls = find(zscore(activations(cls,:)) > 1 );
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls));
            clustmeans{cls} = collmeans(onecls,:); 
            clusttrack{cls} = keeptrack(onecls,:);
            onecls = find(zscore(activations(cls,:)) < -1 );
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls))*-1;
            clustmeans{cls} = [clustmeans{cls};-1*collmeans(onecls,:)]; 
            clusttrack{cls} = [clusttrack{cls};keeptrack(onecls,:)]; 

            if ~isempty(clustmeans{cls})            
                sbplot(row,col,cls)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;            
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',[.6 .6 .6]);
                end;  
                ph = plot(mean(clustmeans{cls},1),'-');set(ph,'color',[.6 .6 .6]);
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(cls),'; ',int2str(size(clustmeans{cls},1)),' members']);
            end;
        end;
        ph=textsc(['Emotion Dim ICA Clusters; Corr Cut: ',num2str(corrcut),'; minIM: ',int2str(minim)],'title');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        
    elseif clustmeth == 'p'
        [pc,eigvec,sv] = runpca(collmeans',nclusts); % use pca to get IM weights
        activations = pc; winv = eigvec;
        figure; msize = 10;polarvec = zeros(1,0);
        row=round(sqrt(size(winv,2))); 
        col=ceil(sqrt(size(winv,2)))+1;  clear clustmeans clusttrack
        for cls = 1:size(winv,2)
            onecls = find(zscore(activations(cls,:)) > 1 );
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls));
            clustmeans{cls} = collmeans(onecls,:); 
            clusttrack{cls} = keeptrack(onecls,:);
            onecls = find(zscore(activations(cls,:)) < -1 );
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls))*-1;
            clustmeans{cls} = [clustmeans{cls};-1*collmeans(onecls,:)]; 
            clusttrack{cls} = [clusttrack{cls};keeptrack(onecls,:)]; 
            %%%  Find combination wts: dims 1 and 2 !! this is only true for PCA decomp, not MD (dim 2 opposite)
            onecls = find(zscore(activations(1,:)) < -.5& zscore(activations(2,:)) > .5 ); % grief/sad/comp
            clustmeans{nclusts+1} = collmeans(onecls,:); 
            clusttrack{nclusts+1} = keeptrack(onecls,:); 
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls));

            %%%  Find combination wts: dims 1 and 2
            onecls = find(zscore(activations(1,:)) < -.5 & zscore(activations(2,:)) < -.5 );% relief/content
            clustmeans{nclusts+2} = collmeans(onecls,:); 
            clusttrack{nclusts+2} = keeptrack(onecls,:); 
            polarvec(end+1:end+length(onecls)) = ones(1,length(onecls));

            if ~isempty(clustmeans{cls})            
                sbplot(row,col,cls)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;            
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',[.6 .6 .6]);
                end;  
                ph = plot(mean(clustmeans{cls},1),'-');set(ph,'color',[.6 .6 .6]);
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(cls),'; ',int2str(size(clustmeans{cls},1)),' members']);
            end;
        end;
        
        sbplot(row,col,cls+1)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;            
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{nclusts+1}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{nclusts+1}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',[.6 .6 .6]);
                end;  
                ph = plot(mean(clustmeans{nclusts+1},1),'-');set(ph,'color',[.6 .6 .6]);
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(nclusts+1),'; ',int2str(size(clustmeans{nclusts+1},1)),' members']);
        
        sbplot(row,col,cls+2)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;            
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{nclusts+2}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{nclusts+2}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',[.6 .6 .6]);
                end;  
                ph = plot(mean(clustmeans{nclusts+2},1),'-');set(ph,'color',[.6 .6 .6]);
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(nclusts+2),'; ',int2str(size(clustmeans{nclusts+2},1)),' members']);
        
        sbplot(row,col,cls+3)
        c1 = 1; c2 = 2; c3 = 3;
        for e = 1:size(winv,1)
            ph=plot3(winv(e,c1),winv(e,c2),winv(e,c3),'.');hold on;
            set(ph,'markersize',25);set(ph,'color',cols(e,:));
            ph = text(winv(e,c1),winv(e,c2),winv(e,c3),emos{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        end;
        zl = get(gca,'zlim');
        for e = 1:size(winv,1)
            pl =plot3([winv(e,c1) winv(e,c1)],[winv(e,c2) winv(e,c2)],[zl(1)  winv(e,c3)]);
            set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
        end;
        set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
        xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);zlabel(['Dim ',int2str(c3)]);
        %ph=textsc(['Emotion Dim PCA Clusters; Corr Cut: ',num2str(corrcut),'; minIM: ',int2str(minim)],'title');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        
        
    elseif clustmeth == 'd'
        %%  use dist/linkage/dendrogram
        
        alldist = pdist(collmeans, 'correlation'); % euc better than seuc
        links = linkage(alldist,'complete');%ward 
        figure;[hnd,idx,perm]=  dendrogram(links,nclusts);% 100 
        figure;  msize = 10;
        cols = jet(15);cols(10,:) = [.9 .9 0];
        row=round(sqrt(max(idx))); 
        col=ceil(sqrt(max(idx)));  clear keepcut keepouts clustmeans clusttrack
        for cls = 1:max(idx)
            onecls = find(idx == cls);
            clustmeans{cls} = collmeans(onecls,:); 
            clusttrack{cls} = keeptrack(onecls,:);
            % determine which members are within ? std of mean
            for iter = 1:5
                if ~isempty(clustmeans{cls})
                    centr = mean(clustmeans{cls},1); clear onedist
                    for mem = 1:size(clustmeans{cls},1)
                        onedist(1,mem) = pdist([centr; clustmeans{cls}(mem,:)], 'correlation');
                    end;
                    cutoff = mean(onedist);
                    
                    keepcut(cls,iter) = cutoff;
                    if cutoff < .65
                        clear outliers
                        %if ~isempty(find(onedist < cutoff))
                        %outliers = find(onedist < cutoff);
                        [x y] = min(onedist);
                        outliers = y;
                        clustmeans{cls}(outliers,:) = [];
                        clusttrack{cls}(outliers,:) = [];  
                        keepouts(cls,iter) = length(outliers);              
                        %end;
                    end;
                end;
            end;        
            if ~isempty(clustmeans{cls})
                
                sbplot(row,col,cls)
                ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;
                
                for em = 1:size(collmeans,2)
                    ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
                    set(ph,'color',cols(em,:));
                    ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2.5);
                    set(ph,'color',cols(em,:));
                end;        
                set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
                title(['Cls ',int2str(cls),'; ',int2str(size(clustmeans{cls},1)),' members']);
            end;
        end;
        ph=textsc(['Emotion Dimension Clusters across all subjects; Corr Cut: ',num2str(corrcut)],'title');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    end;
    
    if plotall == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Plot the actual factors for each subj/dim in a cluster
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for cls = 1:length(clustmeans)
            figure; pl = 1;      pg = 0;    fc = 0;
            row = 6; col = 6; zoom = 1.4;
            q = load('/data/common4/emotion/clustfacs.mat'); 
            cols = jet(15);cols(10,:) = [.9 .9 0];
            faccols = lines(12); 
            for nxx = 1:size(clustmeans{cls},1)
                if pl > row * col
                    pg = pg+1;
                    textsc(['Cluster ',int2str(cls)],'title'); 
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    %str = ['print /data/common4/emotion/EmoDimClusts',int2str(cls),'-',int2str(pg),'.jpg -djpeg']; eval(str);
                    %close
                    figure; pl = 1;   fc=0;  
                end;                
                nx = clusttrack{cls}(nxx,1);
                EEG = pop_loadset('sources1.set', fullpaths{nx});
                s = load([fullpaths{nx},savedat,'.mat']);  
                sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
                wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
                icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
                ws = wts*sph;    activations = ws*icamatall;   clear wts sph ws
                
                maxsz = 25;
                %for fac = 1:length(whichfacs{nx}{clusttrack{cls}(nxx,2)})
                if pl > row * col
                    pg = pg+1;
                    textsc(['Cluster ',int2str(cls)],'title'); 
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    %str = ['print /data/common4/emotion/EmoDimClusts',int2str(cls),'-',int2str(pg),'.jpg -djpeg']; eval(str);
                    %close
                    figure; pl = 1;   fc=0;    
                end;                
                fc = fc+1;
                sbplot(row,col,pl);            
                cmpcols = lines(length(gdcomps{nx}));
                tp = clusttrack{cls}(nxx,2);               
                for cp = 1:length(gdcomps{nx})
                    rcp =cp;
                    if tp < 0
                        ph = plot(s.freqs,activations(abs(tp),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1,'linewidth',1.5); 
                    else
                        ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',1.5); 
                    end;                       
                    hold on; %set(ph,'color',cmpcols(cp,:));                    
                    set(ph,'color',faccols(fc,:));
                end;pl = pl+1;
                set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                title(['Sb:',int2str(nx),' -- IM:',int2str(clusttrack{cls}(nxx,2)) ]);
                %title(['IM:',int2str(whichfacs{nx}{clusttrack{cls}(nxx,2)}(fac)),'(',int2str(fac),') Sb:',int2str(nx)]);
                set(gca,'xlim',[s.freqs(1) s.freqs(end)]); 
                set(gca,'ylim',[-5.5 9.5]);set(gca,'xticklabel',[]);
                set(gca,'box','off');
                set(gca,'xgrid','on');
                set(gca,'ticklength',[.05 .05]);
                %if fac == length(whichfacs{nx}{clusttrack{cls}(nxx,2)})
                %    xlabel('Frequency (Hz)');
                %    ylabel('Relative Power');
                %end; 
                facfac = maxsz - 5;
                %if ~isempty(q.allbigs{nx})
                sbplot(row,col,pl);  
                mydipplot(EEG.dipfit.model(cell2mat(q.allbigs{nx}(clusttrack{cls}(nxx,2)))),'image','mri','gui','off','dipolelength',0,'dipolesize',facfac,'normlen','on','spheres','on','color',{faccols(fc,:)},'projlines','on','projimg','on','coordformat','spherical');                hold on;
                %mydipplot(EEG.dipfit.model(cell2mat(q.allbigs{nx}(abs(whichfacs{nx}{clusttrack{cls}(nxx,2)}(fac))))),'image','mri','gui','off','dipolelength',0,'dipolesize',facfac,'normlen','on','spheres','on','color',{faccols(fc,:)},'projlines','on','projimg','on','coordformat','spherical');                hold on;
                view(60,30); pl = pl+1;
                camzoom(zoom); 
                %else
                %    pl = pl+1;
                %end;
                sbplot(row,col,pl);
                ph = plot([0 length(emos)+1],[0 0],'k-','linewidth',1);set(ph,'color',[.3 .3 .3]);hold on;
                ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'-','linewidth',2);hold on;
                set(ph,'color',[.5 .5 .5]);
                ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'.','markersize',msize);hold on;
                set(ph,'color',[.5 .5 .5]);
                for e = 1:length(emos)
                    ph = plot(e,emomeans{nx}(clusttrack{cls}(nxx,2),e),'k.','markersize',msize*1.5);
                    set(ph,'color',cols(e,:)); hold on;
                end;
               set(gca,'xticklabel',[]); set(gca,'ylim',[-1.5 1.5]);set(gca,'xlim',[0 16]);
                pl = pl+1;title(['Sb:',int2str(nx),' (',fullpaths{nx}(end-4:end-1),')']);
                set(gcf,'color','w');
            end;
            %end;    
            sbplot(row,col,[pl pl+2])
            ph = plot([0 16],[0 0],'k-'); set(ph,'color',[.4 .4 .4]);hold on;        
            ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'-','linewidth',2);hold on;
            set(ph,'color',[.5 .5 .5]);
            ph = plot([1:size(emomeans{nx},2)],mean(clustmeans{cls},1),'.','markersize',msize);hold on;
            set(ph,'color',[.5 .5 .5]);
            for em = 1:size(collmeans,2)
                ph = plot(em,clustmeans{cls}(:,em),'.','markersize',msize);hold on;
                set(ph,'color',cols(em,:));
                ph = plot(em,mean(clustmeans{cls}(:,em),1),'.','markersize',msize*2.5);
                set(ph,'color',cols(em,:));
                ph = text(em,-1.5,emos{em});
                set(ph,'rotation',90);
                set(ph,'color',cols(em,:));             
            end;        
            set(gca,'ylim',[-1.6 1.5]); set(gca,'xlim',[0 16]); set(gca,'xticklabel',[]);
            title(['Cls ',int2str(cls),'; ',int2str(size(clustmeans{cls},1)),' members']);
            
            pg = pg+1;
            textsc(['Cluster ',int2str(cls)],'title'); 
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            str = ['print /data/common4/emotion/EmoDimICAClusts',int2str(cls),'-',int2str(pg),'.eps -depsc']; eval(str);
            close
        end;
    end;    
    
