% clusters components across subjects by 3D dipole locations
% uses bilateral dipole locations, so most have zeros in second half.
%
%          [clustcps, keeptrack ] = ClusterDipoles(setname,fullpaths,gdcomps,nclusts,stdcut,pk);
%
% Inputs:
% setname:  [string or cell array] name of data file from which to retrieve
% dipole locations; If a cell array, the dataset names are different, but
% must correspond to each cell in 'fullpaths'
% fullpaths: full paths of directories where subject directories are found. [string]
% gdcomps: selected components from each subject to cluste. [cell array of vectors]. Note: the length of datpaths and gdcomps must be the same and corresponding indexes must refer to the same subject
% nclusts: number of kmeans clusters to create
% stdcut: number of standard deviations from the mean a component must be to be excluded from a cluster
% pka -- ['p'|'k'|'a'} p for pdist, k for kmeans, 'a' for affinity clustering
% seed -- [x,y,z] coordinates for which to find single cluster of nearby ICs from all subjs
%
% Outputs:
% clustcps: cell array of 'kclusts' cells each containing vectors of components for each subject that below in that cluster. 
% keeptrack: matrix of dimensions (number of components) X 2. The first column is subject number; the second is compoent number.

function [clustcps,ninclusts,outliers,noutliers] = ClusterDipoles(setname,fullpaths,gdcomps,nclusts,stdcut,pka,seed);
    
    
    clustcps=[]; 
    ninclusts=[]; 
    outliers=[]; 
    noutliers=[]; 
    
    allloc = zeros(0,3);pl = 1;clear keeptrack
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
             if iscell(setname)
                EEG = pop_loadset(setname{nx} ,fullpaths{nx});
            else
                EEG = pop_loadset(setname,fullpaths{nx});
            end;
            for cp = 1:length(gdcomps{nx})
                subjwinv = EEG.dipfit.model(gdcomps{nx}(cp)).posxyz;
                if size(subjwinv,1) > 1
                    subjwinv(1,end+1:end+size(subjwinv,2)) = subjwinv(2,:); 
                    subjwinv(2,:) = [];
                else
                    subjwinv(1,end+1:end+3) = [0 0 0];
                end;
                subjwinv = subjwinv(1,1:3);
                allloc(end+1,:) = subjwinv(1,:);
                keeptrack(pl,:) = [nx,gdcomps{nx}(cp)]; pl = pl+1;
            end;
        end;
    end;
    
    if ~exist('seed')
      seed = [];
    end;
    if ~isempty(seed) % cluster based on a single seed location instead
      tmpcomps = cell(1,length(gdcomps));   
      for x=1:size(allloc,1)
        closeby = 0;
        closeby = closeby + (allloc(x,1) > seed(1)-stdcut & allloc(x,1) < seed(1)+stdcut);
        closeby = closeby + (allloc(x,2) > seed(2)-stdcut & allloc(x,2) < seed(2)+stdcut);
        closeby = closeby + (allloc(x,3) > seed(3)-stdcut & allloc(x,3) < seed(3)+stdcut);
        if closeby == 3 % close enough
          tmpcomps{keeptrack(x,1)} = [tmpcomps{keeptrack(x,1)},keeptrack(x,2)];
        end;
      end;
      clustcps{1} = tmpcomps;
    else
      
    if pka == 'k'
        repeat = 1;manyrep = 0;
        for x=1:size(allloc,1)
            allloc(x,:) = allloc(x,:)/std(allloc(x,:));
        end;
        while repeat == 1
            [kout, C,sumd, allsums] = kmeans(allloc,nclusts,'replicates',25);
            
            clear clustcps
            for clust = 1:max(kout)
                cp = find(kout == clust); 
                cp(find(abs(zscore(allsums(cp,clust)))> stdcut)) = [];
                outliers{clust} = find(abs(zscore(allsums(cp,clust)))> stdcut);
                manyout(1,clust) = length(outliers{clust});
                cpidx{clust} = cp;
                relcomps = keeptrack(cp,:); 
                cpoi = cell(1,length(gdcomps));   
                for w = 1:length(cp)
                    cpoi{relcomps(w,1)}(end+1) = relcomps(w,2);
                end;
                clustcps{clust} = cpoi;
            end;
            if isempty(find(manyout > 3))
                repeat = 0;manyrep = manyrep+1;
            end;
        end;
    elseif pka == 'p' % else do pdist       
                      %alldist = pdist(allloc, 'cosine');    
                      %links = linkage(alldist,'complete');
        alldist = pdist(allloc, 'euclidean');    
        links = linkage(alldist,'ward');        
        figure;[hnd,idx,perm]=  dendrogram(links,nclusts);
        
        for cls = 1:nclusts
          incidx = find(idx == cls);
          % find outliers components:---------------
          tmplocs = allloc(incidx,:);
          mnclust = mean(tmplocs); 
          stdclust = std(tmplocs);
          delidx = [];
          if ~isempty(stdcut)
            stdclust = stdclust*stdcut; % mult by input std cut off
            for c = 1:size(tmplocs,1)
              dif = abs(mnclust - tmplocs(c,:));
              if ~isempty(find(dif > stdclust)) % one is over
                delidx = [delidx c];
              end;
            end;
            %  save outliers separately ----
            noutliers(1,cls) = length(delidx);
            delcomps = keeptrack(incidx(delidx),:); 
            stdouts = cell(1,length(gdcomps));
            for nx = 1:length(gdcomps)
              stdouts{nx} = delcomps(find(delcomps(:,1) == nx),2)';
            end;
            outliers{cls} = stdouts; % keep outliers in separate cell array
          end;
          % delete outlier components and save pruned clusters ---------
          incidx(delidx) = []; % delete index out of range
          ninclusts(1,cls) = length(incidx);
          relcomps = keeptrack(incidx,:); 
          newcomps = cell(1,length(gdcomps));
          for nx = 1:length(gdcomps)
            newcomps{nx} = relcomps(find(relcomps(:,1) == nx),2)';
          end;
          clustcps{cls} = newcomps;              
        end; 
    elseif pka == 'a' % else do affinity         
        %addpath('/home/nima/tools/clustering/affinity/') % for apclusterK
        if isempty(nclusts)
            M=size(allloc,1)*size(allloc,1)-size(allloc,1); s=zeros(M,3); % Make ALL N^2-N similarities
            j=1;
            for i=1:size(allloc,1)
                for k=[1:i-1,i+1:size(allloc,1)]
                    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((allloc(i,:)-allloc(k,:)).^2);
                    j=j+1;
                end;
                if ismember(i,[100:100:size(allloc,1)])
                    fprintf('\n%s of %s done..',int2str(i),int2str(size(allloc,1)));
                end;
            end;
            p=median(s(:,3)); % Set preference to median similarity
            if size(s,1) > 3000
                [idx,netsim,dpsim,expref]=apcluster(s,p,'convits',25,'sparse');
            else                
                [idx,netsim,dpsim,expref]=apcluster(s,p,'convits',25);
            end;
        else % otherwise use apclusterK
            alldist = pdist(allloc, 'euclidean');
            alldist=squareform(alldist);
            [idx,netsim,dpsim,expref]=apclusterK(-alldist,nclusts);
        end;
        examplar = unique(idx);
        clustidx = zeros(1,length(examplar));
        for i=1:length(examplar)
            comps = find(idx == examplar(i));
            clustidx(comps) = i;
        end;
        clear clustcps
        for clust = 1:max(clustidx)
            ids = find(clustidx == clust);
            tmpclusts = keeptrack(ids,:);
            tmplocs = allloc(ids,:);
            mnclust = mean(tmplocs); 
            stdclust = std(tmplocs);
            for nx = 1:length(gdcomps)
              onesubj = find(tmpclusts(:,1) == nx);
              clustcps{clust}{nx} = tmpclusts(onesubj,2);
            end;
            
            delidx = [];
            if ~isempty(stdcut)
              stdclust = stdclust*stdcut; % mult by input std cut off
              for c = 1:size(tmplocs,1)
                dif = abs(mnclust - tmplocs(c,:));
                if ~isempty(find(dif > stdclust)) % one is over
                  delidx = [delidx c];
                end;
              end;
              %  save outliers separately ----
              noutliers(1,clust) = length(delidx);
              delcomps = keeptrack(ids(delidx),:); 
              stdouts = cell(1,length(gdcomps));
              for nx = 1:length(gdcomps)
                stdouts{nx} = delcomps(find(delcomps(:,1) == nx),2)';
              end;
              outliers{clust} = stdouts; % keep outliers in separate cell array
            end;
            % delete outlier components and save pruned clusters ---------
            ids(delidx) = []; % delete index out of range
            ninclusts(1,clust) = length(ids);
            tmpclusts = keeptrack(ids,:);
            newcomps = cell(1,length(gdcomps));
            for nx = 1:length(gdcomps)
              onesubj = find(tmpclusts(:,1) == nx);
              newcomps{nx} = tmpclusts(onesubj,2)';
            end;
            clustcps{clust} = newcomps;            
        end;            
    end;
    
    end;