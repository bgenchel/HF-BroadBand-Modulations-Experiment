% takes mean or decile IM weights from specified subjs/IMs and uses ICA or multi-dimensional scaling to find emo space across subjects
%
% [fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,orivecs,numdims,subjlist,useims,meth,mvon,ttl_in,indiv,elimbyrating,permtest)
%
% emomeans -- [cell array] {subj}(IM,emo) or deciles: {subj}(IM,emo,dec)
% orivecs -- [cell array] vector of 1's and -1's that will orient the 
%            values in 'emomeans' (ie, these cell arrays must be the same size.
% subjlist -- [vector] of subject indices to use in analysis
% useims -- [cell array] of IM indices to use from each subject (ideally, these
%           are IMs that were clustered into modulator categories.
% meth -- ['ica', 'mds', 'pca' or 'lda] for ICA, multi-dimensional scaling, PCA and linear discriminant analysis, respectively
% numdims -- [integer] number of dimensions to calculate for between subj solution. [] will choose heuristic or for ICA will pca to rank of matrix.
% mvon -- [1|0] if 1, will make a movie of MD-scaled space from various angles
% ttl_in -- [string] to title the still picture of 3D MD scaling figure
% indiv -- [1 or 0] if 1, will use straight mean weights, if 0 (for across subject decomposition)
%                   will normalize each IM mean weights by dividing by the standard dev of that IM.
% elimbyrating -- ['elim' or 'keepall'] if 'elim', will use subject emotion ratings and eliminate
%                 emotion median values that correspond to poorly embodied emotions. The values 
%                 will then be filled in by the mean of all other IMs in the matrix. Therefore,
%                 this is only appropriate for homogeneous inputs.
% permtest -- ['perm' or 'noperm'] if 'perm', will shuffle emotion dimension before 
%             decomposition for collection of bootstrap limits.
%
% OUTPUT:
% keeptrack -- [ndims x 2 or 3 columns] if emo mean or median: 2 columns, if deciles: 3 columns
%              column 2 is the IM index (not all need to be used), column 3 is the decile of IM
% nsteps -- [number] of steps that ICA (runica) took.

function [fullmd,fullwts,keeptrack,collmeans,nsteps] = EmoSpace(emomeans,orivecs,numdims,subjlist,useims,meth,mvon,ttl_in,indiv,elimbyrating,permtest)

    
    
    r = load('/data/common1/emotion/SubjRatings.mat');% load just in case 'elim'
    % 6.4838 +- 2.4025, (mean +- std)
    
    nsteps = [];
    if ~exist('indiv')
        fprintf('\nPlease specify whether this is a decomposition of a single subject or multiple (indiv = 1 or 0, respectively)\n');
        return;
    end;
    
    cols = jet(15);cols(10,:) = [.9 .9 0]; cols(end+1,:) = [.7 .7 .7]; cols(end+1,:) = [.3 .3 .3];
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited', '  prebase', '  postbase'};

    if length(size(emomeans{1})) > 2
        keeptrack = zeros(0,3);
    else
        keeptrack = zeros(0,2);
    end;
    % take instead IMs that were specified 
    collmeans = []; 
    for nx = 1:length(useims)
        if ~isempty(useims{nx}) & ~isempty(find(ismember(nx,subjlist)))
            if strcmp(elimbyrating,'elim')
                elimemos = find(r.ratings{nx} < 5); % 5 is arbitrary, but I think reasonable
                emomeans{nx}(:,elimemos) = 0; % set to zero
            end;            
            for dm = 1:length(useims{nx})
                if length(size(emomeans{nx})) > 2
                    tmpmeans =  squeeze(emomeans{nx}(abs(useims{nx}(dm)),:,:)); %  emos x deciles 
                    collmeans = [collmeans, tmpmeans]; % emos x IMs*ndeciles
                    keeptrack(end+1:end+size(tmpmeans,2),:) = [repmat(nx,[size(tmpmeans,2) 1]), repmat(useims{nx}(dm),[size(tmpmeans,2) 1]), [1:9]'];
                else   
                    if indiv == 1 % if decomp of a single subject
                        collmeans = [collmeans, emomeans{nx}(abs(useims{nx}(dm)),:)'];% single subject emo space  
                    else % or else you have to normalize
                        %collmeans = [collmeans, emomeans{nx}(abs(useims{nx}(dm)),:)'/std(emomeans{nx}(abs(useims{nx}(dm)),:))];% emos x IMs,/std  (across subject) 
                        collmeans = [collmeans, zscore(emomeans{nx}(abs(useims{nx}(dm)),:))'];% emos x zscore(IMs) (across subject) 
                        if ~isempty(orivecs)% orient by 'onebig' orientation
                            collmeans(:,end) = collmeans(:,end)*orivecs{nx}(dm);
                        end;
                    end;
                    keeptrack(end+1,:) = [nx, useims{nx}(dm)];
                end;
            end;
        end;
    end;
    
    if strcmp(elimbyrating,'elim')
        if indiv ~=1 % if multiple subjects
            for e = 1:size(collmeans,1)
                if ~isempty(find(collmeans(e,:) == 0))
                    othermean = mean(collmeans(e,find(collmeans(e,:)))); % mean of all other values
                    zeroidx = find(collmeans(e,find(collmeans(e,:)==0)));
                    for z = 1:length(zeroidx)
                        collmeans(e,zeroidx(z)) = othermean;
                    end;
                end;
            end;
        else % one subject
            delrows = [];
            for e = 1:size(collmeans,1)
                if ~isempty(find(collmeans(e,:) == 0))
                    delrows = [delrows e]; % delete bad emos
                end;
            end;
            collmeans(delrows,:) = [];
            cols(delrows,:) = []; % fix plotting colors accordingly
            ee=1;
            for e = 1:length(emo2)
                if ~ismember(e,delrows)
                    newemos{ee} = emo2{e};ee=ee+1;
                end;
            end;
            emo2 = newemos;
        end;
    end;
    
    if strcmp(permtest,'perm')
        collmeans = collmeans(randperm(size(collmeans,1)),:);        
    end;
    
    
    if strcmp(meth,'ica')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Calculate ICA decomposition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        if isempty(numdims)
            numdims = rank(collmeans)-1;
        end;
        fdims = numdims;
        %[wts,sph] = binica( collmeans,'pca',numdims,'extended',1,'stop',1e-7,'maxsteps',2000); 
        [wts,sph,compvars,bias,signs,lrates,acts] = runica(collmeans,'pca',numdims,'extended',1,'stop',1e-7,'maxsteps',6000);
        nsteps = length(lrates);
        ws = wts*sph; winv = pinv(ws); acts = ws*collmeans;
        fullmd = winv;  fullwts = acts;
        %fullwts = winv'; fullmd = acts';% for a IM x emotions matrix
    elseif strcmp(meth,'mds')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Calculate MD scaling   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if isempty(numdims)
            numdims = 3;
        end;
        try 
            dd = pdist(collmeans, 'correlation') ;
        catch % sometimes correlation doesn't work
            %dd = pdist(collmeans, 'euclidean') ;
        end            
        % here 'fullwts' is the 'stress'
        % if too few dims will return zeros
        try
            [fullmd,fullwts] = mdscale(dd,numdims,'replicates',25); 
        catch
            try
                if numdims > 1
                    numdims = numdims - 1;
                    [fullmd,fullwts] = mdscale(dd,numdims,'replicates',25); % here 'fullwts' is the 'stress'
                else
                    fullmd = zeros(length(emo2),numdims);fullwts = 0;% last resort
                end;
            catch
                try
                    if numdims > 1
                        numdims = numdims - 1;
                        [fullmd,fullwts] = mdscale(dd,numdims,'replicates',25); % here 'fullwts' is the 'stress'     
                    else
                        fullmd = zeros(length(emo2),numdims);fullwts = 0;% last resort
                    end;
                catch             
                    fullmd = zeros(length(emo2),numdims);fullwts = 0;% last resort
                end;
            end;
        end;
        for x=size(fullmd,2):-1:2
            if isempty(find(fullmd(:,x)))
                fullmd(:,x) = [];
            end;
        end;
        fdims = size(fullmd,2);
        % final number of dims to plot
        
        %eigvals = fullwts; % not returned
        %clear fullwts
        %fullwts = pinv(fullmd)*collmeans; %  recover IM weightings
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(meth,'lda')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Calculate Linear discriminant analsyis %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        [ldadims,P,STATS] = manova1(collmeans,emo2,.01) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(meth,'pdist')
        if isempty(numdims)
            fprintf('/npdist clustering option requires a set number of dimensions to return.\n')
            return;
        end;
        fdims = numdims;collmeans = collmeans';
        alldist = pdist(collmeans, 'correlation'); % euc better than seuc
        links = linkage(alldist,'complete');       
        figure;[hnd,idx,perm]=  dendrogram(links,fdims);
        for cls = 1:max(idx)
            onecls = find(idx == cls);
            fullwts{cls} = collmeans(onecls,:);
            allkeeps{cls} = keeptrack(onecls,:);
            fullmd(:,cls) = mean(fullwts{cls},1);
            fullwts{cls} = fullwts{cls}';
        end;  
        fullwts{cls+1} = collmeans';
        allkeeps{cls+1} = keeptrack;
        keeptrack = allkeeps; collmeans = collmeans';
    elseif strcmp(meth,'pca')
        
        if isempty(numdims)
            numdims = size(collmeans,1);
        end;
        [U,S,fullmd] = svds(collmeans',numdims);% if you scale 'acts', you can't scale the 'eigvec'
        fullwts = (U*S)'; % scale 'activations' for appropriate weighting in decomp      
        fdims = numdims; 
    end;
    if strcmp(permtest,'noperm')
    %figure; 
    % plot still figure:
    if fdims < 3
        % for one point per emotion (mean)
% $$$         c1 = 2; c2=1;
% $$$         for e = 1:size(fullmd,1)
% $$$             if fdims == 1
% $$$                 ph=plot(e,fullmd(e,c1)','.');hold on;
% $$$                 set(ph,'markersize',20);set(ph,'color',cols(e,:));
% $$$             else               
% $$$                 ph=plot(fullmd(e,c1),fullmd(e,c2),'.');hold on;
% $$$                 set(ph,'markersize',20);set(ph,'color',cols(e,:));
% $$$                 %ph=text(fullmd(e,c1),fullmd(e,c2),'N');hold on;
% $$$                 %set(ph,'fontweight','bold'); set(ph,'color',cols(e,:));
% $$$                 ph = text(fullmd(e,c1),fullmd(e,c2),emo2{e});
% $$$                 set(ph,'color',cols(e,:)); 
% $$$             end;
% $$$         end;
% $$$         if fdims == 1
% $$$             xlabel(['Emotions']);ylabel(['Dim ',int2str(c1)]);
% $$$         else               
% $$$             xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);
% $$$         end;
% $$$         title(ttl_in); 
    elseif fdims == 3
        %%%%%%%%%%%
        %%  Plot 3 Dims vs each other:
        %figure; % just 3  dims vs each other
        c1 = 1; c2 = 2; c3 = 3;
        for e = 1:size(fullmd,1)
            ph=plot3(fullmd(e,c1),fullmd(e,c2),fullmd(e,c3),'.');hold on;
            set(ph,'markersize',25); set(ph,'color',cols(e,:));
            ph = text(fullmd(e,c1),fullmd(e,c2),fullmd(e,c3),emo2{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        end;
        zl = get(gca,'zlim');
        for e = 1:size(fullmd,1)
            pl =plot3([fullmd(e,c1) fullmd(e,c1)],[fullmd(e,c2) fullmd(e,c2)],[zl(1)  fullmd(e,c3)]);
            set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
        end;
        set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
        xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);zlabel(['Dim ',int2str(c3)]);
        title(ttl_in); 
    else
        row = round(sqrt(fdims));col = ceil(sqrt(fdims));
        for fd = 1:fdims
            sbplot(row,col,fd)
            for e = 1:size(fullmd,1)
                ph=plot(e,fullmd(e,fd),'.');hold on;
    %ph=plot(e,fullwts{fd}(e,:),'.');hold on;
                set(ph,'markersize',15);set(ph,'color',cols(e,:));
                ph = text(e,fullmd(e,fd),emo2{e});
                set(ph,'color',cols(e,:)); set(ph,'fontsize',11); 
            end;  
            ph = plot([get(gca,'xlim')],[0 0],'k-'); hold on;
            set(gca,'xticklabel',[]);
            title(['Dim ',int2str(fd)]);
        end;   
        textsc(ttl_in,'title');
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    end;
    end;

    if mvon == 1
        ad = 1; figure; ads = [0:2:360]; mvs = [[0:.017:3.1],[3.1:-.017:0]]; clear M
        mov = avifile('EmoSpace.avi','fps',8,'quality',100)
        set(gcf,'color','w');
        for mf = 1:180
            if fdims < 3
                figure; 
                % for one point per emotion (mean)
                for e = 1:size(fullmd,1)
                    ph=plot(fullmd(e,1),fullmd(e,2),'.');hold on;
                    set(ph,'markersize',20);set(ph,'color',cols(e,:));
                    ph = text(fullmd(e,1),fullmd(e,2),emo2{e});
                    set(ph,'color',cols(e,:)); 
                end;
                title(['All Sbj Emo Space (numdims: ',int2str(fdims),')']); 
                xlabel(['Dim 1']);ylabel(['Dim 2']);
            else
                %%%%%%%%%%%
                %%  Plot 3 Dims vs each other:
                %figure; % just 3  dims vs each other
                c1 = 1; c2 = 4; c3 = 14;
                for e = 1:size(fullmd,1)
                    ph=plot3(fullmd(e,c1),fullmd(e,c2),fullmd(e,c3),'.');hold on;
                    set(ph,'markersize',25);                set(ph,'color',cols(e,:));
                    ph = text(fullmd(e,c1),fullmd(e,c2),fullmd(e,c3),emo2{e});
                    set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
                end;
                zl = get(gca,'zlim');
                for e = 1:size(fullmd,1)
                    pl =plot3([fullmd(e,c1) fullmd(e,c1)],[fullmd(e,c2) fullmd(e,c2)],[zl(1)  fullmd(e,c3)]);
                    set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
                end;
                set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
                xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);zlabel(['Dim ',int2str(c3)]);
                %title(['All Sbj Emo Space (numdims: ',int2str(fdims),')']);
                fac = .25;
                set(gca,'xlim',[min(fullmd(:,c1))+min(fullmd(:,c1))*fac max(fullmd(:,c1))+max(fullmd(:,c1))*fac]); 
                set(gca,'ylim',[min(fullmd(:,c2))+min(fullmd(:,c2))*fac max(fullmd(:,c2))+max(fullmd(:,c2))*fac]);  
                set(gca,'zlim',[min(fullmd(:,c3))+min(fullmd(:,c3))*fac max(fullmd(:,c3))+max(fullmd(:,c3))*fac]);
                view(20+ads(mf),32+sin(mvs(mf))*3); 
                M2(mf) = getframe(gcf);
                %mov = addframe(mov,M);
            end;
        end;
        %mov = close(mov);
        movie(M2)
        movie2avi(M2,'/data/common2/emotion/EmoSpace.avi','fps',8,'quality',100)
        save /data/common2/emotion/EmoSpaceMov.mat M
    end;
    %%% 
