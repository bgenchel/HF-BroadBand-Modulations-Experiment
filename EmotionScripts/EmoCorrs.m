function [corrs,clustidxs,keeptrack] = EmoCorrs(emomeans,orivecs,subjlist,useims,indiv,elimbyrating,noemos,corrwhat)

    
    
    w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');% subjective ratings
    r = load('/data/common1/emotion/SubjRatings.mat');% load just in case 'elim'
    % 6.4838 +- 2.4025, (mean +- std)
    
    if length(size(emomeans{1})) > 2
        keeptrack = zeros(0,3);
    else
        keeptrack = zeros(0,2);
    end;
    % take instead IMs that were specified 
    collmeans = []; 
    for nx = 1:length(emomeans)
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
        for e = 1:size(collmeans,1)
            if ~isempty(find(collmeans(e,:) == 0))
                othermean = mean(collmeans(e,find(collmeans(e,:)))); % mean of all other values
                zeroidx = find(collmeans(e,find(collmeans(e,:)==0)));
                for z = 1:length(zeroidx)
                    collmeans(e,zeroidx(z)) = othermean;
                end;
            end;
        end;
    end;

    if ~isempty(noemos) % take out emotions that don't correlate with valence
        collmeans(noemos,:) = [];
        w.emoval(noemos) = [];
        w.emoactiv(noemos) = [];
    end;
    
    if strcmp(corrwhat,'Valence')
        regvec = repmat(zscore(w.emoval),[size(collmeans,2) 1]); % emoval or emoactiv        
    elseif strcmp(corrwhat,'Arousal')
        regvec = repmat(zscore(w.emoactiv),[size(collmeans,2) 1]); % emoval or emoactiv
    end;
    
    [corrs,indx,clustidxs,crs] = matcorr(regvec,collmeans');
