% takes training and test data and performs a classification based on knnclassify()
%
% [mtchpercents,meantestmatchs,subjoneconf,fullconf,subjemeans,keepuseims,nims] = ClassifyIMdata(savedat,fullpaths,subjlist,minsize,testpercent,threshold,method);
%
% max overlap = 2; with 2-sec windows!!
%
% minsize -- [number] minimum size of training data (s.dstrials) (Default: 5)
% threshold -- [number or vector] If a number, will use that number of highest F-score IMs
%              for classification. If a vector, the length of the vector must correspond to 
%              'subjlist' and give, for each subject, the F-score cutoff to be used for 
%              each subject (vector can contain the same number for all subjects).
% method ['euclidean' or 'correlation'] for knnclassify (default: 'euclidean')
%

function [mtchpercents,meantestmatchs,subjoneconf,fullconf,meanconf,subjemeans,keepuseims,nims,subjNdelpoints] = ClassifyIMdata(savedat,fullpaths,subjlist,minsize,testpercent,threshold,method);

        
    
    if ~exist('method')
        method = 'euclidean';
    elseif isempty(method)
        method = 'euclidean';
    end;

     if ~exist('minsize')
        minsize = 5;
    elseif isempty(method)
        minsize = 5;
    end;
    
    % initialize output variables
    mtchpercents = zeros(length(fullpaths),15); 
    meantestmatchs = zeros(length(fullpaths),15);  
    subjoneconf = cell(1,length(fullpaths));             
    fullconf = zeros(15,15); 
    meanconf= zeros(15,15); 
    n4perc = zeros(1,15);
    f = zeros(1,15);
    
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        s = load([fullpaths{nx},savedat,'.mat']);  
        sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
        wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
        ws = wts*sph;    winv = inv(ws); 
        clear wts sph ws 
        speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
        specwts = speceig*winv;  % templates   
        winv = specwts;  clear delpoints emeans  npoints 
        for e = 1:length(s.dstrials) % break up windows into emos...
            ndec = round(s.dstrials(e)*testpercent); 
            epoints = sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e));
            rpoints = [round(length(epoints)/2-ndec/2):round(length(epoints)/2+ndec/2)];
            delpoints{e} = epoints(rpoints);% prediction windows
            npoints(1,e) = length(delpoints{e});
            epoints([rpoints(1)-1,rpoints,rpoints(end)+1]) = [];% take out buffer windows
            for dim = 1:size(winv,2) % for all (new)dims...
                emeans(dim,e) = mean(winv(epoints,dim));
            end;
        end; % makes a 15 dims X 15 emotions matrix
        subjemeans{nx} = emeans;
        subjNdelpoints{nx} = npoints;
        % OR ------------
        useims = [];grp = [];
        for e=1:length(s.dstrials)
            for t = s.keeptrack(e,1):s.keeptrack(e,2)
                grp{t} = s.datset{e};
            end;
        end; clear Fs
        for im = 1:size(winv,2)
            [P,table]=anovan(winv(:,im)',{grp},'display','off');
            Fs(1,im) = table{2,6};
        end;
        tmpF = sort(Fs); 
        if length(threshold) == 1 % set cut off F-score based on function input.
            cutval = tmpF(end-threshold);
        else
            cutval = threshold(nxx);% actual F-score cutoff given
        end;
        useims = find(Fs >= cutval); % use IMs in top 'perctop' percentile
        keepuseims{nx} = useims;
        nims(1,nx) = length(keepuseims{nx});
        %----------------
        if ~isempty(useims)
            clear imcombos colldata scores ntrials  percscores  mtchs svmgrp confmatrix
            for emoidx = 1:length(s.dstrials)
                if s.dstrials(emoidx) > minsize % 130 such that 1 min remains after test data removed
                    clear mtch sgrp  oneepochconf
                    for itr = 1:length(delpoints{emoidx})                    
                        nowtr = delpoints{emoidx}(itr); %use each point individually
                        onetrial = mean(winv(nowtr,useims),1);
                        nearemo = knnclassify(onetrial,emeans(useims,:)',s.datset,1,method);% classify by 1-epoch
                        oneepochconf(1,itr) = find(ismember(s.datset,nearemo));
                        if strcmp(s.datset{emoidx},nearemo) % correct classification?
                            mtch(itr) = 1; 
                        else
                            mtch(itr) = 0;
                        end;
                    end;
                    onetrial = mean(winv(delpoints{emoidx},useims),1);% mean of all 5%
                    nearemo = knnclassify(onetrial,emeans(useims,:)',s.datset,1,method); % classify by full 5%
                    meanepochconf = find(ismember(s.datset,nearemo));
                    meanconf(emoidx,meanepochconf) = [meanconf(emoidx,meanepochconf)+1]; 
                    n4perc(1,emoidx) = n4perc(1,emoidx) + 1;
                    for e = 1:length(s.datset)
                        confmatrix(emoidx,e) = 100*(length(find(oneepochconf == e))/length(find(oneepochconf)));
                        fullconf(emoidx,e) = [fullconf(emoidx,e)+length(find(oneepochconf == e))]; 
                    end
                    f(1,emoidx) = f(1,emoidx) + length(find(oneepochconf));
                    mtchs(emoidx) = sum(mtch)/length(mtch); %
                   if strcmp(s.datset{emoidx},nearemo)
                        alldatmatch(emoidx) = 1; % 1/0 using all test data
                    else
                        alldatmatch(emoidx) = 0;
                    end;
                else
                    mtchs(emoidx) = NaN;
                    alldatmatch(emoidx) = NaN;
                    confmatrix(emoidx,:) = NaN;
                end;
            end;
            mtchpercents(nx,:) = mtchs*100; 
            meantestmatchs(nx,:) = alldatmatch; 
            subjoneconf{nx} = confmatrix;
        else
            mtchpercents(nx,:) = NaN; 
            meantestmatchs(nx,:) = NaN; 
            subjoneconf{nx} = NaN;            
        end;
        fprintf('.');
    end;
    fullconf = 100*(fullconf./repmat(f,[size(fullconf,1) 1])); % percent confusions
    meanconf = 100*(meanconf./repmat(n4perc,[size(meanconf,1) 1])); % percent confusions from mean of test data
