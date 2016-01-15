% takes a {nx}{dim}(facs) cell array and converts to a {em}{nx}(facs) array to get emo/cluster relationships
%
% [faclist allfacs] = Facs2Clusts(complist,fullpaths,emomeans,whichfacs,nxgdims);
%
% 
%
% OUTPUTS:
% faclist -- a {emotion}{subject}(orig factor) cell array to relate subject effects to specific
%            emotions and to then relate to existing clusters...
% allfacs -- a {emotion}{subject}(3D matrix of comps X freqs X facs) to be able to plot the results


function [faclist, allfacs,clustlist,specs] = Facs2Clusts(complist,savedat,fullpaths,emomeans,whichfacs);
    
    
    emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
    % limit to non-network factors (only one whichfacs{nx}{dim})
    clear faclist allfacs
    for nx = 1:length(fullpaths)
        %if ~isempty(nxgdims{nx})
            %keeplist = zeros(1,0);
            s = load([fullpaths{nx},savedat,'.mat']);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);  
            ws = wts*sph;    activations = ws*icamatall;    
            clear wts sph ws

            tpfacs = cell(1,15);
            currfacs = zeros(15,0);
            for dim = 1:size(emomeans{nx},1)
                %if length(whichfacs{nx}{dim}) < 2
                % emocut should be a std above 0 (not zscore)
                    %emocut = max(abs(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim)))*.7;%tried .8,.7
                    %currem = find(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim) > emocut);
                    %currem = find(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim)/std(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim)) > 1); % above 2 stds
                    emocut = 1.5; % 1.5 mean std above zero
                currem = [];
                    currem = find(emomeans{nx}(dim,:) > emocut);
                    if ~isempty(currem)
                    for em = 1:length(currem)
                        currfacs(currem(em),end+1:end+length(whichfacs{nx}{dim})) = whichfacs{nx}{dim};
                        if isempty(tpfacs{currem(em)})
                            g = 1;
                        else
                            g = size(tpfacs{currem(em)},3) + 1;
                        end;                    
                        for fac = 1:length(whichfacs{nx}{dim})
                            if whichfacs{nx}{dim}(fac) > 0
                                for rcp = 1:length(complist{nx})
                                    tpfacs{currem(em)}(rcp,:,g) = activations(abs(whichfacs{nx}{dim}(fac)),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                                end;
                                
                            else
                                for rcp = 1:length(complist{nx})
                                    tpfacs{currem(em)}(rcp,:,g) = activations(abs(whichfacs{nx}{dim}(fac)),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1; 
                                end;
                            end;
                            
                        end;
                    end;
                    end;
                    % do the same for neg-weighted emos:
                    %currem = find(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim) < -emocut);
                    %currem = find(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim)/std(emomeans{nx}(:,abs(maxfacs{nx}{dim}),dim)) < -1);
                    currem = [];
                    currem = find(emomeans{nx}(dim,:) > emocut);
                    if ~isempty(currem)
                    for em = 1:length(currem)
                        currfacs(currem(em),end+1:end+length(whichfacs{nx}{dim})) = -whichfacs{nx}{dim};
                        if isempty(tpfacs{currem(em)})
                            g = 1;
                        else
                            g = size(tpfacs{currem(em)},3) + 1;
                        end;                    
                        for fac = 1:length(whichfacs{nx}{dim})
                            if whichfacs{nx}{dim}(fac) > 0
                                for rcp = 1:length(complist{nx})
                                    tpfacs{currem(em)}(rcp,:,g) = activations(abs(whichfacs{nx}{dim}(fac)),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1;  % opposite because finding neg mean back-proj
                                end;
                                
                            else
                                for rcp = 1:length(complist{nx})
                                    tpfacs{currem(em)}(rcp,:,g) = activations(abs(whichfacs{nx}{dim}(fac)),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                                end;
                            end;
                            
                        end;
                    end;
                    end;
                %end;
            end;
            for em = 1:size(currfacs,1)
                if ~isempty(currfacs(em,find(currfacs(em,:))))
                    faclist{em}{nx} = currfacs(em,find(currfacs(em,:)));
                    allfacs{em}{nx} = tpfacs{em};
                else
                    faclist{em}{nx} = [];
                    allfacs{em}{nx} = [];
                end;                
            end;                
        %end; % to if ~isempty
        fprintf('Subject %s done.\n',int2str(nx));
    end; % to nx = 1
    

    % for each emo, find which clusters are represented in what proportion
    % must compare factor lists from each subject to 'facvec' in each file:
    clear clustlist
    r = load('/data/common4/emotion/CoModAlphaClusts.mat'); % just to get r.frs
    for em = 1:length(emos)
        tpspecspos = cell(1,8);
        tpspecsneg = cell(1,8);
        cls = 1;  % just to initialize matrix:
        tpspecspos{cls} = zeros(0,length(r.frs));
        tpspecsneg{cls} = zeros(0,length(r.frs));
        for nx = 1:length(complist)
            cls = 1; 
            % find members of alpha clusters:    
            r = load('/data/common4/emotion/CoModAlphaClusts.mat');
            for fc = 1:length(r.facvec)
                if find(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx}))                   
                    clustlist{em}{cls}{nx} = faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx})); 
                    
                    cidx = find(r.alphafacs{fc}(:,1) == nx & ismember(r.alphafacs{fc}(:,2),abs(clustlist{em}{cls}{nx})'));
                    multfac = (faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.alphafacs{fc}(cidx,2)))./abs(faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.alphafacs{fc}(cidx,2)))));
                    for memb = 1:length(clustlist{em}{cls}{nx})
                        curridx = find(r.alphafacs{fc}(:,1) == nx & r.alphafacs{fc}(:,2) == abs(clustlist{em}{cls}{nx}(memb)));
                        if multfac(memb) > 0
                            tpspecspos{cls}(end+1,:) = multfac(memb) * r.alphatempls{fc}(curridx,:);
                        else
                            tpspecsneg{cls}(end+1,:) = multfac(memb) * r.alphatempls{fc}(curridx,:);
                        end;
                        
                    end;
                end;
                cls = cls+1;
            end;
            % find members of beta clusters:    
            r = load('/data/common4/emotion/BetaClusters.mat');
            for fc = 1:length(r.facvec)
                if find(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx}))
                    clustlist{em}{cls}{nx} = faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx}));  
                    
                    cidx = find(r.betafacs{fc}(:,1) == nx & ismember(r.betafacs{fc}(:,2),abs(clustlist{em}{cls}{nx})'));
                    multfac = (faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.betafacs{fc}(cidx,2)))./abs(faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.betafacs{fc}(cidx,2)))));
                    
                    for memb = 1:length(clustlist{em}{cls}{nx})
                        curridx = find(r.betafacs{fc}(:,1) == nx & r.betafacs{fc}(:,2) == abs(clustlist{em}{cls}{nx}(memb)));
                        if multfac(memb) > 0
                            tpspecspos{cls}(end+1,:) = multfac(memb) * r.betatempls{fc}(curridx,:);
                        else
                            tpspecsneg{cls}(end+1,:) = multfac(memb) * r.betatempls{fc}(curridx,:);
                        end;
                    end;
                    
                end;
                cls = cls+1;
            end;         
            
            % find members of gamma clusters:    
            r = load('/data/common4/emotion/GammaClusters.mat');
            for fc = 1:length(r.facvec)
                if find(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx}))
                    clustlist{em}{cls}{nx} = faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.facvec{fc}{nx}));                    
                    cidx = find(r.gamafacs{fc}(:,1) == nx & ismember(r.gamafacs{fc}(:,2),abs(clustlist{em}{cls}{nx})'));
                    multfac = (faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.gamafacs{fc}(cidx,2)))./abs(faclist{em}{nx}(ismember(abs(faclist{em}{nx}),r.gamafacs{fc}(cidx,2)))));
                    for memb = 1:length(clustlist{em}{cls}{nx})
                        curridx = find(r.gamafacs{fc}(:,1) == nx & r.gamafacs{fc}(:,2) == abs(clustlist{em}{cls}{nx}(memb)));
                        if multfac(memb) > 0
                        tpspecspos{cls}(end+1,:) = multfac(memb) * r.gamatempls{fc}(curridx,:);
                        else
                        tpspecsneg{cls}(end+1,:) = multfac(memb) * r.gamatempls{fc}(curridx,:);
                        end;
                    end;
                end;
                cls = cls+1;
            end;
            
        end;
        specs{em}{1} = tpspecspos;
        specs{em}{2} = tpspecsneg;
        fprintf('Emotion %s done.\n',emos{em});
    end;
