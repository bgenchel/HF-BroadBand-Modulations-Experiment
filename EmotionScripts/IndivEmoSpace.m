%
%
%
%
%
%
%
% ploton -- [0 | 1] if 1, will plot factors contributing to each emotion. 
%

function [clustlist,specs] = IndivEmoSpace(fullmat,subjlist,fullpaths,savedat,complist,emos,ploton)

    keyboard
    if ~exist('ploton')
        ploton = 0;
    end;    
    
    pcdims = 6; 
    justsubj = unique(subjlist);
    s = load([fullpaths{1},savedat,'.mat']);
        
    clear posemos negemos posfacs negfacs
    for nxx = 1:length(justsubj)
        nx = justsubj(nxx);
        takeout = justsubj;
        takeout(justsubj == justsubj(nxx)) = [];
        tmpmat = fullmat;
        tmpmat(:,find(ismember(subjlist,takeout))) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % md scaling instead %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %clear dd
        %for x = 1:size(fullmat,1)
        %    for y = 1:size(fullmat,1)        
        %dd(x,y) = pdist(fullmat([x y],:),'euclidean');            
        %    end;
        %end;
        %clear dd
        %for x = 1:size(tmpmat,1)
        %    for y = 1:size(tmpmat,1)
        %    dd(x,y) = pdist(tmpmat([x y],:),'euclidean');            
        %    end;
        %end;        
        %[md,STRESS,DISPARITIES] = mdscale(dd,pcdims);
        % This method looks the same as runica, but you don't get the corresponding weights for factors...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ICA instead%%%%%%%%%%%%%%%%%%%%%%%%%%%        
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(tmpmat,'extended',1,'stop',1e-7,'verbose','off','pca',pcdims,'maxsteps',1000);
        ws = weights*sphere; winv = pinv(ws); clear ws weights sphere 
        % find high emos and high facs for each dim
        for dim = 1:size(winv,2)
            posemos{nx}{dim} = find(winv(:,dim) >  max(abs(winv(:,dim))) * .8)';
            negemos{nx}{dim} = find(winv(:,dim) < -max(abs(winv(:,dim))) * .8)';
            pickfac = find(activations(dim,:) > max(abs(activations(dim,:))) * .8);
            for x =  1:length(pickfac)
                for fc = 1:15
                    if ismember(pickfac(x),[fc:15:75])
                        pickfac(x) = fc;
                    end;
                end;
            end;           
            posfacs{nx}{dim} = pickfac;
            pickfac = find(activations(dim,:) < -max(abs(activations(dim,:))) * .8);
            for x =  1:length(pickfac)
                for fc = 1:15
                    if ismember(pickfac(x),[fc:15:75])
                        pickfac(x) = fc;
                    end;
                end;
            end;           
            negfacs{nx}{dim} = pickfac;
        end;        
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now plot for each subject the salient emotions/factors for each dimension
    % Each emotion gets its own figure
    cols = jet(15);cols(10,:) = [.9 .9 0];
    row = round(sqrt(length(justsubj))); col = ceil(sqrt(length(justsubj)));
    for em = 1:length(emos)
        numfacs = 0;
        for nxx = 1:length(justsubj)
            nx = justsubj(nxx);
            tpfacs = zeros(length(complist{nx}),length(s.freqs),0);
            keeplist = zeros(1,0);
            g = 1;
            s = load([fullpaths{nx},savedat,'.mat']);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            %sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
            %wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
            icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);  
            ws = wts*sph;    activations = ws*icamatall;    
            winv = pinv(ws); clear wts sph ws
            for dim = 1:pcdims
                if ismember(posemos{nx}{dim},em)
                    for facs = 1:length(posfacs{nx}{dim})
                        for rcp = 1:length(complist{nx})
                            tpfacs(rcp,:,g) = activations(posfacs{nx}{dim}(facs),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                        end;
                        keeplist(1,end+1) = posfacs{nx}{dim}(facs);
                        g = g + 1;
                        numfacs = numfacs + 1;
                    end;                    
                    for facs = 1:length(negfacs{nx}{dim})
                        for rcp = 1:length(complist{nx})
                            tpfacs(rcp,:,g) = activations(negfacs{nx}{dim}(facs),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1; 
                        end;
                        keeplist(1,end+1) = -negfacs{nx}{dim}(facs);
                        g = g + 1;
                        numfacs = numfacs + 1;
                    end;
                end;
                if ismember(negemos{nx}{dim},em)
                    for facs = 1:length(posfacs{nx}{dim})                
                        for rcp = 1:length(complist{nx})
                            tpfacs(rcp,:,g) = activations(posfacs{nx}{dim}(facs),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1; 
                        end;
                        keeplist(1,end+1) = -posfacs{nx}{dim}(facs);
                        g = g + 1;
                        numfacs = numfacs + 1;
                    end;
                    for facs = 1:length(negfacs{nx}{dim}) 
                        for rcp = 1:length(complist{nx})
                            tpfacs(rcp,:,g) = activations(negfacs{nx}{dim}(facs),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                        end;
                        keeplist(1,end+1) = negfacs{nx}{dim}(facs);
                        g = g + 1;
                        numfacs = numfacs + 1;
                    end;                    
                end;
            end;
            if isempty(find(tpfacs))
                tpfacs = [];
                keeplist = [];
            end;        
            faclist{em}{nx} = keeplist;
            allfacs{em}{nx} = tpfacs;
            fprintf('\n One More SUBJECT Done: %i',nx);
        end;
        if ploton == 1
            figure; pp = 1;   lnwdth = 1.5;
            row = round(sqrt(numfacs));
            col = ceil(sqrt(numfacs));        
            for nx = 1:length(allfacs{em})
                if ~isempty(allfacs{em}{nx})
                    cols = jet(length(complist{nx}));
                    for fac = 1:size(allfacs{em}{nx},3)
                        sbplot(row,col,pp)
                        for cp = 1:size(allfacs{em}{nx},1)
                            ph = plot(s.freqs,allfacs{em}{nx}(cp,:,fac),'linewidth',lnwdth); 
                            hold on; set(ph,'color',cols(cp,:));
                        end;
                        set(gca,'ylim',[-10 10]); set(gca,'xticklabel',[]);
                        plot([10 10],[get(gca,'ylim')],'k-','linewidth',1);
                        title(['S ',int2str(nx),' F ',int2str(faclist{em}{nx}(fac))]);
                        pp = pp+1;
                    end;
                end;
            end;  
            textsc(['All contributing factors to ',emos{em}],'title');axcopy
        end;
        fprintf('\n One More EMOTION Done: %s',emos{em});
    end;
    
    %%%%%%%%%%%%%%%%%%%%%     
               
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
    end;
    

