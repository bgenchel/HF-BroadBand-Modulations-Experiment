% For each subject, find IMs whose weights can differentiate one or more emotions
%
% SigEmoShifts(savedat,fullpaths,subjlist);
%
%
% OUTPUT:
% emosigs -- [cell array] for each subject, number of sig diffs by anova (dim x emos)
% emoPs-- [cell array] significance values for all dims (1 x dim)
%

function [emosigs,emoPs] = SigEmoShifts(savedat,fullpaths,subjlist,emos);
    
    
    for nxx = 1:length(subjlist)
        clear P numsigs  % clear the matrix for diff # of IMs
        nx = subjlist(nxx);
        s = load([fullpaths{nx},savedat,'.mat']);  
        sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
        ws = wts*sph;    winv = pinv(ws);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for dim = 1:size(winv,2) % for all (new)dims...
            clear dimwts szs eqwts 
            for e = 1:length(s.dstrials) % break up windows into emos...
                ewts{e} = winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),dim);
                szs(1,e) = length(ewts{e});
            end;
            for e = 1:length(ewts)
                randidx = randperm(szs(1,e));% choose random windows
                randidx(min(szs)+1:end) = [];% to match min # windows
                eqwts(:,e) = ewts{e}(randidx);
            end;
            [P(1,dim),ANOVATAB,STATS] = anova1(eqwts,emos,'on');close;close;
            comp = multcompare(STATS,'alpha',.001,'ctype','bonferroni'); close;
            comppairs = [];
            for cp = 1:size(comp,1)
                if length(find(comp(cp,[3,5]) == abs(comp(cp,[3,5])))) ~= 1 %(not straddling 0=sig)
                    comppairs = [comppairs;[comp(cp,[1:2])]];
                end;
            end;
            for e = 1:size(eqwts,2)
                numsigs(dim,e) = length(find(comppairs==e));
            end;            
        end; 
        emosigs{nx} = numsigs; % for each dim x emotion, how many sig diffs by anova
        emoPs{nx} = P; % save significance values for all dims
        fprintf('Subject %s of %s total subjs done.\n',int2str(nxx),int2str(length(subjlist)));
    end;
   
    
    
