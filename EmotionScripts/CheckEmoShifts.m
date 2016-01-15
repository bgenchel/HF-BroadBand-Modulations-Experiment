% takes SpecCoModAnal info and checks for distribution differences between two specified emos
%
%
%
%
%
%
%
%

function [diffims,P] = CheckEmoShifts(savedat,fullpaths,subjlist,emopair,alpha);
    
    
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        s = load([fullpaths{nx},savedat,'.mat']);  
        sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0); 
        icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
        ws = wts*sph;   winv = pinv(ws);   clear wts sph ws
        clear distr
        for im = 1:size(winv,2)
            for emo = 1:2
                distr{emo}(:,im) = winv(sum(s.dstrials(1:emopair(emo)-1))+1:sum(s.dstrials(1:emopair(emo))),im);
            end;
        end;
        [H,P(nx,:),CI,STATS] = ttest2(distr{1},distr{2},alpha);
        diffims{nx} = find(P(nx,:) < alpha);
    end;
    
