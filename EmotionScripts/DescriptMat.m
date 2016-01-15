% collects mean, median, variance, skewness, and kurtosis for all factors, all emos
%
% [allsubjmat,alldatmat] = DescriptMat(savedat,fullpaths,subjlist);
% 
% INPUTS:
% savedat -- [string] name of .mat file to load, containing ICA decomp parameters
% fullpaths -- [cell array] of strings with fullpaths to all subjects in expt
% subjlist -- [vector] of subj indexes to use in analysis
% divtrials -- ['on' or []] if 'on', will calculate MI for each set of trials
%              in the variable 'dstrials' in the 'savedat' file
%
% OUTPUTS:
% distchar -- [matrix] emotions x 5 measures (mean,median,var,skew,kurt) x IMs
% alldatmat -- [matrix] same measures as allsubjmat, but for full data 
%              distribution (all emos). [5 measures x IMs x subj] 
% MI -- [cell array] mutual information between weights matrices
% fullxcorr -- dot product of full winv with itself (xcorr = winv' * winv;
% emxcorr -- dot products of 'dstrial' sets with themselves

function [distchar,alldatmat,MI,fullxcorr,emxcorr] = DescriptMat(savedat,fullpaths,subjlist,divtrials);
    
    
    for nxs = 1:length(subjlist)
        nx = subjlist(nxs);
        s = load([fullpaths{nx},savedat,'.mat']);     
        wts = floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numtrials],[],0);
        sph = floatread([fullpaths{nx},savedat,'.sph'],[s.numtrials s.numtrials],[],0);  
        ws = wts*sph;  winv = pinv(ws); % winv is what you need.
        fullxcorr{nx} = winv' * winv;
        
        clear wts sph ws allfacs alltemps
        dc = zeros(15,5,size(winv,2)); % emos x 5 measures x IMs
        alld = zeros(5,size(winv,2)); % 5 measures x IMs
        for e = 1:length(s.dstrials)
            tempmat = winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),:); 
            dc(e,1,:) = mode(tempmat);
            dc(e,2,:) = mean(tempmat);
            dc(e,3,:) = median(tempmat);
            dc(e,4,:) = var(tempmat);
            dc(e,5,:) = skewness(tempmat);
            dc(e,6,:) = kurtosis(tempmat);
            emx(:,:,e) = (winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),:)' * winv(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),:))/s.dstrials(e);
        end;
        emxcorr{nx} = emx;
        alld(1,:) =  mode(winv);
        alld(2,:) =  mean(winv);
        alld(3,:) =  median(winv);
        alld(4,:) =  var(winv);
        alld(5,:) =  skewness(winv);
        alld(6,:) =  kurtosis(winv);
        distchar{nx} = dc;
        alldatmat{nx} = alld;
        
        %calculate mutual information between IMs
        %if ~isempty(divtrials)
        %    for t = 1:length(s.dstrials)
         %      m(:,:,t) = mi_pairs (winv(sum(s.dstrials(1:t-1))+1:sum(s.dstrials(1:t)),:)');
         %   end;
         %   MI{nx} = m;
        %else
            MI{nx} = mi_pairs (winv');
        %end;
        fprintf('\n One More SUBJECT Done: %i',nx);
        clear sph winv wts ws activations icamatall
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
