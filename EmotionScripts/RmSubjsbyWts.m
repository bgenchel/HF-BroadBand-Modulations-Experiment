% takes the 'activations' matrix from 'Emotion Space' decomp and returns subjects without 
% 'highly weighted' (specified by user) factors in the dimensions specified.
%
%
%
%
% INPUTS:
% activations -- matrix of size (dimensions x factors*subjects)
% subjlist -- vector of size 1 x factors*subjects (# columns of activations matrix) indicating
%             for each column which subject is represented
% dims -- [vector] list of dimensions of decomposition to find highly weighted factors for.
% cutwt -- [number] number of standard deviations from the mean to set threshold for rejection.
%
%

function [badsubj] = RmSubjsbyWts(activations,subjlist,dims,cutwt);

    
    for dm = 1:length(dims)
        thresh(dm) = mean(activations(dims(dm),:)) + std(activations(dims(dm),:))*cutwt;
    end;
    
    totsubj = zeros(1,0);
    for dm = 1:length(dims)
        hifacs = find(abs(activations(dims(dm),:)) > thresh(dm));
        hisubj = subjlist(hifacs);
        hisubj = unique(hisubj);
        totsubj(end+1:end+length(hisubj)) = hisubj;
    end;
    totsubj = unique(totsubj);
    allsubj = unique(subjlist);
    
    badsubj = allsubj(find(ismember(allsubj,totsubj) == 0));
    
