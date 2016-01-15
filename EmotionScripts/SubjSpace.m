% takes cluster-inclusion info to cluster subjects
% make a vector for each subject indicating, for each cluster
% whether the subj contributes an IM (1) or not (0)

validsbj = [1:21,23:35];
sspace = zeros(8,length(validsbj));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get stats on all clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common4/emotion/AllCoModAlpha.mat 
% find how many comps/IMs/subjects per cluster
for cls = 1:length(facvec)
    for nxx = 1:length(validsbj)
        nx = validsbj(nxx);
        if ~isempty(facvec{cls}{nx})
          sspace(cls,nxx) = 1;
        end;
    end;
end;
load /data/common4/emotion/AllCoModBeta.mat 
for cls = 1:2
    for nxx = 1:length(validsbj)
        nx = validsbj(nxx);
        if ~isempty(facvec{cls}{nx})
          sspace(cls+3,nxx) = 1;
        end;
    end;
end;
load /data/common4/emotion/AllCoModGama.mat 
for cls = 1:length(facvec)
    for nxx = 1:length(validsbj)
        nx = validsbj(nxx);
        if ~isempty(facvec{cls}{nx})
          sspace(cls+5,nxx) = 1;
        end;
    end;
end;
[pc,eigvec,sv] = runpca(sspace);


[weights,sphere,compvars,bias,signs,lrates,activations]  = runica(sspace,'stop',1e-7,'pca',3);
