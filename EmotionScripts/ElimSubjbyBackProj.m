% takes a full matrix of emotions x factors*subjects and iteratively finds the most conserved
% dimensions by different combinations of subjects included
%
%
%
%
%
%
%
%
%

function [elimsubj] = ElimSubjbyBackProj(fullmat,subjlist)
  
    keyboard
    newsubjlist = subjlist;
    pcdims = 6;
    justsubj = unique(subjlist);
    numfacs = length(find(subjlist == 1));
    [weights,sphere,compvars,bias,signs,lrates,activations] = runica(fullmat,'extended',1,'stop',1e-7,'verbose','on','pca',pcdims,'lrate',.00001);
    ws = weights*sphere; fullwinv = pinv(ws); clear ws weights sphere
    for bpdim = 1:pcdims
        acts = activations;
        allbut = [1:size(acts,1)];  allbut(bpdim) = [];
        acts(allbut,:) = 0;
        
        backproj = fullwinv*acts;
        
        for nx = 1:length(justsubj)
            for fac = 1:numfacs
                allrms(nx,fac) = sqrt(mean(backproj(:,(nx-1)*numfacs+fac).^2));
            end;
        end;
        allmaxes = max(allrms');
        mnrms = mean(allmaxes);
        stdrms = std(allmaxes);
        elimsubj = zeros(1,0);
        for nx = 1:length(justsubj)
            subjmaxes(nx,bpdim) = max(allrms(nx,:));
        end;        
    end;
    mnmaxes = sum(subjmaxes');
    thresh = mean(mnmaxes) -  std(mnmaxes)*1.25;
    elimsubj = justsubj(find(mnmaxes < thresh));
    
