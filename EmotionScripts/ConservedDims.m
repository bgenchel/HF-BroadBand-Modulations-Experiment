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

function [suggestrm,subjprob] = ConservedDims(fullmat,subjlist,iternum)
  
    
    pcdims = 6;
    subnum = 6;
    [weights,sphere,compvars,bias,signs,lrates,activations] = runica(fullmat,'extended',1,'stop',1e-7,'verbose','on','pca',pcdims);
    ws = weights*sphere; fullwinv = pinv(ws); clear ws weights sphere activations
    
    %allkeeps = cell(1,0);
    justsubj = unique(subjlist);
    shuffvec = [1:length(justsubj)];
    
    mn = length(justsubj) - subnum;
    
    %for mn = 1:length(justsubj)/2
        allwinvs = zeros(size(fullwinv,1),size(fullwinv,2),iternum);
        keepdims = zeros(1,0);
        goodsubj = zeros(0,subnum);
        triedpairs = zeros(iternum,subnum);
        corrpairs = zeros(iternum,pcdims);
        for itt = 1:iternum
            randselec = shuffle(shuffvec);
            takeout = randselec(1:mn);
            to = justsubj(takeout);
            tmpmat = fullmat;
            tmpmat(:,find(ismember(subjlist,to))) = [];
             %[weights,sphere] = binica( tmpmat, 'extended',1,'stop',1e-7,'verbose','off');
            [weights,sphere,compvars,bias,signs,lrates,activations] = runica(tmpmat,'extended',1,'stop',1e-7,'verbose','on','pca',pcdims);
            ws = weights*sphere; winv = pinv(ws); clear ws weights sphere activations
            
            allwinvs(:,:,itt) = winv;
            [corr,indx,indy,corrs] = matcorr(allwinvs(:,:,itt)',fullwinv');
            js = justsubj;
            js(takeout) = [];
            triedpairs(itt,:) = js;
            corrpairs(itt,:) = corr';
        end;

        %allkeeps{mn} = unique(keepdims);
    %end;
    
    for nxx = 1:length(justsubj)
        nx = justsubj(nxx);
        subjprob(1,nxx) = length(find(triedpairs == nx));
    end;
    hicorrs = triedpairs(find(corrpairs(:,1) > mean(abs(corrpairs(:,1))) + std(abs(corrpairs(:,1)))),:); % above avg correlations
    for nxx = 1:length(justsubj)
        nx = justsubj(nxx);
        subjprob(2,nxx) = length(find(hicorrs == nx));
    end;
    subjprob(3,:) = subjprob(1,:)./subjprob(2,:);
    figure; bar(subjprob')
    figure; bar(subjprob(3,:)')
    
    rmthresh = mean(subjprob(3,:)) + std(subjprob(3,:));
    suggestrm = justsubj(find(subjprob(3,:)>rmthresh));

    
