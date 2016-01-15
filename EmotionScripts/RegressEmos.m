% takes subjective emotion ratings in valence/arousal and regresses
% with median weights of each IM, each subject
%
%
%
%
%
% OUTPUTS:
% subjvalact -- [cell array] for each subject (cell), gives the correlation coefficient with 
%               (row 1) valence, and (row 2) arousal, for each IM of that subject
% valcorr --  [cell array] for each subject (cell), gives the IM indices that were above 
%              the correlation cutoff for valence. 
% actcorr --  [cell array] for each subject (cell), gives the IM indices that were above 
%              the correlation cutoff for arousal. 
% subjidxval -- [n x 3 matrix] gives, for each entry in valcorr, the 1) subject index, 
%               2) IM index, and 3) the correlation coefficient
% subjidxact -- [n x 3 matrix] gives, for each entry in actcorr, the 1) subject index, 
%               2) IM index, and 3) the correlation coefficient
%

function [subjvalact,valcorr,actcorr,subjidxval,subjidxact] = RegressEmos(emomeans,subjlist,corrcut,useratings);
    
    
    if ~exist('corrcut')
        corrcut = .4;
    end;
    if ~exist('useratings')
        useratings = 4; % less than this rating not used
    end;
    
    w=load('/data/common4/emotion/EmoValence.mat');
    z=load('/data/common2/emotion/SubjRatings.mat');
    
    % JO ratings
    emoactiv = [9, 9, 8, 9, 8, 3.5, 2, 3.5, 5, 3, 2, 3, 6, 8, 9];
    emoval = [1,1,1,2,2.5,1,1,5,9,9,7.5,9,9,9,9];
    w.emoactiv = emoactiv;
    w.emoval = emoval;
    
    
    subjidxval = [];
    subjidxact = [];
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        if ~isempty(useratings)
            delems = find(z.ratings{nx} < useratings);
            emomeans{nx}(:,delems) = [];
            ev = w.emoval; ev = ev-mean(ev);
            ev(delems) = [];
            ea = w.emoactiv;ea = ea - mean(ea);
            ea(delems) = [];
        end;
        % correlation with valence
        mat4corr = repmat(ev,[size(emomeans{nx},1) 1]);
        [corr,indx,indy,corrs] = matcorr(mat4corr,emomeans{nx});
        [val idx] = sort(indy);
        subjvalact{nx}(1,:) = corr(idx)';% put back into abs order
        x = find(abs(subjvalact{nx}(1,:)) > corrcut);
        valcorr{nx} = x;
        subjidxval = [subjidxval;[repmat(nx,[length(x) 1]) x' subjvalact{nx}(1,x)']];
        
        % correlation with arousal
        mat4corr = repmat(ea,[size(emomeans{nx},1) 1]);
        [corr,indx,indy,corrs] = matcorr(mat4corr,emomeans{nx});
        [val idx] = sort(indy);
        subjvalact{nx}(2,:) = corr(idx)';% put back into abs order
        y = find(abs(subjvalact{nx}(2,:)) > corrcut);
        actcorr{nx} = y;    
        subjidxact = [subjidxact;[repmat(nx,[length(y) 1]) y' subjvalact{nx}(2,y)']];
    end;
    % collect highly correlated IMs in each category

    figure; pl = 1;
    row = round(sqrt(length(subjlist))); col = ceil(sqrt(length(subjlist)));
    scols = hsv(length(subjlist));
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        sbplot(row,col,pl);
        for im = 1:size(emomeans{nx},1)
            ph = plot(subjvalact{nx}(2,im),subjvalact{nx}(1,im),'.'); hold on;
            set(ph,'color',scols(pl,:));
        end;
        ph = plot([0 0],[get(gca,'ylim')],'k-');
        ph = plot([get(gca,'xlim')],[0 0],'k-');
        title(['Subj ',int2str(nx)]); pl = pl+1;
    end;
    
