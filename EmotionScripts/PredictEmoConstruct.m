% calls in Emotion data for a single subject and constructs
% a (IM x time*emo) matrix with a corresponding vector with
% emotional valence of each time point
%
% [outmat emovec] = PredictEmoConstruct(datafile,fullpath,emofactor,rmemos);
%
% datafile -- [string] name of IM datafile
% fullpath -- [sting] full data directory location of 'datafile'
% ims -- [vector] IM indices to include in output matrix. [] uses all. 
% emofactor -- [string] 'valence' or 'arousal'
% rmemos -- [cell array] list of emotions to remove from 
%          output matrices. The emotion labels are as follows:
% 'anger','frustration','jealousy','fear','disgust','grief','sad' 
% 'compassion','love','relief','content','awe','happy','joy','excite'
% if labels are not entered exactly, they will not be detected.
%
% OUTPUT:
% outmat -- [matrix] size: (# IMs x # time points)- contains IM weights
%           for each IM over all specified time points (ie, emotions). 
% emovec -- [vector] size: (1 x # time points) contains emotion valence
%           or arousal ratings for each time point



function [outmat emovec] = PredictEmoConstruct(datafile,fullpath,ims,emofactor,rmemos);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import IM data:------------------
    s = load([fullpath,datafile,'.mat']);  % load file info
    sph=floatread([fullpath,datafile,'.sph'],[s.pcs s.pcs],[],0); % load sphere 
    wts=floatread([fullpath,datafile,'.wts'],[s.pcs s.pcs],[],0); % load weights        
    ws = wts*sph;    winv = pinv(ws);  % create inverse weight matrix
    clear wts sph ws     
    speceig = floatread([fullpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);% load pca winv
    winv = speceig*winv;  % pass ICA winv through PCA winv matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ims)
        ims = [1:size(winv,2)]; % default is all ims
    end;
    
    % create valence vector according to target emotions:

    w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');% load emotion ratings
    rmindices = find(ismember(w.emos,rmemos)); %find emotions to remove from matrix
    if strcmp(emofactor,'valence')
        emodata = w.emoval;
    elseif strcmp(emofactor,'arousal')
        emodata = w.emoactiv;
    else
        fprintf('\nNo emofactor has been detected. \nPlease enter either ''valence'' or ''arousal'' as the ''emofactor'' in the function call. ')
        return; 
    end;

    emovec = []; % initialize valence vector
    outmat = []; % build output matrix from winv

    for ds = 1:size(s.keeptrack,1) % ds corresponds to emotions
        if ~ismember(ds,rmindices) % skip row if in remove list
            pnts = [s.keeptrack(ds,1):s.keeptrack(ds,2)];
            emovec = [emovec, repmat(emodata(ds),[1 length(pnts)])];
            outmat = [outmat, winv(pnts,ims)'];
        end;
    end;
