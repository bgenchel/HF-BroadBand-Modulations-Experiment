% use plain ICA to find emotion space for all subjects
%
%
%
%
%
%
% emomeans -- straight means: {nx}(dim,emo), or deciles:  {subj}(IM,emo,dec)
%

function [emomeans,emodeciles,subjpoints] = EmoWeights(savedat,fullpaths,gdcomps,subjlist,nsamples);
    
    if ~exist('nsamples')
        nsamples = [];
    end;
    subjpoints = cell(1,length(subjlist));
    %[emoorders] = FindEmoOrder(fullpaths,emos);
    clear   emomeans
    for nxx = 1:length(subjlist)   
        nx = subjlist(nxx);
        emeans = zeros(15,15);
        s = load([fullpaths{nx},savedat,'.mat']);  
        if isfield(s,'eigfile')
            sph=floatread([fullpaths{nx},savedat,'.sph'],[s.pcs s.pcs],[],0); 
            wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.pcs],[],0);        
        else
            sph=floatread([fullpaths{nx},savedat,'.sph'],[s.numrows s.numrows],[],0); 
            wts=floatread([fullpaths{nx},savedat,'.wts'],[s.pcs s.numrows],[],0); 
        end;
        ws = wts*sph;    winv = pinv(ws); 
        clear wts sph ws 
        if isfield(s,'eigfile')
            speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
            specwts = speceig*winv;  % templates   
            winv = specwts;    
        end;
        newwinv = winv';
        
         epoints = cell(1,size(newwinv,1));
        %%%%%%%%%%  find straight median of each dim/emo %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for dim = 1:size(newwinv,1) % for all (new)dims...
            clear mintrials
            if ~isempty(nsamples) % takes only first 'nsamples' samples into median
                for e = 1:length(s.dstrials) % break up windows into emos...
                    emeans(dim,e) = median(newwinv(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e-1))+nsamples));
                    % save all vals for later  
                    emopoints{e} = newwinv(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e-1))+nsamples);
                    if ~exist('mintrials')
                        mintrials = length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e-1))+nsamples);
                    else
                        mintrials = min(mintrials,length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e-1))+nsamples));
                    end;
                end;                
            else                
                for e = 1:length(s.dstrials) % break up windows into emos...
                    emeans(dim,e) = median(newwinv(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))));
                    % save all vals for later  
                    emopoints{e} = newwinv(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)));
                    if ~exist('mintrials')
                        mintrials = length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)));
                    else
                        mintrials = min(mintrials,length(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))));
                    end;
                end;
            end;
            % downsample for each dim
            for e = 1:length(emopoints)
                rpoints = randperm(length(emopoints{e}));
                epoints{dim} = [epoints{dim},emopoints{e}(rpoints(1:mintrials))'];
            end;
        end; % makes a 15 dims X 15 emotions matrix
        subjpoints{nx} = epoints;
        emomeans{nx} = emeans;
        %%%%%%%%%%  find decile values for each dim/emo %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for dim = 1:size(newwinv,1) % for all (new)dims...
            for e = 1:length(s.dstrials) % break up windows into emos...
                tmpdist = sort(newwinv(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))));                
                for d = 10:10:90
                    emodeciles{nx}(dim,e,d/10) = tmpdist(ceil(length(tmpdist)*(d/100)));
                end;
            end;
        end;
        fprintf('Subject %s of %s total subjs done.\n',int2str(nx),int2str(length(subjlist)));
    end;

    
    
