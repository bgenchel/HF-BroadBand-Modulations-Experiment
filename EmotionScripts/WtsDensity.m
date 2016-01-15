% use plain ICA to find emotion space for all subjects
%
%
%
%
% ims -- [cell array] of two vectors of IMs to analyze weights pairwise
%
%

function [quadprobs,quadwinv,quaddens] = WtsDensity(savedat,fullpaths,ims,btstrap);

    
    s = load([fullpaths,savedat,'.mat']);     
    wts = floatread([fullpaths,savedat,'.wts'],[s.pcs s.numtrials],[],0);
    sph = floatread([fullpaths,savedat,'.sph'],[s.numtrials s.numtrials],[],0);  
    ws = wts*sph;  winv = pinv(ws); % winv is what you need.
    for imm = 1:length(ims{1})
        im1 = ims{1}(imm);
        im2 = ims{2}(imm);
        if im1 < 0
            winv(:,im1) = winv(:,im1)*-1;
        end;
         if im2 < 0
            winv(:,im2) = winv(:,im2)*-1;
        end;
            if ~isempty(btstrap)
                for b = 1:btstrap                        
                    perms(:,:,b) = [winv(randperm(size(winv,1)),im1),winv(randperm(size(winv,1)),im2)];
                    quadvals(b,1) = length(find(perms(:,1,b) < 0 & perms(:,2,b) > 0)')/(size(winv,1)/4);
                    quadvals(b,2) = length(find(perms(:,1,b) > 0 & perms(:,2,b) > 0)')/(size(winv,1)/4);
                    quadvals(b,3) = length(find(perms(:,1,b) < 0 & perms(:,2,b) < 0)')/(size(winv,1)/4); 
                    quadvals(b,4) = length(find(perms(:,1,b) > 0 & perms(:,2,b) < 0)')/(size(winv,1)/4); 
                end;
                quadprobs{imm}(1,:) = max(quadvals,[],1);
                quadprobs{imm}(2,:) = min(quadvals,[],1);
                quadwinv{imm}(1,1) = length(find(winv(:,im1) < 0 &  winv(:,im2) > 0))/(size(winv,1)/4);
                quadwinv{imm}(1,2) = length(find(winv(:,im1) > 0 &  winv(:,im2) > 0))/(size(winv,1)/4);
                quadwinv{imm}(1,3) = length(find(winv(:,im1) < 0 &  winv(:,im2) < 0))/(size(winv,1)/4);
                quadwinv{imm}(1,4) = length(find(winv(:,im1) > 0 &  winv(:,im2) < 0))/(size(winv,1)/4);            
            end;
        for em = 1:length(s.dstrials)
            quaddens{imm}(em,1) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) < 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) > 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
            quaddens{imm}(em,2) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) > 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) > 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
            quaddens{imm}(em,3) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) < 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) < 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
            quaddens{imm}(em,4) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) > 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) < 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);           
        end;
        
    end;
    
    
    
