% takes precomputed correlation measures and compares
% correlations between IMs involving only the specified
% IC
%
%
%
%
%


function [iccorrs,outboot,icothers] = IMcorrsbyIC(corr,bootcorr,subj,ics,ims,allbigs,correl);
    
    
% if ims is not empty, then only return correlations
% between the specified IMs

    iccorrs = []; icothers = []; outboot = [];
    keepcorr = []; othercorr = []; bootstats = [];
    if isempty(ims)
        if length(size(corr)) > 2
            for ic = 1:length(ics)
                for im1 = 1:size(corr,1)
                    if ismember(ics(ic),allbigs{subj}{im1})
                        relcorr = corr(im1,:,:);% all corrs with im1
                        for im2 = 1:length(relcorr)
                            if ismember(ics(ic),allbigs{subj}{im2})
                                keepcorr = [keepcorr relcorr(im2)];
                            else
                                othercorr = [othercorr relcorr(im2)];
                            end;
                        end;
                    end;
                end;
                iccorrs{ic} = keepcorr(find(keepcorr));
                icothers{ic} = othercorr(find(othercorr));
            end;
        end;       

    else
        
        %if length(size(corr)) > 2
            for ic = 1:length(ics)
                bootstats = zeros(size(bootcorr,4),2,0);
                for im1 = 1:length(ims{1})
                    for im2 = 1:length(ims{2})
                        if ims{1}(im1) < 0 
                            flip = -1;
                        else
                            flip = 1;
                        end;
                        if ims{2}(im2) < 0
                            flip = flip * -1; % flip again if both negative
                        end;   

                        keepcorr = [keepcorr squeeze(corr(abs(ims{1}(im1)),abs(ims{2}(im2)),:))];
                        bootstats(:,:,end+1) =  [squeeze(bootcorr(abs(ims{1}(im1)),abs(ims{2}(im2)),1,:)), squeeze(bootcorr(abs(ims{1}(im1)),abs(ims{2}(im2)),2,:))];
                                    
                    
                    end;
                end;
                iccorrs{ic} = keepcorr; keepcorr = []; 
                outboot{ic} = mean(bootstats,3); % mean of collected IM pair bootstats
            end;
        %end;
        
        
    end;
    
