% Calls in ERSP info from each subject for selected components and plots avg
% includes optional masking for binomial probability across subjects
%
% [allcross,times, freqs] = PlotSigAvgCoh(savename,paths,plotcomps,clim,matalpha,binomalpha,cond);
% 
% datafile: name of datafile in respective subject paths (.mat). (string)
% paths: cell array of strings with full paths to datafiles.
% origclust: cell array of vectors with original component lists for each subject (provides index into matrices)
% plotcomps: cell array of component pairs for each subject to include in average.
%            form should be: plotcomps{subject}{IC pair}
%            or in other words: plotcomps{1}{1} = [2 6];% subj 1, cross coh between IC2 and IC6
%                               plotcomps{1}{2} = [5 7];% subj 1, cross coh between IC5 and IC7
%                               plotcomps{2}{1} = [4 10]; % subj 2, cross coh between IC4 and IC10
% clim: [value between 0 and 1] sets color axis for avg coherence plot
% freqscale -- ['log' or 'linear'] for log- or linear-spaced freq bins
% matalpha: p value for individual bootstrap significance for each cross coherence (usually .01).
% binomalpha: if not empty, will mask resulting avg ERSP using binomial
%            probability at the specified alpha level (suggest: 0.000001 (very small)) 
% cond -- [0|1|2|3] '0' means crossf was on single condition, '1' or '2' plots condition 1 or 2, 
%                   respectively, '3' plots 1 minus 2 

function [allcross,times,freqs] = PlotSigAvgCoh(savename,paths,plotcomps,clim,freqscale,matalpha,binomalpha,cond);

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(clim) > 1
        fprintf('Colorlim must be a single value (will scale around zero automatically)');
    end;
    
    allcross = [];    
    for nx = 1:length(plotcomps)
        if ~isempty(plotcomps{nx})
            if strcmp(savename(end-3:end),'.mat')
                s = load([paths{nx},savename]);
            else
                s = load([paths{nx},savename,'.mat']);
            end;     
            if isempty(allcross)
                allcross = zeros(length(s.freqsout),length(s.timesout),0);
            end;

            if cond > 0
                for crs = 1:length(plotcomps{nx})
                    onecross = s.cohercell{find(plotcomps{nx}{crs}(1)==s.complist),find(plotcomps{nx}{crs}(2)==s.complist),cond};
                    oneboot = s.bootmask{find(plotcomps{nx}{crs}(1)==s.complist),find(plotcomps{nx}{crs}(2)==s.complist),cond};                        
                    if cond == 3 
                        onecross(find(onecross < oneboot(:,:,2))) = 0;
                        onecross(find(onecross > oneboot(:,:,1))) = 0;                   
                    else                    
                        oneboot = repmat(oneboot,[1 length(s.timesout)]);
                        onecross(find(onecross <  oneboot)) = 0;           
                    end;
                    allcross(:,:,end+1) = onecross;
                end;
            else
                for crs = 1:length(plotcomps{nx})
                    onecross = s.cohercell{find(plotcomps{nx}{crs}(1)==s.complist),find(plotcomps{nx}{crs}(2)==s.complist)};
                    oneboot = s.bootmask{find(plotcomps{nx}{crs}(1)==s.complist),find(plotcomps{nx}{crs}(2)==s.complist)};          
                    oneboot = repmat(oneboot,[1 length(s.timesout)]);
                    onecross(find(onecross <  oneboot)) = 0;           
                    allcross(:,:,end+1) = onecross;
                end;
            end;             
        end;        
    end;
    [plotmat] = GroupSig(allcross,matalpha,binomalpha,'binom');
    figure;
    if strcmp(freqscale,'linear') % linear spacing
        imagesc(s.timesout, s.freqsout, plotmat,[-clim clim]);hold on;
        set(gca,'ytick',[10:10:s.freqsout(end)]);
    elseif strcmp(freqscale,'log')    % log spacing            
        mylogimagesc(s.timesout, s.freqsout, plotmat,[-clim clim]);hold on;        
    end;
    set(gca,'ydir','norm');hold on;
    set(gca,'ticklength',[.03 .03]);
    plot([0 0],[get(gca,'ylim')],'k-');    
    cbar;
    times = s.timesout;  freqs = s.freqsout;
    
