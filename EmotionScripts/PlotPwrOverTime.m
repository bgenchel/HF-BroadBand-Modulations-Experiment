%
%
% matfile -- [string] name of .mat file containing all time/freq info to plot
% fullpaths -- [cell array] of strings with full directory paths for each subject
% frqlim -- [minHz maxHz] frequency band to integrate over
% tmlim -- [minMS maxMS] time window to integrate over (in ms)
% plotcomps -- [cell array] of component indices for each subject to plot
% mask -- ['on' or 'off']
%


function [pwr] = PlotPwrOverTime(matfile,fullpaths,frqlim,tmlim,plotcomps);
    
    
    s = load ([fullpaths{1},matfile]);
    if isempty(tmlim)
        tmlim = [s.times(1) s.times(end)];
    end;
    tms = find(s.times > tmlim(1) & s.times < tmlim(2));
    pwr = zeros(length(find(s.freqs>frqlim(1)&s.freqs<frqlim(2))),length(tms),0);;
    for nx = 1:length(fullpaths)
        if ~isempty(plotcomps{nx})
            s = load ([fullpaths{nx},matfile]);
            for c = 1:length(plotcomps{nx})
                ic = find(s.complist == plotcomps{nx}(c));
                pwr(:,:,end+1) = s.comp_ersp(find(s.freqs>frqlim(1)&s.freqs<frqlim(2)),tms,ic);
            end;
        end;
        fprintf('.');
    end;
    fprintf('\ndone.\n');
    
    pwr = mean(pwr,1);
    pwr = squeeze(pwr);
    
    figure; 
    for ic = 1:size(pwr,3)
        ph = plot(s.times(tms),pwr,'.');
    end;
    
                
            
    
