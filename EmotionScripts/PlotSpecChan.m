%  Plots spectra of specified components from each subject superimposed
%
% [] = PlotSpecs(datfile,datpath,complist,frqlim)
%
% datfile -- [string] filename for each subject containing spectral data and freqs
%            assumes output from 'CalcSpectra.m'
% datpath -- [string] path for datfile to plot (one subplot per subject)
% complist - [cell array of integers] empty brackets will signify which subjects NOT to plot
% frqlim -- [minfrq maxfrq] If specified, will plot only these freqs, all freqs if []
% Will plot all subplots on one page.
%
function [] = PlotSpecs(datfile,datpath,complist,frqlim);
    s = load([datpath,datfile]);
    figure; pl = 1; count = 0;
    for x = 1:length(s.allspec)
        if ~isempty(complist{x})
            count = count+1;
        end;
    end;
    row = ceil(sqrt(count)); col = round(sqrt(count));
    if isempty(frqlim)
        fr = [1:length(s.freqs)];
    else
        fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(2));
    end;    
    for nx = 1:length(s.allspec)
        if ~isempty(complist{nx})
            plotcomps = [1:size(s.allspec{nx},1)];
            cols = winter(length(plotcomps));
            subplot(row,col,pl)
            for cp = 1:length(plotcomps)
                onespec = s.allspec{nx}(plotcomps(cp),fr);
                ph = plot(s.freqs(fr),onespec,'k-'); hold on;
                set(ph,'color',cols(cp,:)); set(ph,'linewidth',1.5);
            end;
            ph = plot(s.freqs(fr),mean(s.allspec{nx}(:,fr),1),'k-'); 
            set(ph,'linewidth',2);
            set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);
            % find peak of mean power spectrum
            peaksurround = 4;pp=1;
            X = mean(s.allspec{nx}(:,fr),1); mxma = [];
            for n = peaksurround+1:length(X)-peaksurround
                if X(n) > X(n-[1:peaksurround])& X(n) > X(n+[1:peaksurround])%
                    mxma(1,pp) = n; pp = pp+1;                
                end;
            end;
            allmax = s.freqs((fr(1)-1) + mxma)';allmax = round(allmax*100)/100;            
            title(['Subj ',int2str(nx),'; ',num2str(allmax),' Hz']); pl = pl+1;
        end;
    end;
    axcopy
                
        
        
        