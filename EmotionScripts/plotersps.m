% loads an ERSP .mat saved from 'runtimef.m' and plots with optional masking
%
% [] = plotersps(comop_ersp,pathname,maskon,tmrange,frqrange,colorlim,complist,mapset,title,ersp_boot,freqs,times);
% datfile: full path and name of .mat to load with variable comp_ersp, freqs, times and ersp_boot contained. (string)
% pathname: full path where mapset can be found. (string)
% maskon: if 1, then will mask ersp with ersp_boot significance level.
% tmrange: [mintime maxtime] to plot;
% frqrange: [minfreq maxfreq] to plot.
% colorlim: constant color max and min to plot for all components.
% complist: vector of numbers corresponding to components to plot. This algorithm assumes that all the comps in 
%           the comp_ersp are the desired components (no spaces), so this argument is for plot titles.
% mapset: name of dataset to load for associated scalp maps. (string)
% title: desired title to print to top of figure.
% 
function [] = plotersps(datfile,pathname,maskon,tmrange,frqrange,colorlim,complist,mapset,title);
    
    s = load(datfile);
    fr = find(s.freqs > frqrange(1) & s.freqs < frqrange(2));
    tm = find(s.times > tmrange(1) & s.times < tmrange(2));

    if ~isempty(mapset)
        row = ceil(sqrt(size(s.comp_ersp,3)*2))+1;col=ceil(sqrt(size(s.comp_ersp,3)*2))+1;
        if ~iseven(col)
            col = col+1; row = row-1;
        end; 
        EEG = pop_loadset( mapset,pathname,'all' );
        figure; pl = 1;
        for mp = 1:length(complist)
            subplot(row,col,pl)
            topoplot(EEG.icawinv(:,complist(mp)),EEG.chanlocs,'electrodes','off');
            pl = pl+2;
            %title(['Comp: ',int2str(complist(cp))]);
        end;
        pl = 2;
        for cp = 1:length(complist)
            minmask = repmat(s.ersp_boot(fr,1,cp),[1 size(s.comp_ersp,2)]);
            maxmask = repmat(s.ersp_boot(fr,2,cp),[1 size(s.comp_ersp,2)]);
            plotersp = s.comp_ersp(fr,:,cp);
            if maskon == 1
                plotersp(find(plotersp > minmask & plotersp < maxmask)) = 0;    
            end;
            subplot(row, col, pl)
            imagesc(times(tm),freqs(fr),plotersp(:,tm),[-colorlim colorlim]);
            pl = pl+2; 
        end;colorbar;
        ph = textsc(title,'title');set(ph,'fontsize',14);
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    else
        figure; pl = 1;  
        row = ceil(sqrt(size(s.comp_ersp,3)));col=ceil(sqrt(size(s.comp_ersp,3)))+1;
        for cp = 1:length(complist)
            minmask = repmat(s.ersp_boot(fr,1,cp),[1 size(s.comp_ersp,2)]);
            maxmask = repmat(s.ersp_boot(fr,2,cp),[1 size(s.comp_ersp,2)]);
            plotersp = s.comp_ersp(fr,:,cp);
            if maskon == 1
                plotersp(find(plotersp > minmask & plotersp < maxmask)) = 0;    
            end;
            subplot(row, col, cp)
            imagesc(times(tm),freqs(fr),plotersp(:,tm),[-colorlim colorlim]);pl = pl+1;
            %title(['Comp: ',int2str(complist(cp))]);
        end;colorbar;
    end;
    
