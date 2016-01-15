% plots single subject results from 'runtimef.m'
%
% PlotRuntimefOUT(datafile,mapset,pathname,complist,origlist,mask,clim,frqlim,tmlim,ttls,sigdiff);
%
% datafile: full path and file name with .mat extension
% sigdiff: if 1 then assumes two conditions were run, 0 if single condition timef
% mask: if 1 then mask for significance from erspboot
% mapset: dataset name with .set extension from which scalp map data can be extracted
% path: full path where mapset can be found
% clim: color limit to use for all plots. will center the given number around 0 so green is 0
% frqlim: [min max] Frequency limits for plotting
% tmlim: [min max] Time limits to plot (ms)
% ttls: {string(s)}; Cell array of strings for single or two conditions timef plots (1 string per cond)
% complist: components to plot ersp for
% invdiff -- if 1, will plot the inverse (*-1) of the 'difference' image (sigdiff must be 1 as well)
% **  if time warping was used, will look for a vector of latencies in the data variables called 'warpevs'
% requires a field in datafile called 'complist' with all original ICs decomposed.
%
function [freqscale] = PlotRuntimefOUT(datafile,mapset,pathname,plotcomps,mask,clim,frqlim,tmlim,ttls,sigdiff,invdiff,plotparams);
    
row = plotparams(1); col = plotparams(2);pl = plotparams(3);
            
    EEG = pop_loadset(mapset,pathname); 
    s = load ([pathname,datafile]);
    %s = load ([datafile]);
    fr=find(s.freqs >= frqlim(1) & s.freqs<=frqlim(2));  
    tms=find(s.times > tmlim(1) & s.times<tmlim(2));
    if isfield(s,'warpevs')
        cols = lines(length(s.warpevs));
    end;    
    if isempty(plotcomps)
      plotcomps = s.complist;
    end;
    %%%%%%  Check for log frequency spacing  %%%%%%%%%%%%%%%%
    %%%  need to fix this
    if isfield(s,'freqscale') % then this is easy
        freqscale = s.freqscale;
    else
        tst1 = sqrt(s.freqs);
        tst2 = log(s.freqs);
        lfr1 = linspace(sqrt(s.freqs(1)),sqrt(s.freqs(end)),length(s.freqs));
        lfr2 = linspace(log(s.freqs(1)),log(s.freqs(end)),length(s.freqs));

        if tst1 == lfr1
            freqscale = 'quad'; % quadratic spacing
        elseif tst2 == lfr2 
            freqscale = 'log'; % log spacing
        else
            freqscale = 'linear'; % linear spacing
        end;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%
    if sigdiff == 1        
        %figure; row = ceil(length(plotcomps)/2); col = 8;  pl = 1;
        for cp = 1:length(plotcomps)
          if pl > row*col
            figure; pl=1;
          end;
          sbplot(row,col,pl)
            topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.65); pl = pl+1;
            title(['IC ',int2str(plotcomps(cp))]);
            for cond = 1:3
                ep = s.comp_ersp{find(ismember(s.complist,plotcomps(cp))),cond}; % (fr,tms)
                bt = s.ersp_boot{find(plotcomps(cp)==s.complist),cond}; % (fr,2),(fr,tms,2)
                if cond == 3
                    if mask == 1
                        ep(find(ep > bt(:,:,1) & ep < bt(:,:,2))) = 0; % mask for sig
                    end;
                    if invdiff == 1
                        ep = ep * -1;
                    end;
                    sbplot(row,col,pl)
                    if strcmp(freqscale,'quad')
                        quadimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+1;
                    elseif strcmp(freqscale,'log')
                        mylogimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+1;
                    else
                        imagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); hold on;pl = pl+1;
                    end;
                    set(gca,'ydir','norm');plot([0 0],[get(gca,'ylim')],'k-');
                    set(gca,'ticklength',[.05 .05]);
                    if isfield(s,'warpevs')
                        for wp = 1:length(s.warpevs)
                            ph = plot([s.warpevs(wp) s.warpevs(wp)],[get(gca,'ylim')],'k-');
                            set(ph,'color',cols(wp,:));
                        end;
                    end;    
                    if cp < 3
                        title('Difference');
                    end;
                else
                    if mask == 1
                        ep(find(ep > repmat(bt(:,1),[1 length(s.times)]) & ep < repmat(bt(:,2),[1 length(s.times)]))) = 0;
                    end;
                    sbplot(row,col,pl)
                    if strcmp(freqscale,'quad')
                        quadimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+1;
                    elseif strcmp(freqscale,'log')
                        mylogimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+1;
                    else
                        imagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); hold on;pl = pl+1;
                    end;
                    set(gca,'ydir','norm');plot([0 0],[get(gca,'ylim')],'k-');
                    set(gca,'ticklength',[.05 .05]);
                    if isfield(s,'warpevs')
                        for wp = 1:length(s.warpevs)
                            ph = plot([s.warpevs(wp) s.warpevs(wp)],[get(gca,'ylim')],'k-');
                            set(ph,'color',cols(wp,:));
                        end;
                    end;    
                    if cp < 3 & cond == 1
                        title(ttls{1});
                    elseif cp < 3 & cond == 2
                        title(ttls{2});
                    end;
                end;     
            end;
        end;cbar;
        ph = textsc([pathname(end-4:end-1),'/',datafile],'title'); %ph = set(ph,'fontsize',14);
    else        
        figure; row = ceil(sqrt(length(plotcomps)*3)); col = ceil(sqrt(length(plotcomps)*2.5));  pl = 1;
        if col ~= 3 & col ~=6 & col ~= 9 & col ~= 12& col ~= 15
            while col ~= 3 & col ~=6 & col ~= 9 & col ~= 12& col ~= 15
                col = col+1; 
                row = row-1;
            end;
        end;        
        for cp = 1:length(plotcomps)
            sbplot(row,col,pl)
            topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.65); pl = pl+1;
            title(['Comp. ',int2str(plotcomps(cp))]);
            ep = s.comp_ersp(:,:,find(plotcomps(cp)==s.complist)); % (fr,tms)
            bt = s.ersp_boot(:,:,find(plotcomps(cp)==s.complist)); % (fr,2,cp)
            if mask == 1
                ep(find(ep > repmat(bt(:,1),[1 length(s.times)]) & ep < repmat(bt(:,2),[1 length(s.times)]))) = 0;
            end;
            sbplot(row,col,[pl pl+1])            
                   if strcmp(freqscale,'quad')
                        quadimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+2;
                    elseif strcmp(freqscale,'log')
                        mylogimagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+2;
                    else
                        imagesc(s.times(tms),s.freqs(fr),ep(fr,tms),[-clim clim]); 
                        hold on;pl = pl+2;
                    end;
            set(gca,'ydir','norm');plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'ticklength',[.05 .05]);
            if isfield(s,'warpevs')
                for wp = 1:length(s.warpevs)
                    ph = plot([s.warpevs(wp) s.warpevs(wp)],[get(gca,'ylim')],'k-');
                    set(ph,'color',cols(wp,:));
                end;
            end;    
            if pl < col+2
                title(ttls);
            end;
        end;cbar;
        ph = textsc([pathname(end-4:end-1),'/',datafile],'title'); %ph = set(ph,'fontsize',14);
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
    
                
        
