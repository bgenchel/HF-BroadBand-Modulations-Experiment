% plots output of runcrossf() for a single subject (specified components)
%
% [allcross, crossidx]] = PlotCrossCohAng(datset,savename,fullpath,plotcomps,mask,cond)
%
%
% INPUTS:
% datset -- [string] dataset (.set) that contains ICA weights to plot a scalp map
% savename -- [string] name of .mat file given to runcrossf() to save data as in each path
% fullpath -- [string] full directory path where data is saved.
% plotcomps -- [vector] component indices to plot
% mask -- [0|1] 1 will use bootstrap matrix to mask for significance. 0 plots raw coherence.
% freqscale -- ['log' or 'linear'] for log- or linear-spaced freq bins
% cond -- [0|1|2|3] '0' means crossf was on single condition, '1' or '2' plots condition 1 or 2, 
%                   respectively, '3' plots 1 minus 2 
%
% OUTPUT:
% allcross -- [matrix of freqsout x timesout x # of crosses plotted] gives the UNmasked cross 
%             coherence for each pair requested.
% crossidx - [# of crosses x 2 matrix] gives the component pair indices for each coherence
%            in 'allcross' 

function [allcross, crossidx] = PlotCrossCohAng(datset,savename,fullpath,plotcomps,mask,cond)
    
    lim = pi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_loadset( datset,fullpath);
    if strcmp(savename(end-3:end),'.mat')
        s = load([fullpath,savename]);
    else
        s = load([fullpath,savename,'.mat']);
    end;        
    if isfield(s,'freqscale')
        freqscale = s.freqscale;
    else
        freqscale = 'linear'; % default
    end;
    
    figure;
    row = length(plotcomps);
    col = length(plotcomps);
    strows = [col+1:col:(col+1)*(row-1)];
    pl = 2;
    for cp = 2:length(plotcomps)
        sbplot(row,col,pl)
        topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.7);pl = pl+1;
        title(int2str(plotcomps(cp)));
    end;
    for cp = 1:length(strows)
        sbplot(row,col,strows(cp))
        set(gca,'fontsize',5);
        topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs,'electrodes','off','plotrad',.7);pl = pl+1;
        title(int2str(plotcomps(cp)));
    end;
    allcross = zeros(length(s.freqsout),length(s.timesout),0);
    crossidx = zeros(0,2);
    s.timesout = s.timesout/1000;
    if cond > 0
        for cp1 = 1:length(plotcomps)-1
            xx = 0;
            for cp2 = cp1+1:length(plotcomps)
                onecross = s.crsangcell{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist),cond};
                onecrosscoh = s.cohercell{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist),cond};
                if mask == 1
                    oneboot = s.bootmask{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist),cond};
                allcross(:,:,end+1) = onecross;
                crossidx(end+1,:) = [plotcomps(cp1) plotcomps(cp2)];
                    if cond == 3 
                        onecross(find(onecross < oneboot(:,:,2))) = 0;
                        onecross(find(onecross > oneboot(:,:,1))) = 0;                   
                    else                    
                        oneboot = repmat(oneboot,[1 length(s.timesout)]);
                        onecross(find(onecrosscoh <  oneboot)) = 0;           
                    end;
                end;
                sbplot(row,col,strows(cp1)+cp1+xx); xx = xx+1;
                if strcmp(freqscale,'quad')
                    quadimagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;
                elseif  strcmp(freqscale,'log')                
                    mylogimagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;
                else
                    imagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;               
                    set(gca,'ytick',[10:10:s.freqsout(end)]);
                    set(gca,'ydir','norm');hold on;
                end;
                set(gca,'ticklength',[.03 .03]);
                if cp1 ~= cp2-1
                    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
                end;
                plot([0 0],[get(gca,'ylim')],'k-');
            end;
        end;
        cbar;        
    else
        for cp1 = 1:length(plotcomps)-1
            xx = 0;
            for cp2 = cp1+1:length(plotcomps)
                onecross = s.crsangcell{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist)};
                onecrosscoh = s.cohercell{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist)};
                allcross(:,:,end+1) = onecross;
                crossidx(end+1,:) = [plotcomps(cp1) plotcomps(cp2)];
                if mask == 1
                    oneboot = s.bootmask{find(plotcomps(cp1)==s.complist),find(plotcomps(cp2)==s.complist)};
                    oneboot = repmat(oneboot,[1 length(s.timesout)]);
                    onecross(find(onecrosscoh <  oneboot)) = 0;           
                end;
                sbplot(row,col,strows(cp1)+cp1+xx); xx = xx+1;
                if strcmp(freqscale,'quad')
                    quadimagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;
                elseif  strcmp(freqscale,'log')                
                    mylogimagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;
                else
                    imagesc(s.timesout, s.freqsout, onecross,[-lim lim]);hold on;               
                    set(gca,'ytick',[10:10:s.freqsout(end)]);
                    set(gca,'ydir','norm');hold on;
                end;
                set(gca,'ticklength',[.03 .03]);
                if cp1 ~= cp2-1
                    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
                end;
                plot([0 0],[get(gca,'ylim')],'k-');
            end;
        end;
        cbar;       
    end;
    textsc(['All Cross Coherence angles for ',fullpath(end-4:end-1),' Condition ',int2str(cond),' from ',savename],'title');
    axcopy
