% takes cluster information and plots cross coherences in subjects spanning the two input clusters
% Plots the mean scalp maps and coherence (the latter masked by binomial probability across subjects if 'mask' on
%
%     [clim] = PlotClustCoh(datset,savename,fullpaths,clust1,clust2,mask,cond);
%
% INPUTS:
% datset -- [string] .set that contains ICA weights to plot scalp maps
% savename -- [string] name of .mat where cross coherence data is saved
%                      (same as that given to runcrossf())
% fullpaths -- [cell array] for each subject, string indicating full data path
% clust1 -- [cell array] vector of clustered components for each subject
% clust2 -- [cell array] vector of clustered components for each subject
% tmlim -- [minms maxms] in milliseconds
% frqlim -- [min max] in Hz
% clim -- [max coh] uses the specified value for the color scale limits (will make symmetrical
%                   around zero so green = zero always. [] will make heuristic value from data.
% mask -- [0|1] if 1, will mask coherences using bootstrap matrix
% cond -- [0|1|2|3] if 0: single-condition crossf, 1: condition 1, 2: condition 2
%                   3: coherence difference between cond 1 and 2.

function [clim] = PlotClustCoh(datset,savename,fullpaths,clust1,clust2,tmlim,frqlim,clim,mask,cond);
    
    
    howmany = 0;
    for nx = 1:length(clust1)
        if ~isempty(clust1{nx})
            if ~isempty(clust2{nx})
                for cp1 = 1:length(clust1{nx})
                    for cp2 = 1:length(clust2{nx})
                        howmany = howmany + 1;                    
                    end;
                end;
            end;
        end;
    end;
    howmany = howmany*3;
    possibs = [3,6,9,12,15,18,21]; % make columns be multiple of 3
    figure; pl = 1;
    col = round(sqrt(howmany)); 
    while ~ismember(col,possibs)
        col = col+1;
    end;    
    row =ceil(howmany/(col))+1;

    lim = .5; allcrosses = []; 
    pt = 1; clear map1 map2 
    for nx = 1:length(clust1)
        if ~isempty(clust1{nx})
            if ~isempty(clust2{nx})
                EEG = pop_loadset(datset,fullpaths{nx});
                if strcmp(savename(end-3:end),'.mat')
                    s = load([fullpaths{nx},savename]);
                else
                    s = load([fullpaths{nx},savename,'.mat']);
                end;   
                if isfield(s,'comps')
                    s.complist = s.comps;
                    s.bootmask = s.sigdiffs;
                end;
                for cp1 = 1:length(clust1{nx})
                    for cp2 = 1:length(clust2{nx})
                        subjttls(pt,:) = [nx,clust1{nx}(cp1),clust2{nx}(cp2)];
                        sbplot(row,col,pl)
                        [h, grid1(:,:,pt), plotrad, xmesh, ymesh]=topoplot(EEG.icawinv(:,clust1{nx}(cp1)),EEG.chanlocs,'plotrad',.6,'electrodes','off');pl = pl+2;
                        title(['IC1: ',int2str(clust1{nx}(cp1))]);
                        sbplot(row,col,pl)
                        [h, grid2(:,:,pt), plotrad, xmesh, ymesh]=topoplot(EEG.icawinv(:,clust2{nx}(cp2)),EEG.chanlocs,'plotrad',.6,'electrodes','off');pl = pl+1;
                        title(['IC2: ',int2str(clust2{nx}(cp2))]);
                        if cond > 0
                            oneic = find(clust1{nx}(cp1)==s.complist);
                            otheric = find(clust2{nx}(cp2)==s.complist);
                            if oneic < otheric
                                onecross = s.cohercell{oneic,otheric,cond};
                            else
                                onecross = s.cohercell{otheric,oneic,cond};
                            end;
                            if mask == 1
                                if oneic < otheric
                                    oneboot = s.bootmask{oneic,otheric,cond}; 
                                else
                                    oneboot = s.bootmask{otheric,oneic,cond};
                                end;
                                if cond == 3                                     
                                    onecross(find(onecross < oneboot(:,:,2))) = 0;
                                    onecross(find(onecross > oneboot(:,:,1))) = 0;                   
                                else                    
                                    oneboot = repmat(oneboot,[1 length(s.timesout)]);
                                    onecross(find(onecross <  oneboot)) = 0;           
                                end;
                            end;
                        else
                            oneic = find(clust1{nx}(cp1)==s.complist);
                            otheric = find(clust2{nx}(cp2)==s.complist);
                            if oneic < otheric
                                onecross = s.cohercell{oneic,otheric};
                            else
                                onecross = s.cohercell{otheric,oneic};
                            end;
                            if mask == 1
                                if oneic < otheric
                                    oneboot = s.bootmask{oneic,otheric};
                                else
                                    oneboot = s.bootmask{otheric,oneic};
                                end;
                                oneboot = repmat(oneboot,[1 length(s.timesout)]);
                                onecross(find(onecross <  oneboot)) = 0;           
                            end;
                        end;
                        if isempty(allcrosses)
                            allcrosses = onecross; pt = pt+1;
                        else
                            allcrosses(:,:,pt) = onecross;pt = pt+1; 
                        end;
                        
                    end;
                end;
            end;
        end;
    end;
    if ~isempty(tmlim)
        tms = find(s.timesout > tmlim(1) & s.timesout < tmlim(2));
    else
        tms = [1:length(s.timesout)];
    end;
    if ~isempty(frqlim)
        frs = find(s.freqsout > frqlim(1) & s.freqsout < frqlim(2));
    else
        frs = [1:length(s.freqsout)];
    end;    
    
    s.timesout = s.timesout/1000;
    lim = max(max(max(abs(allcrosses(frs,tms,:)))));
    if lim == 0
        lim = 1;
    end;
    pl = 2;
    for cp = 1:size(allcrosses,3)
        sbplot(row,col,pl);
        if strcmp(s.freqscale,'quad')
            quadimagesc(s.timesout(tms), s.freqsout(frs), allcrosses(frs,tms,cp),[-lim lim]);hold on;
        elseif  strcmp(s.freqscale,'log')                
            mylogimagesc(s.timesout(tms), s.freqsout(frs), allcrosses(frs,tms,cp),[-lim lim]);hold on;
        else
            imagesc(s.timesout(tms), s.freqsout, allcrosses(frs,tms,cp),[-lim lim]);hold on;               
            set(gca,'ytick',[10:10:s.freqsout(frs(end))]);
        end;
        set(gca,'ydir','norm');hold on;
        plot([0 0],[get(gca,'ylim')],'k-');
        if pl < row*col-(col-1)
            set(gca,'xticklabel',[]);
        end;
        title(['Subj ',int2str(subjttls(cp,1))]);    pl = pl+3;     
    end;
    cbar;
    
    if mask == 1
        [maskedmat] = GroupSig(allcrosses,s.alpha,.0000001,'binom');
    else
        maskedmat = allcrosses;
    end;    
    x = mean(maskedmat(frs,tms,:),3); 
    if isempty(clim)
        newlim = max(abs(x(:))); 
        if newlim == 0
            newlim = 1;
        end;
    else
        newlim = clim;
    end;
    clim = newlim;
    
    sbplot(row,col,pl-1)
    toporeplot(mean(grid1,3),'plotrad',.6, 'intrad',plotrad);
    sbplot(row,col,pl)
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;
    else
        imagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;               
    end;
    set(gca,'ydir','norm');hold on;
    plot([0 0],[get(gca,'ylim')],'k-');
    title(['Mean Cross Coh']);
    cbar;       
    sbplot(row,col,pl+1)
    toporeplot(mean(grid2,3),'plotrad',.6, 'intrad',plotrad);

    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
    figure;
    if strcmp(s.freqscale,'quad')
        quadimagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;
    elseif  strcmp(s.freqscale,'log')                
        mylogimagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;
    else
        imagesc(s.timesout(tms), s.freqsout(frs), x,[-newlim newlim]);hold on;               
    end;
    set(gca,'ydir','norm');hold on;
    plot([0 0],[get(gca,'ylim')],'k-'); set(gca,'ticklength',[.04 .04]);
    title(['Mean Cross Coh']);
    cbar;       
    
    axcopy
