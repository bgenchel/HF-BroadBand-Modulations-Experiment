% plots results from DecompMat()
%
%
%
%
%
% keeptrack -- [vector] of sorting indices corresponding to number of epochs
%              or, if scalar, uses the corresponding column of saved keeptrack matrix.

function [P PH] = PlotDecompMatWts(wtsphname,writepath,plotfacs,clabels,vertlines,ttl,keeptrack)
    
    
    
    alpha = .01;
    shuffnum = 500;
    
    s = load ([writepath,wtsphname,'.mat']);
    sph=floatread([writepath,wtsphname,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([writepath,wtsphname,'.wts'],[s.pcmat s.numtrials],[],0); 
    icamat=floatread([writepath,wtsphname,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;winv = pinv(ws);    acts = ws*icamat;   
    if length(keeptrack< 2)
        keeptrack = s.keeptrack(:,keeptrack);
    end;
    
    if isempty(plotfacs)
        plotfacs = [1:size(winv,2)];
    end;

    if ~isempty(s.pcmat)
        if length(s.rowmeans) == length(s.freqs)*length(s.times)
            erspeig = floatread([writepath,wtsphname,'ERSPEIG.fdt'],[length(s.freqs)*length(s.times) s.pcmat],[],0);
            erspdat = erspeig*winv(1:s.pcmat,:);
        elseif length(s.rowmeans) == length(s.freqs)
            erspeig = floatread([writepath,wtsphname,'ERSPEIG.fdt'],[length(s.times) inf],[],0);
            windat = erspeig*winv;% WT decomp, no context
        end;
    end;
    if ~isempty(s.pcidx)
        idxeig = floatread([writepath,wtsphname,'IDXEIG.fdt'],[length(s.rowmeanidx) s.pcidx],[],0);
        idxdat = idxeig*winv(s.pcmat+1:end,:);
        %if ~isempty(shuffnum) % plot weighted mean instead
        %    idxdat = idxeig*icamat(s.pcmat+1:end,:);
        %end;
    end;
    possvals = unique(keeptrack);
    figure; 
    row = 3; col = 2; pl = 1;
    for dim = 1:length(plotfacs)
        if pl > row*col
            if ~isempty(shuffnum)
                textsc([wtsphname,'--',ttl,'-- Weighted Means'],'title');
            else
                textsc([wtsphname,'--',ttl,'-- Winv Templates'],'title');
            end;
            figure; pl = 1;
        end;
        %% stats on groups
        if length(possvals) < 3 % ttest if 2
            [H P(1,dim)] = ttest2(acts(dim,find(keeptrack==possvals(1))),acts(dim,find(keeptrack==possvals(2))),alpha);
            if ~isempty(shuffnum)
                % permutation stats
                meansubj = 'on';
                if strcmp(meansubj,'on')
                    clustsubjs = unique(s.keeptrack(:,1));
                    for sbj = 1:length(clustsubjs)
                        subjmeans(1,sbj) = median(acts(dim,find(s.keeptrack(:,1) == clustsubjs(sbj))));
                        newtrack(1,sbj) = unique(s.keeptrack(find(s.keeptrack(:,1)==clustsubjs(sbj)),3));
                    end;
                    for b = 1:shuffnum
                        newidx = shuffle(newtrack);
                        bootdiff(1,b) = mean(subjmeans(find(newidx==possvals(1)))) - mean(subjmeans(find(newidx==possvals(2))));
                    end;                
                else
                    for b = 1:shuffnum
                        newidx = shuffle(keeptrack);
                        bootdiff(1,b) = mean(acts(dim,find(newidx==possvals(1)))) - mean(acts(dim,find(newidx==possvals(2))));
                    end;
                end;
                bootdiff = sort(bootdiff);
                minmask = bootdiff(ceil(shuffnum*alpha));
                maxmask = bootdiff(length(bootdiff) - round(shuffnum*alpha));
                maskedmat = mean(acts(dim,find(keeptrack==possvals(1)))) - mean(acts(dim,find(keeptrack==possvals(2)))); % real diff   
                if maskedmat < minmask | maskedmat > maxmask
                    PH(1,dim) = 1;
                else
                    PH(1,dim) = 0;
                end;            
            end;
        end;         
        sbplot(row,col,pl)
        %[vals sortidx] = sort(keeptrack);
        cols = jet(length(possvals)); cols = [[1 0 0];[0 0 1]];
        plot([0 length(subjmeans)+1],[0 0],'k-'); hold on;
        startx = 1;
        for v = 1:length(possvals)
            
           %ph = plot([startx:startx+length(find(vals==possvals(v)))-1],acts(dim,find(keeptrack==possvals(v))),'k.');hold on;
            %set(ph,'color',cols(v,:)); startx = startx+length(find(keeptrack==possvals(v)));
           ph = plot([startx:startx+length(find(newtrack==possvals(v)))-1],subjmeans(find(newtrack==possvals(v))),'k.');hold on;
            set(ph,'color',cols(v,:)); startx = startx+length(find(newtrack==possvals(v)));
        end;
        set(gca,'xlim',[1 startx]);
        title(['Dim ',int2str(dim)]);
        pl = pl+1;
    end;
