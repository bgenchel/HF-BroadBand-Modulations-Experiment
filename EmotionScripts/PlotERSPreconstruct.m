% Uses TW spectral decomp and reconstructs an ERSP by filtering 
% full data matrix
% plotspecs -- [row,col,place]

function [reersp,place] = PlotERSPreconstruct(savedat,datset,datpath,plotdims,tmlim,frqlim,plotspecs);
    
    
    
    s = load ([datpath,savedat,'.mat']);
    sph=floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);
    dat=floatread([datpath,savedat,'DAT.fdt'],[s.numrows s.numframes],[],0);    
    ws = wts*sph;    winv = pinv(ws); 
    
    speceig = floatread([datpath,s.eigfile],[s.numrows s.pcs],[],0);
    specwts = speceig*winv;  
    winv = specwts;  
    ws = pinv(winv);
    acts = ws*dat;
    
    if isempty(plotdims)
        plotdims = [1:size(winv,2)];
    end;
    
    plotpnts = unique(s.keeptrack(:,3));% will need to define differently later if you try to filter new data
    condidx = [];ntrials= [];
    for c = 1:length(plotpnts):size(s.keeptrack,1)
        condidx = [condidx s.keeptrack(c,1)];
        ntrials = [ntrials s.keeptrack(c,2)];
    end;
    condidx = [condidx size(dat,2)];
    if isfield(s,'datset')
        condttls = s.datset;
    else
        condttls = 'Reconstructed ERSP';
    end;
    % now, filter raw spectral data by this ws for each time point of the ERSP
    % EEG = pop_loadset( datset,datpath);
    % would have to now do spectral decomposition for each time
    % point of interest
    row = plotspecs(1); col = plotspecs(2); place = plotspecs(3);
    for dm = 1:length(plotdims)
        dim = plotdims(dm);
        for cond = 1:length(condidx)-1
            reersp{cond} = zeros(length(s.freqs),length(plotpnts),0);
            for tm = 1:length(plotpnts)
                pl = 1;
                for epoch = condidx(cond)+ntrials(cond)*(tm-1):(condidx(cond)-1)+ntrials(cond)*tm
                    dimwts = ws*dat(:,epoch);
                    reersp{cond}(:,tm,pl) = dimwts(dim).*winv(:,dim);
                    pl=pl+1;
                end;
            end;    
        end;
        %figure;row = round(sqrt(length(condidx)));col = ceil(sqrt(length(condidx))); 
        clim = max(max(abs(mean(reersp{1},3))));
        for cond = 1:length(condidx)-1
            if place > row*col
                figure; place = 1;
            end;
            sbplot(row,col,place);
            if strcmp(s.freqscale,'quad')
                quadimagesc(plotpnts,s.freqs,mean(reersp{cond},3),[-clim clim]);
            elseif strcmp(s.freqscale,'log')
                mylogimagesc(plotpnts',s.freqs,mean(reersp{cond},3),[-clim clim]);
            end;
            title(condttls{cond}); place = place+1;
        end;
        %textsc(['Dimension ',int2str(dim)],'title');
    end;
    
