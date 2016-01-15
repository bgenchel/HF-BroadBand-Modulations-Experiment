% made to plot results of SpecDecomp() with no scalp maps
%
%  PlotSpecDecomp(datpath,complist,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone);
%
% INPUTS:
% datpath: full directory path where dataset can be found
% complist: full list of components/channels in orginal decomposition
% plotcomps: list of components/channels to plot (must be in original 'complist')(default: all)
% factors: list of factor numbers to plot (default: all)
% savedat: filename that float files were saved as (same as input to SpecCoModAnal.m)
% dims = [vector] dimension numbers to plot; [] will plot all
% superpos -- ['y' or 'n'] 'y' will plot a single axis with all comps superimposed.
%             (can only be used for plotting all comps (justone needs to be []) 
% meanvar -- [0 | 1] if 1, will plot mean spectra and trial backprojections, time intensive
% justone: (integer) IM number to plot alone. [] to plot all. Plots mean with spectral variations on top

function PlotSpecDecompTW(datpath,savedat,dims,row,col,pl);
    
    
    
    s = load ([datpath,savedat,'.mat']);
    if ~isfield(s,'pcs')
      s.pcs = s.numrows;
    end;
    sph=floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;   winv = pinv(ws);      
    
    speceig = floatread([datpath,s.eigfile],[s.origdims inf],[],0);
    specwts = speceig*winv;  
    winv = specwts;    
    
    if isempty(dims)
        dims = [1:size(winv,2)];
    end;
    
    % plot the spectral templates:------------------
    lmax = max(winv(:)) + (max(winv(:))*.02); % add 2%
    lmin = min(winv(:)) - (min(winv(:))*.02); % add 2%
    for imm = 1:length(dims)
        im = dims(imm);
        if pl > row*col
            ph=textsc(['Spectral Templates: TW decomposition (',savedat,')'],'title');
            figure; pl=1;
        end;
        sbplot(row,col,pl); pl = pl+1;
        if strcmp(s.freqscale,'log')
            ph = logplot(s.freqs,winv(:,im),2,'k'); hold on;
        elseif strcmp(s.freqscale,'quad')
            ph = quadplot(s.freqs,winv(:,im),2,'k'); hold on;
        else
            ph = plot(s.freqs,winv(:,im),'k-','linewidth',2); hold on;
            set(gca,'xlim',[s.freqs(1) s.freqs(end)]);
        end;
        set(gca,'ylim',[lmin lmax]);
        plot([get(gca,'xlim')],[0 0],'k--');            
        title(['Mode ',int2str(im)]);
    end;
    ph=textsc(['Spectral Templates: TW decomposition (',savedat,')'],'title');
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    

