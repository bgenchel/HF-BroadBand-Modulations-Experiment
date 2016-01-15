% plots the PCA spectral data that is passed to ICA for 
% decompostion in the TW fashion. For SPECTRA only, not
% ERSPs

function PlotTWspecPCs(savedat,datpath);
    
    
    
    s = load([datpath,savedat,'.mat']);  

    if isfield(s,'pcs') % if no context:
        data = floatread([datpath,savedat,'.fdt'],[s.pcs inf],[],0);
        wts = floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);
        sph = floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0);
        ws = wts*sph; winv = inv(ws);
        speceig = floatread([datpath,s.eigfile],[s.numrows s.pcs],[],0);
        specwts = speceig*winv;  
        
    else % there is context, so load then separate
        data = floatread([datpath,savedat,'.fdt'],[s.numrows inf],[],0);
        data(s.pcersp+1:end,:) = []; % get rid of context
        wts = floatread([datpath,savedat,'.wts'],[s.numrows s.numrows],[],0);
        sph = floatread([datpath,savedat,'.sph'],[s.numrows s.numrows],[],0);
        ws = wts*sph; winv = inv(ws);
        speceig = floatread([datpath,s.eigfile],[length(s.mntfpwr) s.pcersp],[],0);
        specwts = speceig*winv(1:s.pcersp,:);  
    end;
    
    if isfield(s,'times') % if it's an ersp
        figure; row = round(sqrt(s.pcersp)); col = ceil(sqrt(s.pcersp));  
        clim = max(speceig(:));
        for dm = 1:s.pcersp
            sbplot(row,col,dm);
            if strcmp(s.freqscale,'quad')
                quadimagesc(s.times,s.freqs,reshape(speceig(:,dm),length(s.freqs),length(s.times)),[-clim clim]);
            elseif strcmp(s.freqscale,'log')
                mylogimagesc(s.times,s.freqs,reshape(speceig(:,dm),length(s.freqs),length(s.times)),[-clim clim]);
            end;
            title(['PC ',int2str(dm)]);
        end;        cbar;
        textsc(['ERSP templates after PCA decomposition'],'title');        
        figure; row = round(sqrt(s.numrows)); col = ceil(sqrt(s.numrows));  
        clim = max(specwts(:))-.4*max(specwts(:));
        for dm = 1:s.numrows
            %clim = max(specwts(:,dm));
            sbplot(row,col,dm);
            if strcmp(s.freqscale,'quad')
                quadimagesc(s.times,s.freqs,reshape(specwts(:,dm),length(s.freqs),length(s.times)),[-clim clim]);
            elseif strcmp(s.freqscale,'log')
                mylogimagesc(s.times,s.freqs,reshape(specwts(:,dm),length(s.freqs),length(s.times)),[-clim clim]);
            end;
            title(['IC ',int2str(dm)]);
        end;        cbar;
        if size(specwts,2) ~= size(speceig,2)
            textsc(['ERSP templates after ICA decomposition (minus context portion)'],'title');        
        else
            textsc(['ERSP templates after ICA decomposition'],'title');        
        end;
    else        
        figure; row = round(sqrt(s.pcs)); col = ceil(sqrt(s.pcs));    
        for dm = 1:s.pcs
            sbplot(row,col,dm);
            if strcmp(s.freqscale,'quad')
                quadplot(s.freqs,speceig(:,dm)/std(speceig(:,dm)),2,'b'); hold on;
                ph = quadplot(s.freqs,specwts(:,dm)/std(specwts(:,dm)),2,'r');
            elseif strcmp(s.freqscale,'log')
                logplot(s.freqs,speceig(:,dm)/std(speceig(:,dm)),2,'b'); hold on;
                ph = logplot(s.freqs,specwts(:,dm)/std(specwts(:,dm)),2,'r');
            else
                plot(s.freqs,speceig(:,dm)/std(speceig(:,dm)),'b','linewidth',2); hold on;
                ph = plot(s.freqs,specwts(:,dm)/std(specwts(:,dm)),'r','linewidth',2);
                
            end;
        end;
        legend({'PCA','ICA'});
        textsc(['Spectral templates after PCA/ICA decomposition'],'title');        
    end;
