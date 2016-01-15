% takes a 2D matrix of EEG data. Calculates power in each epoch (2nd dim) and 
% performs spectral decomposition in the TW fashion
%
%
% INPUTS:
% alldat -- [2D matrix] dimensions: freqs x time*epochs for one IC (designed for long epochs)
% datpath -- [string] data file path where 'savedat.fdt' will be saved
% savedat -- [string] name to save data float file as
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% freqscale -- ['linear','quad' or 'log'] 'quad' makes quadratically spaced freq bins
% pcs -- [integer] number of PCA dimensions to retain in final spectral decomposition 
% srate -- [integer] sampling rate of the input data matrix 'alldat'

function [numtrials,numframes,freqs, pcs, rowmeans] = SpecDecomp(alldat,datpath,wincenters,savedat,frqlim,freqscale,pcs,fac,srate);

    nfreqs = 40;
    hwin = hanning(srate)'; pl = 1; clear keeptrack pwr frqs
    fprintf('\nCalculating pwr for each epoch...\n');
    for epoch = 1:size(alldat,2)
        for win = 1:length(wincenters)
            keeptrack(pl,:) = [epoch,win];
            subwin = alldat(wincenters(win)-(srate/2)+1:wincenters(win)+srate/2,epoch).*hwin';
            %--------------------------------------------      
            [pwr1 frqs] = pwelch(subwin,srate,srate/2,srate*6,srate);
            pwr(pl,:) = pwr1(find(frqs > frqlim(1) & frqs < frqlim(2)))';
            pwr(pl,:) = 10*log10(pwr(pl,:));% convert to log for 'modulator' interpretation
            pl = pl+1;
        end;
        fprintf('.');
    end;   
    %%%% create specified frequency spacing:-------------------
    frqs = frqs(find(frqs > frqlim(1) & frqs < frqlim(2)));
    if strcmp(freqscale, 'quad') % make quadratic freq spacing
        if isempty(nfreqs)
        freqs = linspace(sqrt(frqs(1)), sqrt(frqs(end)), 1.5*(round(frqs(end)-frqs(1))));
        else
            freqs = linspace(sqrt(frqs(1)), sqrt(frqs(end)), nfreqs);
        end;
        freqs = freqs.^2;
        for f = 1:length(freqs)
            ff = frqs - freqs(f);
            [val idx] = min(abs(ff));
            newfreqs(f) = frqs(idx);
            newpwr(:,f) = pwr(:,idx);
        end;
        freqs = newfreqs; 
        pwr = newpwr;
        clear newpwr newfreqs
    elseif strcmp(freqscale,'log')
        %freqs = linspace(log(frqs(1)), log(frqs(end)), 1.5*round(frqs(end)-frqs(1)));
        if isempty(nfreqs)
            freqs = linspace(log(frqs(1)), log(frqs(end)), 1.5*(round(frqs(end)-frqs(1))));
        else
            freqs = linspace(log(frqs(1)), log(frqs(end)), nfreqs);
        end;
        freqs = exp(freqs);
        for f = 1:length(freqs)
            ff = frqs - freqs(f);
            [val idx] = min(abs(ff));
            newfreqs(f) = frqs(idx);
            newpwr(:,f) = pwr(:,idx);
        end;
        freqs = newfreqs; 
        pwr = newpwr;
        clear newpwr newfreqs        
    else
        fr = find(frqs > frqlim(1) & frqs < frqlim(2));
        freqs = frqs(fr); pwr = pwr(:,fr);        
    end;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear alldat
    
    % subtract mean spectrum:----------------------------
    pwr = pwr'; % make a freqs x epochs matrix
    fprintf('\nRemoving mean power spectrum ...\n');
    rowmeans =  squeeze(mean(pwr,2));  % across windows
    pwr = pwr - repmat(rowmeans,[1 size(pwr,2)]);

    fprintf('\nSaving spectral data for ICA...\n');
    
    floatwrite(pwr, [datpath,savedat,'.fdt']);
    
    if isempty(pcs) 
        pcs =  round(sqrt(size(pwr,2)/fac));
        if pcs > size(pwr,1)
            fac = fac*2;
            pcs =  round(sqrt(size(pwr,2)/fac));
            if pcs > size(pwr,1)
                fac = fac*2;
                pcs =  round(sqrt(size(pwr,2)/fac));
            end;
        end;
    end;       
    
    % generate an ICA script 
    ICA_SCRIPT = [datpath,'SpecDecomp.sc'];
    fid = fopen( ICA_SCRIPT, 'w');
    fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
    fprintf(fid, 'chans %d\n', size(pwr,1));
    fprintf(fid, 'frames %d\n',size(pwr,2));
    fprintf(fid, 'WeightsOutFile %s\n', [datpath,savedat,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [datpath,savedat,'.sph']);  
    fprintf(fid, 'sphering on\n');
    fprintf(fid, 'bias on\n');
    fprintf(fid, 'extended 1\n');
    if pcs ~= size(pwr,1)
        fprintf(fid, 'pca %d\n', pcs); 
    end;
    fprintf(fid, 'lrate 1.0e-4\n');
    fprintf(fid, 'blocksize 0\n');
    fprintf(fid, 'stop 1e-07\n');
    fprintf(fid, 'maxsteps 2000\n');
    fprintf(fid, 'posact on\n');
    fprintf(fid, 'annealstep 0.98\n');
    fprintf(fid, 'annealdeg 60\n');
    fprintf(fid, 'momentum 0\n');
    fprintf(fid, 'verbose on\n');
    fclose(fid);
    % run the ICA program
    run_ica_str = [ '/data/common/matlab/ica_linux2.4 < ', ICA_SCRIPT ];
    %run_ica_str = ['/home/jason/binica/ica_linux < ', ICA_SCRIPT]; % Jason's binica
    fprintf('\nRunning ICA in Linux, may take a long time...\n');
    [status, result] = system(run_ica_str);
    %%%%%
    numtrials = size(pwr,1);
    numframes = size(pwr,2);
    str = ['save ',datpath,savedat,'.mat numtrials numframes freqs pcs rowmeans wincenters keeptrack srate freqscale'];eval(str);    
    
    fprintf('\ndone.\n');
