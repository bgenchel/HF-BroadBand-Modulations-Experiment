% takes a 2 or 3D matrix of EEG data. 
% Calculates power in each 1-sec window (overlap = ?? %)
% OR calculates power in each epoch (3rd dim)  
% performs spectral decomposition
%
%
% INPUTS:
% alldat -- [3D matrix] dimensions: channels/components X time X epochs
% datpath -- [string] data file path where 'savedat.fdt' will be saved
% savedat -- [string] name to save data float file as
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% pcs -- [integer] number of PCA dimensions to retain in final spectral decomposition {15}
% pcfac -- [number] number of points per weight to have in final ICA (will PCA down accordingly)
%          If 'pcfac' is > 10, then this variable specifies the number of PCA dimensions
%          to keep (does not scale by the number of data points).
% srate -- [integer] sampling rate of the input data matrix 'alldat'

function [numrows numframes complist freqs pcs pcfac freqscale rowmeans meanpwr keeptrack eigfile] = SpecCoModDecomp(alldat,datpath,savedat,complist,keeptrack,frqlim,freqscale,pcfac,srate,nfreqs,pwrdecomp);


    overlap = 6; % lots of overlap
    windowed = 'on';
    reject = 'off';
    
    if length(size(alldat)) > 2 % if 3d data:
        if strcmp(reject,'on')
            fprintf('\nRejecting noisy epochs by probability estimate...\n');
            [jp rej] = jointprob( alldat, 5, [], 1);
            rej = sum(rej,1);
            alldat(:,:,find(rej)) = [];% remove noisy epochs
        end;
    else % if 2D
        neps = length(srate:srate/2:size(alldat,2)-srate);
        newdat = zeros(size(alldat,1),srate+1,neps);
        e = 1; % index for epochs
        for ep = srate:round(srate/overlap):size(alldat,2)-srate % ?% overlap (see above)
            newdat(:,:,e) = alldat(:,ep-(srate/2):ep+(srate/2));e = e+1;
        end;
        if strcmp(reject,'on')
            fprintf('\nRejecting noisy epochs by probability estimate...\n');
            [jp rej] = jointprob( newdat, 5, [], 1);
            rej = sum(rej,1);
            newdat(:,:,find(rej)) = [];% remove noisy epochs
        end;
        alldat = newdat;
    end;
    if strcmp(freqscale, 'quad') % make quadratic freq spacing
       freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
       freqs = freqs.^2;
    elseif strcmp(freqscale, 'log') % make quadratic freq spacing
       if frqlim(1) < 1
          fprintf('\nlowest frequency must be 1 or greater for log-spaced frequencies'); return;
       end;
       freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
       freqs = exp(freqs);
    else % if linear-spacing
       freqs = linspace(frqlim(1), frqlim(2),nfreqs);
    end;
    if length(pwrdecomp) > 2 %if wavelet, specify freq spacing in advance
        cycles1 = pwrdecomp(1);
        cycles2 = pwrdecomp(2);
        winsize = pwrdecomp(3);
        hwin = hanning(winsize)';
        wavelets = computepsifamilyQodd(freqs,1/srate,cycles1,cycles2,winsize); % 
    else % FFT
        winsize = srate;
        hwin = hanning(winsize)';
    end;             

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\nCalculating pwr for each epoch...\n');
    eplength = size(alldat,2);
    if pwrdecomp == 0 % FFT using pwelch          
        fprintf('\nCalculating pwr for each %s second epoch...\n',int2str(size(alldat,2)/srate));
        hwin = hanning(size(alldat,2))';
        fprintf('\nUsing hanning windowed data...\n');
        % to get number of freqs:
        [p f] = pwelch(alldat(1,:,1),EEG.srate,EEG.srate*3/4,freqs,EEG.srate);
        pwr = zeros(size(alldat,3),length(f),size(alldat,1));
        for cp = 1:length(complist)
            for epoch = 1:size(alldat,3)
                 [pwr(epoch,:,cp) freqs] = pwelch(alldat(cp,:,epoch),EEG.srate,EEG.srate*3/4,freqs,EEG.srate); 
                 pwr(epoch,:,cp) = 10*log10(pwr(epoch,:,cp));% need to convert to log for 'modulator' interpretation
            end;   
            fprintf('\nIC %s done...',int2str(cp));
        end;
    elseif length(pwrdecomp) > 2 % use wavelet decomp
        if size(alldat,2) > winsize
            ndiff = size(alldat,2) - winsize;
            rempnts1 = floor(ndiff/2);
            rempnts2 = ceil(ndiff/2);
            if rempnts2 == 1
                alldat = alldat(:,1:end-rempnts2,:,:);
            elseif rempnts2 > rempnts1
                alldat = alldat(:,rempnts2:end,:);
                alldat = alldat(:,1:end-rempnts2,:);
            else                    
                alldat = alldat(:,rempnts2:end,:);
                alldat = alldat(:,1:end-(rempnts2-1),:);
            end;
        elseif size(alldat,2) < winsize
            fprintf('\nError: Epoch length must be equal or greater to the winsize (%s)\n\n',int2str(winsize));return;
        end;
        pwr = zeros(size(alldat,3),length(freqs),size(alldat,1));
        fprintf('\nCalculating pwr for each %s second epoch...\n',num2str(size(alldat,2)/srate));
        for cp = 1:size(alldat,1)
            pwr1 = wavelets * (squeeze(alldat(cp,:,:)).*repmat(hwin',[1 size(alldat,3)])); 
            pwr1 = pwr1.*conj(pwr1); % pwr1 is now non-complex
            pwr1 = 10*log10(pwr1);% convert to log for 'modulator' interpretation
            pwr(:,:,cp) = pwr1'; clear pwr1
            fprintf('\nIC %s done...',int2str(cp));
        end;
    end; % to decision of decomp method
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % calculate mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cp = 1:size(pwr,3)
        meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
    end;
    
    clear alldat
    % subtract mean spectrum   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nRemoving mean spectral power for each component...\n');
    for cp = 1:size(pwr,3)
        pwr(:,:,cp) = pwr(:,:,cp) - repmat(meanpwr(cp,:),[size(pwr,1) 1]);
    end;    
    pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));
    fprintf('\nRemoving mean for each window (all comps)...\n');
    for tt = 1:size(pwr,1)
        rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
        pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (all comps)
    end; 
    
 
    fprintf('\nSaving original data as %sDAT.fdt.\n',savedat);
    floatwrite(pwr, [datpath,savedat,'DAT.fdt']);
    
    if pcfac > 10 % pcfac specifies number of pcs
        pcs = pcfac;
    else % pcs determined by number of columns and pcfac        
        pcs =  round(sqrt(size(pwr,2)/pcfac)); % determine number of PCA components
        if pcs >= size(pwr,1)
            pcfac = pcfac*2;
            pcs =  round(sqrt(size(pwr,2)/pcfac));
        end;
    end;
    fprintf('\n PCAing to %s dimensions\n',int2str(pcs));%%%%%%%%%%%%%%%%%%%%%%
    
    % eigvec is the winv matrix
    %[eigvec,sv,pc] = pcsquash(pwr,pcs); 
    [U,S,eigvec] = svds(pwr',pcs);% if you scale 'acts', you can't scale the 'eigvec'
    pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp

    eigfile = [savedat,'EIGVEC.fdt'];% save in .mat
    floatwrite(eigvec, [datpath,eigfile]);
    fprintf('\nSaving PCA''d data for ICA as %s.fdt.\n',savedat);%%%%%%%%%%%%%%%%%%%%%%
    floatwrite(pc, [datpath,savedat,'.fdt']);
    numrows = pcs;
    numframes = size(pwr,2);
    clear pc pceigvec  
    
    % generate an ICA script 
    ICA_SCRIPT = [datpath,savedat,'.sc'];
    fid = fopen( ICA_SCRIPT, 'w');
    fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
    fprintf(fid, 'chans %d\n', numrows);
    fprintf(fid, 'frames %d\n',numframes);
    fprintf(fid, 'WeightsOutFile %s\n', [datpath,savedat,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [datpath,savedat,'.sph']);  
    fprintf(fid, 'sphering on\n');
    fprintf(fid, 'bias on\n');
    fprintf(fid, 'extended 1\n');
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
    run_ica_str = [ '/programs/MATLAB/binica/ica_linux < ', ICA_SCRIPT ];
    fprintf('\n\nRun the following line in Linux to complete the ICA decomposition:\n\n');
    fprintf('%s',run_ica_str);
     fprintf('\n\n');
     %[status, result] = system(run_ica_str);
    %%%%%
    numtrials = size(pwr,1);
    numframes = size(pwr,2);
    
    str = ['save ',datpath,savedat,'.mat numrows numframes complist freqs pcs pcfac freqscale rowmeans meanpwr keeptrack eigfile'];eval(str);    
    
    fid = fopen(['/programs/OASIS/VideoCoModCalls.txt'], 'a');
    fprintf(fid, '%s\n', run_ica_str);
    fclose(fid);

    %fprintf('\ndone.\n');
