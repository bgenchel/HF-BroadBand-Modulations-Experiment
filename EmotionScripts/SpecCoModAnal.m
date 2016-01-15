% Calls in data for a subject, runs spectral analysis on short bits of data and then decomposes with ICA
%
% [meanpwr,numtrials,numframes,freqs,keeptrack,rmepochs,dstrials, pcs, rowmeans] = SpecCoModAnal(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath,nfreqs,pwrdecomp,dec);
%
% datset: (cell array of strings) of dataset names of continuous data to use for spectral decomposition
% datpath: (string) fullpath to directory where all datasets are found and to which data will be stored
% complist: (vector) component indices to use for decomposition
% savedat: (string) name of float files in which to save input matrix and output wts and sph matrices
% frqlim -- [minfreq maxfreq]
% freqscale -- 'linear' or 'quad' for quadratic frequency spacing with more bins at lower freqs
% pcfac -- number of points per weight for final decomposition
% overlap -- overlap factor:  4= 75% overlap,2=50% overlap,1=no overlap,.75=25% gap,.5=50% gap
% auxpath -- if not [], uses this directory to store all generated data (.fdt,.wts,.sph)
% nfreqs -- [number] number of output frequency bins (only for 'quad' or 
%           'log'; for 'linear', use frqfac (internal to program))
% pwrdecomp -- [0 or [cycle1 cycle2 winsize]] if 0, will use pwelch to calculate standard
%              FFT. If number of cycles and winsize specified, then will perform wavelet
%              decomp. If cycle2 is larger than cycle1, then a linear scaling of the 
%              number of wavelets will be created between the lowest and highest freqs
%              requested. Default is 0. for wavelets, recommend 6 for low frqs,20-40 for high
%
% OUTPUT:
% meanpwr -- IC x freqs matrix giving removed mean spectrum for each IC
% numtrials -- number of trials or spectral windows used in decomposition
% numframes -- number of data points (columns) in decomposed matrix (freqs*#ICs)
% freqs -- vector of frequencies for power spectra
% keeptrack -- mostly only useful for several dataset decomps. gives start and end windows for each  datset
% dstrials: [vector] Gives total number of epochs for each dataset (sum(dstrials) = numtrials)
% rowmeans -- [matrix] numtrials X 1; mean removed from each row, AFTER removing component means, 
%             but just before ICA
% files automatically saved are:
%    'savedat'.fdt -- raw data submitted to ICA
%    'svaedat'.wts -- unmixing matrix returned by ICA
%    'savedat'.sph -- sphering matrix returned by ICA (W = wts*sph)

function [meanpwr,numtrials,numframes,freqs,keeptrack,rmepochs,dstrials, pcs, rowmeans] = SpecCoModAnal(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath,nfreqs,pwrdecomp,dec);
    
    thresh = 1000000; % threshold for uV auto-rejection (normally 1000)
    
   if ~exist('nfreqs') 
        nfreqs = frqlim(end)-frqlim(1);        
    elseif isempty(nfreqs)
        nfreqs = frqlim(end)-frqlim(1);
    end;
   
    if ~exist('dec')
        dec = [];
    end;
    if ~exist('pwrdecomp')
        pwrdecomp = 0;
    end;
    if pwrdecomp == 0 % FFT
        method = 'FFT';
    else
        method = 'wavelet';
    end;
    EEG = pop_loadset( datset{1},datpath);srate = EEG.srate;% for srate
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
       if freqs(1) == 0
           freqs([1:2]) = []; % eliminate 0 Hz
           cycles1 = freqs([1:2]);
       end;
       wavelets = computepsifamilyQodd(freqs,1/srate,cycles1,cycles2,winsize); %
    else % FFT
       winsize = srate;
       hwin = hanning(winsize)';
       %hwin = hanning(length(subwin))';
    end;
     
    rej = 'on'; % 'on' to reject noisy time points    
    
    if ~exist('auxpath')
        auxpath = [];
    end;
    
    for ds = 1:length(datset)
        EEG = pop_loadset( datset{ds},datpath);srate = EEG.srate;
        if size(EEG.data,3) > 1
            epoched = 'true';
        else
            epoched = 'false';
        end;
        if pwrdecomp == 0 % FFT
            epsize = .5; % one second
        else
            epsize = ((pwrdecomp(3)-1)/EEG.srate)/2 + .05;% need to make just > 1 sec
        end;
        clear dat
        if strcmp(epoched,'true') % if there are events to speak of
            % take windows at specified time points around all events
            % keep track of events and peri-event timing.
        else       
            if ~isempty(EEG.event)
                x = EEG.event(1).type;
            else
                EEG.event(1).type = 'fake';
                EEG.event(1).latency = 1;
                x = EEG.event(1).type;  
            end;
            %%%%%%%%%%  Divide continuous data into ? sec epochs:
            for evtm = epsize/2*EEG.srate:round(EEG.srate/overlap):size(EEG.data,2)-epsize/2*EEG.srate
                % create events to make overlapping 1 sec epochs
                EEG.event(end+1) =  EEG.event(1);% appends events to the end
                EEG.event(end).latency = evtm;
                if ischar(x)
                    EEG.event(end).type = 'fake';% for string event codes  
                else
                    EEG.event(end).type = 1000;% for string event codes 
                end;
            end;
            if ischar(x)
                EEG = pop_epoch( EEG,{'fake'} , [-epsize epsize]);% was-1 1 recently
            else
                EEG = pop_epoch( EEG,{1000} , [-epsize epsize]);
            end;        
            %EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
            EEG = eeg_checkset(EEG);
        end;
        %%%%%%%  Run automatic rejection on newly epoched data:
        %-------------------------------------------------------
        % probability and kurtosis are based on standard deviations, 
        % so iterative rejection is useful
        if strcmp(rej,'on') % if rejection set to 'on'
            fprintf('\nRunning auto-rejection protocol...\n');
            rmep = zeros(1,0);
            alleps = [1:size(EEG.data,3)];
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)] ,-thresh,thresh,EEG.xmin,EEG.xmax,0,0);
            numrej = length(find(EEG.reject.rejthresh));  % count number of epochs marked
            if numrej > 0
                if numrej > size(EEG.data,3)
                    fprintf('\nThreshold rejection has found more than half of your epochs to be artifact. Consider increasing the uV threshold and running again\n');
                    return;
                else
                    rmep(1,end+1:end+length(find(EEG.reject.rejthresh))) = alleps(find(EEG.reject.rejthresh));
                    alleps(find(EEG.reject.rejthresh)) = [];        
                    EEG = pop_rejepoch( EEG,EEG.reject.rejthresh,0); % actually reject high prob epochs
                end;
            end;
            %--------------------------------------------
            startlim = 5; % start probability threshold that will be up'ed if too many trials found
            maxrej = 40;
            if complist == 0 % channels
               EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,startlim,startlim,0,0);% calculate component probabilities
               numrej = length(find(EEG.reject.rejjp));  % count number of epochs marked
               if numrej < maxrej
                  rmep(1,end+1:end+length(find(EEG.reject.rejjp))) = alleps(EEG.reject.rejjp);
                  alleps(EEG.reject.rejjp) = [];
                  EEG = pop_rejepoch( EEG,EEG.reject.rejjp,0); % actually reject high prob epochs
               else
                  fprintf('Re-adjusting probability limits to %s and running again...*************************\n',int2str(startlim));
                  startlim = startlim + .5;
               end;
               repeat = 1; maxiter = 0;
               while repeat == 1 % keep running probability until there are no epochs above threshold
                  if numrej > 0
                     EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,startlim,startlim,0,0);
                     numrej = length(find(EEG.reject.rejjp));
                     if numrej < maxrej
                        rmep(1,end+1:end+length(find(EEG.reject.rejjp))) = alleps(EEG.reject.rejjp);
                        alleps(EEG.reject.rejjp) = [];
                        EEG = pop_rejepoch( EEG,EEG.reject.rejjp,0);
                        fprintf('\nep mat: %s,\tdat mat: %s\n',int2str(length(alleps)),int2str(size(EEG.icaact,3)));
                     else
                        startlim = startlim + .5; EEG.reject.rejjp = [];EEG.reject.rejjpE = [];
                        fprintf('Re-adjusting probability limits to %s and running again...*************************\n',int2str(startlim));
                     end;
                  else
                     if startlim > 5 & maxiter < 8 % don't decrease and startover more than 8 times
                        fprintf('Decreasing probability limits so %s ###########################\n',int2str(startlim));
                        startlim = startlim - .5; numrej = 1; maxiter = maxiter+1; % repeat process back to 5 stds
                     else
                        if maxiter > 8
                           maxrej = 100; % go through last round with a high threshold
                        else
                           repeat = 0;
                        end;
                     end;
                  end;
               end;
               % run kurtosis check
               EEG = pop_rejkurt(EEG,0,complist ,6,6,0,0);
               numrej = length(find(EEG.reject.rejkurt));  % count number of epochs marked
               if numrej > 0
                  rmep(1,end+1:end+length(find(EEG.reject.rejkurt))) = alleps(EEG.reject.rejkurt);
                  alleps(EEG.reject.rejjp) = [];
                  EEG = pop_rejepoch( EEG,EEG.reject.rejkurt,0);
               end;
            else
               EEG = pop_jointprob(EEG,0,complist ,startlim,startlim,0,0);% calculate component probabilities
               numrej = length(find(EEG.reject.icarejjp));  % count number of epochs marked
               if numrej < maxrej
                  rmep(1,end+1:end+length(find(EEG.reject.icarejjp))) = alleps(EEG.reject.icarejjp);
                  alleps(EEG.reject.icarejjp) = [];
                  EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0); % actually reject high prob epochs
               else
                  fprintf('Re-adjusting probability limits to %s and running again...*************************\n',int2str(startlim));
                  startlim = startlim + .5;
               end;
               repeat = 1; maxiter = 0;
               while repeat == 1 % keep running probability until there are no epochs above threshold
                  if numrej > 0
                     EEG = pop_jointprob(EEG,0,complist ,startlim,startlim,0,0);
                     numrej = length(find(EEG.reject.icarejjp));
                     if numrej < maxrej
                        rmep(1,end+1:end+length(find(EEG.reject.icarejjp))) = alleps(EEG.reject.icarejjp);
                        alleps(EEG.reject.icarejjp) = [];
                        EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0);
                        fprintf('\nep mat: %s,\tdat mat: %s\n',int2str(length(alleps)),int2str(size(EEG.icaact,3)));
                     else
                        startlim = startlim + .5; EEG.reject.icarejjp = [];EEG.reject.rejjpE = [];
                        fprintf('Re-adjusting probability limits to %s and running again...*************************\n',int2str(startlim));
                     end;
                  else
                     if startlim > 5 & maxiter < 8 % don't decrease and startover more than 8 times
                        fprintf('Decreasing probability limits to %s ###########################\n',int2str(startlim));
                        startlim = startlim - .5; numrej = 1; maxiter = maxiter+1; % repeat process back to 5 stds
                     else
                        if maxiter > 8
                           maxrej = 100; % go through last round with a high threshold
                        else
                           repeat = 0;
                        end;
                     end;
                  end;
               end;
               % run kurtosis check
               EEG = pop_rejkurt(EEG,0,complist ,6,6,0,0);
               numrej = length(find(EEG.reject.icarejkurt));  % count number of epochs marked
               if numrej > 0
                  rmep(1,end+1:end+length(find(EEG.reject.icarejkurt))) = alleps(EEG.reject.icarejkurt);
                  alleps(EEG.reject.icarejjp) = [];
                  EEG = pop_rejepoch( EEG,EEG.reject.icarejkurt,0);
               end;
            end;
            %--------------------------------------------
            fprintf('\nep mat: %s,\tdat mat: %s\n',int2str(length(alleps)),int2str(size(EEG.data,3)));
               
            rmat(rmep) = 1;
            rmepochs{ds} = rmep; % keep track of which epochs you removed
            EEG.icaact = [];
            EEG = eeg_checkset(EEG);
        else
            rmepochs = [];
        end; % end to rejection step
        
        
        dstrials(1,ds) = size(EEG.data,3);     
        if complist == 0 % channels
           for cp = 1:size(EEG.data,1)
              dat(cp,:,:) = EEG.data(cp,:,:); 
           end;           
        else
           %%%%  Choose only component activations of interest
           for cp = 1:length(complist)
              multfact = sqrt(mean(EEG.icawinv(:,complist(cp)).^2));
              dat(cp,:,:) = EEG.icaact(complist(cp),:,:)*multfact; % added multfact back 12-6-05
           end;
        end;
        if ds == 1
            alldat = dat;            
        else
            alldat(:,:,end+1:end+size(dat,3)) = dat;
        end;
        if ds == 1
            keeptrack(ds,:) = [1 size(dat,3)];
        else
            keeptrack(ds,:) = [keeptrack(ds-1,2)+1 keeptrack(ds-1,2)+size(dat,3)];
        end;
    end;  clear dat
    
    % added for testing:
    if ~isempty(dec)
        randdels = randperm(size(alldat,3));
        randdels = randdels(1:round(size(alldat,3)/dec));
        alldat = alldat(:,:,randdels);
    end;
    alldat = double(alldat);
    %floatwrite(alldat, [datpath,savedat,'RAW.fdt']);% save raw data for subsequent decomps

    %%%%%%%%%%%%%%%%  Begin decompostion of single epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eplength = size(alldat,2);
    if length(pwrdecomp) < 3 & pwrdecomp == 0 % FFT using pwelch          
        fprintf('\nCalculating pwr for each %s second epoch...\n',int2str(size(alldat,2)/srate));
        %hwin = hanning(size(alldat,2))'; fprintf('\nUsing hanning windowed data...\n');
        % to get number of freqs:
        [p f] = pwelch(alldat(1,:,1),EEG.srate,EEG.srate*3/4,freqs,EEG.srate);
        %[p f] = pwelch(alldat(1,:,1).*hwin,EEG.srate,EEG.srate*3/4,EEG.srate*10,EEG.srate);  % old method
        pwr = zeros(size(alldat,3),length(f),size(alldat,1));
        for cp = 1:size(alldat,1)
            for epoch = 1:size(alldat,3)
               % remove 'hwin' multiplication (use all power in epoch win
                [pwr(epoch,:,cp) freqs] = pwelch(alldat(cp,:,epoch),EEG.srate,EEG.srate*3/4,freqs,EEG.srate); 
                %[pwr(epoch,:,cp) frqs] = pwelch(double(alldat(cp,:,epoch).*hwin),EEG.srate,EEG.srate*3/4,EEG.srate*10,EEG.srate);  % old method
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
            elseif rempnts2 == rempnts1        
              alldat = alldat(:,rempnts1+1:end,:);
              alldat = alldat(:,1:end-rempnts2,:);
            end;
        end;
        pwr = zeros(size(alldat,3),length(freqs),size(alldat,1));
        fprintf('\nCalculating pwr for each %s second epoch...\n',num2str(size(alldat,2)/srate));
        for cp = 1:size(alldat,1)
           % remove 'hwin' multiplication (use all power in epoch win)
           pwr1 = wavelets * (squeeze(alldat(cp,:,:))); 
            pwr1 = pwr1.*conj(pwr1); % pwr1 is now non-complex
            pwr1 = 10*log10(pwr1);% convert to log for 'modulator' interpretation
            pwr(:,:,cp) = pwr1'; clear pwr1
            fprintf('\nIC %s done...',int2str(cp));
        end;
    end; % to decision of decomp method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     if pwrdecomp == 0 % only if FFT resample to different freq spacing:   
%         if strcmp(freqscale, 'quad') % make quadratic freq spacing
%             freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
%             freqs = freqs.^2;
%         elseif strcmp(freqscale, 'log') % make quadratic freq spacing
%             freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
%             freqs = exp(freqs);
%         else
%             freqs = linspace(frqlim(1), frqlim(2),nfreqs);
%         end;
%         clear idx
%         for f = 1:length(freqs)
%             ff = frqs - freqs(f);
%             [val idx(f)] = min(abs(ff));
%         end;
%         if length(unique(idx)) < length(idx)
%             fprintf('\nWarning, some frequency bins have been repeated...\n');
%         end;           
%         freqs = frqs(idx);
%         pwr = pwr(:,idx,:); % (epochs x freqs x ICs)
%     end;
    
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
    
    if ~isempty(auxpath)
        datpath = auxpath;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\nSaving spectral data (%s by %s) as %sDAT.fdt\n',int2str(size(pwr,1)),int2str(size(pwr,2)),savedat);%%%%%%%%%%%%%%%%%%%%%%
    floatwrite(pwr, [datpath,savedat,'DAT.fdt']);
    
    if pcfac < 10 % otherwise use as number of pcs
        pcs =  round(sqrt(size(pwr,2)/pcfac));
        if pcs >= size(pwr,1)
            pcfac = pcfac*2;
            pcs =  round(sqrt(size(pwr,2)/pcfac));
        end;
    else
        pcs =  pcfac;
    end;
    % eigvec is the winv matrix
    %[eigvec,sv,pc] = pcsquash(pwr,pcs); 
    [U,S,eigvec] = svds(pwr',pcs);% if you scale 'acts', you can't scale the 'eigvec'
    pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp

    floatwrite(eigvec, [datpath,savedat,'EIGVEC.fdt']);
    eigfile = [savedat,'EIGVEC.fdt'];% save in .mat
    fprintf('\nSaving PCA''d data for ICA as %s.fdt.\n',savedat);%%%%%%%%%%%%%%%%%%%%%%
    floatwrite(pc, [datpath,savedat,'.fdt']);
    numrows = pcs;
    numtrials = size(pwr,1);
    numframes = size(pwr,2);

    
    % generate an ICA script 
%     ICA_SCRIPT = [datpath,'SpecDecomp.sc'];
%     fid = fopen( ICA_SCRIPT, 'w');
%     fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
%     fprintf(fid, 'chans %d\n', pcs);
%     fprintf(fid, 'frames %d\n',numframes);
%     fprintf(fid, 'WeightsOutFile %s\n', [datpath,savedat,'.wts']);
%     fprintf(fid, 'SphereFile %s\n', [datpath,savedat,'.sph']);  
%     fprintf(fid, 'sphering on\n');
%     fprintf(fid, 'bias on\n');
%     fprintf(fid, 'extended 1\n');
%     fprintf(fid, 'lrate 1.0e-4\n');
%     fprintf(fid, 'blocksize 0\n');
%     fprintf(fid, 'stop 1e-07\n');
%     fprintf(fid, 'maxsteps 2000\n');
%     fprintf(fid, 'posact on\n');
%     fprintf(fid, 'annealstep 0.98\n');
%     fprintf(fid, 'annealdeg 60\n');
%     fprintf(fid, 'momentum 0\n');
%     fprintf(fid, 'verbose on\n');
%     fclose(fid);
%     % run the ICA program
%     run_ica_str = [ '/data/common/matlab/ica_linux2.4 < ', ICA_SCRIPT ];
%     fprintf('\nRunning ICA in Linux, may take a long time...\n');
%     [status, result] = system(run_ica_str);
    %%%%%
    [wts sph] = runica(pc,'stop',1e-7,'extended',1,'maxsteps',2000);
    floatwrite(wts, [datpath,savedat,'.wts']);
    floatwrite(sph, [datpath,savedat,'.sph']);
    
    str = ['save ',datpath,savedat,'.mat method meanpwr numrows numframes freqs freqscale keeptrack rmepochs dstrials  pcs rowmeans complist overlap datset eigfile'];eval(str);    
    
    fprintf('\ndone.\n');
