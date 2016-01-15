% Calls in data for a subject, runs FFT analysis on short bits of data and then decomposes with ICA, organizes the matrix as: (freqs*ICs x time windows)
%
% SpecCoModTW(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath);
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
% nfreqs -- [number] number of output frequency bins (only for 'quad' or 'log'; for 'linear', use frqfac (internal to program))
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

function SpecCoModTW(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath,nfreqs);
    
    if ~exist('nfreqs') 
        nfreqs = frqlim(end)-frqlim(1);        
    elseif isempty(nfreqs)
        nfreqs = frqlim(end)-frqlim(1);
    end;
    frqfac = 6; % factor to adjust NUMBER of frequency bins (traditionally 6)
    rej = 'on'; % 'on' to reject noisy time points    
    
    if ~exist('auxpath')
        auxpath = [];
    end;
    
    for ds = 1:length(datset)
        EEG = pop_loadset( datset{ds},datpath);
        clear dat
        x = EEG.event(1).type;
        %%%%%%%%%%  Divide continuous data into 2 sec epochs:
        for evtm = EEG.srate*2:EEG.srate/overlap:size(EEG.data,2)-EEG.srate*2 
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
            EEG = pop_epoch( EEG,{'fake'} , [-1 1]);
        else
            EEG = pop_epoch( EEG,{1000} , [-1 1]);
        end;        
        EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
        EEG = eeg_checkset(EEG);
        %%%%%%%  Run automatice rejection on newly epoched data:
        %-------------------------------------------------------
        % probability and kurtosis are based on standard deviations, 
        % so iterative rejection is useful
        if strcmp(rej,'on') % if rejection set to 'on'
            fprintf('\nRunning auto-rejection protocol...\n');
            rmep = zeros(1,0);
            alleps = [1:size(EEG.icaact,3)];
            EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)] ,-1000,1000,EEG.xmin,EEG.xmax,0,0);
            numrej = length(find(EEG.reject.rejthresh));  % count number of epochs marked
            if numrej > 0
                rmep(1,end+1:end+length(find(EEG.reject.rejthresh))) = alleps(find(EEG.reject.rejthresh));
                alleps(find(EEG.reject.rejthresh)) = [];        
                EEG = pop_rejepoch( EEG,EEG.reject.rejthresh,0); % actually reject high prob epochs
                fprintf('\nRe-baselining after large amplitude artifact removed...\n');
                EEG = pop_rmbase( EEG,[-500 500]);
            end;
            %--------------------------------------------
            
            startlim = 5; % start probability threshold that will be up'ed if too many trials found
            maxrej = 40;
            EEG = pop_jointprob(EEG,0,complist ,startlim,startlim,0,0);% calculate component probabilities
            numrej = length(find(EEG.reject.icarejjp));  % count number of epochs marked
            if numrej < maxrej
                rmep(1,end+1:end+length(find(EEG.reject.icarejjp))) = alleps(EEG.reject.icarejjp);
                alleps(EEG.reject.icarejjp) = [];
                EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0); % actually reject high prob epochs
            else
                fprintf('Re-adjusting probability limits and running again...*************************\n');
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
                        fprintf('Re-adjusting probability limits and running again...*************************\n');
                    end;                    
                else
                    if startlim > 5 & maxiter < 8 % don't decrease and startover more than 8 times
                        fprintf('Decreasing probability limits for final pruning...###########################\n');
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
            %--------------------------------------------
            
            rmat(rmep) = 1;
            rmepochs{ds} = rmep; % keep track of which epochs you removed
            EEG.icaact = [];
            EEG = eeg_checkset(EEG);
        else
            rmepochs = [];
        end; % end to rejection step
        
        
        dstrials(1,ds) = size(EEG.data,3);         
        %%%%  Choose only component activations of interest
        for cp = 1:length(complist)
            multfact = sqrt(mean(EEG.icawinv(:,complist(cp)).^2));
            dat(cp,:,:) = EEG.icaact(complist(cp),:,:)*multfact; % added multfact back 12-6-05
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
    
    
    %%%%%%  Begin decompostion of single epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eplength = size(alldat,2);
    fprintf('\nCalculating pwr for each 1 second epoch...\n');
    %hwin = hanning(size(alldat,2))';
    for cp = 1:length(complist)
        for epoch = 1:size(alldat,3)
            [pwr(epoch,:,cp) frqs] = pwelch(alldat(cp,:,epoch),EEG.srate,EEG.srate/2,EEG.srate*frqfac,EEG.srate); 
            pwr(epoch,:,cp) = 10*log10(pwr(epoch,:,cp));% need to convert to log for 'modulator' interpretation
            
        end;   
    end;
    if strcmp(freqscale, 'quad') % make quadratic freq spacing
        freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
        freqs = freqs.^2;
        for f = 1:length(freqs)
            ff = frqs - freqs(f);
            [val idx] = min(abs(ff));
            newfreqs(f) = frqs(idx);
            newpwr(:,f,:) = pwr(:,idx,:);
        end;
        freqs = newfreqs; 
        pwr = newpwr;
        clear newpwr newfreqs
    elseif strcmp(freqscale, 'log') % make quadratic freq spacing
        freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
        freqs = exp(freqs);
        for f = 1:length(freqs)
            ff = frqs - freqs(f);
            [val idx] = min(abs(ff));
            newfreqs(f) = frqs(idx);
            newpwr(:,f,:) = pwr(:,idx,:);
        end;
        freqs = newfreqs; 
        pwr = newpwr;
        clear newpwr newfreqs
    else
        fr = find(frqs > frqlim(1) & frqs < frqlim(2));
        freqs = frqs(fr); pwr = pwr(:,fr,:);        
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % calculate mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%
    for cp = 1:size(pwr,3)
        meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
    end;
    
    clear alldat
    % subtract mean spectrum   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     fprintf('\nRemoving mean for each component/trial...\n');
% $$$     for cp = 1:size(pwr,3)
% $$$         for ep = 1:size(pwr,1)
% $$$             pwr(ep,:,cp) = pwr(ep,:,cp) - meanpwr(cp,:);
% $$$         end;        
% $$$     end;    
    pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));
    pwr = pwr'; % flip to be freqs*ICs x time windows
    
    fprintf('\nRemoving mean for each freq/IC ...\n');
    for tt = 1:size(pwr,1)
        pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); 
    end; 
    
    if ~isempty(auxpath)
        datpath = auxpath;
    end;
    
    fprintf('\nSaving spectral data for ICA...\n');%%%%%%%%%%%%%
    floatwrite(pwr, [datpath,savedat,'.fdt']);
    
    pcs =  round(sqrt(size(pwr,2)/pcfac));
    if pcs >= size(pwr,1)
        pcfac = pcfac*2;
        pcs =  round(sqrt(size(pwr,2)/pcfac));
    end;
    
    % generate an ICA script 
    ICA_SCRIPT = [datpath,'SpecDecompTW.sc'];
    fid = fopen( ICA_SCRIPT, 'w');
    fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
    fprintf(fid, 'chans %d\n', size(pwr,1));
    fprintf(fid, 'frames %d\n',size(pwr,2));
    fprintf(fid, 'WeightsOutFile %s\n', [datpath,savedat,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [datpath,savedat,'.sph']);  
    fprintf(fid, 'sphering on\n');
    fprintf(fid, 'bias on\n');
    fprintf(fid, 'extended 1\n');
    fprintf(fid, 'pca %d\n', pcs); 
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
    fprintf('\nRunning ICA in Linux, may take a long time...\n');
    [status, result] = system(run_ica_str);
    %%%%%
    numtrials = size(pwr,1);
    numframes = size(pwr,2);
    str = ['save ',datpath,savedat,'.mat meanpwr numtrials numframes freqs freqscale keeptrack rmepochs dstrials pcs complist overlap datset'];eval(str);    
    
    fprintf('\ndone.\n');
