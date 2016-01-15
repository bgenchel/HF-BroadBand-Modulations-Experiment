% takes a float file used for CoMod decomp and sub-samples, then
% re-runs ICA with little or no PCA. Also can be used to rerun 
% comod decomp with other parameters but no subsample. Takes either
% raw data and runs FFT or starts from spectral data.
%
% CoModDecomp(datpath,savedat, newsave,ntrials,pcfac,auxpath,datatype,srate);
%
% datpath: (string) fullpath to directory where all datasets are found and to which data will be stored
% savedat: (string) name of float file containing data matrix ('savedat'RAW.fdt)
% newsave: (string) name of float file to save new data matrix as.
% ntrials -- [number] of rows to retain (at random) from the data matrix
%            [] will skip sub-sample step.
% pcfac -- number of points per weight for final decomposition
% auxpath -- if not [], uses this directory to store all generated data (.fdt,.wts,.sph)
% datatype -- ['raw' or 'spec'] 'raw' indicates raw ICA activations to be
%             submitted to FFT decomp, then ICA. 'spec' skips the FFT step.
% specparams -- [srate freqscale minfrq maxfrq numfreqs] data sampling rate of 'raw' data. 
%               leave [] if submitting 'spec'
%               also indicate 'log', 'quad' or 'linear' frequency spacing, min and max
%               freqs and number of frequency bins to create
% pwrdecomp -- [0 or [cycle1 cycle2 winsize]] if 0, will use pwelch to calculate standard
%              FFT. If number of cycles and winsize specified, then will perform wavelet
%              decomp. If cycle2 is larger than cycle1, then a linear scaling of the 
%              number of wavelets will be created between the lowest and highest freqs
%              requested. Default is 0. for wavelets, recommend 6 for low frqs,20-40 for high

function CoModDecomp(datpath,savedat,newsave,ntrials,pcfac,auxpath,datatype,specparams,pwrdecomp);
    
    
    skipresample = 0; % this is for timefreq decomp that makes quad/log within
    srate = specparams{1};
   freqscale = specparams{2};
   frqlim = [specparams{3} specparams{4}];
   nfreqs = specparams{5};
   
    if ~exist('pwrdecomp')
        pwrdecomp = 0;
    end;
     if length(pwrdecomp) > 2 %if wavelet, specify freq spacing in advance
        if strcmp(freqscale, 'quad') % make quadratic freq spacing
            freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
            freqs = freqs.^2;
        elseif strcmp(freqscale, 'log') % make quadratic freq spacing
            freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
            freqs = exp(freqs);
        else % if linear-spacing
            freqs = linspace(frqlim(1), frqlim(2),nfreqs);
        end;
        cycles1 = pwrdecomp(1);
        cycles2 = pwrdecomp(2);
        winsize = pwrdecomp(3);
        hwin = hanning(winsize)';
        wavelets = computepsifamilyQodd(freqs,1/srate,cycles1,cycles2,winsize); % 
    else % FFT
        winsize = srate;
        hwin = hanning(winsize)';
        %hwin = hanning(length(subwin))';
    end;             
   
    if strcmp(datatype,'raw')
        str = ['load ',datpath,savedat,'.mat'];eval(str);    % for the resave
        clear meanpwr freqs % may want to change these
        s = load([datpath,savedat,'.mat']);
        alldat = floatread([datpath,savedat,'RAW.fdt'],[length(s.complist) srate*2 length(s.rowmeans)],[],0);
        freqscale = specparams{2};    
        
        if ~isempty(ntrials)
            fprintf('\nDecimating number of trials as requested to %s trials.\n',int2str(ntrials));
            randidx = randperm(size(alldat,3));
            randidx = randidx(1:ntrials);
            alldat = alldat(:,:,sort(randidx));             
            for ds = 1:size(s.keeptrack,1)                 
                keepvec(1,s.keeptrack(ds,1):s.keeptrack(ds,2)) = ds;
            end;
            keepvec = keepvec(1,sort(randidx));
            for ds = 1:size(s.keeptrack,1)  
                idxes = find(keepvec == ds);
                newtrack(ds,:) = [idxes(1) idxes(end)];
            end;
            s.keeptrack = newtrack;
        end;
        
    %%%%%%
        %%%%%%%%%%%%%%%%  Begin decompostion of single epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eplength = size(alldat,2);
        fprintf('\nUsing hanning windowed data...\n');
                
        
        if pwrdecomp == 0 % FFT using pwelch
            hwin = hanning(size(alldat,2))';
           % to get number of freqs:
           [p ff] = pwelch(alldat(1,:,1).*hwin,srate,srate*3/4,srate*10,srate);  
            pwr = zeros(size(alldat,3),length(ff),size(alldat,1)); clear p ff
            fprintf('\nCalculating pwr for each %s second epoch...\n',num2str(size(alldat,2)/srate));
            for cp = 1:length(s.complist)
                for epoch = 1:size(alldat,3)
                    [pwr(epoch,:,cp) frqs] = pwelch(alldat(cp,:,epoch).*hwin,srate,srate*3/4,srate*10,srate); 
                    pwr(epoch,:,cp) = 10*log10(pwr(epoch,:,cp));% convert to log  
                end;   
                fprintf('\nIC %s done...',int2str(cp));
            end;
            frqs = frqs';
        elseif length(pwrdecomp) > 2 % use wavelet decomp
            if pwrdecomp(2) == specparams{4}           
                % timefreq from eeglab:
                skipresample = 1;
                frqs = [specparams{3}:.1:specparams{4}];
                winsize = pwrdecomp(3);
                if size(alldat,2) > winsize
                    ndiff = size(alldat,2) - winsize;
                    rempnts1 = floor(ndiff/2);
                    rempnts2 = ceil(ndiff/2);
                    if rempnts2 == 1
                        alldat = alldat(:,1:end-rempnts2,:,:);
                    elseif rempnts2 > rempnts1
                        alldat(:,1:rempnts1,:) = [];
                        alldat(:,end - rempnts1:end,:) = []; % yes, rempnts1
                    elseif rempnts2 == rempnts1 & rempnts1 ~= 0     
                        alldat(:,1:rempnts1,:) = [];
                        alldat(:,end - (rempnts2-1):end,:) = [];
                    end;
                end;
                frqlim(1) = specparams{3}; frqlim(2) = specparams{4};
                freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), specparams{5});
                freqs = freqs.^2;
                pwr = zeros(size(alldat,3),specparams{5},size(alldat,1));
                for cp = 1:size(alldat,1)
                    [tf, freqs, times] = mytimefreq(squeeze(alldat(cp,:,:)),srate,'cycles',[pwrdecomp(1) .5],'wletmethod','dftfilt2','ntimesout',1,'freqs',freqs);% freq scale and nfreqs determined by freqs
                    pwr(:,:,cp) = squeeze(tf)'; clear tf
                end;
            else
                % wavelet analysis made by Rey Ramirez:
                if size(alldat,2) > winsize
                    ndiff = size(alldat,2) - winsize;
                    rempnts1 = floor(ndiff/2);
                    rempnts2 = ceil(ndiff/2);
                    if rempnts2 == 1
                        alldat = alldat(:,1:end-rempnts2,:,:);
                    elseif rempnts2 > rempnts1
                        alldat(:,1:rempnts1,:) = [];
                        alldat(:,end - rempnts1:end,:) = []; % yes, rempnts1
                    elseif rempnts2 == rempnts1 & rempnts1 ~= 0     
                        alldat(:,1:rempnts1,:) = [];
                        alldat(:,end - (rempnts2-1):end,:) = [];
                    end;               
                end;
                pwr = zeros(size(alldat,3),length(frqs),size(alldat,1));
                fprintf('\nCalculating pwr for each %s second epoch...\n',num2str(size(alldat,2)/srate));       
                for cp = 1:size(alldat,1)
                    pwr1 = wavelets * (squeeze(alldat(cp,:,:)).*repmat(hwin',[1 size(alldat,3)])); 
                    pwr1 = pwr1.*conj(pwr1); %  pwr1 is now non-complex
                    pwr1 = 10*log10(pwr1);% convert to log for 'modulator' interpretation
                    pwr(:,:,cp) = pwr1'; clear pwr1
                    fprintf('\nIC %s done...',int2str(cp));
                end;
            end;
        end;
        clear wavelets s 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if skipresample == 0 % always except when using timefreq.m
            if pwrdecomp == 0 % only if FFT resample to different freq spacing:   
                nfreqs = specparams{5};
                frqlim(1) = specparams{3}; frqlim(2) = specparams{4};
                if strcmp(freqscale, 'quad') % make quadratic freq spacing
                    freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
                    freqs = freqs.^2;
                elseif strcmp(freqscale, 'log') % make quadratic freq spacing
                    freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
                    freqs = exp(freqs);
                else
                    freqs = linspace(frqlim(1), frqlim(2),nfreqs);
                end;
                clear idx
                for f = 1:length(freqs)
                    ff = frqs - freqs(f);
                    [val idx(f)] = min(abs(ff));
                end;
                if length(unique(idx)) < length(idx)
                    fprintf('\nWarning, some frequency bins have been repeated...\n');
                end;           
                freqs = frqs(idx);
                pwr = pwr(:,idx,:); % (epochs x freqs x ICs)
            end;
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % calculate mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for cp = 1:size(pwr,3)
            meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
        end;
        
        clear alldat
        % subtract mean spectrum   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nRemoving mean power spectrum for each component...\n');
        for cp = 1:size(pwr,3)
            pwr(:,:,cp) = pwr(:,:,cp) - repmat(meanpwr(cp,:),[size(pwr,1) 1]);
            %for ep = 1:size(pwr,1)
            %    pwr(ep,:,cp) = pwr(ep,:,cp) - meanpwr(cp,:);
            %end;        
        end;    
        pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));
        fprintf('\nRemoving mean for each window (all comps)...\n');
        for tt = 1:size(pwr,1)
            rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
            pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (all comps)
        end; 
        data = pwr; clear pwr
    else % load power data instead
        
        str = ['load ',datpath,savedat,'.mat'];eval(str);    
        s = load([datpath,savedat,'.mat']);
        data = floatread([datpath,savedat,'.fdt'],[s.numtrials s.numframes],[],0);
        floatwrite(data, [datpath,savedat,'DAT.fdt']); % resave as 'DAT'

    end;
    
    fprintf('\nSaving new spectral data as %sDAT.fdt\n',savedat);%%%%%%%%%%%%%%%%%%%%%%
    floatwrite(data, [datpath,newsave,'DAT.fdt']);
    
    % specify pca dimensions:
    if ~isempty(pcfac)
        pcs =  round(sqrt(size(data,2)/pcfac));
        if pcs >= size(data,1)
            pcfac = pcfac*2;
            pcs =  round(sqrt(size(data,2)/pcfac));
        end;
        fprintf('\nPCA''ing to %s\n',int2str(pcs));
    else
        pcs = size(data,1); % no pca
        fprintf('\nNo PCA reduction\n');
    end;
    %%-------------------
    if ~isempty(auxpath)
        datpath = auxpath;
    end;
    
    numframes = size(data,2);

    %[pc,eigvec,sv] = runpca(data,pcs);
    %[eigvec,sv,pc] = pcsquash(data,pcs); % takes longer cuz # pcs is not chosen until end
    [U,S,eigvec] = svds(data',pcs);% if you scale 'acts', you can't scale the 'eigvec'
    pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp
    clear U S data sv 
    eigfile = [newsave,'EIGVEC.fdt'];% save with  data
    floatwrite(eigvec, [datpath,newsave,'EIGVEC.fdt']);
    clear eigvec
    fprintf('\nSaving PCA''d data for ICA as %s.fdt.\n',savedat);%%%%%%%%%%%%%%%%%%%%%%
    floatwrite(pc, [datpath,newsave,'.fdt']);
    clear pc 
    
    numrows = pcs;
    
    % generate an ICA script 
    ICA_SCRIPT = [datpath,'SpecDecomp.sc'];
    fid = fopen( ICA_SCRIPT, 'w');
    fprintf(fid, 'DataFile %s\n', [datpath,newsave,'.fdt']);
    fprintf(fid, 'chans %d\n', pcs);
    fprintf(fid, 'frames %d\n',numframes);
    fprintf(fid, 'WeightsOutFile %s\n', [datpath,newsave,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [datpath,newsave,'.sph']);  
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
    fprintf('\nRunning ICA in Linux, may take a long time...\n');
    fprintf('\n\n');
    fprintf('%s',run_ica_str);
    %[status, result] = system(run_ica_str);
    fprintf('done.\n\n');
    %%%%%

    str = ['save ',datpath,newsave,'.mat meanpwr numrows numframes freqs freqscale keeptrack rmepochs dstrials pcs rowmeans complist overlap datset eigfile'];eval(str);    
    
