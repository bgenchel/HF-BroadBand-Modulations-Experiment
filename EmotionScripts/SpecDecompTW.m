% takes one or more datasets (epoched or continuous) and finds spectral variations
% from the mean across time points
%
% SpecDecompTW(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath,nfreqs,pwrdecomp);
%
% datset: (cell array of strings) of dataset names of continuous data to use for spectral decomposition
% datpath: (string) fullpath to directory where all datasets are found and to which data will be stored
% complist: (vector) component indices to use for decomposition
% savedat: (string) name of float files in which to save input matrix and output wts and sph matrices
% frqlim -- [minfreq maxfreq]
% freqscale -- 'linear' or 'quad' for quadratic frequency spacing with more bins at lower freqs
% pcfac -- number of points per weight for final decomposition
% overlap -- overlap factor:  4= 75% overlap,2=50% overlap,1=no overlap,.75=25% gap,.5=50% gap
%            for pre-epoched datasets, enter a ms value for time point spacing.
% auxpath -- if not [], uses this directory to store all generated data (.fdt,.wts,.sph)
% nfreqs -- [number] number of output frequency bins (only for 'quad' or
%           'log'; for 'linear', use frqfac (internal to program))
% rej = 'on' or 'off'; % 'on' to reject noisy epochs (for continuous of pre-epoced data)
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
%
%
%



function SpecDecompTW(datset,datpath,complist,savedat,frqlim,freqscale,pcfac,overlap,auxpath,nfreqs,rej,pwrdecomp);



typeslash = '\'; % '\' = windows, '/' = linux
sys = 'windows'; %'windows' or 'linux'
continuous = 'false'; % as default
if ~exist('pwrdecomp')
   pwrdecomp = 0;
end;
if ~exist('nfreqs')
   nfreqs = frqlim(end)-frqlim(1);
elseif isempty(nfreqs)
   nfreqs = frqlim(end)-frqlim(1);
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

if ~exist('auxpath')
   auxpath = [];
end;
if iscell(datset) % if not cell, then 'datset' is a RAW data file
   EEG = pop_loadset( datset{1},datpath);
   keeptrack = zeros(0,3);
   for ds = 1:length(datset)
      EEG = pop_loadset( datset{ds},datpath);
      if length(size(EEG.data)) < 3 % continuous data
         continuous = 'true';
         clear dat
         x = EEG.event(1).type;
         %%%%%%%%%%  Divide continuous data into 2 sec epochs:
         for evtm = EEG.srate:EEG.srate/overlap:size(EEG.data,2)-EEG.srate
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
            EEG = pop_epoch( EEG,{'fake'} , [-.5 .5]);% was-1 1 recently
         else
            EEG = pop_epoch( EEG,{1000} , [-.5 .5]);
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
         alleps = [1:size(EEG.icaact,3)];
         EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)] ,-1000,1000,EEG.xmin,EEG.xmax,0,0);
         numrej = length(find(EEG.reject.rejthresh));  % count number of epochs marked
         if numrej > 0
            rmep(1,end+1:end+length(find(EEG.reject.rejthresh))) = alleps(find(EEG.reject.rejthresh));
            alleps(find(EEG.reject.rejthresh)) = [];
            EEG = pop_rejepoch( EEG,EEG.reject.rejthresh,0); % actually reject high prob epochs
            fprintf('\nRe-baselining after large amplitude artifact removed...\n');
            EEG = pop_rmbase( EEG,[EEG.xmin EEG.xmax]);
         end;
         %--------------------------------------------

         startlim = 4.5; % start probability threshold that will be up'ed if too many trials found
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

      if strcmp(continuous,'true') % continuous data
         %%%%  Choose only component activations of interest
         for cp = 1:length(complist)
            multfact = sqrt(mean(EEG.icawinv(:,complist(cp)).^2));
            dat(cp,:,:) = EEG.icaact(complist(cp),:,:)*multfact; % added multfact back 12-6-05
         end;
         if ds == 1
            alldat = dat;
         else
            alldat = cat(3,alldat,dat);
         end;
         keeptrack(end+1,:) = [sum(keeptrack(:,1))+1 sum(keeptrack(:,1))+size(dat,3) 1];

      else % pre-epoched data, choose time point(s)
         timepoints = [EEG.times(1)+EEG.srate/2:EEG.times(end)-EEG.srate/2]; % 
          for tmpnt = 1:length(timepoints) % for all specified time points
            clear dat
            tm = EEG.times- timepoints(tmpnt);
            [val idx] = min(abs(tm)); %
            for cp = 1:length(complist)
               multfact = sqrt(mean(EEG.icawinv(:,complist(cp)).^2));
               if pwrdecomp == 0 % FFT using pwelch
               dat(cp,:,:) = EEG.icaact(complist(cp),idx - EEG.srate/2+1:idx + EEG.srate/2,:)*multfact;
               else % need slightly longer
               dat(cp,:,:) = EEG.icaact(complist(cp),idx - EEG.srate/2+1:idx + EEG.srate/2+1,:)*multfact;
              end;
            end;
            if ds == 1 & tmpnt == 1
               alldat = dat;
            else
               alldat = cat(3,alldat,dat);
            end;
            keeptrack(end+1:end+size(dat,3),:) = [repmat(ds,[size(dat,3) 1]) [1:size(dat,3)]' repmat(timepoints(tmpnt),[size(dat,3) 1])];
         end
      end;
   end;  clear dat
   floatwrite(alldat, [datpath,savedat,'RAW.fdt']);% save raw data for subsequent decomps
else % not datasets, just RAW data
   str = ['load ',datpath,savedat,'.mat meanpwr numrows numcols freqs freqscale keeptrack rmepochs dstrials  pcs rowmeans complist overlap datset eigfile datset srate'];eval(str);
   alldat = floatread([datpath,datset,'RAW.fdt'],[length(complist) srate+1 inf],[],0);
end;
srate = EEG.srate;
EEG = [];ALLEEG=[];
%%%%%%%%%%%%%%%%  Begin decompostion of single epochs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nCalculating pwr for each epoch...\n');
eplength = size(alldat,2);
if pwrdecomp == 0 % FFT using pwelch
   fprintf('\nCalculating pwr for each %s second epoch...\n',int2str(size(alldat,2)/srate));
   hwin = hanning(size(alldat,2))';
   fprintf('\nUsing hanning windowed data...\n');
   % to get number of freqs:
   [p f] = pwelch(alldat(1,:,1),srate,srate*3/4,freqs,srate);
   pwr = zeros(size(alldat,3),length(f),size(alldat,1));
   for cp = 1:length(complist)
      for epoch = 1:size(alldat,3)
         [pwr(epoch,:,cp) freqs] = pwelch(alldat(cp,:,epoch),srate,srate*3/4,freqs,srate);
         pwr(epoch,:,cp) = 10*log10(pwr(epoch,:,cp));% need to convert to log for 'modulator' interpretation
      end;
      fprintf('\nIC %s done...',int2str(cp));
   end;
elseif length(pwrdecomp) > 2 % use wavelet decomp
   cycles1 = pwrdecomp(1);
   cycles2 = pwrdecomp(2);
   winsize = pwrdecomp(3);
   hwin = hanning(winsize)';
   wavelets = computepsifamilyQodd(freqs,1/srate,cycles1,cycles2,winsize); %
   if size(alldat,2) > winsize
      ndiff = size(alldat,2) - winsize;
      rempnts1 = floor(ndiff/2);
      rempnts2 = ceil(ndiff/2);
      if rempnts2 == 1
         alldat = alldat(:,1:end-rempnts2,:,:);
      elseif rempnts2 > rempnts1
         alldat = alldat(:,rempnts2:end,:);
         alldat = alldat(:,1:end-rempnts2,:);
      elseif rempnts1 == rempnts2
         alldat(:,[1:rempnts1,end-(rempnts2-1):end],:) = [];
      end;
   elseif size(alldat,2) < winsize
      fprintf('\nYour epoch length is %s...\n',int2str(size(alldat,2)))
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cp = 1:size(pwr,3)
   meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
end;
clear alldat

pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3))';% comps*freqs x windows

fprintf('\nRemoving mean for each spectral window...\n');
for tt = 1:size(pwr,1)
   rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
   pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (comp/freq)
end;

if ~isempty(auxpath)
   datpath = auxpath;
end;

fprintf('\nSaving spectral data as ...DAT.fdt\n');%%%%%%%%%%%%%%%%%%%%%%
floatwrite(pwr, [datpath,savedat,'DAT.fdt']);
fprintf('\nYour matrix size is %s...\n',int2str(size(pwr)))

if pcfac < 10 % otherwise use as number of pcs
   pcs =  round(sqrt(size(pwr,2)/pcfac));
   if pcs >= size(pwr,1)
      pcfac = pcfac*2;
      pcs =  round(sqrt(size(pwr,2)/pcfac));
   end;
else
   pcs =  pcfac;
end;

fprintf('\nPCA''ing to %s dimensions.\n',int2str(pcs));%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svds(pwr',pcs);% if you scale 'acts', you can't scale the 'eigvec'
pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp of pc
eigvec = V;

% eigvec is the winv matrix
%[eigvec,sv,pc] = pcsquash(pwr,pcs);

floatwrite(eigvec, [datpath,savedat,'EIGVEC.fdt']);
eigfile = [savedat,'EIGVEC.fdt'];% save in .mat
fprintf('\nSaving PCA''d data as for ICA.\n');%%%%%%%%%%%%%%%%%%%%%%
floatwrite(pc, [datpath,savedat,'.fdt']);

numrows = size(pwr,1);
numcols = size(pwr,2);
clear pceigvec

if strcmp(sys,'linux')
   % generate an ICA script
   ICA_SCRIPT = [datpath,'SpecDecomp.sc'];
   fid = fopen( ICA_SCRIPT, 'w');
   fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
   fprintf(fid, 'chans %d\n', pcs);
   fprintf(fid, 'frames %d\n',numcols);
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
   run_ica_str = [ '/data/MATLAB/binica/ica_linux < ', ICA_SCRIPT ];
   fprintf('\nRunning ICA in Linux, may take a long time...\n');
   [status, result] = system(run_ica_str);
else 
   [weights,sphere] = runica(pc,'extended',1,'maxsteps',2000,'stop',1e-7); 
   floatwrite(weights, [datpath,savedat,'.wts']);
   floatwrite(sphere, [datpath,savedat,'.sph']);
   
   % try amica
   %runamica([datpath,savedat,'.fdt'],[datpath,savedat,typeslash],numrows,numcols,'qsub','off','do_reject',0,'max_iter',5000,'do_newton',1,'numrej',3);%
   %mod = loadmodout([datpath,savedat,typeslash]);
   %floatwrite(mod.W, [datpath,savedat,'.wts']);
   %floatwrite(mod.S, [datpath,savedat,'.sph']);
end;
%%%%%
str = ['save ',datpath,savedat,'.mat meanpwr numrows numcols freqs freqscale keeptrack rmepochs dstrials  pcs rowmeans complist overlap datset eigfile datset srate'];eval(str);

fprintf('\ndone.\n');
