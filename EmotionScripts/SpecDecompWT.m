% takes a 3D matrix of EEG IC activations (with winv multiplied). 
% Calculates power in each epoch (2nd dim) and
% performs spectral decomposition in the WT fashion
%
%
% INPUTS:
% alldat -- [3D matrix] dimensions: ICs x time x epochs
% datpath -- [string] data file path where 'savedat.fdt' will be saved
% savedat -- [string] name to save data float file as
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% freqscale -- ['linear','quad' or 'log'] 'quad' makes quadratically spaced freq bins
% overlap -- number of ms between each window center
% times -- [vector of latencies] of length equal to size(alldat,2) giving time stamp in
%          ms of each sample.
% pcfac -- [number] number of points per weight to have in final ICA (will PCA down accordingly)
%          If 'pcfac' is > 10, then this variable specifies the number of PCA dimensions
%          to keep (does not scale by the number of data points).
% srate -- [integer] sampling rate of the input data matrix 'alldat'
% datsetvec -- [vector] of indices with the length = size(alldat,3) saying for each
%              epoch which dataset (trial type) the epoch belongs to. Will be added
%              to the 'keeptrack' output variable.
% complist -- [vector] of component indices corresponding to the rows of alldat
% pwrdecomp -- [0 or [cycle1 cycle2 winsize]] if 0, will use pwelch to calculate standard
%              FFT. If number of cycles and winsize specified, then will perform wavelet
%              decomp. If cycle2 is larger than cycle1, then a linear scaling of the
%              number of wavelets will be created between the lowest and highest freqs
%              requested. Default is 0. for wavelets, recommend 6 for low frqs,20-40 for high

function [numrows,numframes,freqs, pcs,rowmeans] = SpecDecompWT(alldat,datpath,savedat,frqlim,nfreqs,freqscale,timepoints,times,pcfac,srate,datsetvec,complist,pwrdecomp);

typeslash = '/'; % '\' = windows, '/' = linux
sys = 'linux'; %'windows' or 'linux'

if isempty(nfreqs)
   nfreqs = frqlim(end)-frqlim(1);
end;
if ~exist('datsetvec')
   datsetvec = ones(size(alldat,3),1); % assume all one dataset if missing
elseif isempty(datsetvec)
   datsetvec = ones(size(alldat,3),1);
end;

if strcmp(freqscale, 'quad') % make quadratic freq spacing
   freqs = linspace(sqrt(frqlim(1)), sqrt(frqlim(end)), nfreqs);
   freqs = freqs.^2;
elseif strcmp(freqscale, 'log') % make quadratic freq spacing
   freqs = linspace(log(frqlim(1)), log(frqlim(end)), nfreqs);
   freqs = exp(freqs);
else % if linear-spacing
   freqs = linspace(frqlim(1), frqlim(2),nfreqs);
end;
if ~isstr(alldat) % if NOT pre-saved data
  fprintf('\nCollecting data windows from input data...\n');
   %timepoints = [times(1)+500:overlap:times(end)-500]; % 50 ms spacing
   dats = [];
   keeptrack = zeros(0,3);
   for tmpnt = 1:length(timepoints) % for all specified time points
      clear dat
      tm = times- timepoints(tmpnt);
      [val idx] = min(abs(tm)); %
      for cp = 1:length(complist)
         if pwrdecomp == 0 % FFT using pwelch
            dat(cp,:,:) = alldat(cp,idx - srate/2+1:idx + srate/2,:);
         else % need slightly longer
            dat(cp,:,:) = alldat(cp,idx - srate/2+1:idx + srate/2+1,:);
         end;
      end;
      if tmpnt == 1
         dats = dat;
      else
         dats = cat(3,dats,dat);
      end;
      keeptrack(end+1:end+size(dat,3),:) = [datsetvec [1:size(dat,3)]' repmat(timepoints(tmpnt),[size(dat,3) 1])];
   end
   %floatwrite(dats, [datpath,savedat,'RAW.fdt']);% save raw data for subsequent decomps

else % not datasets, just RAW data, already did the windowing
   str = ['load ',datpath,savedat,'.mat'];eval(str);
   dats = floatread([datpath,alldat,'RAW.fdt'],[length(complist) srate+1 inf],[],0);
end;
alldat = dats; clear dats
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
      elseif rempnts2 == rempnts1
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
% calculate and remove mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cp = 1:size(pwr,3)
   meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
   pwr(:,:,cp) = pwr(:,:,cp) - repmat(meanpwr(cp,:),[size(pwr,1) 1]); % remove mean
end;

clear alldat

pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));% windows x comps*freqs

fprintf('\nRemoving mean for each spectral window...\n');
for tt = 1:size(pwr,1)
   rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
   pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (windows)
end;

fprintf('\nSaving spectral data as %s\n',[datpath,savedat,'DAT.fdt']);%%%%%%%%%%%%%%%%%%%%%%
floatwrite(pwr, [datpath,savedat,'DAT.fdt']);
fprintf('\nYour matrix size is %s...\n',int2str(size(pwr)))

if isempty(pcfac)
  percvar = .5;%--------------------------------------------
  fprintf('\nFinding dims accounting for > %s%% of total variance for first 200 dims. \n',int2str(percvar*100));
  [U,S,V] = svds(pwr',200);% if you scale 'acts', you can't scale the 'eigvec'
  S = max(S);
  totvar = sum(S);
  for ndims = 1:length(S)
    x=sum(S(1:ndims));
    if x > totvar*percvar % accounting for ?% of the variance
      break;
    end;
  end;
  pcs = ndims;
  
elseif pcfac < 12 % otherwise use as number of pcs
  pcs =  round(sqrt(size(pwr,2)/pcfac));
  if pcs >= size(pwr,1)
    pcfac = pcfac*2;
    pcs =  round(sqrt(size(pwr,2)/pcfac));
  end;
else
  pcs =  pcfac;
end;

 %--------------------------------------------
  fprintf('\nPCA''ing to %s dimensions.\n',int2str(pcs));%%%%%%%%%%%%%%%%%%%%%%
  [U,S,V] = svds(pwr',pcs);% if you scale 'acts', you can't scale the 'eigvec'
  pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp of pc
  eigvec = V;
  
% eigvec is the winv matrix
%[eigvec,sv,pc] = pcsquash(pwr,pcs);

fprintf('\nSaving PCA''d data as %s for ICA.\n',[datpath,savedat,'.fdt']);%%%%%%%%%%%%%%%%%%%%%%
floatwrite(pc, [datpath,savedat,'.fdt']);

fprintf('\nSaving eigenvectors as %s for reconstruction of trial data.\n',[datpath,savedat,'EIGVEC.fdt']);
floatwrite(eigvec, [datpath,savedat,'EIGVEC.fdt']);
eigfile = [savedat,'EIGVEC.fdt'];% save in .mat

numrows = size(pwr,1);
numframes = size(pwr,2);
clear pceigvec

if strcmp(sys,'linux')
   % generate an ICA script
   ICA_SCRIPT = [datpath,'SpecDecomp.sc'];
   fid = fopen( ICA_SCRIPT, 'w');
   fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
   fprintf(fid, 'chans %d\n', pcs);
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
   fprintf('\nRunning ICA in Linux, may take a long time...\n');
   [status, result] = system(run_ica_str);
else
   [weights,sphere] = runica(pc,'extended',1,'maxsteps',2000,'stop',1e-8);
   floatwrite(weights, [datpath,savedat,'.wts']);
   floatwrite(sphere, [datpath,savedat,'.sph']);

   % try amica
   %runamica([datpath,savedat,'.fdt'],[datpath,savedat,typeslash],numrows,numcols,'qsub','off','do_reject',0,'max_iter',5000,'do_newton',1,'numrej',3);%
   %mod = loadmodout([datpath,savedat,typeslash]);
   %floatwrite(mod.W, [datpath,savedat,'.wts']);
   %floatwrite(mod.S, [datpath,savedat,'.sph']);
end;
%%%%%
str = ['save ',datpath,savedat,'.mat meanpwr numrows numframes freqs freqscale keeptrack  pcs rowmeans complist eigfile srate'];eval(str);

fprintf('\ndone.\n');


