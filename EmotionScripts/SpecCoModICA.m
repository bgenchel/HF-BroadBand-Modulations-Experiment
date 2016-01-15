% takes a 3D matrix of power spectra (epochs,freqs,comps). 
% Calculates power in each epoch (1st dim) and performs spectral decomposition
%
%
% INPUTS:
% alldat -- [3D matrix] dimensions: channels/components X time X epochs
% datpath -- [string] data file path where 'savedat.fdt' will be saved
% savedat -- [string] name to save data float file as
% freqs -- [vector] of frequency bins
% pcs -- [integer] number of PCA dimensions to retain in final spectral decomposition {15}
% srate -- [integer] sampling rate of the input data matrix 'alldat'

function [meanpwr,numtrials,numframes,freqs, pcs, rowmeans] = SpecDecomp(pwr,datpath,savedat,freqs,pcs,submean);

    
  
    percmax = .7;
    if strcmp(submean,'on')
        % subtract mean spectrum
        fprintf('\nRemoving mean for each component/trial...\n');
        for cp = 1:size(pwr,3)
            meanpwr(cp,:) =  squeeze(mean(pwr(:,:,cp),1));   
            for ep = 1:size(pwr,1)
                pwr(ep,:,cp) = pwr(ep,:,cp) - meanpwr(cp,:);
            end;        
        end;    
    else
        meanpwr = [];
    end;
    pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));
    fprintf('\nRemoving mean for each trial (all comps)...\n');
    for tt = 1:size(pwr,1)
        rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
        pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (all comps)
    end;   
    fprintf('\nSaving spectral data for ICA...\n');
    floatwrite(pwr, [datpath,savedat,'.fdt']);
    
    if isempty(pcs)
        fprintf('\n Running PCA to determine number of dims to save...\n',int2str(pcs));
        [pc,pceigvec,sv] = runpca(pwr,30); 
        sv = max(sv); sv = sv/max(sv);
        pcs = length(find(sv > percmax)); % keep dims above specified % of max
        % pc is what you will use (analogoous to act)
        %datmat = pc(find(sv > percmax),:);
    else
        [pc,pceigvec,sv] = runpca(pwr,pcs); 
        %datmat = pc;
    end;    
    fprintf('\n Retaining %s dimensions...\n',int2str(pcs));
    %floatwrite(pceigvec(:,find(sv > percmax)), [datpath,savedat,'ERSPEIG.fdt']);
    clear pc pceigvec  


    
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
    fprintf(fid, 'pca %d\n',pcs);
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
    str = ['save ',datpath,savedat,'.mat numtrials numframes freqs pcs rowmeans meanpwr'];eval(str);    
    
    fprintf('\ndone.\n');
