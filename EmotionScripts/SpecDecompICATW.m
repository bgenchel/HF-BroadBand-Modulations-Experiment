% takes a 2D matrix of spectral power data (freqs x windows). Calculates power in each epoch (2nd dim) and 
% performs spectral decomposition in the TW fashion
%
%
% INPUTS:
% alldat -- [2D matrix] dimensions: freqs x time*epochs for one IC (designed for long epochs)
% datpath -- [string] data file path where 'savedat.fdt' will be saved
% savedat -- [string] name to save data float file as
% freqs -- [vector] corresponds to columns of pwr
% freqscale -- ['linear','quad' or 'log'] 'quad' makes quadratically spaced freq bins
% pcs -- [integer] number of PCA dimensions to retain in final spectral decomposition 

function [numtrials,numframes,freqs, pcs, rowmeans] = SpecDecompICATW(pwr,datpath,savedat,freqs,pcs,subwin);

    
    if strcmp(subwin,'on')
        % subtract mean spectrum:----------------------------
        fprintf('\nRemoving mean power spectrum ...\n');
        rowmeans =  squeeze(mean(pwr,2));  % across windows
        pwr = pwr - repmat(rowmeans,[1 size(pwr,2)]);
    else
        rowmeans = [];
    end;
    fprintf('\nSaving spectral data for ICA...\n');
    
    floatwrite(pwr, [datpath,savedat,'.fdt']);
    
    if isempty(pcs)
        pcs =  round(sqrt(size(pwr,2)/25));
    end;    
    fprintf('\nWill retain %s PCA dimensions.\n',int2str(pcs));
  
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
    str = ['save ',datpath,savedat,'.mat numtrials numframes freqs pcs rowmeans'];eval(str);    
    
    fprintf('\ndone.\n');
