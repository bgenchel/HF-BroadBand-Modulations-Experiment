% takes a 2D matrix and decomposes with ICA. PCA's if necessary
% takes an optional context indicator matrix with same # columns as datmat
%
% [numtrials,numframes,pcmat,pcidx,rowmeans,rowmeanidx,eigfile,cxteigfile] = DecompMat(datmat,idxmat,wtsphname,writepath,pcmat,pcidx,pcfac,freqs,times);
% 
% keeptrack -- matrix with nrows equal to nrows or ncols of datmat, any number of columns
%              to identify the data. Matrix is passed through untouched to .mat file.
% meanpwr -- will be saved in the resulting .mat file. This should be a ICs (or IC pairs for
%            cross spectrum) x freqs matrix of mean spectra taken over all windows and datasets.
%
% pcfac -- [number] of desired points per weight in final ERSP matrix

function [numrows,numframes,pcmat,pcidx,rowmeans,rowmeanidx,eigfile,cxteigfile] = DecompMat(datmat,idxmat,wtsphname,writepath,complist,keeptrack,pcmat,pcidx,pcfac,freqs,times,freqscale,meanpwr);

    
    if isempty(pcfac)
        pcfac = 1;
    end;
    rowmeanidx = [];rowmeanspec=[];cxteigfile = []; eigfile = [];
    rowmeans = mean(datmat,2);
    datmat = datmat - repmat(rowmeans,[1 size(datmat,2)]);
    floatwrite(datmat, [writepath,wtsphname,'DAT.fdt']);
    
   
    fprintf('\n Running PCA on spectral data...\n');
    if isempty(pcmat)
        pcmat =  round(sqrt(size(datmat,2)/pcfac));
    end;
    
    fprintf('\n Retaining %s dimensions...\n',int2str(pcmat));
    % run PCA on input data
    [U,S,pceigvec] = svds(datmat',pcmat);% if you scale 'acts', you can't scale the 'eigvec'
    pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp
    clear U S data sv 
    eigfile = [wtsphname,'MATEIG.fdt'];% save with  data
    floatwrite(pceigvec, [writepath,eigfile]);
    clear pceigvec
    
    %[pc,pceigvec,sv] = runpca(datmat,pcmat); 
    datmat = pc; % take the pcs for decomp
    clear pc 
  
    if ~isempty(idxmat)
        if isempty(pcidx)
            rowmeanidx = mean(idxmat,2);
            idxmat = idxmat - repmat(rowmeanidx,[1 size(idxmat,2)]);
            pcidx =  rank(idxmat)-1;
        end;
        if pcidx < size(idxmat,1)
            fprintf('\n Running PCA on indicator matrix...\n');
            fprintf('\n Retaining %s dimensions...\n',int2str(pcidx));
            [U,S,idxeigvec] = svds(idxmat',pcidx);% 
            pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp
            %[pc,idxeigvec,sv] = runpca(idxmat,pcidx); 
            % pc is what you will use (analogoous to act)
            idxmat = pc;
            cxteigfile = [wtsphname,'IDXEIG.fdt'];
            floatwrite(idxeigvec, [writepath,cxteigfile]);
            clear pc idxeigvec
        end;  
        datmat(end+1:end+size(idxmat,1),:) = idxmat;
    else
        pcidx = [];
    end;
        
    floatwrite(datmat, [writepath,wtsphname,'.fdt']);  
    numrows = size(datmat,1);
    numframes = size(datmat,2);
    
    pcs = size(datmat,1); % combined (or not) full matrix dimensions
    %-----------------------------------------------------
    % generate an ICA script 
    fid = fopen([writepath,wtsphname,'.sc'], 'w');
    %fid = fopen([writepath,'TempICAinstruct.sc'], 'w');
    fprintf(fid, 'DataFile %s\n', [writepath,wtsphname,'.fdt']);
    fprintf(fid, 'chans %s\n', int2str(numrows));
    fprintf(fid, 'frames %s\n',int2str(numframes));
    fprintf(fid, 'WeightsOutFile %s\n', [writepath,wtsphname,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [writepath,wtsphname,'.sph']);  
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
    run_ica_str = ['/data/common/matlab/ica_linux2.4 < ',writepath,wtsphname,'.sc'];
    %run_ica_str = ['/data/common/matlab/ica_linux2.4 < ',writepath,'TempICAinstruct.sc'];
    fprintf('\n Running ICA in Linux. May take a long time...\n');
    
    [status, result] = system(run_ica_str);  % won't see progress this way
                                             % Results get saved to disk 
    [status, result] = system(['\rm ',writepath,wtsphname,'.sc']);  % won't see progress this way
    %[status, result] = system(['\rm ',writepath,'TempICAinstruct.sc']);  % won't see progress this way
    str = ['save ',writepath,wtsphname,'.mat numrows numframes complist pcs pcmat rowmeans rowmeanidx freqs times pcidx eigfile cxteigfile freqscale meanpwr keeptrack pcfac'];eval(str);    
        

    
