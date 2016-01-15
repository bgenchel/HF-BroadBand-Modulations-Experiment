% same as SpecCoModAnal.m but uses newtimef instead to find single-trial (sliding window) spectra
%
%  [numtrials, frameno, freqs, times, keeptrack,meanpwr,dstrials] = SpecCoModTimef(datset,fullpath,complist,pcs,name, frqlims,wave,prat,wsize,timepnts);
%
% INPUTS:
% datset -- [string or cell array of strings] dataset names to run spectral analysis on
% fullpath -- [string] full path where datset is/are found and results will be saved
% complist -- [vector] of component indices to include in analysis
% pcs -- [number] of factors to reduce data to  {15}
% name -- [string] name of file in which to store results (in 'fullpath')
% frqlims -- [minfrq maxfrq] in hz to calculate time/freq values for
% wave -- [0 | number | [number decimal]] 0=FFT, 3=3 wavelets, [3 .5] tapered 3 wavelet window
% prat -- [number] padratio for time/frequency analysis (higher number -> more freqs) {1}
% wsize -- [number] window size for time/frequency analysis {same as sampline rate}
% timepnts -- [integer] number of time points to calculate for each dataset. [] to set by data length (1/sec)


function [numtrials, frameno, freqs, times, keeptrack,meanpwr,dstrials] = SpecCoModTimef(datset,fullpath,complist,pcs,name, frqlims,wave,prat,wsize,timepnts)

    
    keeptrack = zeros(0,4); % comp,datset,start idx,# trials in that dataset
    new = 1;pceigvec = []; addeigvec = [];
    if iscell(datset)
        for n = 1:length(complist)
            for ds = 1:length(datset)
                EEG = pop_loadset( datset{ds},fullpath);
                EEG=eeg_checkset(EEG);
                if isempty(timepnts)
                    tmpnts = size(EEG.icaact,2)/EEG.srate;
                else 
                    tmpnts = timepoints;
                end;
                [tf, freqs, times, itcvals] = timefreq(EEG.icaact(complist(n),:)', EEG.srate,'wavelet',wave, 'tlimits',[0 EEG.xmax*1000],'winsize',wsize,'freqs',frqlims,'padratio',prat, 'ntimesout',tmpnts,'freqscale','log');
                fprintf('\nRemoving mean of each t/f row across trials.\n\n',int2str(n),int2str(length(complist)));
                
                tf = abs(tf).^2;
                tf = 10*log10(tf);
                maxzscore = max(tf(:)); 
                tf(find(tf < -maxzscore)) = -maxzscore; % truncate log dropout points
                if ds == 1
                    onecomp = tf';
                else
                    onecomp(end+1:end+size(tf,2),:) = tf';
                end;                
                if size(keeptrack,1) == 0
                    keeptrack(end+1,:) = [complist(n),ds,1,tmpnts];
                else
                    keeptrack(end+1,:) = [complist(n),ds,keeptrack(end,3)+keeptrack(end,4),tmpnts];
                end;
                if n == 1
                    dstrials(1,ds) = tmpnts;
                end;                
            end;           
            for ff = 1:size(onecomp,2)%take out baseline (mean of each freq)
                meanpwr(n,ff) =  mean(onecomp(:,ff));
                onecomp(:,ff) = onecomp(:,ff) - mean(onecomp(:,ff));
            end;
            if new == 1 
                alltrials = onecomp;new=0;
            else
                alltrials(:,end+1:end+size(onecomp,2)) = onecomp;
            end;            
            fprintf('\nDoneComp: %s of %s\n\n',int2str(n),int2str(length(complist)));
        end;
    else
        EEG = pop_loadset( datset,fullpaths);
        EEG=eeg_checkset(EEG);
        dstrials = size(EEG.icaact,3);
        %%%%%%%%%%%%%%%
        for n = 1:length(complist)
            data = squeeze(EEG.icaact(complist(n),:,:));
            [tf, freqs, times, itcvals] = timefreq(data, EEG.srate,'wavelet',wave, 'tlimits',[EEG.xmin*1000 EEG.xmax*1000],'winsize',wsize,'freqs',frqlims,'padratio',prat, 'ntimesout',timepnts,'freqscale','log');
            fprintf('\nRemoving mean of each t/f row across trials.\n\n',int2str(n),int2str(length(complist)));
            tf = abs(tf).^2;
            tf = 10*log10(tf);
            maxzscore = max(tf(:)); 
            tf(find(tf < -maxzscore)) = -maxzscore; % truncate log dropout points
            
            for ff = 1:size(tf,1)%take out baseline (mean of each freq)
                meanpwr(n,ff) =  mean(tf(ff,:));
                tf(ff,:) = tf(ff,:) - mean(tf(ff,:));
            end;
            if new == 1 
                alltrials = tf';new=0;
            else
                alltrials(:,end+1:end+size(tf,1)) = tf';
            end;            
            if size(keeptrack,1) == 0
                keeptrack(end+1,:) = [complist(n),1,1,timepnts];
            else
                keeptrack(end+1,:) = [complist(n),1,keeptrack(end,4)+keeptrack(end,5),timepnts];
            end;
            fprintf('\nDoneComp: %s of %s\n',int2str(n),int2str(length(complist)));
        end;
        dstrials = timepnts;
        EEG=[];clear tf itcvals 
        clear  basetfdata tf1 tf2 normtf data tfnorm  tfdata   
    end; 

    mntfpwr = mean(alltrials,2);
    for tff = 1:size(alltrials,1)
        alltrials(tff,:) = alltrials(tff,:) - mntfpwr(tff);
    end;
    floatwrite(alltrials, [fullpath,name,'.fdt']);  
    numtrials = size(alltrials,1);
    frameno = size(alltrials,2);
    
    % generate an ICA script 
    fid = fopen([fullpath,'TempICAinstruct.sc'], 'w');
    fprintf(fid, 'DataFile %s\n', [fullpath,name,'.fdt']);
    fprintf(fid, 'chans %s\n', int2str(numtrials));
    fprintf(fid, 'frames %s\n',int2str(frameno));
    fprintf(fid, 'WeightsOutFile %s\n', [fullpath,name,'.wts']);
    fprintf(fid, 'SphereFile %s\n', [fullpath,name,'.sph']);  
    fprintf(fid, 'sphering on\n');
    fprintf(fid, 'bias on\n');
    fprintf(fid, 'extended 1\n');
    fprintf(fid, 'pca %s\n',int2str(pcs)); 
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
    
    run_ica_str = ['/data/common/matlab/ica_linux2.4 < ',fullpath,'TempICAinstruct.sc'];
    fprintf('\n Running ICA in Linux. May take a long time...\n');
    [status, result] = system(run_ica_str);  % won't see progress this way
                                             % Results get saved to disk 
    [status, result] = system(['\rm ',fullpath,'TempICAinstruct.sc']);  % won't see progress this way

