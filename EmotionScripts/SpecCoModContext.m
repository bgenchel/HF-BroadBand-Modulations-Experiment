% Calls in data for a subject, runs spectral analysis on short bits of data and then decomposes with ICA
% * currently for emotion data only!!!*
%
% [meanpwr,numtrials,numframes, freqs, keeptrack,rmepochs,dstrials, pcs,rowmeans,addmat, comment] = SpecCoModContext(datset,datpath,complist,savedat,frqlim,pcs,pcnum,overlap,scalefac,evonly);
%
% datset: (cell array of strings) of dataset names of continuous data to use for spectral decomposition
% datpath: (string) fullpath to directory where all datasets are found and to which data will be stored
% complist: (vector) component indices to use for decomposition
% savedat: (string) name of float files in which to save input matrix and output wts and sph matrices
% pcs -- number of final PCA dimensions to retain before ICA. 
% overlap -- overlap factor:  4= 75% overlap,2=50% overlap,1=no overlap,.75=25% gap,.5=50% gap
% evonly: if not empty, will select data between 'Active' and 'Done' events (NF)
% pcnum -- [integer] number of pca dimensions from spectral data to retain
%
% OUTPUT:
% dstrials: [vector] Gives total number of epochs for each dataset (sum(dstrials) = numtrials)
% comment -- documents the parameters of the decomposition protocol
% rowmeans -- [matrix] numtrials X 1; mean removed from each row, AFTER removing component means, 
%             but just before ICA

function [meanpwr,numtrials,numframes,freqs,keeptrack,rmepochs,dstrials, pcs,rowmeans,addmat, comment] = SpecCoModContext(datset,datpath,complist,savedat,frqlim,pcs,pcnum,overlap,scalefac,evonly);
    
    comment = ['Spectral decomposition of component power in single trials; keeptrack matrix gives the start and end; adds indicator matrix to bottom of spectral modulations for each trial (pca''d specs+questions X trials)'];

    for ds = 1:length(datset)
        rmep = [];
        EEG = pop_loadset( datset{ds},datpath,'all');
        %if ~isempty(evonly)
        %    astart = find(strcmp({EEG.event.type},'Active'));
        %    astop = find(strcmp({EEG.event.type},'Done'));
        %    EEG = pop_select( EEG, 'point',[EEG.event(astart).latency EEG.event(astop).latency] );
        %end; 
        %EEG.event(1:end-1)=[]; 
        clear dat
        for evtm = EEG.srate:EEG.srate/overlap:size(EEG.data,2)-EEG.srate 
            % create events to make overlapping 1 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';% for string event codes    
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-.5 .5]);
        EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
        EEG = eeg_checkset(EEG);
        fprintf('\nRunning auto-rejection protocol...\n');
        EEG = pop_eegthresh(EEG,0, complist ,-300,300,EEG.xmin,EEG.xmax,0,0);
        rmep(end+1:end+length(find(EEG.reject.rejthresh)))=find(EEG.reject.rejthresh);
        EEG = pop_rejkurt(EEG,0,complist ,4,4,0,0);        
        rmep(end+1:end+length(find(EEG.reject.icarejkurt)))=find(EEG.reject.icarejkurt);
        EEG = pop_jointprob(EEG,0,complist ,4,4,0,0);
        rmep(end+1:end+length(find(EEG.reject.icarejjp)))=find(EEG.reject.icarejjp);
        EEG = pop_jointprob(EEG,0,complist ,4,4,0,0);
        rmep(end+1:end+length(find(EEG.reject.icarejjp)))=find(EEG.reject.icarejjp);
        rmat = zeros(1,length(rmep));
        rmat(rmep) = 1;
        EEG = pop_rejepoch( EEG,rmat,0);
        rmep = unique(rmep);  
        rmepochs{ds} = rmep;
        EEG.icaact = [];
        EEG = eeg_checkset(EEG);
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
    fprintf('\nCalculating pwr for each 1 second epoch...\n');
    pwr = zeros(257,size(alldat,3),length(complist));
    for cp = 1:length(complist)
        for epoch = 1:size(alldat,3)
            [pwr(:,epoch,cp) frqs] = pwelch(alldat(cp,:,epoch),256,128,512,EEG.srate); 
            pwr(:,epoch,cp) = 10*log10(pwr(:,epoch,cp));% need to convert to log for 'modulator' interpretation
        end;   
    end;
    fr = find(frqs > frqlim(1) & frqs < frqlim(2));
    freqs = frqs(fr); pwr = pwr(fr,:,:);
    clear alldat
    % calculate mean spectrum from all comps 
    % subtract mean spectrum
    fprintf('\nRemoving mean for each component/trial...\n');
    for cp = 1:size(pwr,3)
        meanpwr(cp,:) = mean(pwr(:,:,cp),2)';
        mnpwr =  mean(pwr(:,:,cp),2);   
        for ep = 1:size(pwr,2)
            pwr(:,ep,cp) = pwr(:,ep,cp) - mnpwr;
        end;        
    end;    
    
    pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));
    fprintf('\nRemoving mean for each trial (all comps)...\n');
    for tt = 1:size(pwr,1)
        rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
        pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (all comps)
    end;  
    
    fprintf('\nSaving spectral data as %sSPECmods.fdt...\n',savedat);
    floatwrite(pwr, [datpath,savedat,'SPECmods.fdt']);

    % pca reduce the data
    [pc,eigvec,sv] = runpca(pwr);
    fprintf('\nSaving PCA eigenvector data as %sEIGVEC.fdt...\n',savedat);
    floatwrite(eigvec(:,1:pcnum), [datpath,savedat,'EIGVEC.fdt']); 
    %divstd = mean(std(pc)); % mean of stds within trials
    %pc = pc/divstd; % make a zscore to get within range of context matrix
    pc = pc*scalefac; % scale up the pc values from spectral decomp
    
    % make the addmat matrix
    addmat = ones(15,sum(dstrials)); addmat = addmat*-1;
    for e = 1:15
        addmat(e,sum(dstrials(1:e-1))+1:sum(dstrials(1:e))) = 1;
    end;     
    addmat = repmat(addmat,[1 length(complist)]);

    combdata = [pc(1:pcnum,:);addmat];
    fprintf('\nSaving combined data as %s.fdt for ICA...\n',savedat);
    floatwrite(combdata, [datpath,savedat,'.fdt']); 
    
    if isempty(pcs)
        pcs =  rank(combdata);
    end;
    % generate an ICA script 
    ICA_SCRIPT = [datpath,'SpecCtxDecomp.sc'];
    fid = fopen( ICA_SCRIPT, 'w');
    fprintf(fid, 'DataFile %s\n', [datpath,savedat,'.fdt']);
    fprintf(fid, 'chans %d\n', size(combdata,1));
    fprintf(fid, 'frames %d\n',size(combdata,2));
    fprintf(fid, 'WeightsOutFile %s\n', [datpath,savedat,'PC',int2str(pcs),'.wts']);
    fprintf(fid, 'SphereFile %s\n', [datpath,savedat,'PC',int2str(pcs),'.sph']);  
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
    numtrials = size(combdata,1);
    numframes = size(combdata,2);
    
    
    fprintf('\ndone.\n');
