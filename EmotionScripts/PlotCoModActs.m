% plots activation timecourses of specified ICs with highlighted windows corresponding to highly weighted specified IMs
%
%
%
%
%
%
%
% stdcut -- [number] of standard deviations from the mean to plot as highlighted
%           to plot positive weights, enter positive std. For neg wts, enter neg std (ie, -1.5)
%           Plots each IM as with a different color background.



function  PlotCoModActs(datset,fullpath,ics,ims,savedat,stdcut);
    
    s = load([fullpath,savedat,'.mat']); 
    sph=floatread([fullpath,savedat,'.sph'],[s.numrows s.numrows],[],0); 
    wts=floatread([fullpath,savedat,'.wts'],[s.pcs s.numrows],[],0); 
    ws = wts*sph;     winv = pinv(ws); clear wts sph ws   
    
    for im = 1:length(ims)
        zwinv = zscore(winv(:,ims(im)));
        if stdcut > 0
            ti = find(zwinv > stdcut);
        else
            ti = find(zwinv < stdcut);
        end;
        trialidx{im} = ti;
    end;
    for ds = 1:length(datset)
        cont = 0;
        if ds == 1
            relidx = [1:s.dstrials(ds)];
        else
            relidx = [sum(s.dstrials(1:ds-1))+1:sum(s.dstrials(1:ds))];
        end;
  
        EEG = pop_loadset( datset{ds},fullpath,'all');
        if length(size(EEG.icaact)) < 3
            for evtm = EEG.srate:EEG.srate/s.overlap:size(EEG.data,2)-EEG.srate 
                % create events to make overlapping 1 sec epochs
                EEG.event(end+1) =  EEG.event(1);% appends events to the end
                EEG.event(end).latency = evtm;
                EEG.event(end).type = 'fake';% for string event codes    
            end;
            EEG = pop_epoch( EEG,{'fake'} , [-.5 .5]);
            EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
            EEG = eeg_checkset(EEG);                
            EEG=pop_select(EEG,'notrial',s.rmepochs{ds});
        end;
        if ds == 1
            colacts = zeros(length(ics),size(EEG.data,2),0);
        end;        
        if length(relidx) ~= size(EEG.icaact,3)                    
            fprintf('\n\nERROR: number of epochs in dataset do not equal that in the trial index vector\n\n'); 
            break;return;
        else
            colacts(:,:,end+1:end+size(EEG.icaact,3)) = EEG.icaact(ics,:,:);
        end;
    end;
    % plot the results
    pl = 1; imcol = hsv(length(trialidx));
    nchn = zeros(1,length(ics));
    for im = 1:length(trialidx)
        for tr = 1:length(trialidx{im})
            rejmat(pl,:) = [(trialidx{im}(tr)-1)*size(colacts,2) (trialidx{im}(tr)-1)*size(colacts,2)+size(colacts,2) imcol(im,1) imcol(im,2) imcol(im,3) nchn];
            pl = pl+1;
        end;
    end;
    
        
    eegplot(colacts,'srate',EEG.srate,'winrej', rejmat,'spacing',12,'winlength',15);
