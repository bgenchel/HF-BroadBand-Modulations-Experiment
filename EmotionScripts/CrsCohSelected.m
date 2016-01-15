% runs crossf on only selected windows/epochs of data
%
% CrsCohSelected(datset,savedat,datpath,rmepochs,complist,pos);
%
% datset -- {cell array} of full dataset names from which to get data for cross coh
% datpath -- {string} of full data directory path where datset can be found
% rmepoch -- {cell array} of indices to remove, specified by output of SpecCoModAnal()
% complist -- [vector] of component indices to find cross coherence between
% pos -- [0|1] if 1, will calculates the positively weighted epochs. if 0, calculates neg epochs
% saves data from each IM in the datpath as 'SelCrsCohIM?.mat'
%
% Author: Julie Onton Aug 18, 2006

function CrsCohSelected(datset,savedat,datpath,rmepochs,complist);
    
    
    s = load([datpath,savedat,'.mat']); 
    sph=floatread([datpath,savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([datpath,savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
    %icamatall = floatread([datpath,savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    %activations = ws*icamatall;    
    winv = pinv(ws); clear wts sph ws
    for ds = 1:length(datset)
        EEG = pop_loadset( datset{ds},datpath,'all');
        clear dat
        for evtm = EEG.srate:EEG.srate/s.overlap:size(EEG.data,2)-EEG.srate 
            % create events to make overlapping 1 sec epochs
            EEG.event(end+1) =  EEG.event(1);% appends events to the end
            EEG.event(end).latency = evtm;
            EEG.event(end).type = 'fake';% for string event codes    
        end;
        EEG = pop_epoch( EEG,{'fake'} , [-.95 1.05]); % need .95 b/c will otherwise skip first epoch
        EEG = pop_rmbase( EEG,[EEG.xmin*1000 EEG.xmax*1000]);
        EEG = eeg_checkset(EEG);
        EEG = pop_rejepoch( EEG,rmepochs{ds},0);
        if size(EEG.data,3) < s.dstrials(ds)
            EEG.data(:,:,end+1) = EEG.data(:,:,end);
            EEG.icaact(:,:,end+1) = EEG.icaact(:,:,end);
        end
        for cp = 1:length(complist)
            multfact = sqrt(mean(EEG.icawinv(:,complist(cp)).^2));
            dat(cp,:,:) = EEG.icaact(complist(cp),:,:)*multfact; % added multfact back 12-6-05
        end;
        if ds == 1
            alldatPRE = dat;            
        else
            alldatPRE(:,:,end+1:end+size(dat,3)) = dat;
        end;  
    end;
    EEG.data = []; EEG.icaact = [];EEG.icawinv = [];
    EEG.icasphere = []; EEG.icaweights = [];
    for im = 1:size(winv,2)
        for pos = 1:2
            imwts = winv(:,im);
            if pos == 1
                imwts = find(imwts > (mean(imwts)+std(imwts)*1.5)); % only positive weights
            else
                imwts = find(imwts < (mean(imwts)-std(imwts)*1.5)); % only negative weights
            end;            
            alldat{pos} = alldatPRE(:,:,imwts);
        end;
        
        cycles = [3 .5];
        alpha = .01;
        padrat = 1;
        wsize = 256;
        frqlim = [0 50];
        type = 'phasecoher';
        
        comps = complist;
        for index1=1:length(complist)-1
            for index2=index1+1:length(complist)
                [coh,mcoh,timesout,freqsout,cohboot,cohang] = newcrossf({alldat{1}(index1,:) ,alldat{2}(index1,:)},{alldat{1}(index2,:) ,alldat{2}(index2,:)}, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000], EEG.srate, cycles, 'alpha', alpha,'padratio', padrat,'winsize',wsize,'newfig','off','type',type,'freqs',frqlim,'freqscale' ,'log','savecoher',0 , 'plotamp' ,'off','plotphase' ,'off' );
                cohercell{index1,index2,1} = coh{1,1};
                cohercell{index1,index2,2} = coh{1,2};
                cohercell{index1,index2,3} = coh{1,3};
                crsangcell{index1,index2,1} = cohang{1,1};
                crsangcell{index1,index2,2} = cohang{1,2};
                crsangcell{index1,index2,3} = cohang{1,3};
                sigdiffs{index1,index2,1} = cohboot{1,1}; 
                sigdiffs{index1,index2,2} = cohboot{1,2};
                sigdiffs{index1,index2,3} = cohboot{1,3};                
                clear coh mcoh cohboot cohang
            end;
            str = ['save ',datpath,'SelCrsCohIM',int2str(im),'.mat cohercell crsangcell sigdiffs freqsout timesout comps alpha']; eval(str)
        end;
        clear cohercell crsangcell sigdiffs
    end;
    
    
    
