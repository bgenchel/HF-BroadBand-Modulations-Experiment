% same as ImportEmoRaw.m except works on data without pressure-sensing button
% imports biosemi emotion data saving button info and NO low pass filter
%
%   ImportEmoRawNew(bdf,elp,savepath);
%
% INPUTS:
% bdf -- [string] name of .bdf file to import into a single .set
% elp -- [string] name of .elp file to merge into bdf datasets
% savepath -- [string] full data path where bdf/elp are saved and new dataset will be saved
% savename -- [string] prefix to name resulting .set files with.
% blocks -- [cell array] of 2-element vectors with [start end] frame numbers
%           to import. Each element of the cell array will be saved as a 
%           separate dataset with indexes 1 thru length(blocks). [] will 
%           attempt to import the entire raw datafile into one .set
%           *Hint: 2500 blocks seems to be the upper limit at sccn (12/9/06)
%

function ImportEmoRawNew(bdf,elp,savepath,savename,blocks);


    
    ALLEEG=[]; EEG=[];

    if isempty(blocks)

        EEG = pop_biosig([savepath,bdf], 'blockrange' , [] ,'ref' ,[261 262]);% Right and left mastoids (EXG5 and EXG6), 271 chan total now
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [savepath(end-4:end-1),' Whole Emotion Session']);
        EEG = pop_select( EEG, 'nochannel',[249:256]);% for 256 montage, blank channels
        EEG = pop_select( EEG, 'nochannel',[257:size(EEG.data,1)]);% for 256 montage, blank channels
        
        EEG=pop_chanedit(EEG,  'load',{[savepath,elp] , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
        
        EEG = pop_select( EEG, 'nochannel',[253:254]);% take out ref channels

        for ev = 1:length(EEG.event)
            EEG.event(ev).type = bitand(255,EEG.event(ev).type);
        end;
        
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = pop_eegfilt( EEG, 1, 0, [], [0]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
        delevs = [];
        for ev = 1:length(EEG.event) 
            if EEG.event(ev).type == 30 
                for evv = ev-1:-1:1
                    if strcmp(EEG.event(evv).type,'prebase_instruct')
                        EEG.event(ev).type = 'prebase'; break
                    elseif  strcmp(EEG.event(evv).type,'postbase_instruct')
                        EEG.event(ev).type = 'postbase'; break
                    end;
                end;
            end;       
            if EEG.event(ev).type == 30
                EEG.event(ev).type = 'prebase';             
            end;          
            if EEG.event(ev).type == 1
                EEG.event(ev).type = 'instruct1';
            end;
            if EEG.event(ev).type == 2
                EEG.event(ev).type = 'prebase_instruct';
            end;
            if EEG.event(ev).type == 3
                EEG.event(ev).type = 'instruct2';
            end;
            if EEG.event(ev).type == 5
                EEG.event(ev).type = 'instruct3';
            end;
            if EEG.event(ev).type == 7
                EEG.event(ev).type = 'instruct4';
            end;
            if EEG.event(ev).type == 25
                EEG.event(ev).type = 'postbase_instruct';
            end;
            if EEG.event(ev).type == 6
                EEG.event(ev).type = 'relax';
            end;
            if EEG.event(ev).type == 8
                EEG.event(ev).type = 'enter';
            end;
            if EEG.event(ev).type == 24
                EEG.event(ev).type = 'exit';
            end;
            if EEG.event(ev).type == 9 
                EEG.event(ev).type = 'awe';
            end;
            if EEG.event(ev).type == 10 
                EEG.event(ev).type = 'frustration';
            end;
            if EEG.event(ev).type == 11 
                EEG.event(ev).type = 'joy';
            end;
            if EEG.event(ev).type == 12 
                EEG.event(ev).type = 'anger';
            end;
            if EEG.event(ev).type == 13 
                EEG.event(ev).type = 'happy';
            end;
            if EEG.event(ev).type == 14 
                EEG.event(ev).type = 'sad';
            end;
            if EEG.event(ev).type == 15 
                EEG.event(ev).type = 'love';
            end;
            if EEG.event(ev).type == 16 
                EEG.event(ev).type = 'fear';
            end;
            if EEG.event(ev).type == 17 
                EEG.event(ev).type = 'compassion';
            end;
            if EEG.event(ev).type == 18 
                EEG.event(ev).type = 'jealousy';
            end;
            if EEG.event(ev).type == 19 
                EEG.event(ev).type = 'content';
            end;
            if EEG.event(ev).type == 20 
                EEG.event(ev).type = 'grief';
            end;
            if EEG.event(ev).type == 21 
                EEG.event(ev).type = 'relief';            
            end;
            if EEG.event(ev).type == 23
                EEG.event(ev).type = 'disgust';
            end;
            if EEG.event(ev).type == 22 
                EEG.event(ev).type = 'excite';
            end;
            if EEG.event(ev).type > 30000
                delevs = [delevs ev];
            end;    
        end;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = pop_editeventvals(EEG, 'delete',delevs);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = pop_saveset( EEG, [savename,'.set'], savepath);  

    else  %%%% for more than one chunk of raw data at a time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ALLEEG=[];EEG=[];
        for bl = 1:length(blocks)
            EEG = pop_biosig([savepath,bdf], 'blockrange' , blocks{bl} ,'ref' ,[261 262]);% Right and left mastoids (EXG5 and EXG6)
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [savepath(end-4:end-1),' Emotion part ',int2str(bl),' out of ',int2str(length(blocks))]);
            
            EEG = pop_select( EEG, 'nochannel',[249:256]);% for 256 montage, blank channels
            EEG = pop_select( EEG, 'nochannel',[257:size(EEG.data,1)]);% for 256 montage, blank channels
            
            EEG=pop_chanedit(EEG,  'load',{[savepath,elp] , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            
            EEG = pop_select( EEG, 'nochannel',[253:254]);% take out ref channels
            
            EEG = pop_eegfilt( EEG, 1, 0, [], [0]); % .5 for most, 1.5 for hf45
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
            
            delevs = [];
            for ev = 1:length(EEG.event) 
                if EEG.event(ev).type == 30 
                    for evv = ev-1:-1:1
                        if strcmp(EEG.event(evv).type,'prebase_instruct')
                            EEG.event(ev).type = 'prebase'; break
                        elseif  strcmp(EEG.event(evv).type,'postbase_instruct')
                            EEG.event(ev).type = 'postbase'; break
                        end;
                    end;
                end;       
                if EEG.event(ev).type == 1
                    EEG.event(ev).type = 'instruct1';
                end;
                if EEG.event(ev).type == 2
                    EEG.event(ev).type = 'prebase_instruct';
                end;
                if EEG.event(ev).type == 3
                    EEG.event(ev).type = 'instruct2';
                end;
                if EEG.event(ev).type == 5
                    EEG.event(ev).type = 'instruct3';
                end;
                if EEG.event(ev).type == 7
                    EEG.event(ev).type = 'instruct4';
                end;
                if EEG.event(ev).type == 25
                    EEG.event(ev).type = 'postbase_instruct';
                end;
                if EEG.event(ev).type == 6
                    EEG.event(ev).type = 'relax';
                end;
                if EEG.event(ev).type == 8
                    EEG.event(ev).type = 'enter';
                end;
                if EEG.event(ev).type == 24
                    EEG.event(ev).type = 'exit';
                end;
                if EEG.event(ev).type == 9 
                    EEG.event(ev).type = 'awe';
                end;
                if EEG.event(ev).type == 10 
                    EEG.event(ev).type = 'frustration';
                end;
                if EEG.event(ev).type == 11 
                    EEG.event(ev).type = 'joy';
                end;
                if EEG.event(ev).type == 12 
                    EEG.event(ev).type = 'anger';
                end;
                if EEG.event(ev).type == 13 
                    EEG.event(ev).type = 'happy';
                end;
                if EEG.event(ev).type == 14 
                    EEG.event(ev).type = 'sad';
                end;
                if EEG.event(ev).type == 15 
                    EEG.event(ev).type = 'love';
                end;
                if EEG.event(ev).type == 16 
                    EEG.event(ev).type = 'fear';
                end;
                if EEG.event(ev).type == 17 
                    EEG.event(ev).type = 'compassion';
                end;
                if EEG.event(ev).type == 18 
                    EEG.event(ev).type = 'jealousy';
                end;
                if EEG.event(ev).type == 19 
                    EEG.event(ev).type = 'content';
                end;
                if EEG.event(ev).type == 20 
                    EEG.event(ev).type = 'grief';
                end;
                if EEG.event(ev).type == 21 
                    EEG.event(ev).type = 'relief';          
                end;
                if EEG.event(ev).type == 23
                    EEG.event(ev).type = 'disgust';
                end;
                if EEG.event(ev).type == 22
                    EEG.event(ev).type = 'excite';
                end;
                if EEG.event(ev).type > 30000
                    delevs = [delevs ev];
                end;    
            end;
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            EEG = pop_editeventvals(EEG, 'delete',delevs);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            EEG = pop_saveset( EEG, [savename,'-',int2str(bl),'.set'], savepath);  % just in case of crash
            ALLEEG=[];EEG=[];
        end;       
    end;

