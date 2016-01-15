% imports biosemi emotion data saving button info and NO low pass filter
% This is specifically for emotion subjects 1:35 who used the pressure-
% sensing button to signal emotion onset
%
% ImportEmoRaw(bdf,elp,savepath,savename,saveto,blocks,getpress);
%
% INPUTS:
% bdf -- [string] name of .bdf file to import into a single .set
% elp -- [string] name of .elp file to merge into bdf datasets
% savepath -- [string] full data path where bdf/elp are saved and new dataset will be saved
% savename -- [string] prefix to name resulting .set files with.
% saveto -- [string] file directory to save .set files to
% blocks -- [cell array] of 2-element vectors with [start end] frame numbers
%           to import. Each element of the cell array will be saved as a 
%           separate dataset with indexes 1 thru length(blocks). [] will 
%           attempt to import the entire raw datafile into one .set
%           *Hint: 2500 blocks seems to be the upper limit at sccn (12/9/06)
% getpress -- [0|1] if 1, will evaluate EEG.button for emotion button presses
%              and will save as events in EEG.event

function [EEG]=ImportEmoRaw(bdf,elp,savepath,savename,saveto,blocks,getpress);

    
    ALLEEG=[]; EEG=[];
    refs = [131 215]; % subj 3 can only be 215 (E3 is bad)
    if strcmp(bdf,'ms82.bdf') % fix subj 1
      refs = [215]; % subj 3 can only be 215 (E3 is bad)
    end;
    if isempty(blocks)
      x=openbdf([savepath,bdf]);
      %EEG = pop_readbdf([savepath,bdf], [] ,257,refs );% subjs 9,35
      EEG = pop_readbdf([savepath,bdf], [] ,265,refs );
        %EEG = pop_biosig([savepath,bdf],'channels' ,[1:x.Head.NS], 'blockrange' , [] ,'ref' ,[131 215]);% E3 and G23
        %EEG = pop_biosig([savepath,bdf],'channels' ,[1:x.Head.NS], 'blockrange' , [] ,'ref' ,[215]);% G23 only for nx=33 (jl83) 
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'setname', [savepath(end-4:end-1),' Whole Emotion Session; No filters']);
        EEG.button = EEG.data(257,:);
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = pop_select( EEG, 'nochannel',[257:264] );
        EEG = pop_select( EEG, 'nochannel',[255:256] );%heart
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
        EEG=pop_chanedit(EEG,  'load',{[savepath,elp] , 'filetype', 'autodetect'}, 'delete',1, 'delete',1, 'delete',1, 'eval', 'chans = pop_chancenter( chans, [],[]);', 'forcelocs',{0, 'X', 'B12'});
        fprintf('\nRemoving reference channels...\n')
        EEG = pop_select( EEG, 'nochannel',refs ); 
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        %EEG = pop_eegfilt( EEG, 1, 0, [], [0]); % reimport 1/25/2010
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
        delevs = [];
        first = 1; % finds first '30' event for prebase and then the second is post base  
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
            if EEG.event(ev).type == 23 % yes, this is the way the tracks were numbered
                EEG.event(ev).type = 'disgust';
            end;
            if EEG.event(ev).type == 22 
                EEG.event(ev).type = 'excite';
            end;
            if EEG.event(ev).type > 30000
                delevs = [delevs ev];
            end;    
        end;
        if strcmp(bdf,'tl81.bdf') % fix subj 1 
          EEG.event(46).type = 'fear';
          EEG.event(56).type = 'jealousy';
          EEG.event(66).type = 'grief';
        elseif strcmp(bdf,'mi83.bdf') % fix subj 2
          EEG.event(129).type = 'fear';
          EEG.event(149).type = 'jealousy';
          EEG.event(177).type = 'grief';            
        elseif strcmp(bdf,'ms82.bdf') % fix subj 3
          EEG.event(47).type = 'fear';
          EEG.event(54).type = 'jealousy';
          EEG.event(68).type = 'grief';
        end;
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = pop_editeventvals(EEG, 'delete',delevs);
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
        %EEG = pop_saveset( EEG, [savename,'.set'], saveto,'savemode','twofiles');  
        if getpress == 1
            [EEG,searchnext] = GetPresses(EEG);
            %EEG = pop_saveset( EEG, [savename,'.set'], saveto,'savemode','twofiles');  
        end;

    else %%%% for more than one chunk of raw data at a time %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for bl = 1:length(blocks)
            ALLEEG=[];EEG=[];
            x=openbdf([savepath,bdf]);
            EEG = pop_readbdf([savepath,bdf], blocks{bl} ,265,[131 215] );
            %EEG = pop_biosig([savepath,bdf],'channels' ,[1:x.Head.NS], 'blockrange' , blocks{bl} ,'ref' ,[131 215]);% E3 and G23
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,  'setname', [savepath(end-4:end-1),' Emotion part ',int2str(bl),' out of ',int2str(length(blocks))]);
            clear x
            if size(EEG.data,1) > 256
                EEG.button = EEG.data(257,:);
                EEG = pop_select( EEG, 'nochannel',[257:264] ); 
            end;
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
            EEG = pop_select( EEG, 'nochannel',[255:256] ); % take out heart because 'types' is ridiculous
            EEG=pop_chanedit(EEG,  'load',{[savepath,elp] , 'filetype', 'autodetect'}, 'changefield',{255, 'type', 'ECG'}, 'changefield',{256, 'type', 'ECG'}, 'forcelocs',{0, 'X', 'B12'}, 'eval', 'chantmp = pop_chancenter( chantmp, [],[255 256] );', 'forcelocs',{0, 'X', 'B12'});
            EEG = pop_select( EEG, 'nochannel',refs ); 
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            EEG = pop_eegfilt( EEG, 1, 0, [], [0]); % .5 may be too low: global offset IC?
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,  'overwrite', 'on');
            
            delevs = [];
            for ev = 1:length(EEG.event) 
                if EEG.event(ev).type == 30 
                    for evv = ev-1:-1:1
                        if strcmp(EEG.event(evv).type, 'prebase_instruct')
                            EEG.event(ev).type = 'prebase'; break
                        elseif  strcmp(EEG.event(evv).type, 'postbase_instruct')
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
            %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            EEG = pop_editeventvals(EEG, 'delete',delevs);
            %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
            EEG = pop_saveset( EEG, [savename,'-',int2str(bl),'.set'], saveto,'savemode','twofiles');  % just in case of crash
            if getpress == 1
                [EEG,searchnext] = GetPresses(EEG);
                EEG = pop_saveset( EEG, [savename,'-',int2str(bl),'.set'], saveto,'savemode','twofiles');  
            end;        
            ALLEEG=[];EEG=[];
        end;       
    end;

