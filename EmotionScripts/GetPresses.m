% to be used from the ImportEmoRaw() function to find pressure-sensing
% button presses and log them as EEG events
% if no press events can be detected for a given emotion, then one initial 
% event is created halfway through the entire emotion period
% the first event ('emotion onset') is type 'press1'
% all other presses ('renewed or continued emotion') is type 'presses'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get pressevents (from EEG.button) for each emo (if they exist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS:
% EEG -- EEGLAB data structure with new events included where button presses detected
% searchnext -- [0 | 1] If 1, then an emotion period was cut off before button detected.



function [EEG,searchnext] = GetPresses(EEG,searchnext);


    emos = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','excite','disgust'};
    emobuts = {'bawe', 'bfrustration','bjoy','banger','bhappy','bsad','blove' ,'bfear','bcompassion','bjealousy','bcontent','bgrief','brelief','bexcite','bdisgust'};
    ALLEEG=[];
      
    if ~exist('searchnext')
      searchnext = 0; % if not defined, indicate that no emotion was cut off before button press 
    elseif isempty(searchnext)
      searchnext = 0; % if not defined, indicate that no emotion was cut off before button press 
    end;
    
    peaksurround = 200;
    chunksize = 7000;
    if searchnext == 1
        fprintf('\nSearching beginning of dataset for leftover button press...\n');
        estart = 1;
        estop = EEG.event(1).latency;
        searchbut = EEG.button(estart:estop);
        searchbut = searchbut - searchbut(1,1);
        
        mindata = min(searchbut);
        maxdata = max(searchbut);
        thresh = (maxdata-mindata)*.25;  %15000

        if thresh-mindata < 1000
            thresh = 90000; % impossibly large because no valid presses
        end;
        X = searchbut; 
        mxma = [];pp=1; 
        for chunk = 1:4000:length(X)-chunksize
            newX = X(1,chunk:chunk+chunksize);
            y = newX-median(newX);
            for n = peaksurround+1:length(y)-peaksurround
                if y(n) > y(n-[1:peaksurround])& y(n) > y(n+[1:peaksurround])&y(n)>thresh
                    mxma(1,pp) = estart + n+(chunk-1); pp = pp+1; 
                end;
            end;
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(mxma) > 1            
            if length(unique(mxma)) < length(mxma)             
                newmxma=zeros(1,0);
                for plc = 1:length(mxma)
                    if length(find(mxma == mxma(plc))) > 1
                        newmxma(1,end+1) = mxma(plc);
                    end;
                end;
                newmxma = unique(newmxma);                
                mxma = newmxma;
            end;               
        end; 
        if ~isempty(mxma)
            EEG.event(end+1).type = 'press1';
            EEG.event(end).latency = mxma(1,1);
            fprintf('\nFound at least one leftover button press...\n');
        else
            fprintf('\nNo leftover button press found...\n');
        end;        
        if length(mxma) > 1
            for newev = 2:length(mxma)
                EEG.event(end+1).type = 'press';
                EEG.event(end).latency = mxma(1,newev);
            end;  
        end; 
    end;
            
    searchnext = 0; % reset to zero after taking care of first button press
    
    for e = 1:length(emos)
        estart = [];
        estop = [];
        relevent = find(strcmp({EEG.event.type},emos{e}));
        if ~isempty(relevent)
            estart = EEG.event(relevent).latency;
            for ev = relevent:length(EEG.event)
                if strcmp(EEG.event(ev).type,'exit')
                    estop = EEG.event(ev).latency; break
                end;
            end; 
        else
            estart = [];
            estop = [];
        end;
        mxma = [];
        if ~isempty(estart) & ~isempty(estop)
            searchbut = EEG.button(estart:estop);
            searchbut = searchbut - searchbut(1,1);
            
            mindata = min(searchbut);
            maxdata = max(searchbut);
            thresh = (maxdata-mindata)*.25;  %15000
            if thresh-mindata < 1000
                thresh = 90000; % impossibly large because no valid presses
            end;
            X = searchbut; 
            mxma = [];pp=1; 
            for chunk = 1:4000:length(X)-chunksize
                newX = X(1,chunk:chunk+chunksize);
                y = newX-median(newX);
                for n = peaksurround+1:length(y)-peaksurround
                    if y(n) > y(n-[1:peaksurround])& y(n) > y(n+[1:peaksurround])&y(n)>thresh
                        mxma(1,pp) = estart + n+(chunk-1); pp = pp+1; 
                    end;
                end;
            end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(mxma) > 1            
                if length(unique(mxma)) < length(mxma)             
                    newmxma=zeros(1,0);
                    for plc = 1:length(mxma)
                        if length(find(mxma == mxma(plc))) > 1
                            newmxma(1,end+1) = mxma(plc);
                        end;
                    end;
                    newmxma = unique(newmxma);                
                    mxma = newmxma;
                end;               
            end; 
            if ~isempty(mxma)
                EEG.event(end+1).type = 'press1';
                EEG.event(end).latency = mxma(1,1);
                fprintf('\nInitial button press detected for emotion: %s\n',emos{e});    
            else
                fprintf('\nNo button presses detected for emotion: %s\n',emos{e});
                tottime = estop - estart; % total time
                EEG.event(end+1).type = 'press1'; 
                EEG.event(end).latency = estart+ tottime/3; % 1/3 way in
                fprintf('\nInserted press1 1/3 of the way through emotion',emos{e});
                
            end;
            
            if length(mxma) > 1
                for newev = 2:length(mxma)
                    EEG.event(end+1).type = 'press';
                    EEG.event(end).latency = mxma(1,newev);
                end;
            else    
                fprintf('\nNo additional button presses detected for emotion: %s\n',emos{e});    
            end; 
        elseif ~isempty(estart) & isempty(estop)
                fprintf('\nEmotion %s cut off: search beginning of next dataset!\n',emos{e});    
            searchnext = 1;            
        end;
    end;
    EEG = eeg_checkset(EEG, 'eventconsistency');        
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
    
