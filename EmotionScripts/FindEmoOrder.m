% find actual order that the emotions were presented for all subj
%
%  [emoorders] = FindEmoOrder(fullpaths,emos)
%
%
%
%
%
%
%

function [emoorders] = FindEmoOrder(fullpaths,emos)


for nx = 1:length(fullpaths)
    if ~isempty(fullpaths{nx})
    EEG = pop_loadset( 'ButtonOnly.set',fullpaths{nx},'all');
    pl = 1;
    for ev = 1:length(EEG.event)
        if ismember(EEG.event(ev).type,emos)
            evord{pl} = emos{strcmp(emos,EEG.event(ev).type)};
            pl = pl+1;
        end;
    end;
    emoorders{nx} = evord;
    end;
end;
