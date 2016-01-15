% This goes through the 'sources.set' from each subject and finds components with dipoles in specified x,y,z quadrants.

% enter CompList up through paths

eeglab
 

frback_x = 0;% neg is back of head, pos is front 
clear newcomps
for nx = 1:length(paths)  
    pl = 1;
    EEG = pop_loadset('sources.set' ,['/data/common2/emotion/',paths{nx}]); 
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cp = 1:length(gdcomps{nx})
        dloc = EEG.dipfit.model(gdcomps{nx}(cp)).posxyz(1,:); % x=front/back, y=side/side, z=up/down
        if dloc(1) > frback_x
            newcomps{nx}(pl) = gdcomps{nx}(cp);pl = pl+1;
        end;
    end;
end;
gdcomps= newcomps;
save /data/common2/emotion/FrontalComps.mat gdcomps
% look at dipoles chosen
% find sources for gdcompsa set
 allbesa =  EEG.dipfit.model(1);
for ss = 1:length(gdcomps)
    if ~isempty(gdcomps{ss})        
        EEG = eeg_retrieve(ALLEEG, ss); CURRENTSET = ss;
        dipsources= EEG.dipfit.model(gdcomps{ss}(1));
        for w = 1:length(gdcomps{ss})
            dipsources(1,w)= EEG.dipfit.model(gdcomps{ss}(w));
        end;           
        allbesa(end+1:end+size(dipsources,2)) = dipsources; 
        dipsources=[];
    end; 
end;
allbesa(1) = [];
figure;
dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',20,'dipolelength',.4,'spheres','off','color',{[1 0 .8]},'pointout','on');view(90,0)
figure;
dipplot(allbesa,'image','mri','gui','off','normlen','on','dipolesize',20,'dipolelength',.4,'spheres','off','color',{[0 1 0]},'pointout','on');view(0,90)


ALLEEG=[];EEG=[];
