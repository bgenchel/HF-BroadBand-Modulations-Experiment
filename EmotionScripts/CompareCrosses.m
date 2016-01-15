% This will line up crosses between comps across emotion conditions
name = {'awe' 'frust' 'joy' 'anger' 'sad' 'surprise' 'happy' 'fear' 'love' 'jealousy'  'compassion' 'content' 'grief' 'relief' 'awe2' 'frust2' 'joy2' 'anger2' 'sad2' 'surprise2' 'happy2' 'fear2' 'love2' 'jealousy2' 'compassion2' 'content2' 'grief2' 'relief2' };%'content2' 
load crossall1.mat
name = {'awe' 'frust' 'joy' 'anger' 'sad' 'surprise' 'happy' 'fear' 'love' 'jealousy'  'compassion' 'content' 'grief' 'relief' 'joy','frust','awe','sad','compassion','surprise','fear','happy','jealousy','love','anger','grief','content','relief' };%'content2' 

paths = {'/data/common1/emotion/ap80/ersps/awe/','/data/common1/emotion/ap80/ersps/frust/','/data/common1/emotion/ap80/ersps/joy/','/data/common1/emotion/ap80/ersps/anger/','/data/common1/emotion/ap80/ersps/sad/','/data/common1/emotion/ap80/ersps/surprise/','/data/common1/emotion/ap80/ersps/happy/','/data/common1/emotion/ap80/ersps/fear/','/data/common1/emotion/ap80/ersps/love/','/data/common1/emotion/ap80/ersps/jealousy/','/data/common1/emotion/ap80/ersps/compassion/','/data/common1/emotion/ap80/ersps/content/','/data/common1/emotion/ap80/ersps/grief/','/data/common1/emotion/ap80/ersps/relief/','/data/common1/emotion/ap80/ersps/awe/','/data/common1/emotion/ap80/ersps/frust/','/data/common1/emotion/ap80/ersps/joy/','/data/common1/emotion/ap80/ersps/anger/','/data/common1/emotion/ap80/ersps/sad/','/data/common1/emotion/ap80/ersps/surprise/','/data/common1/emotion/ap80/ersps/happy/','/data/common1/emotion/ap80/ersps/fear/','/data/common1/emotion/ap80/ersps/love/','/data/common1/emotion/ap80/ersps/jealousy/','/data/common1/emotion/ap80/ersps/compassion/','/data/common1/emotion/ap80/ersps/content/','/data/common1/emotion/ap80/ersps/grief/','/data/common1/emotion/ap80/ersps/relief/'};
paths = {'/data/common1/emotion/ap80/ersps/awe/','/data/common1/emotion/ap80/ersps/frust/','/data/common1/emotion/ap80/ersps/joy/','/data/common1/emotion/ap80/ersps/anger/','/data/common1/emotion/ap80/ersps/sad/','/data/common1/emotion/ap80/ersps/surprise/','/data/common1/emotion/ap80/ersps/happy/','/data/common1/emotion/ap80/ersps/fear/','/data/common1/emotion/ap80/ersps/love/','/data/common1/emotion/ap80/ersps/jealousy/','/data/common1/emotion/ap80/ersps/compassion/','/data/common1/emotion/ap80/ersps/content/','/data/common1/emotion/ap80/ersps/grief/','/data/common1/emotion/ap80/ersps/relief/','/data/common1/emotion/ap80/ersps/joy/','/data/common1/emotion/ap80/ersps/frust/','/data/common1/emotion/ap80/ersps/awe/','/data/common1/emotion/ap80/ersps/sad/','/data/common1/emotion/ap80/ersps/compassion/','/data/common1/emotion/ap80/ersps/surprise/','/data/common1/emotion/ap80/ersps/fear/','/data/common1/emotion/ap80/ersps/happy/','/data/common1/emotion/ap80/ersps/jealousy/','/data/common1/emotion/ap80/ersps/love/','/data/common1/emotion/ap80/ersps/anger/','/data/common1/emotion/ap80/ersps/grief/','/data/common1/emotion/ap80/ersps/content/','/data/common1/emotion/ap80/ersps/relief/'};
% all comps before
comps = [7:12,14,17:20,22,24,26]; % ap80
for wh = 4; % which component pair to see crosses for
    wh2 = 19; % 
    tl = int2str(comps(wh));
    fr = find (freqsout<30);
    tm = find (timesout>-1000 & timesout<2000);
    lim = 1;
    figure; pl = 1;
    for p = 15:length(paths)
        cd (paths{p})
        if p < 15
            load crossall1.mat
        else
            load crossall2.mat
        end;
        subplot(4,4,pl)
        onecross = cohercell{1,comps(wh)};
        bootmat1 = crsbootcell {1,comps(wh)};   %chooses a 3D array from cell array
        minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
        minmask = repmat (minboot,[1,200,1]); %makes (63,200,71) matrix of minboot
        onecross(find(onecross <= minmask)) = 0;%finds and zeros cohmats less than minmask
        imagesc(timesout(tm),freqsout(fr),onecross(fr,tm,wh2),[-lim lim]);
        title(name{p});
        pl=pl+1;
        if p == 14
            pl = 2; end;
    end;
    textsc(tl,'title');
    axcopy
end;
colorbar
% all comps after   
comps = [7:12,14,17:20,22,24,26]; % ap80
 wh = 19; % which component pair to see crosses for
 for  wh2 = 11:14; % 
    tl = int2str(comps(wh2));
    fr = find (freqsout<30);
    tm = find (timesout>-1000 & timesout<2000);
    lim = .8;
    figure; pl = 1;
    for p = 1:length(paths)*2
        cd (paths{p})
        if p < 15
            load crossall1.mat
        else
            load crossall2.mat
        end;
        subplot(5,6,pl)
        onecross = cohercell{1,comps(wh)};
        bootmat1 = crsbootcell {1,comps(wh)};   %chooses a 3D array from cell array
        minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
        minmask = repmat (minboot,[1,200,1]); %makes (63,200,71) matrix of minboot
        onecross(find(onecross <= minmask)) = 0;%finds and zeros cohmats less than minmask
        imagesc(timesout(tm),freqsout(fr),onecross(fr,tm,comps(wh2)),[-lim lim]);
        title(name{p});
        pl=pl+2;
        if p == 14
            pl = 2; end;
    end;
    textsc(tl,'title');
    axcopy
end;

    




