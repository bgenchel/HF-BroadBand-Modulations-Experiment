%emos = {'awepress.set' 'awe2press.set' 'frustpress.set' 'frust2press.set' 'joypress.set' 'joy2press.set' 'angerpress.set' 'anger2press.set' 'sadpress.set' 'sad2press.set' 'surprisepress.set' 'surprise2press.set' 'happypress.set' 'happy2press.set' 'fearpress.set' 'fear2press.set' 'lovepress.set' 'love2press.set' 'jealousypress.set' 'jealousy2press.set' 'compassionpress.set' 'compassion2press.set' 'contentpress.set' 'content2press.set' 'griefpress.set' 'grief2press.set' 'reliefpress.set' 'relief2press.set'};
emos = {'awe2press.set' 'frust2press.set'  'joy2press.set' 'anger2press.set'  'sad2press.set' 'surprise2press.set' 'happy2press.set'  'fear2press.set' 'love2press.set' 'jealousy2press.set' 'compassion2press.set' 'content2press.set' 'grief2press.set' 'relief2press.set'};
paths = {'/data/common1/emotion/ap80/ersps/awe/','/data/common1/emotion/ap80/ersps/frust/','/data/common1/emotion/ap80/ersps/joy/','/data/common1/emotion/ap80/ersps/anger/','/data/common1/emotion/ap80/ersps/sad/','/data/common1/emotion/ap80/ersps/surprise/','/data/common1/emotion/ap80/ersps/happy/','/data/common1/emotion/ap80/ersps/fear/','/data/common1/emotion/ap80/ersps/love/','/data/common1/emotion/ap80/ersps/jealousy/','/data/common1/emotion/ap80/ersps/compassion/','/data/common1/emotion/ap80/ersps/content/','/data/common1/emotion/ap80/ersps/grief/','/data/common1/emotion/ap80/ersps/relief/'};
cd /data/common1/emotion/ap80/
load wts
load sph
load aplocs
a = [7:12,14,17:20,22,24,26];
for e = 1:length(emos)
    EEG = pop_loadset( emos{e},'/data/common1/emotion/ap80/');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_select( EEG, 'nochannel',72);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'overwrite', 'on');
    EEG.chanlocs = chan;    
    % now recalculate ica.act
    EEG.icaweights = wts;
    EEG.icasphere = sph;
    EEG.icaact = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
    coher=zeros(65,200,26);
    crsboot=zeros(65,2,26);
    crsangle=zeros(65,200,26);
    cohercell = cell (1,26);
    crsbootcell = cell (1,26);
    crsanglecell = cell (1,26);
    %-- always list lowest to highest comp numbers-------------
    %done =
    figure;
    for index1=1:length(a)-1
        for index2=index1+1: length(a)
            [coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher, alltfX, alltfY] = crossf(EEG.icaact(a(index1),:),EEG.icaact(a(index2),:), EEG.pnts, [-1000 2000], EEG.srate, [3 .5], 'alpha', 0.02,'padratio', 4, 'topovec', EEG.icawinv(:, [a(index1),a(index2)])', 'elocs', EEG.chanlocs,'winsize',256);
            coher(:,:,a(index2)) = coh;
            cohercell{1,a(index1)} = coher;
            crsboot(:,:,a(index2)) = cohboot;
            crsbootcell{1,a(index1)} = crsboot;
            crsangle(:,:,a(index2)) = cohangles;
            crsanglecell{1,a(index1)}= crsangle;
            clf
        end
    end;
    comment = 'ap80 Emotion pilot. Emotion: second time thru (with-out scenerios). epoched on emotion button presses. sig at .02; padratio=4;winsize=256';
    cd (paths{e})
    save crossall1.mat cohercell crsbootcell crsanglecell timesout freqsout mcoh comment;
    close
ALLEEG = pop_delset( ALLEEG, [1] );
end;



lim = .8;
compsrun = a;
%compsrun = [7:12,14,17:20,22,24,26];
f = find (freqsout<30);
t = find (timesout>-1000 & timesout<2000);
figure;
p = ((length(compsrun)-1)/2)*length(compsrun);
for w=1:length(compsrun)-1
a = compsrun(w);
b = compsrun(w+1:end);
Diff = zeros(26,200,71);
Diff = cohercell{1,a};
bootmat1 = crsbootcell {1,a};   %chooses a 3D array from cell array
minboot = bootmat1(:,2,:);   %chooses upper boot threshold,makes a (63,1,71)
minmask = repmat (minboot,[1,200,1]); %makes (63,200,71) matrix of minboot
Diff(find(Diff <= minmask)) = 0;%finds and zeros cohmats less than minmask
for m=w+1:length(b)+w
subplot(length(compsrun)-1,length(compsrun)-1,((length(compsrun)-1)*(w-1))+(m-1))
imagesc(timesout(t),freqsout(f),Diff(f,t,b(m-w)),[-lim lim]);
hold on
set(gca,'ticklength',[0.05 0.05]); 
set(gca,'tickdir','out');      
set(gca,'xtick',[0,200,400,600,800,1000]); 
set(gca,'xticklabel',[]);    
set(gca,'fontsize',5);    
set(gca,'ytick',[5:5:30]);% 
set(gca,'yticklabel',[]); 
plot([0 0],[0 30],'k-');
title(int2str(b(m-w)))
end;
end;
set(gca,'xticklabel',[0,2,4,6,8,10]);    
set(gca,'yticklabel',[5:5:30]);    
colorbar
%textsc(tl,'title');
axcopy
