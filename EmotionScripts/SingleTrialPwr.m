% Finds power in single trials during Emotion expt to get error bars on power spectra

eeglab

emoset = {'awe.set', 'frustration.set','joy.set','anger.set','sad.set','happy.set','fear.set','love.set','jealousy.set','compassion.set','emabarrass.set','content.set','grief.set','relief.set'};

subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % jo74 emot
gdcomps = {subj1 };

paths = {'/data/common1/emotion/ap80/','/data/common1/emotion/rr80/','/data/common1/emotion/jo74/'}
%%%%%%%%%%%%%%  Component Spectra   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sph=floatread('sph252-160.sph',[252 252]); 
wts=floatread('wts252-160.wts',[160 252]); 

for nx = 1:length(emos)
    EEG = pop_loadset( datset{nx},paths{nx});
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.icaweights=wts;
    EEG.icasphere=sph;EEG.icawinv=[];
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    for cmp = 1:length(subj1)
        for tr = 1:size(EEG.icaact,3)
            [Pxx,F] = psd(EEG.icaact(subj1(cmp),:,tr),256,EEG.srate,256);% much smoother
            cspec(subj1(cmp),:,tr) = log10(Pxx(:,1))';
        end;
        
    end;
    allemopwr{1,nx} = cspec;
    ALLEEG=[];EEG=[];
end;
comment = 'jo74 psd power spectra in good comps (5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46) during each emotion (the different cell arrays) in this order: awe, frustration,joy,anger,sad,happy,fear,love ,jealousy,compassion,emabarrass,content,grief,relief; starts at 10 sec to end of emotin(space bar)';
freqs = F;
save /data/common1/emotion/jo74/sngltrpwr.mat allemopwr freqs comment
save /home/julie/EmoTF/sngltrpwr.mat allemopwr freqs comment

           
           
% visualize results with error bars
load /data/common1/emotion/jo74/sngltrpwr.mat
           
color = {'k' 'b' 'c' 'g' 'm' 'r' 'k' 'b' 'c' 'g' 'm' 'r' 'k' 'b'};
figure;
%for c = 1:size(oneemo,1)
    for nx=1:length(emos)
    oneemo = allemopwr{nx};
        for f = 1:size(oneemo,2)
            sdemo(subj1(c),f) = std(oneemo(subj1(c),f,:))/sqrt(size(oneemo,3));
        end;
errorbar(freqs,mean(oneemo(subj1(c),:,:),3),-sdemo(subj1(c),:),sdemo(subj1(c),:),color{nx});hold on;
    end;
set(gca,'xlim',[2 30]);
    %end;

          
% plot error manually 
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief','rest'};
%plotemos = [1:10,12:15];
plotemos = [1,3,6,8,10,12,14,15];
plotemos = [2,4,5,7,9,13,15];
posset = [1,3,6,8,10,12,14];
negset = [2,4,5,7,9,13];
colneg = winter(6);
colpos = autumn(7);
figure;n=1; p=1;  pl = 1; %(for error bar spacing)
for c = 1:length(subj1)
    n=1; p=1;
    subplot(4,4,c)
    for nx=1:length(plotemos)
        nx = plotemos(nx);
        oneemo = allemopwr{nx};
        for f = 1:size(oneemo,2)
            sdemo(subj1(c),f) = std(oneemo(subj1(c),f,:))/sqrt(size(oneemo,3));
        end;
        ph=plot(freqs,mean(oneemo(subj1(c),:,:),3));hold on;
        if ismember(nx,negset)
            set(ph,'color',colneg(n,:));
            currcol = colneg(n,:);n=n+1;
        elseif ismember(nx,posset)
            set(ph,'color',colpos(p,:));
            currcol = colpos(p,:);p=p+1;
        else
            set(ph,'color','k');
            currcol = 'k';%p=p+1;            
        end;
        
        set(ph,'linewidth',2);    
        for ff = 1:45
            xf = ff*pl;
            ph=plot([freqs(xf)  freqs(xf)],[mean(oneemo(subj1(c),xf,:),3)-sdemo(subj1(c),xf) mean(oneemo(subj1(c),xf,:),3)+sdemo(subj1(c),xf)]);
            set(ph,'color',currcol);
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(c),xf,:),3)+sdemo(subj1(c),xf) mean(oneemo(subj1(c),xf,:),3)+sdemo(subj1(c),xf)]);
            set(ph,'color',currcol); 
            set(ph,'linewidth',1);    
             ph = plot([freqs(xf)-.1  freqs(xf)+.1],[mean(oneemo(subj1(c),xf,:),3)-sdemo(subj1(c),xf) mean(oneemo(subj1(c),xf,:),3)-sdemo(subj1(c),xf)]);
            set(ph,'color',currcol); 
            set(ph,'linewidth',1);    
        end;
    end;
    set(gca,'xlim',[2 45]);
    set(gca,'fontsize',16);
    set(gca,'box','off');
    title(int2str(subj1(c)));
end;
axcopy
figure;
for c = 1:length(posset)
ph=barh(c,1);hold on;
set(ph, 'FaceColor',colpos(c,1:3));
end;
set(gca,'xticklabel',[]);
set(gca,'ytick',[1:7]);
set(gca,'yticklabel', emos(posset));
figure;
for c = 1:length(negset)
ph=barh(c,1);hold on;
set(ph, 'FaceColor',colneg(c,1:3));
end;
set(gca,'xticklabel',[]);
set(gca,'ytick',[1:6]);
set(gca,'yticklabel', emos(negset));

% plot power as 2D image instead 
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46];
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief','rest'};
%plotemos = [15,1,3,6,8,10,12,14,2,4,5,7,9,13];
plotemos = [1,3,6,8,10,12,14,15];
plotemos = [2,4,5,7,9,13,15];
posset = [1,3,6,8,10,12,14];
negset = [2,4,5,7,9,13];
colneg = winter(6);
colpos = autumn(7);
fr = find(freqs < 50);

for c = 1:length(subj1)
figure; pl=1;clear plotemo
    for nx=1:length(plotemos)
        nx = plotemos(nx);
        oneemo = allemopwr{nx};
        plotemo(:,pl) = mean(oneemo(subj1(c),fr,:),3)';pl = pl+1;
    end;    
        ph=imagesc([1:length(plotemos)],freqs(fr),plotemo,[-2 2]);hold on;
        set(gca,'ydir','norm');
        set(gca,'xtick',[1:length(plotemos)]);
        set(gca,'xticklabel', emos(plotemos));
        title(int2str(subj1(c)));
end;
axcopy
% plot area under the curve (sum) of power in one comp vs another

figure;
comp1 = 6; 
comp2 = 7;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,1)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 7; 
comp2 = 23;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,2)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 10; 
comp2 = 12;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,3)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 13; 
comp2 = 19;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,4)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 24; 
comp2 = 25;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,5)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 11; 
comp2 = 23;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,6)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 6; 
comp2 = 18;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,7)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 16; 
comp2 = 18;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,8)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);
comp1 = 16; 
comp2 = 28;
fr = find(freqs > 2 & freqs<30);
subplot(3,3,9)
for nx=1:length(emos)
    oneemo = allemopwr{nx};
    c1 = oneemo(comp1,fr,:);
    c1 = mean(c1,3);
    c2 = oneemo(comp2,fr,:);
    c2 = mean(c2,3);
    c1 = sum(c1);
    c2 = sum(c2);
    ph = plot(c1,c2,'b*'); hold on;
    set(ph,'markersize',5);
    text(c1+.1,c2,emos{nx});    
end;        
title(['Comp ',int2str(comp1),' vs ','Comp ',int2str(comp2)]);

subj1 = [5,6,7,10,11,12,13,15,16,18,19,22,23,24,25,28,36,40,46];  % jo74 emot
figure;
for cmp = 1:length(subj1)
       psd(EEG.icaact(subj1(cmp),:),256,EEG.srate,256);hold on;% much smoother
end;
set(gca,'xlim',[0 50]);