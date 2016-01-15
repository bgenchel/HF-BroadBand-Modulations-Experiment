% displays certain comps for all emotions, corrected for button press speed
% use PlotPresses.m to get avg button press trajectories
eeglab
subj1 = [6,12,18,22,23,24];
sph=floatread('/data/common1/emotion/jo74/sph252-160.sph',[252 252]); 
wts=floatread('/data/common1/emotion/jo74/wts252-160.wts',[160 252]); 
load /data/common1/emotion/jo74/ersps/EpochedErspLateEyes.mat
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
emotimes = [2000,800,750,650,2000,1000,800,1600,750,2000,900,1400,1600,2000];
EEG = pop_loadset( 'sources.set', '/data/common1/emotion/jo74/');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);  
EEG.icaweights=wts;
EEG.icasphere=sph;EEG.icawinv=[];
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fr = find(freqs< 20);
emoset = [2,3,4,5,6,10,11,12,14];
for kk = 1:length(subj1)
    figure; row = 4;col = 5;
    subplot(row,col,6);
    dipplot(EEG.dipfit.model(subj1(kk)),'image','mri','gui','off','normlen','on','dipolesize',25,'pointout','on','dipolelength',0,'projimg','on','spheres','on','color',{'r'},'projcol',[1 0 .2],'projlines','on'); 
    view(55,20);
    subplot(row,col,1);
    topoplot(EEG.icawinv(:,subj1(kk)),EEG.chanlocs,'electrodes','off','shrink','off','plotrad',.5);hold on;
    camzoom(1.5);
    title(int2str(subj1(kk)));
    pl = 2;
    for n = 1:length(emoset)
        subplot(row,col,pl)
        oneersp = erspcell{emoset(n)};
        plotersp = oneersp(:,:,subj1(kk));
        oneboot = ebootcell{emoset(n)};
        minmask = oneboot(1,:,:);
        maxmask = oneboot(2,:,:);
        minmask = repmat(minmask,[200 1 1]);
        maxmask = repmat(maxmask,[200 1 1]);
        onemin = minmask(:,:,subj1(kk))';
        onemax = maxmask(:,:,subj1(kk))';
        %plotersp(find(plotersp>onemin & plotersp < onemax))=0;
        imagesc(times(tm),freqs(fr),plotersp,[-lim lim]); hold on;
        if n == length(emoset)
            colorbar;
        end;        
        plot([0 0],[50 0],'k-');
        ph =  title(emos{emoset(n)});
        set(ph,'fontsize',14);
        set(gca,'fontsize',14);
        set(gca,'ydir','norm');
        set(gca,'xlim',[-400 emotimes(emoset(n))]);
         set(gca,'xtick',[-400:200:2000]); 
        set(gca,'xticklabel',{[] [] 0 [] [] [] .8 [] [] [] 1.6 [] []});
        subplot(row,col,pl+5)
        yl = max(allavg(emoset(n),:)); yl = yl+250000;
        ph=plot(EEG.times,allavg(emoset(n),:),'k');hold on;
        set(ph,'linewidth',3);
        set(gca,'xlim',[-400 emotimes(emoset(n))]);
        set(gca,'ylim',[-10000 2700000]);
        plot([0 0],[get(gca,'ylim')],'k-');
        set(gca,'xtick',[-400:200:2000]); 
        set(gca,'xticklabel',{-.4 [] 0 [] .4 [] .8 [] 1.2 [] 1.6 [] 2});
        set(gca,'yticklabel',[]);
        set(gca,'box','off');
        ph =  title(emos{emoset(n)});
        set(ph,'fontsize',14);        
        pl = pl+1;
        if n == 4 
            pl = pl+5;
        end; 
    end;
set(gcf,'color','w');
end;

        set(gca,'xticklabel',{[] [] 0 [] [] [] .8 [] [] [] 1.6 [] []});
