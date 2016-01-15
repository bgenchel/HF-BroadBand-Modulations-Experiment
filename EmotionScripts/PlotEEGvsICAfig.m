% plots EEG and ICA acts from same time window to show separation (for fig)

eeglab
paths = {'/ap82/Sternberg/','/cj82/Sternberg/','/ds76/Sternberg/','/ec81/Sternberg/','/jo74/Sternberg/','/ke70/Sternberg/','/km81/Sternberg/','/mk79/Sternberg/','/nf68/Sternberg/','/tp62/Sternberg/','/ds80/Sternberg/','/kb77/Sternberg/','/cz84/Sternberg/','/gm84/Sternberg/','/ts79/Sternberg/','/ny84/Sternberg/','/ft84/Sternberg/','/gv84/Sternberg/','/ka83/Sternberg/','/cy82/Sternberg/','/jb84/Sternberg/','/rd81/Sternberg/','/km81/Sternberg2/','/jo74/Shortst/','/bt78/','/as78/'};

nx = 1;

EEG = pop_loadset( 'stern.set',['/data/common1/stern/eeg/',paths{nx}]);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
start = 100; finish = 3000;

row = 40; col = 1; start = 24000; finish = 26000;
figure; pl=1;
for d = 41:71
subplot(row,col,pl)
ph = plot(EEG.data(d,start:finish),'b');pl = pl+1;
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
axis('off');
end;
for d = 1:26
subplot(row,col,pl)
ph = plot(EEG.icaact(d,start:finish),'r');pl = pl+1;
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
axis('off');
end;

chans = [2,5,10,26,30,33,37,39,44,49,59,65,70]; % pick specifics for demonstration
comps = [2,3,8,9,10,11,13,31,33];
row = length(comps)+length(chans); col = 1;
figure; pl=1;
for d = 1:length(chans)
subplot(row,col,pl)
ph = plot(EEG.data(chans(d),start:finish),'b');pl = pl+1;
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'xlim',[1 finish-start+1]);
axis('off');
end;
for d = 1:length(comps)
subplot(row,col,pl)
ph = plot(EEG.icaact(comps(d),start:finish),'r');pl = pl+1;
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'xlim',[1 finish-start+1]);
axis('off');
end;set(gcf,'color','w');
 x={EEG.chanlocs(chans).labels};
  Columns 1 through 8

    'REYE'    'FP1'    'AFZ'    'FC5'    'FC2'    'FT8'    'C3'    'CZ'

  Columns 9 through 13

    'TP9'    'CPZ'    'PZ'    'PO7'    'PO10'

pop_topoplot(EEG,0, comps , '',[length(comps) 1] ,0, 'electrodes', 'off', 'plotrad',0.5, 'masksurf', 'on');
