% plot, as an interpolated scalp map, the max pvaf at each chan or pvaf of good comps at each chan

eeglab
load /data/common4/RewTwoback/GoodComps.mat 
datset = 'AllFBepochs.set';
nx=1;
artcomps{1} = [1,2,3,9,10,11,12,14,16,17,21,22,23,24];

EEG = pop_loadset(datset, fullpaths{nx});
[pvaf,pvafswhole,vars] = eeg_pvaf(EEG,gdcomps{nx}); % whole data pvaf
[pvaf,pvafsnoart,vars] = eeg_pvaf(EEG,gdcomps{nx}, 'artcomps',artcomps{nx});% - artifacts 
[pvaf,pvafsallart,vars] = eeg_pvaf(EEG,artcomps{nx});% - artifacts 
allpvafs{1} = pvafswhole;
allpvafs{2} = pvafsnoart;
allpvafs{3} = pvafsallart;
ttls{1} = 'Dipolar Components | Whole Data';
ttls{2} = 'Dipolar Components | Non-Artifact Data';
ttls{3} = 'Artifact Components | Whole Data';
,'cbar',0
row = 1; col = 3;
clear M2
for mf = 1:180
    figure;ads = [0:2:360]; mvs = [[0:.017:3.1],[3.1:-.017:0]];
    set(gcf,'position',[10 302 1000 300]);
    for pl = 1:3
        sbplot(row,col,pl)
        if pl == 3
            headplot(allpvafs{pl},'/data/common4/RewTwoback/bt73/bt73HeadNew.spl','electrodes','off','maplimits' ,[0 100]);
        view(20+ads(mf),32+sin(mvs(mf))*3);
        title(ttls{pl});
cbar('vert',[1:64],[0,100]);
        else
            headplot(allpvafs{pl},'/data/common4/RewTwoback/bt73/bt73HeadNew.spl','electrodes','off','maplimits' ,[0 100]);
        view(20+ads(mf),32+sin(mvs(mf))*3);
         title(ttls{pl});
       end;            
    end;
    M2(mf) = getframe(gcf);             
    close
end;
    movie2avi(M2,'/data/common4/RewTwoback/bt73/fullpvaf.avi','fps',8,'quality',100)

figure;topoplot(pvafsall,EEG.chanlocs,'electrodes','off','plotrad',.75);cbar
clear pvafs
for cmp = 1:length(gdcomps{nx})
    cp = gdcomps{nx}(cmp);
    [pvaf,pvafs(cmp,:),vars] = eeg_pvaf(EEG,cp);% w
end;
for ch = 1:size(pvafs,2)
    [maxpv(1,ch) y] = max(pvafs(:,ch));
end;
%figure;topoplot(maxpv,EEG.chanlocs,'electrodes','off','plotrad',.75);cbar
figure;m=1; clear M2
for mf = 1:180
    figure;ads = [0:2:360]; mvs = [[0:.017:3.1],[3.1:-.017:0]];
    headplot(pvafsall,'/data/common4/RewTwoback/bt73/bt73HeadNew.spl','electrodes','off','maplimits' ,[0 100],'cbar',0);
    view(20+ads(mf),32+sin(mvs(mf))*3); 
    %M2(mf) = getframe;             
    M2(mf) = getframe(gcf);             
    close
end;
%[pvaf,pvafsall,vars] = eeg_pvaf(EEG,gdcomps{nx});
%figure;topoplot((maxpv./pvafsall)*100,EEG.chanlocs,'electrodes','off','plotrad',.75);cbar

    end;
end;
, 'title', ['Subj ',int2str(nx),': max pvaf at all channels']
    movie2avi(M2,'/data/common4/RewTwoback/bt73/chanpvaf.avi','fps',8,'quality',100)
eyec = EEG.icawinv(1,:);
frontch = EEG.data(96,:);
eyeproj = eyec(96)*EEG.icaact(1,:);
eegplot([eyeproj;frontch]);

            [pvaf,pvafs,vars] = eeg_pvaf(EEG,gdcomps{nx});
         figure;topoplot(pvafsall,EEG.chanlocs,'electrodes','off','plotrad',.75);cbar
