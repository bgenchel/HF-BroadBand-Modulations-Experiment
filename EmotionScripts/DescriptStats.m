The scores are in the order of the following. NaNs were unscored emotions. The first two scores were overall geunuineness and intensity.

overallexperience genuine	intense	awe	joy	frustration	anger	sadness	happy	content	love	fear	compassion	jealousy	grief	relief	excitement	disgust
cats = {'overall','genuine?','intense?','awe','joy','frustration','anger','sadness','happiness','contentment','love','fear','compassion','jealousy','grief','relief','excitement','disgust'};

jw77 = [6 ,7,5,9,9,9,7,9,7,4,1,9,8,6,6,5,7,6];
rr83 = [8 ,8,7,6,7,8,9,8,9,8,8,8,7,9,7,6,6,8];
jl83 = [8 ,7,6,6,7,4,7,9,5,7,8,6,4,5,9,7,3,8];
hs83 = [7 ,7,6,5,7,9,7,8,9,8,8,7,3,7,9,9,6,5];
ef76 = [4 ,2,2,9,2,3,3,3,2,2,2,2,3,2,2,2,2,2];
as82= [7 ,7,7,6,4,7,6,9,6,7,7,6,7,4,6,8,7,9];
ab75 = [9 ,8,3,7,6,3,2,3,6,6,7,4,6,4,4,8,8.5,3];
ps82 = [8 ,6,7,7,9,9,9,9,8,9,8,7,6,9,9,7,9,8]; % 10->9
%ps82 = [8 ,6,7,7,9,10,10,9,8,10,8,7,6,10,10,7,9,8];
mm78 = [6 ,7,6,7,5,6,8,7,8,6,5,7,6,6,5,6,6,7];
es76 =[6 ,6,6,6,NaN,8,7,8,NaN,NaN,9,NaN,7,1,NaN,NaN,8,6];
ts79 =[5 ,5,6,5,5,7,8,9,6,6,6,6,3,6,5,5,6,NaN];
kc66=[4 ,6.5,5.5,8,6,8,8,4,7,7,8,6,6,3,4,7,8,7];
cj82=[9 ,7,7,7,9,7,6,9,8,5,9,9,5,6,8,6,7,8];
mr722=[5 ,5,5,6,4,4,3,5,6,5,7,4,6,3,6,4,4,6];
jw84=[6 ,8,6,7,6,8,6,8,6,7,8,8,7,5,8,6,8,7];
an82=[7 ,9,8,7,7,9,7,8,7,7,8,7,9,7,9,8,5,7];
md85=[7 ,7,7,6,7,8,8,8,7,6,8,9,6,6,9,7,8,7];
dg75=[4 ,6,4,6,9,4,7,4,9,8,5,3,4,3,3,4,4,4];
eb79=[9 ,8,8,9,9,9,9,NaN,9,9,9,9,9,7,9,9,9,9];
dn86=[8 ,7,7,7,6,8,8,6,7,8,7,9,8,9,5,7,8,9];
kl80=[9 ,9,8,8,8,8,9,9,9,8,9,9,9,9,7,8,9,9];
jo82=[8 ,7,7,8,8,7,7,6,8,8,9,5,9,9,5,9,8,9];
ar81=[9 ,8,8,7,7,7,8,8,8,8,8,6,6,8,8,8,8,8];
mr71=[8 ,8,6,7,NaN,6,8,6,5,6,8,6,6,1,7,5,8,7];
dk74=[4 ,8,6,7,7,9,9,5,7,8,7,6,6,4,5,7,7,8];
mr72=[6.5 ,7,7,8,7,5,5,6,7,8,8,2,8,1,7,7,6,6];
sg75=[5 ,6,4,6,5,7,7,8,5,8,6,1,7,5,6,6,8,7];
sr81=[6 ,7,6,7,8,8,9,6,7,8,7,7,8,5,8,6,7,7];
an70=[8 ,7,8,8,7,7,7,7,7,7,8,8,7,7,6,7,7,7];
tv81=[7 ,7,5,5,7,8,7,9,8,NaN,6,3,8,8,5,7,8,6];
kw78=[9 ,6,7,8.5,9,7,6,6,9,9,7,8,5,2,4,1,4,9];
js75=[8 ,8,7,6,7,6,7,7,9,8,8,9,7,5,9,6,7,7];
ms82=[8 ,6,5,5,6,4,4,9,5,6,3,5,5,7,8,5,6,7];
mi83=[9 ,7,5,7,8,7,7,8,8,8,8,3,7,8,6,7,NaN,5];
tl81=[9 ,7.5,8,8,8,8,8,8,8,8,8,8,8,1,8,1,8,8];


allsubj = {tl81,mi83,ms82,js75,kw78,jo82,kl80,ar81,eb79,dg75,an82,jw84,tv81,sr81,an70,sg75,mr72,dk74,dn86,mr71,md85,mr722,cj82,kc66,ts79,es76,mm78,ps82,ab75,as82,ef76,hs83,jl83,rr83,jw77};% same order as 'fullpaths'
%allsubj{31} = NaN(1,17); % bad subj
%allsubj{22} =  NaN(1,17); % repeat subj
% go through and normalize to z score
% do own zscore to not count NaNs
for w = 1:length(allsubj)
    if ~isempty(find(~isnan(allsubj{w})))
        zsubjs{w} = (allsubj{w}-mean(allsubj{w}(find(~isnan(allsubj{w})))))/std(allsubj{w}(find(~isnan(allsubj{w}))));
    end;
end;
clear plotmat
nowsubjs = zsubjs;
%nowsubjs = allsubj;
delsubjs = [];
for w = 1:length(nowsubjs)
    for q = 1:length(nowsubjs{w})
        plotmat(w,q) =nowsubjs{w}(q);        
    end;
    if isempty(nowsubjs{w})
            delsubjs = [delsubjs w];
    else
        if  ~isempty(find(isnan(plotmat(w,:)))) | ~isempty(find(plotmat(w,:) == 0))
            delsubjs = [delsubjs w];
        end;
    end;
end;
%plotmat(delsubjs,:) = []; % delete subjs with NaN or 0's
cats = {'overall','genuine?','intense?','awe','joy','frustration','anger','sad','happy','content','love','fear','compassion','jealousy','grief','relief','excite','disgust'};
pmat = zeros(size(plotmat,1),0);
for e = 1:length(emos)
    pmat(:,end+1) = plotmat(:,find(ismember(cats,emos{e})));
end
zrats = pmat;
fullrats=pmat;
save /data/common1/emotion/SubjRatingsMatrix.mat fullrats zrats emos plotmat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATS on emotion ratings:---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /data/common1/emotion/SubjRatingsMatrix.mat fullrats zrats emos plotmat

zrats([22,31],:)=[];
fullrats([22,31],:)=[];

[P,anovtable,STATS] = anova1(zrats,emos);
P = 3.8855e-06;
figure;compare = multcompare(STATS,'alpha',.001);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
str = ['print /home/julie/Meetings/sfn2009/EmoRatingsANOVA.eps -depsc -adobe -painters'];eval(str)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run ICA on the subject scores:---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotmat(delsubjs,:) = [];
[weights,sphere] = runica(plotmat(:,4:18),'pca',2);
winv = pinv(weights*sphere);acts = (weights*sphere)*plotmat(:,4:18);
toplot = winv;

emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sadness','  compassion','  love','  relief','  contentment','  awe','  happiness','  joy','  excitement'}; % for plotting purposes
cols = lines(25);
figure;
c1 = 2; c2=1;
for e = 1:size(toplot,1)
    ph=plot(toplot(e,c1),toplot(e,c2),'.');hold on;
    set(ph,'markersize',20);set(ph,'color',cols(e,:));
    %ph = text(toplot(e,c1),toplot(e,c2),emo2{e});
    %set(ph,'color',cols(e,:)); 
end;              
xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);
title('Activations (15 x subjs decomp mat)');
title('Activations (subjs x 15 decomp mat)');
title('Winv (subjs x 15 decomp mat)');

figure;% 3D
c1 = 1; c2 = 2; c3 = 3;
for e = 1:size(toplot,1)
    ph=plot3(toplot(e,c1),toplot(e,c2),toplot(e,c3),'.');hold on;
    set(ph,'markersize',25); set(ph,'color',cols(e,:));
    %ph = text(toplot(e,c1),toplot(e,c2),toplot(e,c3),emo2{e});
    %set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
end;
zl = get(gca,'zlim');
for e = 1:size(toplot,1)
    pl =plot3([toplot(e,c1) toplot(e,c1)],[toplot(e,c2) toplot(e,c2)],[zl(1)  toplot(e,c3)]);
    set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
end;
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['Dim ',int2str(c1)]);ylabel(['Dim ',int2str(c2)]);zlabel(['Dim ',int2str(c3)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR chart of subjective emotion ratings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear pm
for q = 1:size(plotmat,2)
    pm(q) = mean(plotmat(find(abs(plotmat(:,q))>0),q));
    pmstd(q) = std(plotmat(find(abs(plotmat(:,q))>0),q));
end;
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sadness','compassion','love','relief','contentment','awe','happiness','joy','excitement','overall','how genuine','how intense'}; 
cats = {'overall','how genuine','how intense','awe','joy','frustration','anger','sadness','happiness','contentment','love','fear','compassion','jealousy','grief','relief','excitement','disgust'};
% adjust color scale to be valence scale:----
w=load('/data/common1/emotion/BehavRatings/RatingTally.mat');
evals = round(w.emoval*100);
evals = evals-min(evals);evals = evals+1;
cols = jet(max(evals));
cols(end+1,:) = [.75 .75 .75];cols(end+1,:) = [.5 .5 .5];cols(end+1,:) = [.35 .35 .35];
evals(end+1:end+3) = [size(cols,1)-2,size(cols,1)-1,size(cols,1) ];
%---------------------------------------------------------------------
figure;
for q = 1:length(pm)
    e = find(ismember(emos,cats{q}));
    ph = bar(e,pm(q));    
    set(ph,'facecolor',cols(evals(e),:));hold on;
    ph = plot([e e],[pm(q)-pmstd(q) pm(q)+pmstd(q)],'k-');
    ph=text(e,.2,emos{e});
    set(ph,'rotation',90);    set(ph,'fontsize',14);
end;
set(gca,'xtick',[]);set(gca,'xticklabel',[]);set(gca,'xlim',[0 19]);
ylabel('Zscore emotion intensity (original range: 1-9)');
print /home/julie/Manuscripts/Gamma/Frontiers/ZscoreRatings.jpg -djpeg
print /home/julie/Manuscripts/Gamma/Frontiers/ZscoreRatings.eps -depsc
ylabel('Emotion intensity ratings (range: 1-9)');
print /home/julie/Manuscripts/Gamma/Frontiers/ActualRatings.jpg -djpeg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINE plot of emotion ratings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sadness','compassion','love','relief','contentment','awe','happiness','joy','excitement','overall','genuine?','intense?'}; 
cats = {'overall','genuine?','intense?','awe','joy','frustration','anger','sadness','happiness','contentment','love','fear','compassion','jealousy','grief','relief','excitement','disgust'};
figure;
for s = 1:size(plotmat,1)
    ph = plot(plotmat(s,:),'k-'); hold on;
end;
for q = 1:size(plotmat,2)
    e = find(ismember(emos,cats{q}));set(gca,'xlim',[0 19]);
    ph = plot(e,plotmat(:,e),'.','color',cols(evals(e),:),'markersize',25); hold on;
end;
set(gca,'xtick',[1:1:18]); set(gca,'xticklabel',emos);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
ylabel('Zscore of emotion intensity (original range: 1-9)');
title('Subject ratings across emotions');
print /home/julie/Manuscripts/Gamma/Frontiers/SubjectRatingsLine.jpg -djpeg
print /home/julie/Manuscripts/Gamma/Frontiers/SubjectRatingsLine.eps -depsc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCATTER plot of general genuine vs intensity:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; % scatter plot    
plotvals1 = []; plotvals2 = [];
for nx = 1:size(plotmat,1)
    if ~isnan(plotmat(nx,3)) | ~isnan(plotmat(nx,2))
        ph = plot(plotmat(nx,3),plotmat(nx,2),'k.','markersize',23);hold on;
        set(ph,'color',[.4 .4 .4]);
        plotvals1 = [plotvals1,plotmat(nx,3)];
        plotvals2 = [plotvals2,plotmat(nx,2)];
    end;
end;
ph = plot(mean(plotvals1),mean(plotvals2),'r.','markersize',33);hold on;
text(mean(plotvals1),mean(plotvals2),'  Mean');
plot([0 0],[get(gca,'ylim')],'k-');plot([get(gca,'xlim')],[0 0],'k-');
xlabel('Zscore of how intense');ylabel('Zscore of how genuine');
 print /home/julie/Manuscripts/Gamma/Frontiers/IntesityVsGenuiness.jpg -djpeg
 print /home/julie/Manuscripts/Gamma/Frontiers/IntesityVsGenuinessRaw.jpg -djpeg


%%  Plot 3 Dims vs each other:    
[x,y,z]=sphere(20);x=x/30;y=y/30; z=z/30;
l=sqrt(x.*x+y.*y+z.*z);
normals = reshape([x./l y./l z./l],[21 21 3]);
figure; % just 3  dims vs each other     
c1 = 3; c2 = 2; c3 = 1;cols = hsv(size(plotmat,1));
for nx = 1:size(plotmat,1)
    colorarray = repmat(reshape([.25 .25 .25] , 1,1,3), [size(z,1) size(z,2) 1]);
    %colorarray = repmat(reshape(cols(nx,:) , 1,1,3), [size(z,1) size(z,2) 1]);
    ph=surf(plotmat(nx,c1)+x,plotmat(nx,c2)+y,plotmat(nx,c3)+z,colorarray, 'EdgeColor','none', 'VertexNormals', normals,'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);hold on;
    %ph=plot3(plotmat(nx,c1),plotmat(nx,c2),plotmat(nx,c3),'.');hold on;
    %set(ph,'markersize',25); set(ph,'color',cols(nx,:));
end;
zl = get(gca,'zlim');
for nx = 1:size(plotmat,1)
    [cx,cy,cz] = cylinder([1 1],200);cx=cx/100;cy=cy/100; cz=cz/100;
    cx(1,:) = cx(1,:) + plotmat(nx,c1);
    cy(1,:) = cy(1,:) + plotmat(nx,c2);
    cz(1,:) = cz(1,:) + plotmat(nx,c3);
    cx(2,:) = cx(2,:) + plotmat(nx,c1);
    cy(2,:) = cy(2,:) + plotmat(nx,c2);
    cz(2,:) = cz(2,:) + zl(1);
    colorarray = repmat(reshape([.25 .25 .25] , 1,1,3), [size(cz,1) size(cz,2) 1]);
    %colorarray = repmat(reshape(cols(nx,:) , 1,1,3), [size(cz,1) size(cz,2) 1]);
    ph=surf(cx,cy,cz,colorarray,'EdgeColor','none','backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.15);hold on;
    %pl =plot3([plotmat(nx,c1) plotmat(nx,c1)],[plotmat(nx,c2) plotmat(nx,c2)],[zl(1)  plotmat(nx,c3)]);
    %set(pl,'color',cols(nx,:)); set(pl,'linewidth',2)             
end;
%view(-19,32);
lighting phong; material shiny;
camlight right;camlight left; 
light;light;light;lightangle(-19,-32); 
set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
xlabel(['How Intense?']);ylabel(['How Genuine?']);zlabel(['Overall Experience?']);
 print /home/julie/Manuscripts/Gamma/Frontiers/IntesityVsGenuiness3D.jpg -djpeg


% replace Nan with 0 for computation

tl81=[7.5,8,8,8,8,8,8,8,8,8,8,8,1,8,1,8,8];
mi83=[7,5,7,8,7,7,8,8,8,8,3,7,8,6,7,0,5];
ms82=[6,5,5,6,4,4,9,5,6,3,5,5,7,8,5,6,7];
js75=[8,7,6,7,6,7,7,9,8,8,9,7,5,9,6,7,7];
kw78=[6,7,8.5,9,7,6,6,9,9,7,8,5,2,4,1,4,9];
jo82=[7,7,8,8,7,7,6,8,8,9,5,9,9,5,9,8,9];
kl80=[9,8,8,8,8,9,9,9,8,9,9,9,9,7,8,9,9];
ar81=[8,8,7,7,7,8,8,8,8,8,6,6,8,8,8,8,8];
eb79=[8,8,9,9,9,9,0,9,9,9,9,9,7,9,9,9,9];
dg75=[6,4,6,9,4,7,4,9,8,5,3,4,3,3,4,4,4];
an82=[9,8,7,7,9,7,8,7,7,8,7,9,7,9,8,5,7];
jw84=[8,6,7,6,8,6,8,6,7,8,8,7,5,8,6,8,7];
tv81=[7,5,5,7,8,7,9,8,0,6,3,8,8,5,7,8,6];
sr81=[7,6,7,8,8,9,6,7,8,7,7,8,5,8,6,7,7];
an70=[7,8,8,7,7,7,7,7,7,8,8,7,7,6,7,7,7];
sg75=[6,4,6,5,7,7,8,5,8,6,1,7,5,6,6,8,7];
mr72=[7,7,8,7,5,5,6,7,8,8,2,8,1,7,7,6,6];
dk74=[8,6,7,7,9,9,5,7,8,7,6,6,4,5,7,7,8];
dn86=[7,7,7,6,8,8,6,7,8,7,9,8,9,5,7,8,9];
mr71=[8,6,7,0,6,8,6,5,6,8,6,6,1,7,5,8,7];
md85=[7,7,6,7,8,8,8,7,6,8,9,6,6,9,7,8,7];
mr722=[5,5,6,4,4,3,5,6,5,7,4,6,3,6,4,4,6];
cj82=[7,7,7,9,7,6,9,8,5,9,9,5,6,8,6,7,8];
kc66=[6.5,5.5,8,6,8,8,4,7,7,8,6,6,3,4,7,8,7];
ts79 =[5,6,5,5,7,8,9,6,6,6,6,3,6,5,5,6,0];
es76 =[6,6,6,0,8,7,8,0,0,10,0,7,1,0,0,8,6];
jw77 = [7,5,9,9,9,7,9,7,4,1,9,8,6,6,5,7,6];
rr83 = [8,7,6,7,8,9,8,9,8,8,8,7,9,7,6,6,8];
jl83 = [7,6,6,7,4,7,9,5,7,8,6,4,5,9,7,3,8];
hs83 = [7,6,5,7,9,7,8,9,8,8,7,3,7,9,9,6,5];
ef76 = [2,2,9,2,3,3,3,2,2,2,2,3,2,2,2,2,2];
as82= [7,7,6,4,7,6,9,6,7,7,6,7,4,6,8,7,9];
ab75 = [8,3,7,6,3,2,3,6,6,7,4,6,4,4,8,8.5,3];
ps82 = [6,7,7,9,10,10,9,8,10,8,7,6,10,10,7,9,8];
mm78 = [7,6,7,5,6,8,7,8,6,5,7,6,6,5,6,6,7];

allsubj = {tl81,mi83,ms82,js75,kw78,jo82,kl80,ar81,eb79,dg75,an82,jw84,tv81,sr81,an70,sg75,mr72,dk74,dn86,mr71,md85,mr722,cj82,kc66,ts79,es76};
for w = 1:length(allsubj)
newsubj{w} = allsubj{w}/max(allsubj{w}); % divides all scores by max score to make internally consistent
end;

% Same without first two overall scores, 
cats = {'awe','joy','frustration','anger','sad','happy','content','love','fear','compassion','jealousy','grief','relief','excite','disgust'};

tl81=[8,8,8,8,8,8,8,8,8,8,1,8,1,8,8];
mi83=[7,8,7,7,8,8,8,8,3,7,8,6,7,0,5];
ms82=[5,6,4,4,9,5,6,3,5,5,7,8,5,6,7];
js75=[6,7,6,7,7,9,8,8,9,7,5,9,6,7,7];
kw78=[8.5,9,7,6,6,9,9,7,8,5,2,4,1,4,9];
jo82=[8,8,7,7,6,8,8,9,5,9,9,5,9,8,9];
kl80=[8,8,8,9,9,9,8,9,9,9,9,7,8,9,9];
ar81=[7,7,7,8,8,8,8,8,6,6,8,8,8,8,8];
eb79=[9,9,9,9,0,9,9,9,9,9,7,9,9,9,9];
dg75=[6,9,4,7,4,9,8,5,3,4,3,3,4,4,4];
an82=[7,7,9,7,8,7,7,8,7,9,7,9,8,5,7];
jw84=[7,6,8,6,8,6,7,8,8,7,5,8,6,8,7];
tv81=[5,7,8,7,9,8,0,6,3,8,8,5,7,8,6];
sr81=[7,8,8,9,6,7,8,7,7,8,5,8,6,7,7];
an70=[8,7,7,7,7,7,7,8,8,7,7,6,7,7,7];
sg75=[6,5,7,7,8,5,8,6,1,7,5,6,6,8,7];
mr72=[8,7,5,5,6,7,8,8,2,8,1,7,7,6,6];
dk74=[7,7,9,9,5,7,8,7,6,6,4,5,7,7,8];
dn86=[7,6,8,8,6,7,8,7,9,8,9,5,7,8,9];
mr71=[7,0,6,8,6,5,6,8,6,6,1,7,5,8,7];
md85=[6,7,8,8,8,7,6,8,9,6,6,9,7,8,7];
mr722=[6,4,4,3,5,6,5,7,4,6,3,6,4,4,6];
cj82=[7,9,7,6,9,8,5,9,9,5,6,8,6,7,8];
kc66=[8,6,8,8,4,7,7,8,6,6,3,4,7,8,7];
ts79 =[5,5,7,8,9,6,6,6,6,3,6,5,5,6,0];
es76 =[6,0,8,7,8,0,0,10,0,7,1,0,0,8,6];
jw77 = [9,9,9,7,9,7,4,1,9,8,6,6,5,7,6];
rr83 = [6,7,8,9,8,9,8,8,8,7,9,7,6,6,8];
jl83 = [6,7,4,7,9,5,7,8,6,4,5,9,7,3,8];
hs83 = [5,7,9,7,8,9,8,8,7,3,7,9,9,6,5];
ef76 = [9,2,3,3,3,2,2,2,2,3,2,2,2,2,2];
as82= [6,4,7,6,9,6,7,7,6,7,4,6,8,7,9];
ab75 = [7,6,3,2,3,6,6,7,4,6,4,4,8,8.5,3];
ps82 = [7,9,10,10,9,8,10,8,7,6,10,10,7,9,8];
mm78 = [7,5,6,8,7,8,6,5,7,6,6,5,6,6,7];

allsubj = {tl81,mi83,ms82,js75,kw78,jo82,kl80,ar81,eb79,dg75,an82,jw84,tv81,sr81,an70,sg75,mr72,dk74,dn86,mr71,md85,mr722,cj82,kc66,ts79,es76,mm78,ps82,ab75,as82,ef76,hs83,jl83,rr83,jw77};


for w = 1:length(allsubj)
newsubj{w} = allsubj{w}/max(allsubj{w}); % divides all scores by max score to make internally consistent
end;
% Same without first two overall scores, low scores to drop set to 0

tl81=[8,8,8,8,8,8,8,8,8,8,0,8,0,8,8];
mi83=[7,8,7,7,8,8,8,8,0,7,8,6,7,0,5];
ms82=[5,6,4,4,9,5,6,0,5,5,7,8,5,6,7];
js75=[6,7,6,7,7,9,8,8,9,7,0,9,6,7,7]; % zero out jealousy
%js75=[6,7,6,7,7,9,8,8,9,7,5,9,6,7,7];
kw78=[8.5,9,7,6,6,9,9,7,8,5,0,4,0,4,9];
jo82=[8,8,7,7,6,8,8,9,5,9,9,5,9,8,9];
kl80=[8,8,8,9,9,9,8,9,9,9,9,7,8,9,9];
ar81=[7,7,7,8,8,8,8,8,6,6,8,8,8,8,8];
eb79=[9,9,9,9,9,9,9,9,9,9,7,9,9,9,9];
dg75=[6,9,4,7,4,9,8,5,3,4,3,3,4,4,4];
an82=[7,7,9,7,8,7,7,8,7,9,7,9,8,5,7];
jw84=[7,6,8,6,8,6,7,8,8,7,5,8,6,8,7];
tv81=[5,7,8,7,9,8,0,6,0,8,8,5,7,8,6];
sr81=[7,8,8,9,6,7,8,7,7,8,5,8,6,7,7];
an70=[8,7,7,7,7,7,7,8,8,7,7,6,7,7,7];
sg75=[6,5,7,7,8,5,8,6,0,7,5,6,6,8,7];
mr72=[8,7,5,5,6,7,8,8,2,8,0,7,7,6,6];
dk74=[7,7,9,9,5,7,8,7,6,6,0,5,7,7,8];
dn86=[7,6,8,8,6,7,8,7,9,8,9,5,7,8,9];
mr71=[7,0,6,8,6,5,6,8,6,6,0,7,5,8,7];
md85=[6,7,8,8,8,7,6,8,9,6,6,9,7,8,7];
mr722=[6,4,4,0,5,6,5,7,4,6,0,6,4,4,6];
cj82=[7,9,7,6,9,8,5,9,9,5,6,8,6,7,8];
kc66=[8,6,8,8,0,7,7,8,6,6,0,0,7,8,7];
ts79 =[5,5,7,8,9,6,6,6,6,0,6,5,5,6,0];
es76 =[6,7,8,7,8,7,7,10,8,7,0,7,7,8,6]; % substituted some values for blank ratings
jw77 = [9,9,9,7,9,7,4,0,9,8,6,6,5,7,6];
rr83 = [6,7,8,9,8,9,8,8,8,7,9,7,6,6,8];
jl83 = [6,7,4,7,9,5,7,8,6,4,5,9,7,0,8];
hs83 = [5,7,9,7,8,9,8,8,7,0,7,9,9,6,5];
ef76 = [9,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
as82= [6,4,7,6,9,6,7,7,6,7,4,6,8,7,9];
ab75 = [7,6,0,0,0,6,6,7,4,6,4,4,8,8.5,0];
ps82 = [7,9,10,10,9,8,10,8,7,6,10,10,7,9,8];
mm78 = [7,5,6,8,7,8,6,5,7,6,6,5,6,6,7];

allsubj = {tl81,mi83,ms82,js75,kw78,jo82,kl80,ar81,eb79,dg75,an82,jw84,tv81,sr81,an70,sg75,mr72,dk74,dn86,mr71,md85,mr722,cj82,kc66,ts79,es76,mm78,ab75,hs83,ps82,as82,ef76,jl83,rr83,jw77};

clear ratings % put in 'emos' order 
for e = 1:length(emos)
    for nx = 1:length(allsubj)
        ratings{nx}(e) = allsubj{nx}(find(strcmp(cats,emos{e})));
    end;
end;

% reorder emo scores
cats = {'awe','joy','frustration','anger','sad','happy','content','love','fear','compassion','jealousy','grief','relief','excite','disgust'}; % this is the order of the scores
% pick ONE of the following :   **********************8
realorder = {'awe', 'frustration','joy','anger','happy','sad','love' ,'fear','compassion','jealousy','content','grief','relief','disgust','excite'}; % this is the order of data (spectra)
realorder = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % this is the reordered order in terms of valence/activity
%   *********************************
for w = 1:length(allsubj)
for q = 1:length(allsubj{w})
    newsubj{w}(q) = allsubj{w}(find(ismember(cats,realorder{q})));
end;
end;% I checked and this works

for w = 1:length(allsubj)
    newsubj{w} = newsubj{w}/max(newsubj{w}); % divides all scores by max score to make internally consistent
end;
