% rating of each emotion on 3 scales, Then plot in 3D space

% rate [active-passive, valence, strong-weak]; on a scale of -1 1
aw = [-.9,.9,.5];
fr = [.9,-.9,.7];
jy = [.9,.9,.8];
an = [.9,-.9,.8];
sd = [-1,-.9,-.7];
ha = [.3,.8,.8];
fe = [.9,-.5,-.4];
lv = [-.5,.8,.8];
jl = [.8,-.9,.1];
cp = [-.7,.3,-.2];
em = [.1,-.1,-.6];
ct = [.1,.7,.6];
gr = [-.8,-.5,-.7];
rl = [-.9,.9,.9];

allrates = {aw,fr,jy,an,sd,ha,fe,lv,jl,cp,em,ct,gr,rl};

emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','embarrass','content','grief','relief'};

% plot the 3D graph

cols = jet(14);
figure;
for x = 1:length(allrates)
    plot3(allrates{x}(1),allrates{x}(2),allrates{x}(3)); hold on;
    ph = text(allrates{x}(1),allrates{x}(2),allrates{x}(3),emos{x});
    col = [1 1 1] - allrates{x};
    col = col/max(col);
    set(ph,'color',col);
    set(ph,'fontsize',17);
end;
set(gca,'xgrid','on');
set(gca,'ygrid','on');
set(gca,'zgrid','on');

