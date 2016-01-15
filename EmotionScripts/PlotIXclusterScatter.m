function PlotIXclusterScatter(stem,fullpath,icidx);

ixpairs = [];p=1; clear plotwts
for ic1 = 1:size(icidx{1},1)-1
  scatwts = [];
  for ic2 = ic1+1:size(icidx{1},1)
    name = [stem,int2str(icidx{1}(ic1,1))];  
    s = load([fullpath,name,'.mat']);  
    wts = floatread([fullpath,name,'.wts'],[s.numrows s.numrows],[],0);
    sph = floatread([fullpath,name,'.sph'],[s.numrows s.numrows],[],0);     
    ws = wts*sph; acts = ws*dat;  
    scatwts(1,:) = acts(abs(icidx{1}(ic1,2)),:);
    wtsori = (icidx{1}(ic1,2)./abs(icidx{1}(ic1,2)))'; % 1s and -1s according to orientation
    scatwts(1,:) = scatwts(1,:)*wtsori;
    
    name = [stem,int2str(icidx{1}(ic2,1))];  
    s = load([fullpath,name,'.mat']);  
    wts = floatread([fullpath,name,'.wts'],[s.numrows s.numrows],[],0);
    sph = floatread([fullpath,name,'.sph'],[s.numrows s.numrows],[],0);     
    ws = wts*sph; acts = ws*dat;  
    scatwts(2,:) = acts(abs(icidx{1}(ic2,2)),:);
    wtsori = (icidx{1}(ic2,2)./abs(icidx{1}(ic2,2)))'; % 1s and -1s according to orientation
    scatwts(2,:) = scatwts(2,:)*wtsori;
    ixpairs = [ixpairs; [icidx{1}(ic1,1),icidx{1}(ic1,2),icidx{1}(ic2,1),icidx{1}(ic2,2)]];
    plotwts{p} = scatwts*wtsori;p=p+1;
    end;
end;

figure; pl=1;
row = round(sqrt(length(plotwts)));
col = ceil(sqrt(length(plotwts)));
mrkers = {'b.','r+','go','m^','cs','kd','yx','bh','r*','gv','m<','c.','k>','yp'};
for ic = 1:length(plotwts)
  sbplot(row,col,pl);pl=pl+1;
plot(plotwts{ic}(1,:),plotwts{ic}(2,:),mrkers{ic});hold on;
end;
