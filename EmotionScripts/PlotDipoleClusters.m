% Plots clusters of ICs as color-coded dipoles in MRI images (side, rear, top and oblique angles)
%
% PlotDipoleClusters(datset,fullpaths,origlist,clustcps,clustvec,row,col,place,ttl,viewnum,colvec,centroid);
%
% datset -- [string] name of dataset containing dipoles for all subjects
% fullpaths -- [cell array of strings] full paths for subject directories where datset is found
% origlist -- [cell array of numbers] list of 'good components' for each subject
% clustcps -- [cell array of numbers] subset of origlist indicating components in current cluster
%             clustcps{clust}{subj}
% clustvec -- [vector of numbers] list of clusters to plot in same head space
% row -- [number] number of rows in figure (by default plots 4 views, so 2 rows, 2 col is standard)
% col -- [number] number of columns in figure
% place -- [number] subplot position to start plotting
% ttl -- [string] figure title
% viewnum -- [vector] list of views to plot: 1=top, 2=side, 3=rear, 4 is an oblique view;
%                     length(viewnum) gives the number of subplots that will be produced and the
%                     values within the vector tell the orientation and order of views
% colvec -- [vector or matrix] if 1 x 3 vector of RGB values, this will plot all dipoles as the
%           same color. ex. [1 0 0] is red, [0 0 1] is blue, [0 1 0] is green.
%           If a matrix, should be n x 3, with the number of rows equal to the number
%           of clusters to be plotted and the columns should be RGB values for each.
%           If [], will plot clusters as 'jet' colorscale from the first to the last cluster
%           requested (therefore an alternate way to control dipole color is to input a specific
%           order of clusters).
%             [] will assign colors from hsv color scale.
% centroid -- ['only', 'add', 'off'] 'only' will plot only cluster centroids, 'add' will superimpose
%             centroids over cluster dipoles, 'off' will skip centroid plotting
%

function [centrstr1,centrstr2] = PlotDipoleClusters(datset,fullpaths,origlist,clustcps,clustvec,row,col,place,ttl,viewnum,colvec,centroid);


centrstr1 = [];
centrstr2 = [];


if isempty(clustvec)
   clustvec = [1:length(clustcps)];
end;

if ~exist('onecolor')
   onecolor = [];
end;
if ~exist('centroid')
   centroid = 'off';
end;
%col = 2;
%row = 2;

for nx = 1:length(origlist)
   dipsources = []; 
   if ~isempty(origlist{nx})
      if iscell(datset)
         EEG = pop_loadset(datset{nx} ,fullpaths{nx});
      else
         EEG = pop_loadset(datset,fullpaths{nx});
      end;
      if isfield(EEG.dipfit.model,'diffmap')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');
      end;
      if isfield(EEG.dipfit.model,'active')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');
      end;
      if isfield(EEG.dipfit.model,'select')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');
      end;
      dipsources.posxyz = EEG.dipfit.model(origlist{nx}(1)).posxyz;
      dipsources.momxyz = EEG.dipfit.model(origlist{nx}(1)).momxyz;
      dipsources.rv = EEG.dipfit.model(origlist{nx}(1)).rv;p=1;

      for w = 1:length(origlist{nx})
         dipsources(1,p).posxyz = EEG.dipfit.model(origlist{nx}(w)).posxyz;
         dipsources(1,p).momxyz = EEG.dipfit.model(origlist{nx}(w)).momxyz;
         dipsources(1,p).rv = EEG.dipfit.model(origlist{nx}(w)).rv;
         p=p+1;
      end;
      allbesa1{nx} = dipsources; new = 0;
   end;
end;
% adjust color matrix for dipoles:---------------
if isempty(colvec)
   cols = hsv(length(clustvec));
else
   cols = colvec;
end;
if ~isempty(onecolor)
   cols = onecolor;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%
new = 1;     pp=1;  bic = 1;
centrstr1 = struct('posxyz',[0 0 0],'momxyz',[0 0 0],'rv',0);
centrstr2 = struct('posxyz',[0 0 0],'momxyz',[0 0 0],'rv',0);
for clst =1:length(clustvec)
   clust = clustvec(clst);
   centr1 = [];
   centr2 = [];
   for nx = 1:length(clustcps{clust})
      if ~isempty(clustcps{clust}{nx})
         for k = 1:length(clustcps{clust}{nx})
            if new == 1
               allbesa = allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k)));
               centr1 = [centr1; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(1,:)];
               if size(allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz,1) > 1& allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,1) ~= 0 % actual values, not zero
                  if allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,2) > 0 % on the wrong side, switch with centr1
                     centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                     centr2(end,2) = centr2(end,2)*-1; centr1(end,2) = centr1(end,2)*-1;
                  else
                     centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                  end;
               end;
               new = 0;
            else
               allbesa(1,end+1) = allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k)));
               centr1 = [centr1; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(1,:)];
               if size(allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz,1) > 1 & allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,1) ~= 0 % actual values, not zero
                  if allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,2) > 0 % on the wrong side, switch with centr1
                     centr2 = [centr2; allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                     centr2(end,2) = centr2(end,2)*-1; centr1(end,2) = centr1(end,2)*-1;
                  else
                     centr2 = [centr2;allbesa1{nx}(find(origlist{nx} == clustcps{clust}{nx}(k))).posxyz(2,:)];
                  end;
               end;
               % $$$                         if length(allbesa) > 1
               % $$$                             if allbesa(end-1).posxyz(1,3) < -40
               % $$$                                 allbesa(end-1) = [];
               % $$$                                 centr1(end-1,:) = [];
               % $$$                                 %centr2(end-1,:) = [];
               % $$$                             end;
               % $$$                         end;
               % $$$                         if allbesa(end).posxyz(1,3) < -40
               % $$$                             allbesa(end) = [];
               % $$$                             centr1(end,:) = [];
               % $$$                             %centr2(end,:) = [];
               % $$$                         end;
            end;
            colset{pp} = cols(clst,:); pp = pp+1;
         end;
      end;
   end;
   if length(allbesa) > 1
      centr1 = mean(centr1,1);
      centr2 = mean(centr2,1);
      centrstr1(clst).posxyz = centr1;
      centrstr1(clst).momxyz = allbesa(end).momxyz(1,:); % arbitrary?
      centrstr1(clst).rv = 2; % arbitrary?
      centcols{clst} = cols(clst,:);
      centcols2{clst} = cols(clst,:)/2;
      if ~isempty(find(centr2))
         centrstr2(bic).posxyz = centr2;
         centrstr2(bic).momxyz = allbesa(end).momxyz(1,:);
         centrstr2(bic).rv = 2;
         secondcol1{bic} = cols(clst,:);
         secondcol2{bic} = cols(clst,:)/2;
         bic = bic + 1; % separate count for bilaterals
      end;
   end;
end;

if length(allbesa) > 1
   for sbpt = 1:length(viewnum)
      if viewnum(sbpt) < 4
         prjimg = 'off';
      else
         prjimg = 'on';
      end;
      sbplot(row,col,place)
      if strcmp(centroid,'only')
         dipplot(centrstr1,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);
         if ~isempty(find(centrstr2(1).posxyz)) % only if there were bilaterals
            dipplot(centrstr2,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',secondcol1,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0); camzoom(.8)
         else
            camzoom(1)
         end;
      elseif strcmp(centroid,'add')
         dipplot(allbesa,'image','mri','gui','off','dipolelength',0,'dipolesize',25,'normlen','on','spheres','on','color',colset,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);
         dipplot(centrstr1,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',centcols2,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);  camzoom(.8)
         if ~isempty(find(centrstr2(1).posxyz)) % only if there were bilaterals
            dipplot(centrstr2,'image','mri','gui','off','dipolelength',0,'dipolesize',40,'normlen','on','spheres','on','color',secondcol2,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);camzoom(.8)
         else
            camzoom(1)
         end;
      else % centroid off
         dipplot(allbesa,'image','mri','gui','off','dipolelength',0,'dipolesize',25,'normlen','on','spheres','on','color',colset,'projlines','off','projimg',prjimg,'coordformat',EEG.dipfit.coordformat);hold on; view(90,0);  camzoom(1)
      end;
      if viewnum(sbpt) == 3
         view(0,0)
      elseif viewnum(sbpt) == 1
         view(0,90)
      elseif viewnum(sbpt) == 4
         view(63,22);
      end;
      place = place+1;
   end;
   if ~isempty(ttl)
      ph = text(-75,-75,125,ttl); set(ph,'color','r');
      %ph = title(ttl); set(ph,'color','r');
   end;
   set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
end;

