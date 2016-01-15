% plots Spectral Co-mod (clusters?) as dipoles color-coded by template strength/orientation
%
%
% [angs,dists,newdips,realdips,countbil] = PlotCoModasDipoles(complist,justcomps,paths,datset,row,col,place,zoom,pairwise,viewnum,modwts,justwts,soloblack,btstrap);
%
%
%
%
%
% viewnum -- [vector] list of views to plot: 1=top, 2=side, 3=rear; length(viewnum) gives the number
%                     of subplots that will be produced and the values within the vector tell the
%                     orientation and order of views
%
% soloblack -- if not [], will plot non-comodulated comps in black (enter 1 for non-[])
% btstrap -- [integer, cell array, or []], if integer, will collect specified number
%            of connection angles from randomly paired dipoles within the plotted
%            dipoles with comodulation connections. A cell array of vectors containing,
%            for all subjects corresponding to 'complist', 'good ICs' that were used in
%            all analyses for each subject.
%
%
%
%
% OUTPUTS:
% angs -- if btstrap is off: [1 x # pairs] vector of connection angles
%         if btstrap is on: cell array: {1} is a [1 x # pairs] vector of connection angles
%                                       {2} is a vector of all connection angles collected by
%                                           randomizing dipole pairs for specified # of iterations
%
%

function [angs,dists,newdips,realdips,countbil] = PlotCoModasDipoles(complist,justcomps,paths,datset,row,col,place,zoom,pairwise,viewnum,modwts,justwts,soloblack,btstrap);

shuffnum = 400;
bilatcol = 'y';
solocol = [.35 0 .8]; % usually black: [0 0 0] for non-comodulated ICs or purple:[.35 0 .8]
monocol = []; % specify a monocolor (plots all the same color)
linecol = 'g';
coresz = 25;  % size of non-comod spheres; 26 n0ormally ****
regsize = 25; % size of comod spheres; 25 normally  ******
linsz = 1.5; % width of comod lines; 1.5 normally    ****
% set linsz = 0 for no comod lines!

angs = []; newdips = []; realdips = []; dists = []; countbil = [];
newdips = zeros(0,6);  realdips = zeros(0,6);
% create hot/cold color scheme
% $$$     colormtmp = hot;
% $$$     colorm = colormtmp;
% $$$     colorm(:,1) =  colormtmp(:,2);
% $$$     colorm(:,2) =  colormtmp(:,3);
% $$$     colorm(:,3) =  colormtmp(:,1);
% $$$     colorm = [colorm; colormtmp(end:-1:1,:) ];
% $$$     colorm = squeeze(permute(colorm, [1 3 2]));
colorm = jet(128); % to make it jet instead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pairwise == 1 % indicates that cells for each subject should be paired
   list1 = complist{1};
   list2 = complist{2};
   if ~isempty(modwts)
      for nx = 1:length(list1)
         wtscell1{nx} = round(100*modwts{1}{nx})/100;
         wtscell2{nx} = round(100*modwts{2}{nx})/100;
      end;
   else
      for nx = 1:length(list1)
         wtscell1{nx} = ones(1,length(list1{nx}));
         wtscell2{nx} = ones(1,length(list1{nx}));
      end;
   end;
else  % otherwise need to construct pairwise sets of comps for all subjects
   % create pairwise combinations between all:
   list1 = cell(1,length(complist));
   list2 = cell(1,length(complist));
   wtscell1 = cell(1,length(complist));
   wtscell2 = cell(1,length(complist));
   for nx = 1:length(complist)
      for sets = 1:length(complist{nx})
         if length(complist{nx}{sets}) > 1
            tmplist = zeros(1,0);
            tmplist2 = zeros(1,0);
            newwts1 = zeros(1,0);
            newwts2 = zeros(1,0);
            for cp = 1:length(complist{nx}{sets})
               tls = complist{nx}{sets};
               tls = tls(cp+1:end);
               tmplist(end+1:end+length(tls)) = complist{nx}{sets}(cp);
               tmplist2(end+1:end+length(tls)) = tls;
               if ~isempty(modwts)
                  newwts1(end+1:end+length(tls)) = round(100*modwts{nx}{sets}(cp))/100;
                  newwts2(end+1:end+length(tls)) = round(100*modwts{nx}{sets}(cp+1:end))/100;
               else
                  newwts1(end+1:end+length(tls)) = 1;
                  newwts2(end+1:end+length(tls)) = 1;
               end;
            end;
            list1{nx} = [list1{nx} tmplist];
            list2{nx} = [list2{nx} tmplist2];
            wtscell1{nx} = [wtscell1{nx} newwts1];
            wtscell2{nx} = [wtscell2{nx} newwts2];
         else
         end;
      end;
   end;
end;
EEG=[]; allbesa1 = []; allbesa2=[]; allbesa3=[];allmodwts1=[];allmodwts2=[]; plotjcwts = [];
new=1;new2 = 1 ; new3 = 1 ;  count = 0;
for nx = 1:length(list1)
   if ~isempty(list1{nx})
      EEG = pop_loadset(datset ,paths{nx});
      if isfield(EEG.dipfit.model,'diffmap')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');
      end;
      if isfield(EEG.dipfit.model,'active')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');
      end;
      if isfield(EEG.dipfit.model,'select')
         EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');
      end;
      dipsources = [];modwts1 = [];
      dipsources.posxyz = EEG.dipfit.model(list1{nx}(1)).posxyz;
      dipsources.momxyz = EEG.dipfit.model(list1{nx}(1)).momxyz;
      dipsources.rv = EEG.dipfit.model(list1{nx}(1)).rv;p=1;

      for w = 1:length(list1{nx})
         dipsources(1,p).posxyz = EEG.dipfit.model(list1{nx}(w)).posxyz;
         dipsources(1,p).momxyz = EEG.dipfit.model(list1{nx}(w)).momxyz;
         dipsources(1,p).rv = EEG.dipfit.model(list1{nx}(w)).rv;
         if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
            count = count+1;
         end;
         modwts1(1,p) = wtscell1{nx}(p);
         p=p+1;
      end;
      if new == 1
         allbesa1 = dipsources; new = 0;
         allmodwts1 = modwts1;
      else
         allbesa1(end+1:end+size(dipsources,2)) = dipsources;
         allmodwts1 = [allmodwts1 modwts1];
      end;   dipsources = [];   modwts2=[];
      dipsources.posxyz = EEG.dipfit.model(list2{nx}(1)).posxyz;
      dipsources.momxyz = EEG.dipfit.model(list2{nx}(1)).momxyz;
      dipsources.rv = EEG.dipfit.model(list2{nx}(1)).rv;p=1;
      for w = 1:length(list2{nx})
         dipsources(1,p).posxyz = EEG.dipfit.model(list2{nx}(w)).posxyz;
         dipsources(1,p).momxyz = EEG.dipfit.model(list2{nx}(w)).momxyz;
         dipsources(1,p).rv = EEG.dipfit.model(list2{nx}(w)).rv;
         if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
            count = count+1;
         end;
         modwts2(1,p) = wtscell2{nx}(p);
         p=p+1;
      end;
      if new2 == 1
         allbesa2 = dipsources; new2 = 0;
         allmodwts2 = modwts2;
      else
         allbesa2(end+1:end+size(dipsources,2)) = dipsources;
         allmodwts2 = [allmodwts2 modwts2];
      end;
   end;
   countbil{1} = count;
   if ~isempty(justcomps)
      if ~isempty(justcomps{nx})
         EEG = pop_loadset(datset ,paths{nx});
         dipsources = [];  jwts = [];
         for p = 1:length(justcomps{nx})
            dipsources(1,p).posxyz = EEG.dipfit.model(justcomps{nx}(p)).posxyz;
            dipsources(1,p).momxyz = EEG.dipfit.model(justcomps{nx}(p)).momxyz;
            dipsources(1,p).rv = EEG.dipfit.model(justcomps{nx}(p)).rv;
            jwts(1,p) = round(100*justwts{nx}(p))/100;
         end;
         if new3 == 1
            allbesa3 = dipsources; new3 = 0;
            plotjcwts = jwts;
         else
            allbesa3(end+1:end+size(dipsources,2)) = dipsources;
            plotjcwts = [plotjcwts jwts];
         end;
      end;
   end;
end;
if ~isempty(monocol)
   allab = allbesa1;
   allab(end+1:end+size(allbesa2,2)) = allbesa2;
   plotmodwts = [allmodwts1,allmodwts2];;
   % assign color cell to all approx wts
   mid = round(size(colorm,1)/2);
   % pos wts:
   for x = 1:length(plotmodwts)
      if plotmodwts(x) > 0
         colcellMod{x} = monocol;
      end;
   end;
   % neg wts:
   for x = 1:length(plotmodwts)
      if plotmodwts(x) < 0
         colcellMod{x} = monocol;
      end;
   end;
   if ~isempty(plotjcwts) % assign 'solo' IC colors
      % pos justcomps wts:
      for x = 1:length(plotjcwts)
         if plotjcwts(x) > 0
            if isempty(soloblack)
               colcelljc{x} = monocol;
            else
               colcelljc{x} = solocol;
            end;
         end;
      end;
      % neg justcomps wts:
      for x = 1:length(plotjcwts)
         if plotjcwts(x) < 0
            if isempty(soloblack)
               colcelljc{x} = monocol;
            else
               colcelljc{x} = solocol;
            end;
         end;
      end;
   end;
else
   allab = allbesa1;
   allab(end+1:end+size(allbesa2,2)) = allbesa2;
   plotmodwts = [allmodwts1,allmodwts2];;
   % assign color cell to all approx wts
   mid = round(size(colorm,1)/2);
   % pos wts:
   for x = 1:length(plotmodwts)
      if plotmodwts(x) > 0
         colcellMod{x} = colorm(mid+round(plotmodwts(x) * mid),:);
      end;
   end;
   % neg wts:
   for x = 1:length(plotmodwts)
      if plotmodwts(x) < 0
         colcellMod{x} = colorm((mid+1)+round(plotmodwts(x) * mid),:);
      end;
   end;
   if ~isempty(plotjcwts)
      % pos justcomps wts:
      for x = 1:length(plotjcwts)
         if plotjcwts(x) > 0
            if isempty(soloblack)
               colcelljc{x} = colorm(mid+round(plotjcwts(x) * mid),:);
            else
               colcelljc{x} = solocol;
            end;
         end;
      end;
      % neg justcomps wts:
      for x = 1:length(plotjcwts)
         if plotjcwts(x) < 0
            if isempty(soloblack)
               colcelljc{x} = colorm((mid+1)+round(plotjcwts(x) * mid),:);
            else
               colcelljc{x} = solocol;
            end;
         end;
      end;
   end;
end;
%figure;
for sbpt = 1:length(viewnum)
   sbplot(row,col,place);        place = place + 1;
   if ~isempty(allbesa1)
      [sources realX realY realZ XE YE ZE] = dipplot(allab ,'dipolelength',0,'dipolesize',regsize,'normlen','on','gui','off','image','mri','spheres','on','color',colcellMod,'coordformat',EEG.dipfit.coordformat);
      htmp = get(gca,'children');
      pl=1;clear  totnum  % Find the total number of sources plotted
      for idx = 1:length(htmp)
         ctmp = get(htmp(idx),'userdata');
         if isstruct(ctmp)
            x = ctmp.name;
            totnum(pl,:) = x;pl = pl+1;
         end;
      end;
      numdips = length(totnum);
      totnum = unique(totnum);
      totsource = size(totnum,1);  %***

      pl=1;clear newstruct  % make new structure to name paired components alike
      for idx = 1:length(htmp)-1
         ctmp = get(htmp(idx),'userdata');
         if isstruct(ctmp)
            if ctmp.name <= totsource/2
               ctmp.cnum = ctmp.name;
            else
               ctmp.cnum = ctmp.name-(totsource/2);
            end;
            %newstruct((numdips+1)-pl) = ctmp;% reorder for backwards compatibility
            newstruct(pl) = ctmp;
            pl = pl+1;
         end;
      end;

      ff={newstruct.cnum};  % next decide where halfway is, taking bilaterals into account
      ff=cell2mat(ff); clear allcoo
      haf = find(ff==1);
      if length(haf) == 4
         haf = haf(2);
      elseif length(haf) == 3 & haf(2) == haf(3)-1
         haf = haf(1);
      elseif length(haf) == 3 & haf(1) == haf(2)-1
         haf = haf(2);
      elseif length(haf) == 2
         haf = haf(1);
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      pl=1;bilines=[];plc=1;allcoo=[];
      for cc = 1:totsource/2
         coords2 = []; coords3 = []; tpbi = []; clear coords   %%%% Do first of the dipole pair
         ids = find(ff(1:haf)==cc);
         for ns = 1:length(ids)
            tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
            if length(ids) > 1 & ns == 1
               tpbi(1,:) =  tp;
            end;
            if ns == 1
               coords(1,:) = tp;
            elseif ns == 2
               tpbi(2,:) =  tp;
               coords2(1,:) = tp;
            end;
         end;
         %%%% Do second of the dipole pair
         ids = find(ff(haf+1:end)==cc); ids = ids+haf;
         for ns = 1:length(ids)
            tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
            if length(ids) > 1 & ns == 1
               tpbi(1,:) =  tp;
            end;
            if ns == 1
               coords(2,:) = tp;
               if ~isempty(coords2)
                  coords2(2,:) = tp;
               end;
            elseif ns == 2
               tpbi(2,:) =  tp;
               if  ~isempty(coords2)
                  coords3 = coords2;
                  coords3(2,:) = tp;
               end;
            end;
         end;
         allcoo{plc} = coords; plc=plc+1;
         if ~isempty(coords2) % check for one-sided bilateral
            allcoo{plc} = coords2;plc=plc+1;
         end;
         if ~isempty(coords3)% check for two--sided bilateral
            allcoo{plc} = coords3;plc=plc+1;
         end;
         if ~isempty(tpbi)  % check for bilateral dipoles
            bilines{pl} = tpbi; pl = pl+1;
         end;
      end;
      cnt = 1;
      for cc = 1:size(allcoo,2)
         tp = allcoo{cc};
         if linsz > 0
            ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color',linecol);
            hold on;set(ph,'linewidth',linsz);
         end;
         dists(1,cnt) = sqrt((tp(1,1)-tp(2,1))^2 + (tp(1,2)-tp(2,2))^2 + (tp(1,3)-tp(2,3))^2);

         if viewnum(sbpt) == 1
            % find leftmost dipole
            if tp(1,1) < tp(2,1)
               newdip = tp(2,:) - tp(1,:); % subtract from dip 2 to move vector to origin
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(2,:) - tp(1,:);
               newdips(end,4:6) = tp(1,:) - tp(2,:);
            else
               newdip = tp(1,:) - tp(2,:);
               realdips(end+1,1:3) = tp(2,:);
               realdips(end,4:6) = tp(1,:);
               newdips(end+1,1:3) = tp(1,:) - tp(2,:);
               newdips(end,4:6) = tp(2,:) - tp(1,:);
            end;
            angs(1,cnt) = atan2(newdip(1,2),newdip(1,1));
            cnt = cnt+1;
         elseif viewnum(sbpt) == 2
            if tp(1,2) < tp(2,2)
               newdip = tp(2,:) - tp(1,:);
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(2,:) - tp(1,:);
               newdips(end,4:6) = tp(1,:) - tp(2,:);
            else
               newdip = tp(1,:) - tp(2,:);
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(1,:) - tp(2,:);
               newdips(end,4:6) = tp(2,:) - tp(1,:);
            end;
            angs(1,cnt) = atan2(newdip(1,3),newdip(1,2));
            cnt = cnt+1;
         end;
      end;
      if ~isempty(bilines)
         for cc = 1:size(bilines,2)
            tp = bilines{cc};
            if linsz > 0
               ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color',bilatcol);hold on;
               set(ph,'linestyle','--'); set(ph,'linewidth',linsz);
            end;
         end;
      end;
   end;
   if ~isempty(allbesa3) % solo ICs
      [solosources realX realY realZ XE YE ZE] = mydipplot(allbesa3 ,'dipolelength',0,'dipolesize',coresz,'normlen','on','gui','off','image','mri','spheres','on','color',colcelljc,'coordformat',EEG.dipfit.coordformat);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      pl=1;bilines=[];
      for cc = 1:length(solosources)
         tp ={solosources(cc).eleccoord}; tp = tp{1};
         if size(tp,1) > 1
            bilines{pl} = tp; pl = pl+1;
         end;
      end;
      if ~isempty(bilines)
         for cc = 1:size(bilines,2)
            tp = bilines{cc};
            if linsz > 0
               ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color',bilatcol);hold on;
               set(ph,'linestyle','--'); set(ph,'linewidth',linsz);
            end;
         end;
      end;
   end;
   set(gcf,'color','w');%
   if viewnum(sbpt) == 3
      view(0,0)
      %camzoom(zoom);
   elseif viewnum(sbpt) == 2
      view(90,0)
      %camzoom(zoom);
   elseif viewnum(sbpt) == 4
      view(63,22)
      %camzoom(zoom);
   end;
   if ~isempty(allbesa3)
      if isempty(allbesa1)
         if sbpt == 1
            camzoom(zoom*1.1);
         else
            camzoom(zoom*1.3);
         end;
      else
         if sbpt == 1
            camzoom(zoom*.75); % if you have both comod and not
         else
            camzoom(zoom*.9);
         end;
      end;
   else
      if sbpt == 1
         camzoom(zoom*1.1);
      else
         camzoom(zoom*1.3);
      end;
   end;
end;
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to collect random connection angles among only cluster dipoles:
if ~isempty(btstrap) & ~iscell(btstrap) & length(viewnum) == 1
   fprintf('\nCollecting %s randomized measures of angle and distance between dipole pairs\n',int2str(btstrap));
   fprintf('Output variables will be changed to cell arrays with bootstrap values in 2nd cell\n\n');


   angs2 = angs; clear angs % remake angs output to be a cell array
   angs{1} = angs2; clear angs2 % first cell is real data
   dists2 = dists; clear dists
   dists{1} = dists2; clear dists2
   % since within cluster, # of bilateral dipoles same, so leave as 1x1 cell

   figure;
   for bt = 1:btstrap % specified number of iterations
      clf
      allab = shuffle(allab);
      sbpt = 1;
      if ~isempty(allab)
         [sources btX(bt,:) btY(bt,:) btZ(bt,:) XE YE ZE] = dipplot(allab ,'dipolelength',0,'dipolesize',15,'normlen','on','gui','off','image','mri','spheres','on','color',colcellMod,'coordformat','spherical');
         htmp = get(gca,'children');
         pl=1;clear  totnum  % Find the total number of sources plotted
         for idx = 1:length(htmp)
            ctmp = get(htmp(idx),'userdata');
            if isstruct(ctmp)
               x = ctmp.name;
               totnum(pl,:) = x;pl = pl+1;
            end;
         end;
         numdips = length(totnum);
         totnum = unique(totnum);
         totsource = size(totnum,1);  %***

         pl=1;clear newstruct  % make new structure to name paired components alike
         for idx = 1:length(htmp)-1
            ctmp = get(htmp(idx),'userdata');
            if isstruct(ctmp)
               if ctmp.name <= totsource/2
                  ctmp.cnum = ctmp.name;
               else
                  ctmp.cnum = ctmp.name-(totsource/2);
               end;
               newstruct(pl) = ctmp;
               pl = pl+1;
            end;
         end;

         ff={newstruct.cnum};  % next decide where halfway is, taking bilaterals into account
         ff=cell2mat(ff); clear allcoo
         haf = find(ff==1);
         if length(haf) == 4
            haf = haf(2);
         elseif length(haf) == 3 & haf(2) == haf(3)-1
            haf = haf(1);
         elseif length(haf) == 3 & haf(1) == haf(2)-1
            haf = haf(2);
         elseif length(haf) == 2
            haf = haf(1);
         end;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%
         pl=1;bilines=[];plc=1;allcoo=[];
         for cc = 1:totsource/2
            coords2 = []; coords3 = []; tpbi = []; clear coords   %%%% Do first of the dipole pair
            ids = find(ff(1:haf)==cc);
            for ns = 1:length(ids)
               tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
               if length(ids) > 1 & ns == 1
                  tpbi(1,:) =  tp;
               end;
               if ns == 1
                  coords(1,:) = tp;
               elseif ns == 2
                  tpbi(2,:) =  tp;
                  coords2(1,:) = tp;
               end;
            end;
            %%%% Do second of the dipole pair
            ids = find(ff(haf+1:end)==cc); ids = ids+haf;
            for ns = 1:length(ids)
               tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
               if length(ids) > 1 & ns == 1
                  tpbi(1,:) =  tp;
               end;
               if ns == 1
                  coords(2,:) = tp;
                  if ~isempty(coords2)
                     coords2(2,:) = tp;
                  end;
               elseif ns == 2
                  tpbi(2,:) =  tp;
                  if  ~isempty(coords2)
                     coords3 = coords2;
                     coords3(2,:) = tp;
                  end;
               end;
            end;
            allcoo{plc} = coords; plc=plc+1;
            if ~isempty(coords2) % check for one-sided bilateral
               allcoo{plc} = coords2;plc=plc+1;
            end;
            if ~isempty(coords3)% check for two--sided bilateral
               allcoo{plc} = coords3;plc=plc+1;
            end;
            if ~isempty(tpbi)  % check for bilateral dipoles
               bilines{pl} = tpbi; pl = pl+1;
            end;
         end;

         cnt = 1;
         for cc = 1:size(allcoo,2)
            tp = allcoo{cc};
            ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
            set(ph,'linewidth',1.5);
            distsbt(bt,cnt) = sqrt((tp(1,1)-tp(2,1))^2 + (tp(1,2)-tp(2,2))^2 + (tp(1,3)-tp(2,3))^2);
            if viewnum(sbpt) == 1
               % find leftmost dipole
               if tp(1,1) < tp(2,1)
                  newdip = tp(2,:) - tp(1,:); % subtract from dip 2 to move vector to origin
                  realdips(end+1,1:3) = tp(1,:);
                  realdips(end,4:6) = tp(2,:);
                  newdips(end+1,1:3) = tp(2,:) - tp(1,:);
                  newdips(end,4:6) = tp(1,:) - tp(2,:);
               else
                  newdip = tp(1,:) - tp(2,:);
                  realdips(end+1,1:3) = tp(2,:);
                  realdips(end,4:6) = tp(1,:);
                  newdips(end+1,1:3) = tp(1,:) - tp(2,:);
                  newdips(end,4:6) = tp(2,:) - tp(1,:);
               end;
               angsbt(bt,cnt) = atan2(newdip(1,2),newdip(1,1));
               cnt = cnt+1;
            elseif viewnum(sbpt) == 2
               if tp(1,2) < tp(2,2)
                  newdip = tp(2,:) - tp(1,:);
                  realdips(end+1,1:3) = tp(1,:);
                  realdips(end,4:6) = tp(2,:);
                  newdips(end+1,1:3) = tp(2,:) - tp(1,:);
                  newdips(end,4:6) = tp(1,:) - tp(2,:);
               else
                  newdip = tp(1,:) - tp(2,:);
                  realdips(end+1,1:3) = tp(1,:);
                  realdips(end,4:6) = tp(2,:);
                  newdips(end+1,1:3) = tp(1,:) - tp(2,:);
                  newdips(end,4:6) = tp(2,:) - tp(1,:);
               end;
               angsbt(bt,cnt) = atan2(newdip(1,3),newdip(1,2));
               cnt = cnt+1;
            end;
         end;
         if viewnum(sbpt) == 3
            view(0,0)
         elseif viewnum(sbpt) == 2
            view(90,0)
         elseif viewnum(sbpt) == 4
            view(63,22)
         end;
      end;
   end;
   % collect bootstrap values into cell array
   angs{2} = angsbt; % second cell contains bootstrap angle values
   dists{2} = distsbt;
   close
   %%%%%%%**********************************************************************************

elseif ~isempty(btstrap) & iscell(btstrap) & length(viewnum) == 1
   % do bootstrapping by shuffling connections among all ICs within subject
   % for angle analysis
   gdcomps = btstrap; % in this case, btstrap is a cell array of IC vectors
   angs2 = angs; clear angs % remake angs output to be a cell array
   angs{1} = angs2; clear angs2 % first cell is real data
   dists2 = dists; clear dists
   dists{1} = dists2; clear dists2

   % first load all gdcomps from all relevant subjects so you only have to load datasets once:
   EEG=[]; allbesa = []; new=1;
   for nx = 1:length(gdcomps)
      if ~isempty(complist{nx})
         EEG = pop_loadset(datset ,paths{nx});
         if isfield(EEG.dipfit.model,'diffmap')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');
         end;
         if isfield(EEG.dipfit.model,'active')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');
         end;
         if isfield(EEG.dipfit.model,'select')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');
         end;
         dipsources = [];
         dipsources.posxyz = EEG.dipfit.model(gdcomps{nx}(1)).posxyz;
         dipsources.momxyz = EEG.dipfit.model(gdcomps{nx}(1)).momxyz;
         dipsources.rv = EEG.dipfit.model(gdcomps{nx}(1)).rv;p=1;
         for w = 1:length(gdcomps{nx})
            dipsources(1,p).posxyz = EEG.dipfit.model(gdcomps{nx}(w)).posxyz;
            dipsources(1,p).momxyz = EEG.dipfit.model(gdcomps{nx}(w)).momxyz;
            dipsources(1,p).rv = EEG.dipfit.model(gdcomps{nx}(w)).rv;
            p=p+1;
         end;
      end;
      subjdips{nx} = dipsources; dipsources = [];
   end;
   %%%%%%%%%%%%%%%%%%%%%% Now do iterations of new comp pairs:
   figure;
   for bt = 1:400 % hard code to 200 iterations
      clf;  clear btlist
      for nx = 1:length(complist)
         if ~isempty(complist{nx})
            for im = 1:length(complist{nx})
               newics = shuffle(gdcomps{nx});
               newics = newics(1:length(complist{nx}{im}));
               btlist{nx}{im} = newics;
            end;
         end;
      end;
      list1 =[]; list2 = [];
      for nx = 1:length(btlist)
         for sets = 1:length(btlist{nx})
            if length(btlist{nx}{sets}) > 1
               tmplist = zeros(1,0);
               tmplist2 = zeros(1,0);
               for cp = 1:length(btlist{nx}{sets})
                  tls = btlist{nx}{sets};
                  tls = tls(cp+1:end);
                  tmplist(end+1:end+length(tls)) = btlist{nx}{sets}(cp);
                  tmplist2(end+1:end+length(tls)) = tls;
               end;
               list1{nx}{sets} = tmplist;
               list2{nx}{sets} = tmplist2;
            else
               list1{nx}{sets} = [];
               list2{nx}{sets} = [];
               list1{nx} = [];
               list2{nx} = [];
            end;
         end;
      end;

      EEG=[]; allbesa1 = []; allbesa2=[];
      new=1;new2 = 1 ;  count=0;
      for nx = 1:length(list1)
         if ~isempty(list1{nx})
            if ~isempty(list1{nx}{1})
               for sets = 1:length(list1{nx})
                  dipsources = [];
                  dipsources.posxyz = subjdips{nx}(1).posxyz;
                  dipsources.momxyz = subjdips{nx}(1).momxyz;
                  dipsources.rv = subjdips{nx}(1).rv;
                  for w = 1:length(list1{nx}{sets})
                     currcp = find(list1{nx}{sets}(w) == gdcomps{nx});
                     dipsources(1,w).posxyz = subjdips{nx}(currcp).posxyz;
                     dipsources(1,w).momxyz = subjdips{nx}(currcp).momxyz;
                     dipsources(1,w).rv = subjdips{nx}(currcp).rv;
                     if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                        count = count+1;
                     end;
                  end;
                  if new == 1
                     allbesa1 = dipsources; new = 0;
                  else
                     allbesa1(end+1:end+size(dipsources,2)) = dipsources;
                  end;
                  dipsources = [];
                  dipsources.posxyz = subjdips{nx}(1).posxyz;
                  dipsources.momxyz = subjdips{nx}(1).momxyz;
                  dipsources.rv = subjdips{nx}(1).rv;
                  for w = 1:length(list2{nx}{sets})
                     currcp = find(list2{nx}{sets}(w) == gdcomps{nx});
                     dipsources(1,w).posxyz = subjdips{nx}(currcp).posxyz;
                     dipsources(1,w).momxyz = subjdips{nx}(currcp).momxyz;
                     dipsources(1,w).rv = subjdips{nx}(currcp).rv;
                     if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                        count = count+1;
                     end;
                  end;
                  if new2 == 1
                     allbesa2 = dipsources; new2 = 0;
                  else
                     allbesa2(end+1:end+size(dipsources,2)) = dipsources;
                  end;
               end;
            end;
         end;
      end;
      countbil{2}(bt) = count;
      allab = allbesa1;
      allab(end+1:end+size(allbesa2,2)) = allbesa2;

      [sources realX realY realZ XE YE ZE] = dipplot(allab ,'dipolelength',0,'dipolesize',regsize,'normlen','on','gui','off','image','mri','spheres','on','color',colcellMod,'coordformat',EEG.dipfit.coordformat);
      htmp = get(gca,'children');
      pl=1;clear  totnum  % Find the total number of sources plotted
      for idx = 1:length(htmp)
         ctmp = get(htmp(idx),'userdata');
         if isstruct(ctmp)
            x = ctmp.name;
            totnum(pl,:) = x;pl = pl+1;
         end;
      end;
      numdips = length(totnum);
      totnum = unique(totnum);
      totsource = size(totnum,1);  %***

      pl=1;clear newstruct  % make new structure to name paired components alike
      for idx = 1:length(htmp)-1
         ctmp = get(htmp(idx),'userdata');
         if isstruct(ctmp)
            if ctmp.name <= totsource/2
               ctmp.cnum = ctmp.name;
            else
               ctmp.cnum = ctmp.name-(totsource/2);
            end;
            %newstruct((numdips+1)-pl) = ctmp;% reorder for backwards compatibility
            newstruct(pl) = ctmp;
            pl = pl+1;
         end;
      end;

      ff={newstruct.cnum};  % next decide where halfway is, taking bilaterals into account
      ff=cell2mat(ff); clear allcoo
      haf = find(ff==1);
      if length(haf) == 4
         haf = haf(2);
      elseif length(haf) == 3 & haf(2) == haf(3)-1
         haf = haf(1);
      elseif length(haf) == 3 & haf(1) == haf(2)-1
         haf = haf(2);
      elseif length(haf) == 2
         haf = haf(1);
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      pl=1;bilines=[];plc=1;allcoo=[];
      for cc = 1:totsource/2
         coords2 = []; coords3 = []; tpbi = []; clear coords   %%%% Do first of the dipole pair
         ids = find(ff(1:haf)==cc);
         for ns = 1:length(ids)
            tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
            if length(ids) > 1 & ns == 1
               tpbi(1,:) =  tp;
            end;
            if ns == 1
               coords(1,:) = tp;
            elseif ns == 2
               tpbi(2,:) =  tp;
               coords2(1,:) = tp;
            end;
         end;
         %%%% Do second of the dipole pair
         ids = find(ff(haf+1:end)==cc); ids = ids+haf;
         for ns = 1:length(ids)
            tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
            if length(ids) > 1 & ns == 1
               tpbi(1,:) =  tp;
            end;
            if ns == 1
               coords(2,:) = tp;
               if ~isempty(coords2)
                  coords2(2,:) = tp;
               end;
            elseif ns == 2
               tpbi(2,:) =  tp;
               if  ~isempty(coords2)
                  coords3 = coords2;
                  coords3(2,:) = tp;
               end;
            end;
         end;
         allcoo{plc} = coords; plc=plc+1;
         if ~isempty(coords2) % check for one-sided bilateral
            allcoo{plc} = coords2;plc=plc+1;
         end;
         if ~isempty(coords3)% check for two--sided bilateral
            allcoo{plc} = coords3;plc=plc+1;
         end;
         if ~isempty(tpbi)  % check for bilateral dipoles
            bilines{pl} = tpbi; pl = pl+1;
         end;
      end;
      cnt = 1;
      for cc = 1:size(allcoo,2)
         tp = allcoo{cc};
         ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
         set(ph,'linewidth',1.5);
         distsbt(bt,cnt) = sqrt((tp(1,1)-tp(2,1))^2 + (tp(1,2)-tp(2,2))^2 + (tp(1,3)-tp(2,3))^2);

         if viewnum(sbpt) == 1
            % find leftmost dipole
            if tp(1,1) < tp(2,1)
               newdip = tp(2,:) - tp(1,:); % subtract from dip 2 to move vector to origin
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(2,:) - tp(1,:);
               newdips(end,4:6) = tp(1,:) - tp(2,:);
            else
               newdip = tp(1,:) - tp(2,:);
               realdips(end+1,1:3) = tp(2,:);
               realdips(end,4:6) = tp(1,:);
               newdips(end+1,1:3) = tp(1,:) - tp(2,:);
               newdips(end,4:6) = tp(2,:) - tp(1,:);
            end;
            angsbt(bt,cnt) = atan2(newdip(1,2),newdip(1,1));
            cnt = cnt+1;
         elseif viewnum(sbpt) == 2
            if tp(1,2) < tp(2,2)
               newdip = tp(2,:) - tp(1,:);
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(2,:) - tp(1,:);
               newdips(end,4:6) = tp(1,:) - tp(2,:);
            else
               newdip = tp(1,:) - tp(2,:);
               realdips(end+1,1:3) = tp(1,:);
               realdips(end,4:6) = tp(2,:);
               newdips(end+1,1:3) = tp(1,:) - tp(2,:);
               newdips(end,4:6) = tp(2,:) - tp(1,:);
            end;
            angsbt(bt,cnt) = atan2(newdip(1,3),newdip(1,2));
            cnt = cnt+1;
         end;
      end;
      % collect bootstrap values into cell array
      angs{2} = angsbt; % second cell contains bootstrap angle values
      dists{2} = distsbt;
   end;
   close
end;
