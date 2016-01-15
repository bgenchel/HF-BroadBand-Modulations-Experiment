% Plots relative power at a particular frequency between ICs either within
% subject or between clusters of ICs
%
%
% complist -- [cell array of indices *or* cell array of cell arrays from
% all subjects]
% paths -- [cell array of strings] directory, for each subject, where
% 'datset' and 'savedat' can be found.
% savedat -- [string] name of .fdt file that saved with relative power
% values in MeanPwrDiffs.m
% freq -- [number] frequency to plot. Function will choose the closest
% frequency from the frequency bins saved in 'savedat'
% figplace -- [row col subplot]
%
% PlotCoModasDipoles(complist,justcomps,paths,datset,row,col,place,zoom,pairwise,viewnum,modwts,justwts,soloblack,btstrap);

function PlotRelativePwr(complist,paths,datset,savedat,freq,viewnum,figplace);

if ~exist('figplace')
    figplace = [2 2 1];
elseif isempty(figplace)
    figplace = [2 2 1];
end;
row = figplace(1); col = figplace(2); place = figplace(3);
if ~iscell(complist) % if only one subject
    for cc = 1:length(complist)
        c{cc}{1} = complist(cc);
    end;
    complist = c; clear c
end;

list1 = cell(1,length(complist{1}));
list2= cell(1,length(complist{1}));
reldiffs= cell(1,length(complist{1}));
for clust1 = 1:length(complist)-1
    for clust2 = clust1+1:length(complist) % for all cluster pairs
        for nx = 1:length(complist{clust1}) % find all intra-subj pairs
            s = load([paths{nx},savedat,'.mat']);  % load header info
            [val frbin] = min(abs(s.freqs - freq));% closest freq bin
            diffs = floatread([paths{nx},savedat,'.fdt'],[length(s.complist)-1 length(s.complist) length(s.freqs)],[],0);
            if ~isempty(complist{clust2}{nx}) % if subj is in both clusters
                for ic1 = 1:length(complist{clust1}) % for all ics in cluster 1
                    list1{nx} = [list1{nx},complist{clust1}{nx}(ic1)];
                    for ic2 = 1:length(complist{clust2})
                        if complist{clust1}{nx}(ic1) > complist{clust2}{nx}(ic2) % wrong order
                            flip = -1; % this ensures list1 > list2 directionality
                        else
                            flip = 1;
                        end;
                        list2{nx} = [list2{nx},complist{clust2}{nx}(ic2)];
                        idx1 = find(ismember(s.complist,complist{clust1}{nx}(ic1)));
                        idx2 = find(ismember(s.complist,complist{clust2}{nx}(ic2)));
                        if isempty(idx1) || isempty(idx2)
                            fprintf('\n\nERROR: One or more requested ICs were not found in the power differences array\n');
                            return;
                        end;
                        reldiffs{nx} = [reldiffs{nx}, flip*diffs(idx1,idx2,frbin)];
                    end;
                end;
            end;
        end;
    end;
end;

% collect dipole information for all ICs
EEG=[]; allbesa1 = []; allbesa2=[]; plotmodwts=[];
new=1;new2 = 1 ; count = 0;
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
        dipsources.rv = EEG.dipfit.model(list1{nx}(1)).rv;
        for w = 1:length(list1{nx})
            dipsources(1,w).posxyz = EEG.dipfit.model(list1{nx}(w)).posxyz;
            dipsources(1,w).momxyz = EEG.dipfit.model(list1{nx}(w)).momxyz;
            dipsources(1,w).rv = EEG.dipfit.model(list1{nx}(w)).rv;
            if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                count = count+1;
            end;
            modwts1(1,w) = reldiffs{nx}(w);
        end;
        if new == 1
            allbesa1 = dipsources; new = 0;
            plotmodwts = modwts1;
        else
            allbesa1(end+1:end+size(dipsources,2)) = dipsources;
            plotmodwts = [plotmodwts modwts1];
        end;
        dipsources = [];
        dipsources.posxyz = EEG.dipfit.model(list2{nx}(1)).posxyz;
        dipsources.momxyz = EEG.dipfit.model(list2{nx}(1)).momxyz;
        dipsources.rv = EEG.dipfit.model(list2{nx}(1)).rv;
        for w = 1:length(list2{nx})
            dipsources(1,w).posxyz = EEG.dipfit.model(list2{nx}(w)).posxyz;
            dipsources(1,w).momxyz = EEG.dipfit.model(list2{nx}(w)).momxyz;
            dipsources(1,w).rv = EEG.dipfit.model(list2{nx}(w)).rv;
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

allab = allbesa1;
allab(end+1:end+size(allbesa2,2)) = allbesa2;


bilatcol = 'y';
solocol = [.35 0 .8]; % usually black: [0 0 0] for non-comodulated ICs or purple:[.35 0 .8]
monocol = {'k'}; % specify a monocolor (plots all the same color)
linecol = 'g';
coresz = 25;  % size of non-comod spheres; 26 n0ormally ****
regsize = 25; % size of comod spheres; 25 normally  ******
linsz = 1.5; % width of comod lines; 1.5 normally    ****
newdips = []; realdips = [];  countbil = [];
newdips = zeros(0,6);  realdips = zeros(0,6);  
for sbpt = 1:length(viewnum)
    sbplot(row,col,place);        place = place + 1;
    if ~isempty(allbesa1)
        [sources realX realY realZ XE YE ZE] = dipplot(allab ,'dipolelength',0,'dipolesize',regsize,'normlen','on','gui','off','image','mri','spheres','on','color',monocol,'coordformat',EEG.dipfit.coordformat);
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
        %cnt = 1;
        cols = jet(200);
        normwts = round(100*(plotmodwts/max(plotmodwts) - .01));
        for cc = 1:size(allcoo,2)
            tp = allcoo{cc};
            if linsz > 0
                midpnt = [(tp(1,1)+tp(2,1))/2,(tp(1,2)+tp(2,2))/2,(tp(1,3)+tp(2,3))/2];
                colspot = 101 - abs(normwts(cc));
                if normwts(cc) > 0
                    ph=line([tp(1,1) midpnt(1,1)],[tp(1,2) midpnt(1,2)],[tp(1,3) midpnt(1,3)],'color',cols(colspot,:),'linewidth',2);
                    ph=line([tp(2,1) midpnt(1,1)],[tp(2,2) midpnt(1,2)],[tp(2,3) midpnt(1,3)],'color',cols(201-colspot,:),'linewidth',3);
                else
                    ph=line([tp(1,1) midpnt(1,1)],[tp(1,2) midpnt(1,2)],[tp(1,3) midpnt(1,3)],'color',cols(201-colspot,:),'linewidth',2);
                    ph=line([tp(2,1) midpnt(1,1)],[tp(2,2) midpnt(1,2)],[tp(2,3) midpnt(1,3)],'color',cols(colspot,:),'linewidth',3);
                end;
                %ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color',linecol);
                hold on;%set(ph,'linewidth',linsz);
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
    if viewnum(sbpt) == 3
        view(0,0)
    elseif viewnum(sbpt) == 2
        view(90,0)
    elseif viewnum(sbpt) == 4
        view(63,22)
    end;
end;
axcopy
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
