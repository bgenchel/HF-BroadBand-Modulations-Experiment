% Plots 2 sets of corresponding dipole spheres and connects them with lines (bilateral dipoles have dashed)
% Accommodates pairwise combinations or multiple connections per dipole (see complist)
%
%  PlotCrossLinesWts(complist,paths,datset,corecomps,wtsmat,wtlims,frqmat,frqlims,row,col,place,zoom,pairwise);    
%
% INPUTS:
% complist -- Cell array of integers containing interconnected components from all subjects OR
%             Cell array containing list1 and list2 which are paired component lists for all
%             subjects (ie, list1 and list2 are of cells of size: 1 X #subjs)
% paths -- Cell array of stings corresponding to full path to all subjects (must correspond to 
%          complist order)
% datset -- [string] dataset to open where dipole locations for each subject can be found.
% corecomps -- Cell array (one cell per subj) of integers pertaining to components to plot in a 
%               different color to indicate a core set of components (ie, strongest pwr modulation)
%              [] to not plot corecomps any differently from other comps
% wtsmat -- cell array (one cell per subj) of decimals corresponding to the normalized dot
%           product with the highest RMS component in the 'corecomps'. Line color and width
%           will be scaled by this.
% wtlims -- if not empty, should be two number vector indicating max/min wts for line width/color scaling.
% frqlims -- if not empty, should be two number vector indicating max/min frqs to form colorbar.
% place -- [number] indicates the  
% zoom -- [decimal] magnification factor, 1 is fine for <4 sbplots, up to 2.5 for many sbplots
% pairwise -- [1|0] 1 indicates that 'complist' gives two lists per subject to plot lines between
%                   Otherwise, it will look at all lists for each subject as separate inter-connected sets
% orihist -- if not [], will calculate the spherical histogram of all the lines connecting dipoles.
% viewnum -- [2|3] plots either 2 or 3 different views of the brain
%

function [angs] = PlotCrossLinesWts(complist,paths,datset,corecomps,wtsmat,wtlims,frqmat,frqlims,row,col,place,zoom,pairwise,orihist,viewnum);
    
    angs = [];
    allbesa1 = [];
    if ~isempty(frqmat) & ~isempty(wtsmat)
       	error('Both frqmat and wtsmat cannot be plotted at once. One of these cell arrays must be [].');
    end;    
    if pairwise == 1 % indicates that cells for each subject should be paired
        for nx = 1:length(complist)
            if ~isempty(complist{nx}{1})
                list1{nx}{1} = complist{nx}{1};
                list2{nx}{1} = complist{nx}{2};
                wts1{nx}{1} = wtsmat{nx}{1};
                wts2{nx}{1} = wtsmat{nx}{2};                
            end;
        end;
    else% otherwise need to construct pairwise sets of comps for all subjects
        % create pairwise combinations between all:
        for nx = 1:length(complist)
            if ~isempty(complist{nx})
                for sets = 1:length(complist{nx})
                    if length(complist{nx}{sets}) > 1
                        tmplist = zeros(1,0);
                        tmplist2 = zeros(1,0);
                        tmpwts1 = zeros(1,0);
                        tmpwts2 = zeros(1,0);
                        tmpfrqs1 = zeros(1,0);
                        tmpfrqs2 = zeros(1,0);
                        for cp = 1:length(complist{nx}{sets})
                            tls = complist{nx}{sets};
                            tls = tls(cp+1:end);
                            tmplist(end+1:end+length(tls)) = complist{nx}{sets}(cp);
                            tmplist2(end+1:end+length(tls)) = tls;
                            list1{nx}{sets} = tmplist;
                            list2{nx}{sets} = tmplist2;
                            if ~isempty(wtsmat)
                                wttls = wtsmat{nx}{sets};
                                wttls = wttls(cp+1:end);
                                tmpwts1(end+1:end+length(wttls)) = wtsmat{nx}{sets}(cp);
                                tmpwts2(end+1:end+length(wttls)) = wttls;
                                wts1{nx}{sets} = tmpwts1;
                                wts2{nx}{sets} = tmpwts2;
                            end;
                            if ~isempty(frqmat)                                
                                tmpfrqs1(end+1:end+length(tls)) = repmat(frqmat{nx}{sets},[1 length(tls)]);
                                tmpfrqs2(end+1:end+length(tls)) = repmat(frqmat{nx}{sets},[1 length(tls)]);;
                                frqs1{nx}{sets} = tmpfrqs1;
                                frqs2{nx}{sets} = tmpfrqs2;
                            end;
                        end;
                    else
                        list1{nx}{sets} = [];
                        list2{nx}{sets} = [];
                        if ~isempty(wtsmat)
                            wts1{nx}{sets} = [];
                            wts2{nx}{sets} = [];
                        end;
                        
                        if ~isempty(frqmat)
                            frqs1{nx}{sets} = [];
                            frqs2{nx}{sets} = [];
                        end;                            
                    end;
                end;                    
            end; 
        end;
    end;
   
    
    new=1;new2 = 1 ; new3=1;concatwts1 = zeros(1,0);concatwts2 = zeros(1,0);concatwts3 = zeros(1,0);
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
            for sets = 1:length(list1{nx})
                if ~isempty(list1{nx}{sets})
                    dipsources = [];
                    dipsources.posxyz = EEG.dipfit.model(list1{nx}{sets}(1)).posxyz;
                    dipsources.momxyz = EEG.dipfit.model(list1{nx}{sets}(1)).momxyz;
                    dipsources.rv = EEG.dipfit.model(list1{nx}{sets}(1)).rv;p=1;                
                    for w = 1:length(list1{nx}{sets})                
                        dipsources(1,p).posxyz = EEG.dipfit.model(list1{nx}{sets}(w)).posxyz;
                        dipsources(1,p).momxyz = EEG.dipfit.model(list1{nx}{sets}(w)).momxyz;
                        dipsources(1,p).rv = EEG.dipfit.model(list1{nx}{sets}(w)).rv; 
                        if ~isempty(wtsmat)
                            concatwts1(1,end+1) = wts1{nx}{sets}(w);% compile list of weights 4 lines
                        elseif ~isempty(frqmat)
                            concatwts1(1,end+1) = frqs1{nx}{sets}(w);  
                        end;
                        
                        p=p+1;
                    end;
                    if new == 1
                        allbesa1 = dipsources; new = 0;
                    else
                        allbesa1(end+1:end+size(dipsources,2)) = dipsources; 
                    end;   dipsources = [];            
                    dipsources.posxyz = EEG.dipfit.model(list2{nx}{sets}(1)).posxyz;
                    dipsources.momxyz = EEG.dipfit.model(list2{nx}{sets}(1)).momxyz;
                    dipsources.rv = EEG.dipfit.model(list2{nx}{sets}(1)).rv;p=1;                
                    for w = 1:length(list2{nx}{sets})                
                        dipsources(1,p).posxyz = EEG.dipfit.model(list2{nx}{sets}(w)).posxyz;
                        dipsources(1,p).momxyz = EEG.dipfit.model(list2{nx}{sets}(w)).momxyz;
                        dipsources(1,p).rv = EEG.dipfit.model(list2{nx}{sets}(w)).rv;  
                        if ~isempty(wtsmat)
                            concatwts2(1,end+1) = wts2{nx}{sets}(w);% compile list of weights 4 lines
                        elseif ~isempty(frqmat)
                            concatwts2(1,end+1) = frqs2{nx}{sets}(w);  
                        end;
                        p=p+1;
                    end;
                    if new2 == 1
                        allbesa2 = dipsources; new2 = 0;
                    else
                        allbesa2(end+1:end+size(dipsources,2)) = dipsources; 
                    end;   
                end;
                dipsources = [];  
                if ~isempty(corecomps)
                if ~isempty(corecomps{nx})
                    for p = 1:length(corecomps{nx})
                        dipsources(1,p).posxyz = EEG.dipfit.model(corecomps{nx}(p)).posxyz;
                        dipsources(1,p).momxyz = EEG.dipfit.model(corecomps{nx}(p)).momxyz;
                        dipsources(1,p).rv = EEG.dipfit.model(corecomps{nx}(p)).rv; 
                    end;
                    if ~isempty(frqmat)
                        concatwts3(1,end+1) = frqmat{nx}{sets};  
                    end;
                    if new3 == 1
                        allbesa3 = dipsources; new3 = 0;
                    else
                        allbesa3(end+1:end+size(dipsources,2)) = dipsources; 
                    end; 
                end;
                else allbesa3 = [];
                end;
            end;
        end;
    end;
    
    linecol = hot(10); % for each .1 step a new color
    linecol([1,2,9],:) = [];
    linecol([2,5,6],:) = [];
    coolcol = winter(10);
    coolcol([1,2,9],:) = [];
    coolcol([2,5,6],:) = [];
    
    linewide = [2 1.75 1.5 1];

    concatwts = [concatwts1 concatwts2]; clear colcell colcellcore
    if ~isempty(frqmat)
        concatwts = concatwts*10;
        concatwts3 = concatwts3*10;
        if isempty(frqlims)% if weighting cxns instead of dipoles
            frqlims(1) = floor(min(concatwts));
            frqlims(2) = ceil(max(concatwts));
            frqcols = jet(length([frqlims(1):1:frqlims(2)]));
            for f = 1:length(concatwts)
                colcell{f} = frqcols(round(concatwts(f))-(frqlims(1)-1),:);
            end; 
            for f = 1:length(concatwts3)
                colcellcore{f} = frqcols(round(concatwts3(f))-(frqlims(1)-1),:);
            end;            
        else
            frqlims = frqlims*10;
            frqcols = jet(length([frqlims(1):1:frqlims(2)]));
            for f = 1:length(concatwts)
                idxn = round(concatwts(f))-(frqlims(1)-1);
                if idxn < 0
                    colcell{f} = [.1 .1 .1];
                    %idxn = 1;
                elseif idxn > size(frqcols,1)
                    idxn = size(frqcols,1);
                    colcell{f} = frqcols(idxn,:);
                else
                    colcell{f} = frqcols(idxn,:);
                end;                
            end;             
            for f = 1:length(concatwts3)
                idxn = round(concatwts3(f))-(frqlims(1)-1);
                if idxn < 0
                    %idxn = 1;
                colcellcore{f} = [.1 .1 .1];
                elseif idxn > size(frqcols,1)
                    idxn = size(frqcols,1);
                    colcellcore{f} = frqcols(idxn,:);
                else
                    colcellcore{f} = frqcols(idxn,:);
                end;      
            end;            
            frqlims = frqlims/10;
        end;            
        concatwts = concatwts/10;
        concatwts3 = concatwts3/10;
    end;
    
    
    
    if ~isempty(allbesa1)
        allab = allbesa1;
        allab(end+1:end+size(allbesa2,2)) = allbesa2;
    else
        allab = [];
    end;
    %figure;
    coresz = 33; %24
    othersz = 25; %18
    for sbpt = 1:length(viewnum)
        sbplot(row,col,place); allcoo = [];
        place = place + 1;
        if ~isempty(allab)
            if ~isempty(frqmat)
                dipplot(allab ,'dipolelength',0,'dipolesize',coresz,'normlen','on','gui','off','image','mri','spheres','on','color',colcell,'coordformat','spherical');
            else
                dipplot(allab ,'dipolelength',0,'dipolesize',othersz,'normlen','on','gui','off','image','mri','spheres','on','color',{[.6 .4 .8]},'coordformat','spherical');
            end;        
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
            ff=cell2mat(ff); 
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
            pl=1;bilines=[];plc=1;allcoo=[]; nxt = 1;clear wtsvec
            
            for cc = 1:totsource/2
                collwts = zeros(1,0);
                coords2 = []; coords3 = []; tpbi = []; clear coords   %%%% Do first of the dipole pair
                ids = find(ff(1:haf)==cc);
                for ns = 1:length(ids)            
                    tp ={newstruct(ids(ns)).eleccoord}; tp = tp{1};
                    if length(ids) > 1 & ns == 1
                        tpbi(1,:) =  tp;  
                    end;
                    if ns == 1
                        coords(1,:) = tp;
                        collwts(1,end+1) = concatwts(1,nxt);
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
                        collwts(1,end+1) = concatwts(1,length(concatwts)/2 +nxt);
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
                nxt = nxt+1;
                allcoo{plc} = coords; 
                wtsvec{plc} = sum(collwts);plc=plc+1;
                if ~isempty(coords2) % check for one-sided bilateral
                    allcoo{plc} = coords2;
                    wtsvec{plc} = sum(collwts);plc=plc+1;
                end;
                if ~isempty(coords3)% check for two--sided bilateral
                    allcoo{plc} = coords3;
                    wtsvec{plc} = sum(collwts);plc=plc+1;
                end;               
                if ~isempty(tpbi)  % check for bilateral dipoles
                    bilines{pl} = tpbi; pl = pl+1;
                end;
            end;
            findmax = cell2mat(wtsvec);
            if isempty(wtlims)
                maxx = max(findmax);
                minn = min(findmax);
            else
                maxx = wtlims(2);
                minn = wtlims(1);
            end; 
            
            
            incr = (maxx-minn)/(length(linewide)+1);
            if ~isempty(wtsmat)
                cnt = 1;
                for cc = 1:size(allcoo,2)
                    if wtsvec{cc} < minn+incr
                        typ = 4; % was opposite (1:4), don't know why...
                    elseif wtsvec{cc} > minn+incr & wtsvec{cc} < minn+incr+incr
                        typ = 3;
                    elseif wtsvec{cc} > minn+incr+incr & wtsvec{cc} < maxx - incr
                        typ = 2;
                    elseif wtsvec{cc} > maxx - incr 
                        typ = 1;
                    end;            
                    
                    tp = allcoo{cc};    
                    ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
                    set(ph,'linewidth',linewide(typ));
                    if sbpt == 1
                        if ~isempty(orihist)
                            minu = [0 0 0] - tp(1,:);
                            newdip = tp(2,:) + minu;
                            angs(1,cnt) = atan(newdip(1,2)/newdip(1,1));
                            angs(2,cnt) = atan(newdip(1,2)/newdip(1,3));
                            %angs(1,cnt) = (atan(newdip(1,2)/newdip(1,1)))*(180/pi);
                            %angs(2,cnt) = (atan(newdip(1,2)/newdip(1,3)))*(180/pi);
                            %allspots(cnt,:) = [cos(ang1),sin(ang1),cos(ang2)]; 
                            cnt = cnt+1;
                        else
                            allspots = [];
                        end;
                  end;
                    if ~isempty(wtsmat)
                        set(ph,'color',linecol(typ,:));
                    else
                        set(ph,'color',[.6 0 .7]);
                    end;                
                end;
            end;
            if ~isempty(bilines)
                for cc = 1:size(bilines,2)
                    tp = bilines{cc};    
                    if ~isempty(wtsmat)
                        ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','c');hold on;
                    else
                        ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
                        
                    end;
                    set(ph,'linestyle','--'); set(ph,'linewidth',1);
                end;
            end;
            hold on;
        end;
        if ~isempty(allbesa3)
            if ~isempty(frqmat)
                dipplot(allbesa3 ,'dipolelength',0,'dipolesize',coresz+10,'normlen','on','gui','off','image','mri','spheres','on','color',colcellcore,'coordformat','spherical');       
            else
                dipplot(allbesa3 ,'dipolelength',0,'dipolesize',coresz,'normlen','on','gui','off','image','mri','spheres','on','color',{[.6 0 .9]},'coordformat','spherical');
            end;
        end;
        set(gcf,'color','w');%
        if viewnum(sbpt) == 1
            camzoom(zoom*.85);
            view(0,90)
        else
            camzoom(zoom);
        end;
        
        if viewnum(sbpt) == 3
            view(0,0); 
            if ~isempty(frqmat)
                cbar('vert',1:64,[frqlims(1) frqlims(2)],8);
            end;
        elseif viewnum(sbpt) == 4         
            view(63,22);
        elseif viewnum(sbpt) == 2            
            view(90,0)
            if viewnum(sbpt) == length(viewnum)
                if ~isempty(frqmat)
                    cbar('vert',[1:64],[frqlims(1) frqlims(2)],8);
                end;
                %set(gca,'yticklabel',{frqlims(1) [] 
            end;               
        end;   
        if isempty(allcoo)
            camzoom(1.45)
        end;        
    end; % to sbpt (subplots)       
    
