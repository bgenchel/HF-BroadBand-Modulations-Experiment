% Plots 2 sets of corresponding dipole spheres and connects them with lines (bilateral dipoles have dashed)
% Accommodates pairwise combinations or multiple connections per dipole (see complist)
%
%  PlotCrossLines(
%
% INPUTS:
% complist -- Cell array of integers containing interconnected components from all subjects OR
%             Cell array containing list1 and list2 which are paired component lists for all
%             subjects (ie, list1 and list2 are of cells of size: 1 X #subjs)
% justcomps -- cell array of indexes for all components without comod lines
% paths -- Cell array of stings corresponding to full path to all subjects (must correspond to 
%          complist order)
% datset -- [string] dataset to open where dipole locations for each subject can be found.
% place -- [number] indicates the 
% zoom -- [decimal] magnification factor, 1 is fine for <4 sbplots, up to 2.5 for many sbplots
% pairwise -- [1|0] 1 indicates that 'complist' gives two lists per subject to plot lines between
%                   Otherwise, it will look at all lists for each subject as separate inter-connected sets
% viewnum -- [2|3] plots either 2 or 3 different views of the brain
%

function PlotCrossLines(complist,justcomps,paths,datset,row,col,place,zoom,pairwise,viewnum);
    
    
    if pairwise == 1 % indicates that cells for each subject should be paired
        list1 = complist{1};
        list2 = complist{2};
    else  % otherwise need to construct pairwise sets of comps for all subjects
        % create pairwise combinations between all:
        for nx = 1:length(complist)
            %if ~isempty(complist{nx})
                for sets = 1:length(complist{nx})
                    if length(complist{nx}{sets}) > 1
                        tmplist = zeros(1,0);
                        tmplist2 = zeros(1,0);
                        for cp = 1:length(complist{nx}{sets})
                            tls = complist{nx}{sets};
                            tls = tls(cp+1:end);
                            tmplist(end+1:end+length(tls)) = complist{nx}{sets}(cp);
                            tmplist2(end+1:end+length(tls)) = tls;
                        end;
                        list1{nx}{sets} = tmplist;
                        list2{nx}{sets} = tmplist2;
                    else
                        list1{nx}{sets} = [];
                        list2{nx}{sets} = [];
                    end;                    
                end; 
            %end;
        end;
    end;    
    EEG=[]; allbesa1 = []; allbesa2=[]; allbesa3=[];
    new=1;new2 = 1 ; new3 = 1 ; 
    for nx = 1:length(list1)
        if ~isempty(list1{nx}) 
            if ~isempty(list1{nx}{1}) 
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
                            p=p+1;
                        end;
                        if new2 == 1
                            allbesa2 = dipsources; new2 = 0;
                        else
                            allbesa2(end+1:end+size(dipsources,2)) = dipsources; 
                        end;   
                    end; 
                end;  
            end;EEG=[];
        end;
        if ~isempty(justcomps{nx})
            EEG = pop_loadset(datset ,paths{nx}); 
            dipsources = [];          
            for p = 1:length(justcomps{nx})
                dipsources(1,p).posxyz = EEG.dipfit.model(justcomps{nx}(p)).posxyz;
                dipsources(1,p).momxyz = EEG.dipfit.model(justcomps{nx}(p)).momxyz;
                dipsources(1,p).rv = EEG.dipfit.model(justcomps{nx}(p)).rv; 
            end;
            if new3 == 1
                allbesa3 = dipsources; new3 = 0;
            else
                allbesa3(end+1:end+size(dipsources,2)) = dipsources; 
            end; 
        end;
    end;
        
    coresz = 26;
    allab = allbesa1;
    allab(end+1:end+size(allbesa2,2)) = allbesa2;
    %figure;
    for sbpt = 1:viewnum
        sbplot(row,col,place)
        place = place + 1;
        if ~isempty(allbesa1)
        dipplot(allab ,'dipolelength',0,'dipolesize',25,'normlen','on','gui','off','image','mri','spheres','on','color',{'r'},'coordformat','spherical');
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

        for cc = 1:size(allcoo,2)
            tp = allcoo{cc};    
            ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','g');hold on;
            set(ph,'linewidth',1);
        end;
        if ~isempty(bilines)
            for cc = 1:size(bilines,2)
                tp = bilines{cc};    
                ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
                set(ph,'linestyle','--'); set(ph,'linewidth',2);
            end;
        end;
        end;
        if ~isempty(allbesa3)
            %dipplot(allbesa3 ,'dipolelength',0,'dipolesize',coresz,'normlen','on','gui','off','image','mri','spheres','on','color',{[.6 0 .9]},'coordformat','spherical');% purple
            dipplot(allbesa3 ,'dipolelength',0,'dipolesize',coresz,'normlen','on','gui','off','image','mri','spheres','on','color',{'b'},'coordformat','spherical');
        end;
        set(gcf,'color','w');%
        if sbpt == 3
            view(0,0)
        elseif sbpt == 2
            view(90,0)
        elseif viewnum(sbpt) == 4
            view(63,22)
        end;      
        if ~isempty(allbesa3)
            if isempty(allbesa1)
                if sbpt == 1
                    camzoom(zoom*.99);
                else
                    camzoom(zoom*1.22);
                end;                
            else                
                if sbpt == 1
                    camzoom(zoom*.69);
                else
                    camzoom(zoom*.85);
                end;
            end;
        else
            if sbpt == 1
                camzoom(zoom*.99);
            else
                camzoom(zoom*1.22);
            end;
        end;            
        %if isempty(allcoo)
        %    camzoom(1.45)
        %end;
    end;
    axcopy
  set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
  
