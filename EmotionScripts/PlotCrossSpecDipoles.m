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
% primary -- ['on' | 'off'] if 'on', will only plot lines between primary 
%                           cross spectral modulation (not comodulation within IM)
%
%
%
%

function PlotCrossSpecDipoles(complist,paths,datset,row,col,place,zoom,viewnum,primary);
    
    
    if ~exist('primary')
        primary = 'on'; % default is to only plot primary cross coherence
    end;
    corecolor = {'r'};
    corelinecol = 'g';
    comodcolor = {'b'};
    linecol = 'c';
    regsize = 25; % size of comod spheres; 25 normally  ******
    linsz = 1.5; % width of comod lines; 1.5 normally    ****
                 % set linsz = 0 for no comod lines!
    
    angs = []; newdips = []; realdips = []; dists = []; countbil = [];
    newdips = zeros(0,6);  realdips = zeros(0,6);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create pairwise combinations between all:
    list1 = cell(1,length(complist)); % for comods
    list2 = cell(1,length(complist));
    corecross1 = cell(1,length(complist)); % for solo cross spec
    corecross2 = cell(1,length(complist));
    
    for nx = 1:length(complist)
        if ~isempty(complist{nx})
        if ~isempty(complist{nx}{1})
            list1{nx} = cell(1,length(complist{nx}));
            list2{nx} = cell(1,length(complist{nx}));
            corecross1{nx} = cell(1,length(complist{nx}));
            corecross2{nx} = cell(1,length(complist{nx}));
            for sets = 1:length(complist{nx})
                if size(complist{nx}{sets},2) > 1 % comod
                    corecross1{nx}{sets} = complist{nx}{sets}(1,:);
                    corecross2{nx}{sets} = complist{nx}{sets}(2,:);
                    tls = unique(complist{nx}{sets})';
                    for cp = 1:length(tls)
                        tptl = tls(:,cp+1:end);
                        list1{nx}{sets} = [list1{nx}{sets} repmat(tls(cp),[1 length(tptl)])];
                        list2{nx}{sets} = [list2{nx}{sets} tptl];
                    end;
                else % 'solo' cross specs
                    corecross1{nx}{sets} = complist{nx}{sets}(1,:);
                    corecross2{nx}{sets} = complist{nx}{sets}(2,:);
                    list1{nx}{sets} = [];
                    list2{nx}{sets} = [];
                end;                    
            end; 
        end; 
        end; 
    end;
    EEG=[]; allbesa1 = []; allbesa2=[]; allbesa3=[];
    new=1;new2 = 1 ; new3 = 1 ;  count = 0;
    for nx = 1:length(list1)
        if ~isempty(list1{nx}) 
            if ~isempty(list1{nx}{1}) |~isempty(corecross1{nx}{1})
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
                            if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                                count = count+1;
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
                            if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                                count = count+1;
                            end;
                            p=p+1;
                        end;
                        if new2 == 1
                            allbesa2 = dipsources; new2 = 0;
                        else
                            allbesa2(end+1:end+size(dipsources,2)) = dipsources; 
                        end;   
                    end; 
                end;  
                countbil{1} = count;
                for sets = 1:length(corecross1{nx})
                    dipsources = [];    
                    dipsources.posxyz = EEG.dipfit.model(corecross1{nx}{sets}(1)).posxyz;
                    dipsources.momxyz = EEG.dipfit.model(corecross1{nx}{sets}(1)).momxyz;
                    dipsources.rv = EEG.dipfit.model(corecross1{nx}{sets}(1)).rv;p=1;                
                    for w = 1:length(corecross1{nx}{sets})                
                        dipsources(1,p).posxyz = EEG.dipfit.model(corecross1{nx}{sets}(w)).posxyz;
                        dipsources(1,p).momxyz = EEG.dipfit.model(corecross1{nx}{sets}(w)).momxyz;
                        dipsources(1,p).rv = EEG.dipfit.model(corecross1{nx}{sets}(w)).rv;  
                        if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                            count = count+1;
                        end;
                        p=p+1;
                    end;
                    if ~exist('corebesa1')
                        corebesa1 = dipsources; 
                    else
                        corebesa1(end+1:end+size(dipsources,2)) = dipsources; 
                    end; 
                    dipsources = [];    
                    dipsources.posxyz = EEG.dipfit.model(corecross2{nx}{sets}(1)).posxyz;
                    dipsources.momxyz = EEG.dipfit.model(corecross2{nx}{sets}(1)).momxyz;
                    dipsources.rv = EEG.dipfit.model(corecross2{nx}{sets}(1)).rv;p=1;                
                    for w = 1:length(corecross2{nx}{sets})                
                        dipsources(1,p).posxyz = EEG.dipfit.model(corecross2{nx}{sets}(w)).posxyz;
                        dipsources(1,p).momxyz = EEG.dipfit.model(corecross2{nx}{sets}(w)).momxyz;
                        dipsources(1,p).rv = EEG.dipfit.model(corecross2{nx}{sets}(w)).rv;  
                        if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                            count = count+1;
                        end;
                        p=p+1;
                    end;
                    if ~exist('corebesa2')
                        corebesa2 = dipsources; 
                    else
                        corebesa2(end+1:end+size(dipsources,2)) = dipsources; 
                    end;                     
                end;
            end;
        end;
    end;
    
    
    %%%%  Do actual plotting:---------------------------
    allab = allbesa1;
    allab(end+1:end+size(allbesa2,2)) = allbesa2;
    coreab = corebesa1;
    coreab(end+1:end+size(corebesa2,2)) = corebesa2;
    %figure;
    for sbpt = 1:length(viewnum)
        sbplot(row,col,place);        place = place + 1;
        if ~isempty(corebesa1)
            [sources realX realY realZ XE YE ZE] = dipplot(coreab ,'dipolelength',0,'dipolesize',regsize,'normlen','on','gui','off','image','mri','spheres','on','color',corecolor,'coordformat','spherical');
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
                    ph=line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color',corelinecol);
                    hold on;set(ph,'linewidth',linsz);% 
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
                        ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
                        set(ph,'linestyle','--'); set(ph,'linewidth',linsz);
                    end;
                end;
            end;
        end;
        if ~isempty(allbesa1) & strcmp(primary,'off')
            [sources realX realY realZ XE YE ZE] = dipplot(allab ,'dipolelength',0,'dipolesize',regsize,'normlen','on','gui','off','image','mri','spheres','on','color',comodcolor,'coordformat','spherical');
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
                        ph = line([tp(1,1) tp(2,1)],[tp(1,2) tp(2,2)],[tp(1,3) tp(2,3)],'color','m');hold on;
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
        if ~isempty(corebesa1)
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
        %if isempty(allcoo)
        %    camzoom(1.45)
        %end;
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
