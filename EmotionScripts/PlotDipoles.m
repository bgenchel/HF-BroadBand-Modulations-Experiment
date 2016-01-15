% Plots selected dipoles from multiple subjects
%
% [dens,minmask,maxmask,mri] = PlotDipoles(datset, paths,gdcomps,wtcell,btstrap,dipargs,densargs,color,clim,norm,gaussblur);
%
% INPUTS:
% datset: (string) name of dataset where dipole locations can be 
%         found for each subj
% paths: (cell array of strings) full directory paths where datset 
%         can be found.
% gdcomps: (cell array of integers) selected components to include 
%           in plot
% wtcell -- cell array of number vectors corresponding to each gdcomps 
%            entry to weight by in plot ([] if none)
% btstrap -- [cell array] of vectors including all 'good components' 
%            for each subject represented in 'gdcomps' above. This list
%            should be larger than 'gdcomps' if un-weighted, but must be the 
%            same list as 'gdcomps' if weighted because in the latter case masking 
%            will be done by shuffling 'wtcell' weights.  [] will skip masking.
% dipargs: [cell array of strings] for dipole plotting (not density), gives
%          the dipplot arguments. 
%          ex.  dipargs = {'image','mri','gui','off','dipolelength',0,...
%          'normlen','on','spheres','on','projlines','off','projimg','off',...
%          'view', [0 0 1]};
% densargs -- [cell array] if not empty, will use the specified 
%             arguments in dipoledensity which correspond to 
%             parameters for mri3dplot.m
%           example: densargs = {'mrislices',[63:-8:-18],'mriview','top',...
%          'geom',[4,3],'cmax',.0015};
% color -- [cell array (dipoles) or string (density)] for dipole plotting, 
%          specify by a cell array of string(s); 
%        ex. {'r'} the desired color for spheres. 'color' is redefined for
%        dipoles that are weighted by 'wtcell', in this case enter [] for 'color'.
%        for dipole density plots: specify 'yred' or 'bred' for 
%        yellow to red or blue to red color scale.  (for density 
%        differences (ie, when positive and negative 'wtcell' included),
%        'color' must be 'yred').
% clim -- [number] maximum value for dipole weightings (dipole plotting only)
% norm -- ['on',[], vector ] For DENSITY plots: 'on' will normalize dipole density 
%         by number of dipoles. If 'on' and weighted by 'wtcell', will normalize 
%         voxel-by-voxel with unweighted density (same procedure while collecting 
%         bootstrap). [] will not normalize density values. For DIPOLE plots: 
%         enter a vector of z-plane coordinates to plot sequentially (similar to 
%         density plot slices). [] will plot standard dipole plot. 
%         Ex [63:-8:-25]
% gaussblur -- [number above 0.1] For most purposes, 10 is a good value for the 
%              gaussian blur of each dipole. But for some applications, the user 
%              may choose to increase or decrease this value. Default: 10.
%
% OUTPUT:
% Returns the 3D density matrix in 'density'
% Plot is either of the specified dipoles as 
% spheres, or of the group dipole density in units of distributed dipoles/cm3.
% minmask and maxmask are only returned if 'btstrap' is not empty.
%
% Author: Julie Onton

function [dens,minmask,maxmask,mri] = PlotDipoles(datset, paths,gdcomps,wtcell,btstrap,dipargs,densargs,color,clim,norm,gaussblur);
    
    mri = [];
    shuffnum = 400;% hard-coded to be ??? iterations
    pval = .01; % for permutation statistics masking
    distmeasure = 'alldistance'; % options are 'alldistance' and 'distance' (latter takes closest dipole from
                                 % each subject, 'alldistance' treats all dipoles equally                   
    if ~iscell(dipargs) & ~isempty(dipargs)
        fprintf('\nPlotDipoles has been modified. Please see ''help PlotDipoles'' for new input arguments.\n');
        return;
    end;
    
    minmask = []; maxmask = []; dens=[];% initialize in case plotting dipoles
    if ~exist('btstrap')
        btstrap = [];
    end;
    if ~exist('clim')
        clim = []; % this is for color-coding for dipole weightings
    end;
    if ~exist('norm')
        norm = 'off'; % no voxel-by-voxel norming (for wtd densities)
    elseif strcmp(norm,'dips')
        norm = 'on';
    end;
    if ~exist('gaussblur')
        gaussblur = 10;
    elseif isempty(gaussblur)
        gaussblur = 10;
    end;
    
    if isempty(color)
        if isempty(densargs)
            color = {'r'};
        else
            color = 'bred';
        end;
    end; 
    
    new=1;
    if ~isempty(densargs) % density plot
        colset = zeros(1,length(gdcomps));pwts= zeros(1,0);
        count = 0;
        for nx = 1:length(gdcomps)
            if ~isempty(gdcomps{nx})        
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
                
                dipsources.posxyz = EEG.dipfit.model(gdcomps{nx}(1)).posxyz;
                dipsources.momxyz = EEG.dipfit.model(gdcomps{nx}(1)).momxyz;
                dipsources.rv = EEG.dipfit.model(gdcomps{nx}(1)).rv; 
                for w = 1:length(gdcomps{nx})
                    dipsources(1,w).posxyz = EEG.dipfit.model(gdcomps{nx}(w)).posxyz;
                    dipsources(1,w).momxyz = EEG.dipfit.model(gdcomps{nx}(w)).momxyz;
                    dipsources(1,w).rv = EEG.dipfit.model(gdcomps{nx}(w)).rv;
                    if size(dipsources(1,w).posxyz,1) > 1 & dipsources(1,w).posxyz(2,1) ~= 0
                        count = count+2;
                    else
                        count = count +1;
                    end;
                    if ~isempty(wtcell)
                        pwts(1,end+1) = wtcell{nx}(w);
                    else
                        pwts(1,end+1) = 1;
                    end;                    
                end;   
                if new == 1
                    allbesa = dipsources; new = 0;
                else
                    allbesa(end+1:end+size(dipsources,2)) = dipsources; 
                end;  
                colset(nx) = length(dipsources);
                dipsources = [];
            end; 
        end;
        pl = 0;clear subjidxmat
        for nx = 1:length(gdcomps)
            if colset(nx) ~= 0
                for nn = 1:colset(nx)
                    subjidxmat(pl+nn) = nx;
                end;
                pl = pl+nn;
            end;
        end;
        
        if ~isempty(btstrap) % detour for bootstrapping:---------------------------
            for nx = 1:length(btstrap)
                if ~isempty(btstrap{nx})        
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
                    
                    dipsources.posxyz = EEG.dipfit.model(btstrap{nx}(1)).posxyz;
                    dipsources.momxyz = EEG.dipfit.model(btstrap{nx}(1)).momxyz;
                    dipsources.rv = EEG.dipfit.model(btstrap{nx}(1)).rv;
                    for w = 1:length(btstrap{nx})
                        dipsources(1,w).posxyz = EEG.dipfit.model(btstrap{nx}(w)).posxyz;
                        dipsources(1,w).momxyz = EEG.dipfit.model(btstrap{nx}(w)).momxyz;
                        dipsources(1,w).rv = EEG.dipfit.model(btstrap{nx}(w)).rv;
                    end;   
                    subjdips{nx} = dipsources;
                    dipsources = [];
                end; 
            end;
            EEG.data =[]; EEG.icasphere = []; EEG.icaweights = [];
            EEG.icawinv = [];EEG.icaact=[];
            if ~isempty(wtcell) % if weighted, normalize by unweighted density
                optdipplot = {allbesa,'normlen','on','gui','off','image','mri','coordformat',EEG.dipfit.coordformat};
                [btunwtdens,btdensneg mri]= mydipoledensity( optdipplot, 'subjind',subjidxmat,'methodparam',gaussblur,'method',distmeasure,'plotargs',densargs); 
            end
            for bt = 1:shuffnum 
                new = 1; 
                btcount = 0;
                for nx = 1:length(gdcomps)
                    if ~isempty(gdcomps{nx})
                        % btstrap list MUST be >= gdcomps list
                        newics = btstrap{nx}(randperm(length(btstrap{nx})));
                        newics = newics([1:length(gdcomps{nx})]);
                        % permute weights as well (this way gdcomps can = btstrap)
                        newwts = pwts(randperm(length(pwts)));
                        
                        dipsources(1,1).posxyz = subjdips{nx}(1).posxyz;
                        dipsources(1,1).momxyz = subjdips{nx}(1).momxyz;
                        dipsources(1,1).rv = subjdips{nx}(1).rv;  
                        for ic = 1:length(newics)
                            curric = find(btstrap{nx} == newics(ic));
                            dipsources(1,ic).posxyz = subjdips{nx}(curric).posxyz;
                            dipsources(1,ic).momxyz = subjdips{nx}(curric).momxyz;
                            dipsources(1,ic).rv = subjdips{nx}(curric).rv; 
                            if size(dipsources(1,ic).posxyz,1) > 1 & dipsources(1,ic).posxyz(2,1) ~= 0
                                btcount = btcount + 2;
                            else
                                btcount = btcount + 1;
                            end;
                        end;
                        if new == 1
                            allboot = dipsources;new = 0;
                        else
                            allboot = [allboot dipsources];
                        end;
                        dipsources = [];
                    end;
                end;
                optdipplot = {allboot,'normlen','on','gui','off','image','mri','coordformat',EEG.dipfit.coordformat};
                [btdens,btdensneg, mri]= mydipoledensity( optdipplot, 'subjind',subjidxmat,'methodparam',gaussblur,'weight',newwts,'plotargs',densargs); 
                if strcmp(norm,'on') % if dip normalization on
                    if isempty(wtcell) % unweighted density: normalize by # dipoles
                        btdens = btdens/btcount; % # dipoles including duals
                    else % normalize by unweighted density
                        btdens = btdens./btunwtdens;% will divide by zeros
                        btdens(find(isnan(btdens))) = 0;
                    end;
                end; 
                if ~isempty(btdensneg) % if pos and neg weights (take difference)
                    btdens = btdens + btdensneg; % pos-neg wtd density
                end;
                allbtdens(:,:,:,bt) = btdens;    
                clear btdens
                if bt==50|bt==100|bt==150|bt==200|bt==250|bt==300|bt==350|bt==400|bt==450|bt==500|bt==550
                    fprintf('\nBOOTSTRAP ITERATION %s COMPLETED.\n\n',int2str(bt));
                end;
            end;
            clear btdens btdensneg btunwtdens mri allboot 
            %allbtdens = allbtdens/std(abs(allbtdens(:))); % added for normalization below
            [allvals y] = sort(allbtdens,4,'ascend');
            minmask = zeros(size(allbtdens,1),size(allbtdens,2),size(allbtdens,3));
            maxmask = zeros(size(allbtdens,1),size(allbtdens,2),size(allbtdens,3));
            minmask = allvals(:,:,:,ceil(shuffnum*pval));
            maxmask = allvals(:,:,:,size(allvals,4) - round(shuffnum*pval));
            % find median bootstrap value for all voxels
            %medbtdens = median(allbtdens,4); % for ratio version
        end;   % end bootstrapping ---------------------------------------------
               %-----------------------------------------------------------------------
        clear allvals
        % Collect density measures:---------------------------------------     
        
        optdipplot = {allbesa,'normlen','on','gui','off','image','mri','coordformat',EEG.dipfit.coordformat};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find density of input dipoles:------------------------------
        [dens,densneg, mri]= mydipoledensity( optdipplot, 'subjind',subjidxmat,'methodparam',gaussblur,'weight',pwts,'method',distmeasure,'plotargs',densargs); 
        %-------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(densneg) % if pos and neg weights (take difference)
            dens = dens + densneg; % pos-neg wtd density
        end;       
        if strcmp(norm,'on') % if dip normalization on
            if pwts(1) == 1 & pwts(end) == 1 % unweighted density: normalize by # dipoles
                dens = dens/count; % # dipoles including duals
            else % normalize by unweighted density
                [unwtdens,densneg, mri]= mydipoledensity( optdipplot, 'subjind',subjidxmat,'methodparam',gaussblur,'weight',[],'method',distmeasure,'plotargs',densargs);  
                dens = dens./unwtdens;% will divide by zeros
                if strcmp(color,'yred')
                    dens(find(isnan(dens))) = 0;
                end;
            end;
        end; 
        % Do actual plotting of dipole density (masked or unmasked):-----------------
        if ~isempty(btstrap)
            %%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % do real divby median boostrap ratio
            %fprintf('\nDividing by median of bootstrap values involves dividing by zero, the NaNs produced will be changed back to zero later.\n');
            %densrat = dens./medbtdens;    
            %minmask = minmask./medbtdens; 
            %maxmask = maxmask./medbtdens; 
            %densrat = log(densrat);
            %minmask = log(minmask);
            %maxmask = log(maxmask);
            %densrat(isnan(densrat)) = 0;
            %minmask(isnan(minmask)) = 0;
            %maxmask(isnan(maxmask)) = 0;
            %densrat(find(dens>minmask & dens<maxmask)) = 0; 
            % because of negative density with weighted density, use normalized
            % limits instead:
            %normdens = dens/std(abs(dens(:)));
            %dens(find(normdens>minmask & normdens<maxmask)) = 0; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dens(find(dens>minmask & dens<maxmask)) = 0; 

            if ~isempty(find(dens)) % if not all masked out
                if ~isempty(densargs)
                    if strcmp(color,'bred')
                        mri3dplot(dens, mri,densargs{:});
                    elseif strcmp(color,'yred')                        
                        mymriplot(dens, mri,densargs{:});
                    end;
                else                    
                    if strcmp(color,'bred')
                        mri3dplot(dens, mri);
                    elseif strcmp(color,'yred')                        
                        mymriplot(dens, mri);
                    end;
                end;   
            end;
        else                
            if ~isempty(densargs)
                if strcmp(color,'bred')
                    mri3dplot(dens, mri,densargs{:});
                elseif strcmp(color,'yred')                        
                    mymriplot(dens, mri,densargs{:});
                end;
            else
                if strcmp(color,'bred')
                    mri3dplot(dens, mri);
                elseif strcmp(color,'yred')                        
                    mymriplot(dens, mri);
                end;
            end;                                   
        end;            
        density = 0;
        
    else % plot dipoles instead (densargs = [];)
        weights = []; clear allbesa
        for nx = 1:length(gdcomps)
            if ~isempty(gdcomps{nx})        
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
                dipsources.posxyz = EEG.dipfit.model(gdcomps{nx}(1)).posxyz;
                dipsources.momxyz = EEG.dipfit.model(gdcomps{nx}(1)).momxyz;
                dipsources.rv = EEG.dipfit.model(gdcomps{nx}(1)).rv;
                for w = 1:length(gdcomps{nx})
                    dipsources(1,w).posxyz = EEG.dipfit.model(gdcomps{nx}(w)).posxyz;
                    dipsources(1,w).momxyz = EEG.dipfit.model(gdcomps{nx}(w)).momxyz;
                    dipsources(1,w).rv = EEG.dipfit.model(gdcomps{nx}(w)).rv;
                end;           
                if ~exist('allbesa')
                    allbesa = dipsources;
                else
                    allbesa(end+1:end+size(dipsources,2)) = dipsources; 
                end;
                dipsources = []; 
                if ~isempty(wtcell)
                    weights = [weights wtcell{nx}];
                end;
            end; 
        end;
        
        if ~isempty(weights) % all colorscales will be blue to red, green is near zero (but not exactly!!)
            clear color
            if max(abs(weights)) < 1
                weights =weights*10;
                if ~isempty(clim)
                    clim = clim*10;
                end;
            end;
            %weights = round(weights*100); 
            %if ~isempty(clim)
            %   clim = round(clim*100);
            %end;
            if max(abs(weights)) < 100
                weights = weights*10; 
                if ~isempty(clim)
                    clim = clim*10;
                end;
            end;
            if max(abs(weights)) < 100 % again (need to be in 100 range)
                weights = weights*10; 
                if ~isempty(clim)
                    clim = clim*10;
                end;
            end;
            if isempty(clim) % create clim
                [clim idx] = max(abs(weights));
                if min(weights) < 0 % neg weights
                    weights = (weights + clim); 
                    colscale = jet(ceil(clim*2));
                else
                    clim = clim - min(weights);
                    weights = (weights - min(weights)) + 1; % start weights at 1
                    colscale = jet(ceil(clim));
                end;
            else % clim specified
                if min(weights) < 0 % neg weights
                    weights = (weights + clim) ; 
                    colscale = jet(ceil(clim*2));
                else
                    clim = clim - min(weights);
                    weights = (weights - min(weights)) + 1; % start weights at 1
                    colscale = jet(ceil(clim));
                end;
            end;
            weights = round(weights);
            weights(find(weights > size(colscale,1))) = size(colscale,1);% last color (if saturated)
            weights(find(weights < 0)) = 1;% saturate blue if still neg
            for d = 1:length(allbesa)
                color{1,d} = colscale(weights(d),:);
            end;
        else % if not weighted: constant color
            if isempty('color')
                color = {'r'};
               color = repmat(color,[1 length(allbesa)]) ;  
            end;
        end;        

        if isempty(dipargs)
            dipargs = {'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat',EEG.dipfit.coordformat};
        else
            dipargs = [dipargs,'coordformat',EEG.dipfit.coordformat];
        end;
        if ~isempty(norm) % if plotting  just a slice of dipoles at a time
            whichview = find(strcmp(dipargs,'view'))+1;
            if dipargs{whichview}(3) == 1  % top
              xyz = 3;
            elseif dipargs{whichview}(1) == 1 % side
              xyz = 1;
            elseif dipargs{whichview}(2) == -1 % rear
              xyz = 2;
            end;
            figure; pl = 1;
            [sources X Y Z XE YE ZE] = mydipplot(allbesa,dipargs{:},'color',color); % to get talcoord
            clf
            % get slices:
            for ss = 1:length(sources)
                zs(ss) = sources(ss).eleccoord(1,xyz);
            end;
            row = round(sqrt(length(norm))); col = ceil(sqrt(length(norm)));% for # of subplots
             %row = 1; col = 5;% for # of subplots
            intvl = abs(norm(2) - norm(1)); % distance between slices
            norm = sort(norm); % always go lowest to highest slice
            for sl = 1:length(norm) % for each slice
              if sl == 1 % at the ends get all stragglers
                nowplot = find(zs > norm(sl)-intvl*2 & zs <= norm(sl)+intvl/2);
              elseif sl == length(norm)
                nowplot = find(zs >= norm(sl)+intvl/2 & zs < norm(sl)+intvl*2);
              else % normally just take slice
                nowplot = find(zs >= norm(sl)-intvl/2 & zs <= norm(sl)+intvl/2);
              end;
              if ~isempty(nowplot)
                newbesa = allbesa(nowplot);
                newcolor = color(nowplot);
                sbplot(row,col,pl);pl = pl+1;
                try % if weirdly placed dipoles for axis tight
                  [sources X Y Z XE YE ZE] = dipplot(newbesa,dipargs{:},'color',newcolor,'axistight','on'); 
                catch
                  nowplot(find(cell2mat({newbesa.rv}) > .5)) = [];
                end;
                if xyz == 2  % good for row=3; col=4;
                  camzoom(1.4); % 1.1
                elseif xyz == 3 % top
                  camzoom(1.2); % 1.1
                elseif xyz == 1 % side
                  camzoom(1.1); % 1.1
                end;                   
              else
                pl = pl+1;
              end;
            end;
        else
          [sources X Y Z XE YE ZE] = mydipplot(allbesa,dipargs{:},'color',color); %,'color',color
        end;
        dens = sources; % make sources the output
    end;        
    %set(gcf,'PaperOrientation','portrait'); set(gcf,'PaperPosition',[0.25 0.25 8 10.5]); 
    set(gcf,'PaperPosition',[0 0 7 1.1]);
    set(gcf,'Position',[520 680 560 420]) ;
    set(gcf,'Papersize',[7 1.1]); 
