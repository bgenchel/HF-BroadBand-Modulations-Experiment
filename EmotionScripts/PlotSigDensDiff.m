% takes component arrays across subjects and the min- and maxmask computed by SigDensityDiff() and plot sig diffs
%
% [densdiff] = PlotSigDensDiff(datset1,datset2,fullpaths1,fullpaths2,complist1,complist2,weights1,weights2,minmask,maxmask,densargs,norm);
%
% INPUTS:
% datset1 -- [string] name of .set where dipole information can be found for complist1
% datset2 -- [string] name of .set where dipole information can be found for complist2 
%             Note: datset1 and datset2 can be the same.
% fullpaths1 -- [cell array] of strings with full paths of directories where 'datset1' is found for all subjects
% fullpaths2 -- [cell array] of strings with full paths of directories where 'datset2' is found for all subjects
% complist1 -- [cell array] of component indices for density 1
% complist2 -- [cell array] of component indices for density 2
% minmask -- [matrix] 3d matrix of minimum significance threshold returned by SigDensityDiff()
% maxmask -- [matrix] 3d matrix of maximum significance threshold returned by SigDensityDiff()
% densargs -- [cell array] of parameters for dipole density difference plot.
%             ex.: densargs = {'mrislices',[63:-8:-25],'mriview','top','geom',[4,4],'cmax',.001};
% norm -- ['subj','dips', 'off', or cell array of 2 cell arrays] if 'subj', will normlalize dipole density 
%         (dip/cm^3) by dividing by the number of subjects in the group. If
%         'dips', will normlalize dipole density by dividing by the number of
%         dipoles in the group. If a cell array of arrays, this signifies a normalization
%         based on a larger IC list from each subject. This, in effect, will return the 
%         ratio, per voxel, of observed dipole density and the 'total' density. The 'total'
%         dipole list is given in each sub-cell array and should be analogous to 'complist1'
%         and 'complist2'.
%
% Author: Julie Onton

function [densdiff] = PlotSigDensDiff(datset1,datset2,fullpaths1,fullpaths2,complist1,complist2,weights1,weights2,minmask,maxmask,densargs,norm);
    
    
    gaussblur = 10;
    new=1;nsubj1 = 0;pwts1 = [];
    for nx = 1:length(complist1)
        if ~isempty(complist1{nx})        
            EEG = pop_loadset(datset1 ,fullpaths1{nx});
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            
            mod1 = EEG.dipfit.coordformat;
            dipsources.posxyz = EEG.dipfit.model(complist1{nx}(1)).posxyz;
            dipsources.momxyz = EEG.dipfit.model(complist1{nx}(1)).momxyz;
            dipsources.rv = EEG.dipfit.model(complist1{nx}(1)).rv;
            for w = 1:length(complist1{nx})
                dipsources(1,w).posxyz = EEG.dipfit.model(complist1{nx}(w)).posxyz;
                dipsources(1,w).momxyz = EEG.dipfit.model(complist1{nx}(w)).momxyz;
                dipsources(1,w).rv = EEG.dipfit.model(complist1{nx}(w)).rv;
                if ~isempty(weights1)
                    pwts1 = [pwts1 weights1{nx}(w)];
                else
                    pwts1 = [pwts1 1];
                end;                    
            end;   
            if new == 1
                allbesa1 = dipsources; new = 0;
            else
                allbesa1(end+1:end+size(dipsources,2)) = dipsources; 
            end;  
            dipsources = [];
            nsubj1 = nsubj1 + 1;
        end; 
    end;
    new=1;nsubj2 = 0; pwts2 = [];
    for nx = 1:length(complist2)
        if ~isempty(complist2{nx})        
            EEG = pop_loadset(datset2 ,fullpaths2{nx});
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            
            mod2 = EEG.dipfit.coordformat;
            dipsources.posxyz = EEG.dipfit.model(complist2{nx}(1)).posxyz;
            dipsources.momxyz = EEG.dipfit.model(complist2{nx}(1)).momxyz;
            dipsources.rv = EEG.dipfit.model(complist2{nx}(1)).rv;
            for w = 1:length(complist2{nx})
                dipsources(1,w).posxyz = EEG.dipfit.model(complist2{nx}(w)).posxyz;
                dipsources(1,w).momxyz = EEG.dipfit.model(complist2{nx}(w)).momxyz;
                dipsources(1,w).rv = EEG.dipfit.model(complist2{nx}(w)).rv;
                if ~isempty(weights2)
                    pwts2 = [pwts2 weights2{nx}(w)];
                else
                    pwts2 = [pwts2 1];
                end;                    
            end;   
            if new == 1
                allbesa2 = dipsources; new = 0;
            else
                allbesa2(end+1:end+size(dipsources,2)) = dipsources; 
            end;  
            dipsources = [];
            nsubj2 = nsubj2 + 1;
        end; 
    end;

    if ~isempty(weights1)
        optdipplot = {allbesa1,'normlen','on','gui','off','image','mri','coordformat',mod1};
        [unwtdens1]= mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur); 
        optdipplot = {allbesa2,'normlen','on','gui','off','image','mri','coordformat',mod2};
        [unwtdens2] = mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur); 
        mndens = (unwtdens1+unwtdens2)/2;% mean of unwted densities    
        unwted = 0;
    else
        unwted = 1;
    end;
    optdipplot = {allbesa1,'normlen','on','gui','off','image','mri','coordformat',mod1};
    figure;[density1]= mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur,'weight',pwts1);     
    optdipplot = {allbesa2,'normlen','on','gui','off','image','mri','coordformat',mod2};
    figure;[density2,unwt2, mri]= mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur,'weight',pwts2); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct density normalization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iscell(norm) % if normalization by full complist
      fulllist1 = norm{1};
      fulllist2 = norm{2};
      new=1;
      for nx = 1:length(fulllist1)
        if ~isempty(fulllist1{nx})        
          EEG = pop_loadset(datset1 ,fullpaths1{nx});
          if isfield(EEG.dipfit.model,'diffmap')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
          end;
          if isfield(EEG.dipfit.model,'active')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
          end;
          if isfield(EEG.dipfit.model,'select')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
          end;            
          mod1 = EEG.dipfit.coordformat;
          dipsources.posxyz = EEG.dipfit.model(fulllist1{nx}(1)).posxyz;
          dipsources.momxyz = EEG.dipfit.model(fulllist1{nx}(1)).momxyz;
          dipsources.rv = EEG.dipfit.model(fulllist1{nx}(1)).rv;
          for w = 1:length(fulllist1{nx})
            dipsources(1,w).posxyz = EEG.dipfit.model(fulllist1{nx}(w)).posxyz;
            dipsources(1,w).momxyz = EEG.dipfit.model(fulllist1{nx}(w)).momxyz;
            dipsources(1,w).rv = EEG.dipfit.model(fulllist1{nx}(w)).rv;                   
          end;   
          if new == 1
            allbesa1 = dipsources; new = 0;
          else
            allbesa1(end+1:end+size(dipsources,2)) = dipsources; 
          end;  
          dipsources = [];
        end; 
      end;
      new=1; 
      for nx = 1:length(fulllist2)
        if ~isempty(fulllist2{nx})        
          EEG = pop_loadset(datset2 ,fullpaths2{nx});
          if isfield(EEG.dipfit.model,'diffmap')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
          end;
          if isfield(EEG.dipfit.model,'active')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
          end;
          if isfield(EEG.dipfit.model,'select')
            EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
          end;            
          mod2 = EEG.dipfit.coordformat;
          dipsources.posxyz = EEG.dipfit.model(fulllist2{nx}(1)).posxyz;
          dipsources.momxyz = EEG.dipfit.model(fulllist2{nx}(1)).momxyz;
          dipsources.rv = EEG.dipfit.model(fulllist2{nx}(1)).rv;
          for w = 1:length(fulllist2{nx})
            dipsources(1,w).posxyz = EEG.dipfit.model(fulllist2{nx}(w)).posxyz;
            dipsources(1,w).momxyz = EEG.dipfit.model(fulllist2{nx}(w)).momxyz;
            dipsources(1,w).rv = EEG.dipfit.model(fulllist2{nx}(w)).rv;                  
          end;   
          if new == 1
            allbesa2 = dipsources; new = 0;
          else
            allbesa2(end+1:end+size(dipsources,2)) = dipsources; 
          end;  
          dipsources = [];
        end; 
      end;
      optdipplot = {allbesa1,'normlen','on','gui','off','image','mri','coordformat',mod1};
      figure;[fulldensity1]= mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur);     
      optdipplot = {allbesa2,'normlen','on','gui','off','image','mri','coordformat',mod2};
      figure;[fulldensity2]= mydipoledensity( optdipplot, 'method','alldistance','methodparam',gaussblur); 
      density1 = density1./fulldensity1;
      density2 = density2./fulldensity2;
      
    else % otherwise do basic normalization:
      if unwted == 0 % if weighted
                     % divide by non-weighted density, then multiply by mean of non-weighted densities
        density1 = density1./unwtdens1;density1=density1.*mndens;
        density2 = density2./unwtdens2;density2=density2.*mndens;
      elseif strcmp(norm,'dips')  & unwted == 1 % if normalization by # dipoles(unweighted density only)
        density1 = density1/length(allbesa1);
        density2 = density2/length(allbesa2);
      elseif strcmp(norm,'subj')  & unwted == 1 % if normalization by # subjs(unweighted density only)
        density1 = density1/nsubj1;
        density2 = density2/nsubj2;       
      end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform difference and masking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    densdiff = density1 - density2;
    if ~isempty(minmask) % mask it if you got it.
      densdiff(find(densdiff>minmask & densdiff<maxmask)) = 0;
    end;
    % get rid of NaNs from divide by zero:
    nns = isnan(densdiff);
    densdiff(nns) = 0;
    
    if ~isempty(densargs)
        mymriplot(densdiff, mri,densargs{:});
    else
        mymriplot(densdiff, mri);
    end;
    
    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 8 10.5]);
