% takes a spectral co-modulation decomp file and finds percent variance accounted for of specified factors
%
% [pcared pv] = SpecCoModPvaf(savedat,fullpath,factors,channels,oforig,permut);
%
%
% INPUTS:
% savedat -- [string] name of datafile containing all relevant info to call in wts and original data
% fullpath -- [string] full data path where 'savedat' can be found
% factors -- [vector] of modulation factor indices to find pvaf for; 
%                     [] for all factors. Pass a 0 to skip factor by factor pv.
% channels -- [integer or vector] of input dimensions for which to test pvaf; 
%                                 [] for all
% oforig -- [0 | 1] 1 to return pvaf of the factor(s) relative to original data, 
%                   0 relative to the PCA-reduced data used for ICA decomposition
% permut -- if 1, will permute the activations; if 2, will permute winv for bootstrapping; input 0 to skip permutation step.
%
% OUPUTS:
% pcared -- [percent] variance accounted for by the PCA-reduced data given to ICA
% pv -- [vector] of percents corresponding to the number of total factors in decomp
%                0 represents a factor that was not calculated.
%                (if oforig == 1: sum of pv should equal pcared)
% won't work for context decomps

function [pcared pv] = SpecCoModPvaf(savedat,fullpath,factors,channels,oforig,permut);
    
        
    if ~exist('permut')
        permut = 0;
    elseif isempty(permut)
        permut = 0;
    end;
    
    s = load([fullpath,savedat,'.mat']); 

    sph=floatread([fullpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpath,savedat,'.wts'],[s.pcs s.pcs],[],0);        
    icamat = floatread([fullpath,savedat,'DAT.fdt'],[length(s.rowmeans) s.numframes],[],0);    
    pcamat = floatread([fullpath,savedat,'.fdt'],[s.pcs s.numframes],[],0);    
    ws = wts*sph;    activations = ws*pcamat;    winv = pinv(ws); 
    clear wts sph ws pcamat
    speceig = floatread([fullpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);
    specwts = speceig*winv;  
    winv = specwts; clear speceig specwts   
    
    if permut == 1
        activations = shuffle(activations,2);
    elseif permut == 2
        winv = shuffle(winv,1);
    end;    
    
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    
    if ~isempty(channels)
        newdat = winv(channels,:)*activations;
        icamat = icamat(channels,:);
    else        
        newdat = winv*activations;
        channels = [1:size(winv,1)];
    end;    
    
    % to calculate pvaf of PCA-reduced data:-----------
    pcared = 1 - (var(icamat' - newdat')/var(icamat'));
    pcared = pcared*100; 
    
    if factors ~= 0
        for ffs = 1:length(factors)
            fac = factors(ffs);
            newmindat = winv(channels,fac)*activations(fac,:); % back proj
            if oforig == 1
                pv(1,fac) = 1 - (var(icamat' - newmindat')/var(icamat'));% calculate pvaf of orig data
            else
                pv(1,fac) = 1 - (var(newdat' - newmindat')/var(newdat'));% calculate pvaf of PCA-reduced data
            end;      
        end;    
        pv = pv*100;
    else
        pv = [];
    end;
    
