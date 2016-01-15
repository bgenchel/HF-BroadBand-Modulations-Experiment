% takes two Newtimef files containing 2 conditions and plots a tile of all differences (straight) of masked
% Plots mean ERSPs for a cluster of components across subjects
%
%
% datfile1 [string] name of first data file with 3 condition ERSP data
% datfile2 [string] name of second data file with 3 condtion ERSP data
% dipset -- [string] name of dataset to load for dipole locations
% fullpaths -- [cell array] cell array of strings with all subject data paths
% origlist -- [cell array] cell array of number vectors corresponding to components run in newtimef
% clustcomps -- [cell array] subset of origlist to plot in avg ERSP plots
% freqlims -- [minfr maxfr]
% timelims -- [mintm maxtm]
% clim -- [real] maximum value for color scale on avg ERSP plots
% shuffnum -- [positive number] Number of repititions of shuffle permutation to perform. alpha = .01
% dfttls -- [cell array of strings] names for two data files for plot titles (will concatenate with ttls)
% ttls -- [cell array of strings] names for conditions (ex, {'Corr','Wrong','Diff'})
% dipson -- [0|1] 1 to plot dipoles, 0 to not plot dipoles
% dfdiff -- [0|1] 1 to plot within file difference by bootstrap instead of the binominal probability

function PlotERSPcondDiffs(datfile1,datfile2,dipset,dipcolor,fullpaths,origlist,clustcomps,freqlims,timelims,clim,shuffnum,maskalpha,dfttls,ttls,dipson,dfdiff);
    
    
    datfiles = {datfile1,datfile2};
    if length(clim) > 1
        fprintf('Colorlim must be a single value (will scale around zero automatically)');
    end;
    s = load ([fullpaths{1},datfiles{1}]);
    fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
    tms=find(s.times > timelims(1) & s.times<timelims(2));
    for df = 1:length(datfiles)
        filename = datfiles{df};
        for cond = 1:3
            eallersps = zeros(length(fr),length(tms),0);
            for nx = 1:length(fullpaths)
                if ~isempty(clustcomps{nx})
                    s = load ([fullpaths{nx},filename]);
                    fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
                    tms=find(s.times > timelims(1) & s.times<timelims(2));
                    
                    for k=1:length(clustcomps{nx})
                        eformask = s.ersp_boot{find(origlist{nx} == clustcomps{nx}(k)),cond};% fr X 2
                        if cond ~= 3
                            eminmask = eformask(fr,1);
                            emaxmask = eformask(fr,2);
                            eminmask = repmat(eminmask,[1 length(tms)]);
                            emaxmask = repmat(emaxmask,[1 length(tms)]);
                        else
                            eminmask = eformask(fr,tms,1);
                            emaxmask = eformask(fr,tms,2);
                        end;
                        ersp = s.comp_ersp{find(origlist{nx} == clustcomps{nx}(k)),cond};
                        ersp = ersp(fr,tms);
                        ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;          
                        eallersps(:,:,end+1) = ersp; 
                    end;
                end;                
            end;
            condersps{df,cond} = eallersps;                
            fprintf('\nDatafile %s, Condition %s done.',int2str(df),int2str(cond));          
            %eallersps = condersps{df,cond};
            [plotmat,minmask,maxmask,shuffnum] = GroupSig(eallersps,.01,binomalpha,'permt');
            plotersps{df,cond} = plotmat;
        end;      
    end;
    
    fprintf('\nPerforming %s permutation tests as requested...',int2str(shuffnum));   
    if dfdiff == 1
    for df = 1:size(condersps,1)
       % start with shuffle within datafile
        tmpmat = mean(condersps{df,1},3) - mean(condersps{df,2},3);
        randidx = [1:size(condersps{df,1},3)*2];
        diffcond = zeros(size(condersps{df,1},1),size(condersps{df,1},2),shuffnum);
        combersps = condersps{df,1};
        combersps(:,:,end+1:end+size(condersps{df,2},3)) = condersps{df,2};
        for rep = 1:shuffnum
            randidx = shuffle(randidx);
            cond1 = combersps(:,:,randidx(1,1:size(condersps{df,1},3)));
            cond2 = combersps(:,:,randidx(1,size(condersps{df,2},3)+1:end));
            diffcond(:,:,rep) = mean(cond1,3) - mean(cond2,3);
        end;
        % make mask
        minmask = zeros(size(diffcond,1),size(diffcond,2));
        maxmask = zeros(size(diffcond,1),size(diffcond,2));
        for ff = 1:size(diffcond,1)
            for tt= 1:size(diffcond,2)
                distr = squeeze(diffcond(ff,tt,:)); distr = sort(distr);
                minmask(ff,tt) = distr(ceil(shuffnum*maskalpha));
                maxmask(ff,tt) = distr(length(distr) - round(shuffnum*maskalpha));
            end;
        end;
        tmpmat = reshape(tmpmat,1,size(tmpmat,1)*size(tmpmat,2));
        minmask = reshape(minmask,1,size(minmask,1)*size(minmask,2));
        maxmask = reshape(maxmask,1,size(maxmask,1)*size(maxmask,2));
        tmpmat(find(tmpmat > minmask & tmpmat < maxmask)) = 0;        
        plotersps{df,3} = reshape(tmpmat,size(condersps{df,1},1),size(condersps{df,1},2));
        fprintf('.');          
    end;
    end;
    % then find diffs across files, but same trial type  
    for xcond = 1:3
        tmpmat = mean(condersps{1,xcond},3) - mean(condersps{2,xcond},3);
        combersps = condersps{1,xcond};
        combersps(:,:,end+1:end+size(condersps{2,xcond},3)) = condersps{2,xcond};
        randidx = [1:size(condersps{df,1},3)*2];
        diffcond = zeros(size(condersps{1,1},1),size(condersps{2,1},2),shuffnum);
        for rep = 1:shuffnum
            randidx = shuffle(randidx);
            cond1 = combersps(:,:,randidx(1,1:size(condersps{1,xcond},3)));
            cond2 = combersps(:,:,randidx(1,size(condersps{2,xcond},3)+1:end));
            diffcond(:,:,rep) = mean(cond1,3) - mean(cond2,3);
        end;
        % make mask
        minmask = zeros(size(diffcond,1),size(diffcond,2));
        maxmask = zeros(size(diffcond,1),size(diffcond,2));
        for ff = 1:size(diffcond,1)
            for tt= 1:size(diffcond,2)
                distr = squeeze(diffcond(ff,tt,:)); distr = sort(distr);
                minmask(ff,tt) = distr(ceil(shuffnum*maskalpha));
                maxmask(ff,tt) = distr(length(distr) - round(shuffnum*maskalpha));
            end;
        end;
        tmpmat = reshape(tmpmat,1,size(tmpmat,1)*size(tmpmat,2));
        minmask = reshape(minmask,1,size(minmask,1)*size(minmask,2));
        maxmask = reshape(maxmask,1,size(maxmask,1)*size(maxmask,2));
        tmpmat(find(tmpmat > minmask & tmpmat < maxmask)) = 0;        
        plotersps{3,xcond} = reshape(tmpmat,size(condersps{1,xcond},1),size(condersps{1,xcond},2));
        fprintf('.');          
    end;
    
    if dipson == 1
    % find dipole locations
    new=1;
    for nx = 1:length(clustcomps)
        if ~isempty(clustcomps{nx})        
            dipsources = []; 
            EEG = pop_loadset(dipset ,fullpaths{nx});
            cformat = EEG.dipfit.coordformat;
            if isfield(EEG.dipfit.model,'diffmap')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'diffmap');      
            end;
            if isfield(EEG.dipfit.model,'active')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'active');      
            end;
            if isfield(EEG.dipfit.model,'select')
                EEG.dipfit.model = rmfield(EEG.dipfit.model,'select');      
            end;
            dipsources.posxyz = EEG.dipfit.model(clustcomps{nx}(1)).posxyz;
            dipsources.momxyz = EEG.dipfit.model(clustcomps{nx}(1)).momxyz;
            dipsources.rv = EEG.dipfit.model(clustcomps{nx}(1)).rv;
            for w = 1:length(clustcomps{nx})
                dipsources(1,w).posxyz = EEG.dipfit.model(clustcomps{nx}(w)).posxyz;
                dipsources(1,w).momxyz = EEG.dipfit.model(clustcomps{nx}(w)).momxyz;
                dipsources(1,w).rv = EEG.dipfit.model(clustcomps{nx}(w)).rv;
            end;           
            if new == 1
                allbesa = dipsources;new = 0;
            else
                allbesa(end+1:end+size(dipsources,2)) = dipsources; 
            end;
            dipsources = []; 
        end; 
    end;
    end;
    if dipson == 1
        row = 4; else row = 3;
    end;
    figure;    
    col = 3; pl = 1;
    for df = 1:size(plotersps,1)
        for cond = 1:size(plotersps,2)
            sbplot(row,col,pl)
            if cond == 3
                imagesc(s.times(tms),s.freqs(fr), plotersps{df,cond},[-clim clim]);
                pl = pl+1;
            else
                imagesc(s.times(tms),s.freqs(fr), plotersps{df,cond},[-clim clim]);
                pl = pl+1;
            end;               
            set(gca,'ydir','norm'); hold on;
            plot([0 0],[get(gca,'ylim')],'k-');
            if df < 3
                title([dfttls{df},' ',ttls{cond}]);
            else
                title(['Diff ',ttls{cond}]);
            end;   
            if df < 3
                set(gca,'xticklabel',[]);
            end;            
        end;
    end; 
    cbar;
    if dipson == 1
    sbplot(row,col,pl)
    mydipplot(allbesa,'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','color',{dipcolor},'projlines','off','projimg','off','coordformat',cformat);hold on;view(0,90); camzoom(1.5)
    sbplot(row,col,pl+1)
    mydipplot(allbesa,'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','color',{dipcolor},'projlines','off','projimg','off','coordformat',cformat);hold on;view(90,0);camzoom(1.5)
    sbplot(row,col,pl+2)
    mydipplot(allbesa,'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','color',{dipcolor},'projlines','off','projimg','off','coordformat',cformat);hold on;view(0,0);camzoom(1.5)
    
    textsc([datfiles{1},' and ',datfiles{2},' Interactions; masked by group permutations (',int2str(shuffnum),')'],'title');
    end;
    textsc(['Cross Condition ERSPs; Diferences masked by shuffle permutation with ',int2str(shuffnum),' reps at p < ',num2str(maskalpha)],'title');
    set(gcf,'color','w');
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
