% takes two newtimef files and plots a tile of all differences 
% Plots mean ERSPs for a cluster of components across subjects
%
%
% datfile1 [string] name of first data file with ERSP data
% datfile2 [string] name of second data file with ERSP data
% dipset -- [string] name of dataset to load for dipole locations
% fullpaths -- [cell array] cell array of strings with all subject data paths
% clustcomps -- [cell array] subset of origlist to plot in avg ERSP plots
% freqlims -- [minfr maxfr]
% timelims -- [mintm maxtm]
% clim -- [real] maximum value for color scale on avg ERSP plots
% alphas [1 x 2 vector] p values for 1) binomial probability and 2) shuffle permutation. ex: [.00001 .01] 
% shuffnum -- [positive number] Number of repititions of shuffle permutation to perform. alpha = .01
% dfttls -- [cell array of strings] names for two data files for plot titles (will concatenate with ttls)
% ttls -- [cell array of strings] names for conditions (ex, {'Corr','Wrong','Diff'})
% dipcolor -- [vector] 1 x 3 vector of RGB color for dipoles

function PlotERSPdiffs(datfile1,datfile2,dipset,fullpaths,clustcomps,freqlims,timelims,clim,shuffnum,alphas,dfttls,ttls,dipcolor);
    
    
    
    datfiles = {datfile1,datfile2};
    if length(clim) > 1
        fprintf('Colorlim must be a single value (will scale around zero automatically)');
    end;
    s = load ([fullpaths{1},datfiles{1}]);
    fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
    tms=find(s.times > timelims(1) & s.times<timelims(2));
    fprintf('\nLoading ERSP data...\n');

    if iscell(s.comp_ersp) % then plot whole tile of diffs
        for df = 1:length(datfiles)
            filename = datfiles{df};
            for cond = 1:3
                eallersps = zeros(length(fr),length(tms),0);
                for nx = 1:length(clustcomps)
                    if ~isempty(clustcomps{nx})
                        s = load ([fullpaths{nx},filename]);
                        origlist{nx} = s.complist;
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
                [plotersps{df,cond}] = GroupSig(eallersps,s.alpha,alphas(1),'binom');
            end;      
        end;
    else
        for df = 1:length(datfiles)
            filename = datfiles{df};
            eallersps = zeros(length(fr),length(tms),0);
            for nx = 1:length(clustcomps)
                if ~isempty(clustcomps{nx})
                    s = load ([fullpaths{nx},filename]);
                    fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
                    tms=find(s.times > timelims(1) & s.times<timelims(2));
                    origlist{nx} = s.complist;                    
                    for k=1:length(clustcomps{nx})
                        eformask = s.ersp_boot(fr,:,find(origlist{nx} == clustcomps{nx}(k)));% fr X 2
                        eminmask = repmat(eformask(:,1),[1 length(tms)]);
                        emaxmask= repmat(eformask(:,2),[1 length(tms)]);
                        ersp = s.comp_ersp(fr,tms,find(origlist{nx} == clustcomps{nx}(k)));
                        ersp(find(ersp > eminmask & ersp < emaxmask)) = 0;          
                        eallersps(:,:,end+1) = ersp; 
                    end;
                end;                
            end;
            condersps{df,1} = eallersps;                
            fprintf('\nDatafile %s done.',int2str(df));          
            [plotersps{df,1}] = GroupSig(eallersps,s.alpha,alphas(1),'binom');
        end;    
    end;
    fprintf('\nPerforming %s permutation tests as requested...\n',int2str(shuffnum));   
    if size(condersps,2) > 1
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
                    minmask(ff,tt) = distr(ceil(shuffnum*alphas(2)));
                    maxmask(ff,tt) = distr(length(distr) - round(shuffnum*alphas(2)));
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
    for xcond = 1:size(condersps,2)
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
                minmask(ff,tt) = distr(ceil(shuffnum*alphas(2)));
                maxmask(ff,tt) = distr(length(distr) - round(shuffnum*alphas(2)));
            end;
        end;
        tmpmat = reshape(tmpmat,1,size(tmpmat,1)*size(tmpmat,2));
        minmask = reshape(minmask,1,size(minmask,1)*size(minmask,2));
        maxmask = reshape(maxmask,1,size(maxmask,1)*size(maxmask,2));
        tmpmat(find(tmpmat > minmask & tmpmat < maxmask)) = 0;        
        plotersps{3,xcond} = reshape(tmpmat,size(condersps{1,xcond},1),size(condersps{1,xcond},2));
        fprintf('.');          
    end;
    
    if ~isempty(dipset)
        row = size(condersps,2) + 2;
    else
        row = size(condersps,2) + 1;
    end;
    figure;    
    col = size(condersps,1) + 1; pl = 1;
    for df = 1:size(plotersps,1)
        for cond = 1:size(plotersps,2)
            sbplot(row,col,pl); pl = pl+1;
            if strcmp(s.freqscale,'quad')
                quadimagesc(s.times(tms),s.freqs(fr), plotersps{df,cond},[-clim clim]);
            elseif  strcmp(s.freqscale,'log')
                mylogimagesc(s.times(tms),s.freqs(fr), plotersps{df,cond},[-clim clim]);
            end;
            
            set(gca,'ydir','norm'); hold on;
            plot([0 0],[get(gca,'ylim')],'k-');
            if df < 3
                title([dfttls{df},' ',ttls{cond}]);
            else
                title(['Diff ',ttls{cond}]);
            end;   
        end;
    end; 
    cbar;
    if ~isempty(dipset)
        if cond == 1 % if only one condition
            pl = pl+3; % skip a line
        end;
        clustcps{1} = clustcomps;
        ttl = '';
        viewnum=[1,2,3]; % top,sag,rear
        PlotDipoleClusters(dipset,fullpaths,origlist,clustcps,1,row,col,pl,ttl,viewnum,dipcolor);
        pl = pl+length(viewnum);
    end;
    
    textsc([datfiles{1},' and ',datfiles{2},' Interactions; masked by group permutations (p<',num2str(alphas(2)),')'],'title');
    set(gcf,'color','w');
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
