


function PlotClusterIMs(savedat,fullpaths,matchsets,formatch,subjclusts,subjbigs,mnchange,plotopt,frqlim);
    

    s = load([fullpaths{1},[savedat,'Day',int2str(1)],'.mat']); % load variables (any subj)
    fr = find(s.freqs>=frqlim(1) & s.freqs<=frqlim(2));

    styles = {'-','--',':','-.','-','--'}; % for clusters with more than one ic from a subj/day

    for nx = 1:length(fullpaths)
        clsims = cell(1,length(subjclusts{nx})); 
        figure;pl = 1; 
        if strcmp(mnchange,'on') 
            row = 6; col = 5;
        else
            row = 9; col = 9;
        end;
        ttlsoff =  zeros(1,4); ttlsoff1 =  zeros(1,4); ttlsoff2 = zeros(1,4); mnmap = zeros(67,67,0);
        for cls = 1:length(subjclusts{nx})
            if pl > row*col
                textsc(['Subject ',int2str(nx)],'title');
                figure; pl = 1;
                ttlsoff =  zeros(1,4);ttlsoff1 =  zeros(1,4); ttlsoff2 = zeros(1,4);
            end;
            nmems = numel(subjclusts{nx}{cls}{1})+numel(subjclusts{nx}{cls}{2})+numel(subjclusts{nx}{cls}{3})+numel(subjclusts{nx}{cls}{4});
            nosig = 0;
            if nmems > 2
                sbplot(row,col,pl);  pl = pl + 1;
                [mnmap(:,:,end+1)] = MeanClusterMap(matchsets,{fullpaths{nx},fullpaths{nx},fullpaths{nx},fullpaths{nx}}, subjclusts{nx}{cls});
                ph = title(['Clust ',int2str(cls),' (',int2str(nmems),' ICs)']);
                for day = 1:length(subjclusts{nx}{cls})
                    ci = [];
                    ics = subjclusts{nx}{cls}{day};% ICs in the subj cluster across sessions
                    ims = formatch{nx}{day}; % IMs w/high condition diff
                    for im = 1:length(ims)
                        if ~isempty(find(ismember(ics,subjbigs{day}{nx}{ims(im)})))
                            ci = [ci ims(im)];
                        end;
                    end;
                    clsims{cls}{day} = ci; % all sig IMs from all clustered ICs for that day
                    if ~isempty(ics) && ~isempty(clsims{cls}{day})
                        s = load([fullpaths{nx},[savedat,'Day',int2str(day)],'.mat']);
                        sph=floatread([fullpaths{nx},[savedat,'Day',int2str(day)],'.sph'],[s.pcs s.pcs],[],0);
                        wts=floatread([fullpaths{nx},[savedat,'Day',int2str(day)],'.wts'],[s.pcs s.pcs],[],0);
                        icamatall = floatread([fullpaths{nx},[savedat,'Day',int2str(day)],'.fdt'],[s.pcs s.numframes],[],0);
                        ws = wts*sph;    acts = ws*icamatall;    winv = pinv(ws);
                        clear wts sph ws icamatall
                        speceig = floatread([fullpaths{nx},s.eigfile],[length(s.rowmeans) s.pcs],[],0);
                        specwts = speceig*winv;
                        winv = specwts;
                        if pl > row*col
                            textsc(['Subject ',subjs{nx},' (green = PRE, red = POST)'],'title'); figure; pl=1;
                        end;
                        if strcmp(mnchange,'on') % plot back-proj IMs as change from mean spec
                            sbplot(row,col,pl);
                            for c = 1:length(ics)
                                tmpls1 = zeros(0,length(s.freqs));
                                tmpls2 = zeros(0,length(s.freqs));
                                rcp = find(ismember(s.complist,ics(c)));
                                backproj1 = winv(1:sum(s.dstrials(1:2)),ims) * acts(ims,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                                backproj1 = s.meanpwr(find(ismember(s.complist,ics(c))),:)+mean(backproj1,1);
                                tmpls1(end+1,:) = backproj1 ;
                                backproj2 = winv(sum(s.dstrials(1:2))+1:end,ims) * acts(ims,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                                backproj2 = s.meanpwr(find(ismember(s.complist,ics(c))),:)+mean(backproj2,1);
                                tmpls2(end+1,:) = backproj2 ;
                                ph = quadplot(s.freqs(fr),tmpls1(:,fr),2,'g'); hold on;
                                set(ph,'linestyle',styles{c});
                                ph = quadplot(s.freqs(fr),tmpls2(:,fr),2,'r');
                                set(ph,'linestyle',styles{c});set(gca,'xgrid','on');
                                if pl <= (row-1)*col
                                    set(gca,'xticklabel',[]);
                                end;
                                if ttlsoff(day) == 0
                                    ph=title(['Day',int2str(day)]);
                                    ttlsoff(day) = 1; pl = pl+1;
                                end;
                            end;
                        else
                            condwts= zeros(0,2);
                            tmpls1 = zeros(0,length(s.freqs));
                            tmpls2 = zeros(0,length(s.freqs));
                            for c = 1:length(ics)                                
                                for i = 1:length(clsims{cls}{day})
                                    if ismember(ics(c),subjbigs{day}{nx}{clsims{cls}{day}(i)})% only if IC is involved in current IM
                                        nowim = clsims{cls}{day}(i);
                                        condwts(end+1,1) = median(winv(1:sum(s.dstrials(1:2)),nowim));
                                        condwts(end,2) = median(winv(sum(s.dstrials(1:2))+1:end,nowim));
                                        rcp = find(ismember(s.complist,ics(c)));
                                        backproj1 = winv(1:sum(s.dstrials(1:2)),nowim) * acts(nowim,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                                        backproj1 = mean(backproj1,1); tmpls1(end+1,:) = backproj1 ;
                                        backproj2 = winv(sum(s.dstrials(1:2))+1:end,nowim) * acts(nowim,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                                        backproj2 = mean(backproj2,1);tmpls2(end+1,:) = backproj2 ;
                                    end;
                                end;
                            end;
                            sbplot(row,col,pl); ph = quadplot(s.freqs(fr),tmpls1(:,fr),2);
                            pl = pl + 1;hold on;    ph = plot([get(gca,'xlim')],[0 0],'k-');
                            set(gca,'xgrid','on');
                            if pl <= (row-1)*col
                                set(gca,'xticklabel',[]);
                            end;
                            if ttlsoff1(day) == 0
                                ph=title(['Day',int2str(day),' Pre']);
                                set(ph,'color','g');ttlsoff1(day) = 1;
                            end;
                            sbplot(row,col,pl); ph = quadplot(s.freqs(fr),tmpls2(:,fr),2);
                            pl = pl + 1;hold on;    ph = plot([get(gca,'xlim')],[0 0],'k-');
                            if pl <= (row-1)*col
                                set(gca,'xticklabel',[]);
                            end;
                            if ttlsoff2(day) == 0
                                ph=title(['Day',int2str(day),' Post']);
                                set(ph,'color','r');ttlsoff2(day) = 1;
                            end;
                            set(gca,'xgrid','on');
                        end;                    
                    else
                        nosig = nosig+1;
                        if strcmp(mnchange,'on')
                            pl = pl+1;
                        else
                            pl = pl+2;
                        end;
                    end;
                end;
                if nosig == length(subjclusts{nx}{cls}) % if only one day represented
                    set(ph,'visible','off');cla; 
                    if strcmp(mnchange,'on')
                        pl=pl-5;
                    else
                        pl = pl-9;
                    end;
                    mnmap(:,:,end)=[];
                else
                    %textsc(['Subject ',subjs{nx},' (green = PRE, red = POST)'],'title');
                end;
            end;
        end;
        textsc(['Subject ',int2str(nx),' (green = PRE, red = POST)'],'title');
        EEG = pop_loadset( 'filename',matchsets{day} , 'filepath',fullpaths{nx}); 
        isnotnans = find(~isnan(mnmap(:,:,1)));
        plotmaps{nx} = mnmap;
        mnmap = reshape(mnmap,[size(mnmap,1)*size(mnmap,2),size(mnmap,3)]);
        mnmap = mnmap(isnotnans,:);
        clustermaps{nx} = mnmap; clusterchans{nx} = EEG.chanlocs;
    end;
