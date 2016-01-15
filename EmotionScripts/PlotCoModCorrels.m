% takes 2D or 3D from selected subjects, as well as cluster info
% to plot correlation relationships between clusters
%
%
%
%
%
% clustlabels -- [cell array] of strings referring to the clusters that will
%                be compared pairwise in the output plots.
% condlabels -- [cell array] of strings referring to the 3rd dimension of input
%               matrix, which are different conditions. ([] if 2D input)
% plottype -- ['corr',varcorr','mi'] plots correlations, variance correlations, 
%             and mutual information, respectively.
%

function [savecorrs,keeppairs,newlabels,tth,comppairs] = PlotCoModCorrels(corr,bootstats,subjlist,clustlabels,condlabels,allvec,plottype);
    
    
    anovalpha = .00001; % anova1 significance (across conditions)
    ttestalpha = .01; % ttest sig ( vs bootstrap values)
    
    if ~exist('plottype')
        plottype = 'corr';
    elseif isempty(plottype)
        plottype = 'corr';
    end;
    keeppairs = []; savecorrs = []; newlabels = [];tth = [];comppairs = [];
    if length(size(corr{subjlist(1)})) < 3 % if 2D (not divided by condition)
        figure; row =3; col = 4; pl = 1;
        for cls1 = 1:length(allvec)-1
            for cls2 = cls1+1:length(allvec)
                currcorrs = []; lim = zeros(0,2);% for bootstrap vals
                for nx = 1:length(allvec{cls1})
                    if ~isempty(allvec{cls1}{nx}) & ~isempty(allvec{cls2}{nx})
                        cls1ims = unique(allvec{cls1}{nx}');
                        cls2ims = unique(allvec{cls2}{nx}');
                        intsc = intersect(cls1ims,cls2ims);
                        if ~isempty(intsc)
                            for in = 1:length(intsc)
                                cls1ims(find(cls1ims == intsc(in))) = [];
                                cls2ims(find(cls2ims == intsc(in))) = [];
                            end;
                        end;
                        if ~isempty(cls1ims) & ~isempty(cls2ims)
                            for im1 = 1:length(cls1ims)
                                for im2 = 1:length(cls2ims)
                                    if cls1ims(im1) < cls2ims(im2)
                                        if strcmp(plottype,'corr')% take abs value of corr:
                                            currcorrs = [currcorrs corr{nx}(cls1ims(im1),cls2ims(im2))];
                                            if ~isempty(bootstats) & cls2ims(im2)~=cls1ims(im1)
                                                lim(end+1,:) = [abs(bootstats{nx}(cls1ims(im1),cls2ims(im2),1)) abs(bootstats{nx}(cls1ims(im1),cls2ims(im2),2))];
                                            end;
                                        elseif strcmp(plottype,'varcorr') | strcmp(plottype,'mi')
                                            currcorrs = [currcorrs corr{nx}(cls1ims(im1),cls2ims(im2))];
                                            if ~isempty(bootstats) & cls2ims(im2)~=cls1ims(im1)
                                                lim(end+1,:) = [bootstats{nx}(cls1ims(im1),cls2ims(im2),1) bootstats{nx}(cls1ims(im1),cls2ims(im2),2)];
                                            end;                                            
                                        end;                                        
                                        
                                    else
                                        if strcmp(plottype,'corr')
                                            currcorrs = [currcorrs corr{nx}(cls2ims(im2),cls1ims(im1))];
                                            if ~isempty(bootstats) & cls2ims(im2)~=cls1ims(im1)
                                                lim(end+1,:) = [bootstats{nx}(cls2ims(im2),cls1ims(im1),1) bootstats{nx}(cls2ims(im2),cls1ims(im1),2)];
                                            end;
                                        elseif strcmp(plottype,'varcorr') | strcmp(plottype,'mi')
                                            currcorrs = [currcorrs corr{nx}(cls2ims(im2),cls1ims(im1))];
                                            if ~isempty(bootstats) & cls2ims(im2)~=cls1ims(im1)
                                                lim(end+1,:) = [bootstats{nx}(cls2ims(im2),cls1ims(im1),1) bootstats{nx}(cls2ims(im2),cls1ims(im1),2)];
                                                
                                            end;
                                        end;
                                    end; 
                                end;
                            end;
                        end;
                    end;
                end;
                lim = max(abs(lim),[],1); lim(1) = lim(1)*-1;
                if pl>row*col
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    textsc('Histograms of correlations between clusters of IMs within subject (red/blue=bootstrap limits)','title');
                    figure; pl=1;
                end;
                sbplot(row,col,pl)
                hist(currcorrs,40);hold on;
                if ~strcmp(plottype,'mi')
                    set(gca,'xlim',[-.75 .75]); % set viewing limits for correlation
                end;
                plot([lim(1) lim(1)],[get(gca,'ylim')],'b-','linewidth',1.5);
                plot([lim(2) lim(2)],[get(gca,'ylim')],'r-','linewidth',1.5);
                plot([0 0],[get(gca,'ylim')],'g-','linewidth',1.5);
                ph = plot([median(currcorrs) median(currcorrs)],[get(gca,'ylim')],'y-','linewidth',2);
                set(ph,'color',[1 .5 0]);
                title([clustlabels{cls1},'/',clustlabels{cls2}]);pl = pl+1;
            end;           
        end;
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        textsc('Histograms of correlations between clusters of IMs within subject (red/blue=bootstrap limits)','title');
    else % if 3D (and divided by expt'al condition)
% $$$         if length(size(corr{subjlist(1)})) > 2 % if 3D
% $$$             corrmeans = [];
% $$$             corrstds = [];
% $$$             for nxx = 1:length(subjlist)
% $$$                 nx = subjlist(nxx);
% $$$                 corrmeans = [corrmeans, unique(reshape(mean(abs(corr{nx}),3),1,size(corr{nx},1)*size(corr{nx},2)))];
% $$$                 corrstds =  [corrstds,  unique(reshape(std(abs(corr{nx}),0,3),1,size(corr{nx},1)*size(corr{nx},2)))];
% $$$             end;
% $$$         end;
        cc = 1;
        for cls1 = 1:length(allvec)-1
            for cls2 = cls1+1:length(allvec)  
                currcorrs = [];lim = zeros(0,15,2);
                for nxx = 1:length(subjlist)
                    nx = subjlist(nxx);
                    cls1ims = unique(allvec{cls1}{nx}');
                    cls2ims = unique(allvec{cls2}{nx}');
                    intsc = intersect(cls1ims,cls2ims);
                    if ~isempty(intsc)
                        for in = 1:length(intsc)
                            cls1ims(find(cls1ims == intsc(in))) = [];
                            cls2ims(find(cls2ims == intsc(in))) = [];
                        end;
                    end;
                    if ~isempty(cls1ims) & ~isempty(cls2ims)
                        for im1 = 1:length(cls1ims)
                            for im2 = 1:length(cls2ims)
                                if cls1ims(im1) < cls2ims(im2)
                                    if strcmp(plottype,'corr')% for correlation
                                        currcorrs = [currcorrs; abs(squeeze(corr{nx}(cls1ims(im1),cls2ims(im2),:)))'];
                                        if ~isempty(bootstats)
                                            lim(end+1,:,:) = [abs(squeeze(bootstats{nx}(cls1ims(im1),cls2ims(im2),1,:))), abs(squeeze(bootstats{nx}(cls1ims(im1),cls2ims(im2),2,:)))];
                                        end;
                                    elseif strcmp(plottype,'varcorr') | strcmp(plottype,'mi')
                                        currcorrs = [currcorrs; squeeze(corr{nx}(cls1ims(im1),cls2ims(im2),:))'];
                                        if ~isempty(bootstats)
                                            lim(end+1,:,:) = [squeeze(bootstats{nx}(cls1ims(im1),cls2ims(im2),1,:)), squeeze(bootstats{nx}(cls1ims(im1),cls2ims(im2),2,:))];
                                        end;
                                    end;
                                else
                                    if strcmp(plottype,'corr')% for correlation
                                        currcorrs = [currcorrs; abs(squeeze(corr{nx}(cls2ims(im2),cls1ims(im1),:)))'];
                                        if ~isempty(bootstats)
                                            lim(end+1,:,:) = [abs(squeeze(bootstats{nx}(cls2ims(im2),cls1ims(im1),1,:))), abs(squeeze(bootstats{nx}(cls2ims(im2),cls1ims(im1),2,:)))];
                                        end;
                                    elseif strcmp(plottype,'varcorr') | strcmp(plottype,'mi')
                                        currcorrs = [currcorrs; squeeze(corr{nx}(cls2ims(im2),cls1ims(im1),:))'];
                                        if ~isempty(bootstats)
                                            lim(end+1,:,:) = [squeeze(bootstats{nx}(cls2ims(im2),cls1ims(im1),1,:)), squeeze(bootstats{nx}(cls2ims(im2),cls1ims(im1),2,:))];% 4th Dim is condition
                                        end;
                                    end;
                                end; 
                            end;
                        end;
                    end;
                end;  
                
                newlabels{cc} = [clustlabels{cls1},'/',clustlabels{cls2}];
                savecorrs{cc} = currcorrs;
                if ~isempty(bootstats) & ~isempty(savecorrs{cc})
                    if strcmp(plottype,'corr')% for correlation
                        saveboots{cc} = mean(lim,3); % because I'm only looking now at abs(correlations)
                        for c = 1:size(lim,2)
                            nsigs(1,c) = length(find(savecorrs{cc}(:,c) > saveboots{cc}(:,c)));
                        end;
                        tth{cc} = find(nsigs >= (round(size(savecorrs{cc},1))/4)*3);
                        tth{cc} = find(nsigs == size(savecorrs{cc},1) );  % require all                   
                    else      
                        saveboots{cc} = lim;
                        for c = 1:size(lim,2)
                            nsigs(1,c) = length(find(savecorrs{cc}(:,c) > saveboots{cc}(:,c,1) & savecorrs{cc}(:,c) < saveboots{cc}(:,c,2)));
                        end;
                         %tth{cc} = find(nsigs <= round(size(savecorrs{cc},1))/4);                        
                         tth{cc} = find(nsigs == size(savecorrs{cc},1) ); % require all                     
                    end;
                    %[hh ttp{cc}] = ttest(savecorrs{cc},saveboots{cc},ttestalpha); % see which are non-random
                    %tth{cc} = find(hh); % get condition idx
                end;
                if ~isempty(savecorrs{cc})
                    [P,ANOVATAB,STATS] = anova1(savecorrs{cc},condlabels,'off'); %close;
                    comp = multcompare(STATS,'alpha',anovalpha,'ctype','bonferroni'); %close;
                    comppairs = [];
                    for cp = 1:size(comp,1)
                        if length(find(comp(cp,[3,5]) == abs(comp(cp,[3,5])))) ~= 1 %(not straddling 0=sig)
                            comppairs = [comppairs;[comp(cp,[1:2])]];
                        end;
                    end;
                    keeppairs{cc} = comppairs;
                    tagemos{cc} = unique(comppairs);
                else
                    keeppairs{cc} = [];
                    tagemos{cc} = [];
                end;
                cc = cc+1;
            end;
        end;
        fprintf('\nPlotting inter-cluster measures...\n');
        cols = jet(size(corr{subjlist(1)},3));
        figure; row = 4; col = 4; pl = 1;
        for clss = 1:length(savecorrs)
            if ~isempty(savecorrs{clss})
                if pl > row*col
                    if strcmp(plottype,'corr')% for correlation
                        textsc(['Abs correlations bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
                    elseif strcmp(plottype,'varcorr')
                        textsc(['Variance correl bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
                    elseif strcmp(plottype,'mi')
                        textsc(['Mutual information bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
                    end;
                    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
                    figure; pl=1;
                end;
                sbplot(row,col,pl); pl = pl + 1;
                for ds = 1:size(savecorrs{clss},2)
                    ph = plot(ds,savecorrs{clss}(:,ds),'.','markersize',4); hold on;
                    set(ph,'color',cols(ds,:));                
                    ph = plot(ds,mean(savecorrs{clss}(:,ds)),'k.','markersize',10);
                    if ~isempty(find(ismember(ds,tagemos{clss})))
                        ph = plot(ds,[.9 .9],'*','markersize',8);
                        set(ph,'color',[.7 0 .9]);
                    end;
                    if ~isempty(bootstats)
                        if ~isempty(find(ismember(ds,tth{clss})))
                            ph = plot(ds,[.8 .8],'m*','markersize',8);
                        end;
                    end;
                end;
                ph = plot(mean(savecorrs{clss},1),'k-');
                plot([get(gca,'xlim')],[0 0],'k-'); %set(gca,'ylim',[0 1]);% for correlation assessment
                ph = plot(squeeze(median(saveboots{clss},1)),'r--'); % plot permut value means
                set(gca,'xlim',[0 size(savecorrs{clss},2)+1]);
                set(gca,'xticklabel',[]);
                title(newlabels{clss});
            end;
        end;        
        if strcmp(plottype,'corr')% for correlation
            textsc(['Abs correlations bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
        elseif strcmp(plottype,'varcorr')
            textsc(['Variance correl bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
        elseif strcmp(plottype,'mi')
            textsc(['Mutual information bet IM wts across conds (x-axis); purple*=diff bet conds by ANOVA, p<',num2str(anovalpha),'; magenta*: Diff by permut test(p<.01)'],'title');
        end;
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    end;
