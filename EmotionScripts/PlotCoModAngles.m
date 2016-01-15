% plot histograms of 2D angles of comodulation between dipoles
%
%
%
%
%
%
%
%
%
% btstrap -- [empty, number or cell array] if not empty, will collect
%             the specified number of randomized connections between
%             the plotted dipoles (number), or will collect samples 
%             of randomized connections between all dipoles specified
%             in the IC vectors of the cell array (each cell should 
%             correspond to cells in 'gdcomps'.

function [angles,bilats] = PlotCoModAngles(clustims,gdcomps,fullpaths,bigwts,allbigs,orivec,btstrap);
    

    clear comods justcomps wtsmat1 jcwts
    for cls = 1:length(clustims)
        for nx = 1:length(gdcomps)
            pl = 1; usedim = []; jc = zeros(1,0);jw = zeros(1,0);
            if ~isempty(find(clustims{cls}(:,1) == nx))
                subtemps = clustims{cls}(find(clustims{cls}(:,1) == nx),:);
                for im = 1:size(subtemps,1)
                    currim = subtemps(im,2);
                    if length(find(subtemps(:,2) == currim)) > 1 & ~ismember(currim,usedim)
                        comods{cls}{nx}{pl} = subtemps(find(subtemps(:,2) == currim),3)';
                        ex = find(ismember(allbigs{nx}{currim},subtemps(find(subtemps(:,2) == currim),3)'));
                        
                        wtsmat1{cls}{nx}{pl} = bigwts{nx}{currim}(ex).*orivec{nx}{currim}(ex);
                        usedim = [usedim currim];
                        pl = pl+1;
                    elseif length(find(subtemps(:,2) == currim))==1 & ~ismember(currim,usedim)
                        jc(end+1) = subtemps(find(subtemps(:,2) == currim),3);
                        ex = find(ismember(allbigs{nx}{currim},subtemps(find(subtemps(:,2) == currim),3)'));
                        jw(end+1) = abs(bigwts{nx}{currim}(ex).*orivec{nx}{currim}(ex));
                    end;
                end;
                justcomps{cls}{nx} = jc;
                jcwts{cls}{nx} = jw;
            else
                comods{cls}{nx}= [];
                wtsmat1{cls}{nx} = [];
                %comods{cls}{nx}{pl} = [];
                %wtsmat1{cls}{nx}{pl} = [];
                jcwts{cls}{nx} = [];
                justcomps{cls}{nx} = [];
            end;
        end;
    end;
    
    figure; % need to edit PlotCoModasDipoles to plot every other subplot
    row=length(comods); 
    row=3;

    zoom=1; 

    
    if iscell(btstrap)
        col = 5;
    else
        col = 4; % # bilateral dipole irrelevant for this one
    end;
    
    kurtos = 0; % 1 to plot kurtosis instead of distance measures
    
    pl = 1;
    for clust = 1:length(comods)
        
        viewnum = [1];
        [angles,dists,newdips,realdips,countbil]  = PlotCoModasDipoles(comods{clust},justcomps{clust},fullpaths,'sources1.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust},1,btstrap);
        title(['Clust ',int2str(clust)]);
        %cla;  
        
        pl = pl+1;
        if iscell(angles)
        degs = (angles{1}*180)/pi;
        degs = (angles{2}*180)/pi;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Calculate histograms with angles
        % distance = sqrt((x2 - x1)^2 + (y2-y1)^2+ (z2-z1)^2);
        gry = .7; 
        [hh bins] = hist(angles{1}(1,:),16);
        negs = find(bins < 0);
        pos = length(bins) - length(negs);
        usebins = zeros(1,0); useh = zeros(1,0);
        usebins(end+1:end+pos) = bins(negs(end)+1:end);
        useh(end+1:end+pos) = hh(negs(end)+1:end);
        
        usebins(end+1:end+length(negs)) = pi + bins(1:negs(end)) ;
        useh(end+1:end+length(negs)) = hh(1:negs(end));
        
        usebins(end+1:end+pos) = pi + bins(negs(end)+1:end);
        useh(end+1:end+pos) = hh(negs(end)+1:end);    
        
        usebins(end+1:end+length(negs)) = bins(1:negs(end));
        useh(end+1:end+length(negs)) = hh(1:negs(end)); 
        
        usebins(end+1) = usebins(1);
        useh(end+1) = useh(1);
        %%%  project to x & y
        clear x y
        for g = 1:length(hh)
            x(g) = hh(g) * cos(bins(g));
            y(g) = hh(g) * sin(bins(g));
        end;
        % compute variance in each axis
        realvarx = var(x);
        realvary = var(y);
        % compute kurtosis in each axis
        realkurtx = kurt(x);
        realkurty = kurt(y);
        
        %%%%%%%%%%%%%%%%%%%%%%%%*******************************
        %%%%%%%%%%%%%%%%%%%%%%%%*******************************
        % calculate bootstrap limits:
        clear hh numangs
        for bt = 1:size(angles{2},1)
            [hh(bt,:) bins] = hist(angles{2}(bt,find(angles{2}(bt,:)~=0)),bins);
            numangs(bt) = length(find(angles{2}(bt,:)~=0));
        end;
        clear x y
        %%%  project to x & y
        for bt = 1:size(hh,1)
            for g = 1:size(hh,2)
                x(bt,g) = hh(bt,g) * cos(bins(g));
                y(bt,g) = hh(bt,g) * sin(bins(g));
            end;
        end;
        % compute variance in each axis
        for bt = 1:size(hh,1)
            btvarx(bt) = var(x(bt,:));
            btvary(bt) = var(y(bt,:));
        end;
        % compute kurtosis in each axis
        for bt = 1:size(hh,1)
            btkurtx(bt) = kurt(x(bt,:));
            btkurty(bt) = kurt(y(bt,:));
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        negs = find(bins < 0);
        pos = length(bins) - length(negs);
        for bt = 1:size(angles{2},1)
            useh2 = zeros(1,0);
            useh2(end+1:end+pos) = hh(bt,negs(end)+1:end);
            useh2(end+1:end+length(negs)) = hh(bt,1:negs(end));
            useh2(end+1:end+pos) = hh(bt,negs(end)+1:end);    
            useh2(end+1:end+length(negs)) = hh(bt,1:negs(end));            
            useh2(end+1) = useh2(1);            
            allvals(bt,:) = useh2;            
        end;
        for bns = 1:size(allvals,2)
            distr = allvals(:,bns); distr = sort(distr);
            btvals(1,bns) = distr(1 + round(length(distr) * .01));
            btvals(2,bns) = distr(end-round(length(distr) * .01));
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot the results    %%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%
        %%  Polar plot of connection angles with filled bootstrap area in green
        sbplot(row,col,pl);
        % divide by number of real angles used
        if max(useh/length(angles{1})) > max(btvals(2,:)/mean(numangs))
            ph = polar(usebins,useh/length(angles{1}),'r-'); set(ph,'linewidth',3);hold on;
            hold on;        
            % divide by mean number of angles
            ph = polar(usebins,btvals(2,:)/mean(numangs),'g-'); 
            xval1 = get(ph,'xdata'); yval1 = get(ph,'ydata');set(ph,'linewidth',.01);
            ph = polar(usebins,btvals(1,:)/mean(numangs),'g-');   
            xval2 = get(ph,'xdata'); yval2 = get(ph,'ydata');set(ph,'linewidth',.01);
        else
            % divide by mean number of angles
            ph = polar(usebins,btvals(2,:)/mean(numangs),'g-');  
            hold on;        
            xval1 = get(ph,'xdata'); yval1 = get(ph,'ydata');set(ph,'linewidth',.01);
            ph = polar(usebins,btvals(1,:)/mean(numangs),'g-');   
            xval2 = get(ph,'xdata'); yval2 = get(ph,'ydata');set(ph,'linewidth',.01);
        end;
        f1 = [xval1,xval2(end:-1:1)];
        f2 = [yval1,yval2(end:-1:1)];
        ph = fill(f1,f2,'g');hold on;
        set(ph,'edgealpha',0)
        set(ph,'facealpha',.3)
        ph = polar(usebins,useh/length(angles{1}),'r-'); set(ph,'linewidth',3);hold on;
        pl = pl+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Plot variance in x vs y axes:
        sbplot(row,col,pl); 
        plot(btvarx,btvary,'g.','markersize',15);
        hold on;
        ph =plot(realvarx,realvary,'r*');
        set(ph,'markersize',12);
        set(gca,'ylim',[0 max([btvary realvary])+5]);
        set(gca,'xlim',[0 max([btvarx realvarx])+5]);
        xlabel('Variance in X projection');
        ylabel('Variance in Y projection');
        legend('Bootstrap values','Measured value','Location','southwest');
        % find projection through bt values and project all vals down
        %[pc,eigvec,sv] = runpca([btvarx;btvary])
        pl = pl+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT  how many bilateral ICs in real and bootstrap distributions
        if length(countbil) > 1
            [bih bibin] = hist(countbil{2},20);
            sbplot(row,col,pl)
            bar(bibin,bih,'g');hold on;
            bar(countbil{1},1,2,'r');
            legend('Bootstrap Values','Measured Value');
            title('Numbers of bilateral ICs');
            pl = pl+1;
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%
        % plot connection distance (3D) OR kurtosis
        if clust > 1
        %if kurtos ==0
            for bt = 1:size(dists{2},1)
                [hhb(bt,:) bins] = hist(dists{2}(bt,find(dists{2}(bt,:) ~=0)),16);
            end;
            [hh bins] = hist(dists{1},bins);
            for bns = 1:size(hhb,2)
                distr = hhb(:,bns); distr = sort(distr);
                btdist(1,bns) = distr(round(length(distr) * .01));
                btdist(2,bns) = distr(end-(round(length(distr) * .01)+1));
            end;
            
            sbplot(row,col,pl)
            f1 = [bins,bins(end:-1:1)];
            f2 = [btdist(1,:),btdist(2,end:-1:1)];
            ph = fill(f1,f2,'g');hold on;
            set(ph,'edgealpha',0)
            set(ph,'facealpha',.3)

            plot(bins,hh,'r-','linewidth',2);
            plot(bins,hh,'r.','markersize',24);
            
            set(gca,'xlim',[0 round(bins(end)+5)]);
            legend('Bootstrap Limits','Measured Values');
            
            xlabel('Equiv. Dipole Distance (mm)');
            ylabel('Number of Observations');
            set(gca,'box','off');
        else
            sbplot(row,col,pl); 
            plot(btkurtx,btkurty,'g.','markersize',15);
            hold on;
            ph =plot(realkurtx,realkurty,'r*');
            set(ph,'markersize',12);
            set(gca,'ylim',[min([btkurty realkurty])-.1 0]);
            xlabel('Kurtosis in X projection');
            ylabel('Kurtosis in Y projection');
            legend('Bootstrap values','Measured value','Location','northeast');
            
            % find projection through bt values and project all vals down
            %[pc,eigvec,sv] = runpca([btkurtx;btkurty],1)            
        end;          
        pl = pl+1;
        end;
        bilats{clust} = countbil;
    end;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    
