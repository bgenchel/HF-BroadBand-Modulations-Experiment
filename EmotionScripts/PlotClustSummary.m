% uses output of ClustERSP.m and images templates with dipole locations
%
%  PlotClustSummary(datset,paths,activations,kmeanswts,clustcps,conds,inclist,clust,freqs,times);
%
%
%
% conds -- {strings} cell array of strings with condition titles
% inclist -- [integers] list of dimensions included in kmeans clustering
% clust -- integer of the SINGLE cluster to plot. Otherwise will plot all on separate figures

function PlotClustSummary(datset,paths,activations,kmeanswts,clustcps,conds,inclist,clust,freqs,times);
    
    if isempty(clust)
        for comp = 1:length(clustcps)
            figure; row = round(sqrt(length(conds))); col =  ceil(sqrt(length(conds)));pl=1;
            for cond = 1:length(conds) 
                clear pic
                for template = 1:length(inclist)
                    pic(template,:) = activations(inclist(template),(cond-1)*length(freqs)*length(times)+1:cond*length(freqs)*length(times))*kmeanswts(comp,template);
                end;
                pic = mean(pic,1);pic = reshape(pic,length(freqs),length(times));
                lim = max(abs(pic(:)));                sbplot(row,col,pl); 
                ph = imagesc(times,freqs,pic,[-lim lim]); pl = pl+1;hold on;
                set(gca,'ydir','norm');plot([0 0],[get(gca,'ylim')],'k-');
                title(conds{cond});
            end; 
            sbplot(row,col,pl); 
            PlotDipoles('sources.set',paths, clustcps{comp},[],'on','on','off','r','off');
            set(gcf,'color','w');
            manysubj = 0; manycomp = 0;
            for nx = 1:length(clustcps{comp})
                if ~isempty(clustcps{comp}{nx})
                    manysubj = manysubj + 1;
                    manycomp = manycomp + length(clustcps{comp}{nx});
                end;
            end;
            ph = textsc(['Cluster ',int2str(comp),'; ',int2str(manysubj),' Subjects; ',int2str(manycomp),' Components'],'title');
            set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        end;            
    else    
        figure; row = round(sqrt(length(conds))); col =  ceil(sqrt(length(conds)));pl=1;
        for cond = 1:length(conds) 
            clear pic
            for template = 1:length(inclist)
                pic(template,:) = activations(inclist(template),(cond-1)*length(freqs)*length(times)+1:cond*length(freqs)*length(times))*kmeanswts(clust,template);
            end;
            pic = mean(pic,1);pic = reshape(pic,length(freqs),length(times));
            lim = max(abs(pic(:)));            sbplot(row,col,pl); 
            ph = imagesc(times,freqs,pic,[-lim lim]); pl = pl+1;hold on;
            set(gca,'ydir','norm');plot([0 0],[get(gca,'ylim')],'k-');
            title(conds{cond});
        end;            
        sbplot(row,col,pl); 
        PlotDipoles('sources.set',paths, clustcps{clust},[],'on','on','off','r','off');    
    set(gcf,'color','w');
            manysubj = 0; manycomp = 0;
            for nx = 1:length(clustcps{clust})
                if ~isempty(clustcps{clust}{nx})
                    manysubj = manysubj + 1;
                    manycomp = manycomp + length(clustcps{clust}{nx});
                end;
            end;
            ph = textsc(['Cluster ',int2str(clust),'; ',int2str(manysubj),' Subjects; ',int2str(manycomp),' Components'],'title');
            set(ph,'fontsize',14);
set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    end;
