% create a density plot of comodulated pairs to see mean angles of cxn
%
%
%
%
%
%
%

function CoModDensity(clustims,gdcomps,fullpaths,bigwts,allbigs,orivec)
    

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
                comods{cls}{nx}{pl} = [];
                wtsmat1{cls}{nx}{pl} = [];
                jcwts{cls}{nx} = [];
                justcomps{cls}{nx} = [];
            end;
        end;
    end;
    keyboard
    figure; % need to edit PlotCoModasDipoles to plot every other subplot
    row=length(comods); 
    zoom=1; viewnum=2;col=viewnum+2; pl = 1;
    for clust = 1:length(comods)
        [angles,newdips,realdips]  = PlotCoModasDipoles(comods{clust},justcomps{clust},fullpaths,'sources1.set',row,col,pl,zoom,0,viewnum,wtsmat1{clust},jcwts{clust});
        pl = pl+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  plot newdips as 3D contour
        % plot 2d: x vs y (col 1 (left/right) and 2 (ant-post); 5 and 6)
        % distance = sqrt((x2 - x1)^2 + (y2-y1)^2);
% $$$     sbplot(row,col,pl)
% $$$     for prs = 1:size(realdips,1)
% $$$         dd = sqrt((realdips(prs,4) - realdips(prs,1))^2 + (realdips(prs,5)-realdips(prs,2))^2);
% $$$         x(prs) = cos(angles(1,prs)) * dd;%takes radians, dd = hyp length
% $$$         y(prs) = sin(angles(1,prs)) * dd;        
% $$$         plot([0 x],[0 y],'r-','linewidth',3);hold on; 
% $$$         text(x, y,[' ',int2str(prs)]);         
% $$$     end;
% $$$     pl = pl+1;
        sbplot(row,col,pl)
        for prs = 1:size(realdips,1)
            ph = plot([realdips(prs,1) realdips(prs,4)],[realdips(prs,2) realdips(prs,5)],'r-'); hold on;
            ph = plot(realdips(prs,1), realdips(prs,2),'b.'); hold on;
            %text(realdips(prs,1), realdips(prs,2),[' ',int2str(prs)]); 
            ph = plot(realdips(prs,4), realdips(prs,5),'b.'); hold on;
        end;
        pl = pl+2;
% $$$     sbplot(row,col,pl)
% $$$     for prs = 1:size(newdips,1)
% $$$         ph = plot([newdips(prs,1) newdips(prs,4)],[newdips(prs,2) newdips(prs,5)],'r-'); hold on;
% $$$         ph = plot(newdips(prs,1), newdips(prs,2),'b.'); hold on;
% $$$         ph = plot(newdips(prs,4), newdips(prs,5),'b.'); hold on;
% $$$         text(newdips(prs,1), newdips(prs,2),[' ',int2str(prs)]); 
% $$$     end;
        
        sbplot(row,col,pl)
        for prs = 1:size(realdips,1)
            ph = plot([realdips(prs,2) realdips(prs,5)],[realdips(prs,3) realdips(prs,6)],'r-'); hold on;
            ph = plot(realdips(prs,2), realdips(prs,3),'b.'); hold on;
            %text(realdips(prs,1), realdips(prs,2),[' ',int2str(prs)]); 
            ph = plot(realdips(prs,5), realdips(prs,6),'b.'); hold on;
        end;
        pl = pl+1;
        
        
    end;
