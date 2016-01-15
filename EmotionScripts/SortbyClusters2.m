% takes the IMs found in Emo Space and quantifies for each spectral cluster:
%
%
% [clsmems,nims] = SortbyClusters(subjfacs,row,col,pl)
%
%
% subjfacs -- [matrix] of nfacs x subj,im indexes for a given emo space dim
% row -- number of rows in the current figure
% col -- number of columns in the current figure
% pl -- sub-plot number to plot into
% if row is [], will pop a new figure
% OUTPUTS:
% clsmems -- [vector] giving the number of IMs from 'subjfacs' in each 
%                     cluster (Delta,Theta,Alpha,Beta,Gamma) 
% nims -- [vector] giving the number of total IMs in each cluster 
%                 (Delta,Theta,Alpha,Beta,Gamma)

function [clsmems,maxmems,allidx2,iclistp,iclistn,templatesn,templatesp,maxy] = SortbyClusters2(subjfacs,row,col,pl,maxy)
    
    
    maxnsubj = max(subjfacs(:,1));
    
    strs = {
    'load /data/common2/emotion/DeltaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/ThetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/AlphaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/BetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/GammaClust.mat finaltempls finalidx finalmean freqs'};
    
    clear emovec
    for nx = 1:max(subjfacs(:,1))
        emovec{nx} = subjfacs(find(subjfacs(:,1) == nx),:); 
    end;
    
    allidx2 = cell(length(strs),4); pp = 1; ss = 1;
    for s = 1:length(strs) % spectral clusters
        eval(strs{s}); 
        facvec = cell(1,length(finalidx));
        for cls = 1:length(finalidx)
            maxmems(1,ss) = size(finalidx{cls},1);ss = ss+1;
            facvec{cls} = cell(1,maxnsubj);
            for nx = 1:max(finalidx{cls}(:,1))
                facvec{cls}{nx} = unique(finalidx{cls}(find(finalidx{cls}(:,1) == nx),2));
            end;
        end;        
        
        tmpidx = []; 
        numpos = 0;numneg = 0;
        for cls = 1:length(facvec)% spectral sub-category
            numpos = 0;numneg = 0;
            tmpidx2 = [];
            for nx = 1:length(emovec)
                % these will find only how many IMs were involved (regardless of how many ICs)
                numpos = numpos+length(find(ismember(emovec{nx}(:,2),facvec{cls}{nx})));
                numneg = numneg+length(find(ismember(-1*emovec{nx}(:,2),facvec{cls}{nx})));
                clustidxs = finalidx{cls}(find(finalidx{cls}(:,1)==nx),:);
                tmat = finaltempls{cls}(find(finalidx{cls}(:,1)==nx),:);
                ilp = []; iln = [];ttmatn = [];ttmatp = [];
                for im = 1:size(emovec{nx},1)
                    % this will find all instances of the im (all ICs)
                    subidxs = clustidxs(find(abs(clustidxs(:,2))==abs(emovec{nx}(im,2))),:);
                    ttmat = tmat(find(abs(clustidxs(:,2))==abs(emovec{nx}(im,2))),:);
                    if ~isempty(subidxs)
                        for ic = 1:size(subidxs,1)
                            if emovec{nx}(im,2) < 0 & subidxs(ic,2) > 0
                                iln = [iln subidxs(ic,3)'];
                                ttmatn = [ttmatn; ttmat(ic,:)*-1];
                            elseif emovec{nx}(im,2) > 0 & subidxs(ic,2) < 0
                                iln = [iln subidxs(ic,3)'];                       
                                ttmatn = [ttmatn; ttmat(ic,:)*-1];
                            elseif emovec{nx}(im,2) > 0 & subidxs(ic,2) > 0
                                ilp = [ilp subidxs(ic,3)'];                       
                                ttmatp = [ttmatp; ttmat(ic,:)];
                            elseif emovec{nx}(im,2) < 0 & subidxs(ic,2) < 0
                                ilp = [ilp subidxs(ic,3)'];                       
                                ttmatp = [ttmatp; ttmat(ic,:)];
                            end;
                            if size(emovec{nx},2) == 3 % for the decile approach
                                subidxs(:,4) = emovec{nx}(im,3);
                            end;
                        end;
                    end;
                    %tmpidx = [tmpidx;subidxs];
                    tmpidx2 = [tmpidx2;subidxs];
                end;
                ilnxp{nx} =  ilp;
                ilnxn{nx} =  iln;
                tmplsn{nx} = ttmatn;
                tmplsp{nx} = ttmatp;
            end;          
            iclistp{pp} =  ilnxp;
            iclistn{pp} =  ilnxn;
            templatesn{pp} = tmplsn;
            templatesp{pp} = tmplsp;
            allidx2{pp} = tmpidx2;
            %allidx{pp} = tmpidx;
            clsmems(1,pp) =  numpos;
            clsmems(2,pp) =  numneg*-1;pp = pp+1;            
        end;
        if isempty(maxy)
            maxy = max(abs(clsmems(:)))+4;
        end;
    end;
    cols = jet(length(clsmems));
    if isempty(row)
        figure; 
        %ph = plot([0 length(clsmems)+1],[100 100],'r-'); hold on;
        for s = 1:length(clsmems)
            ph = bar([s-.2 s+.2],clsmems(:,s)/maxmems(s)); hold on;
            set(ph,'facecolor',cols(s,:));
            set(ph,'edgecolor',cols(s,:));
        end;
        set(gca,'xlim',[0 length(clsmems)+1]); %set(gca,'ylim',[-maxy maxy]);
        set(gca,'xtick',[1:length(clsmems)]);
        set(gca,'xticklabel',{'Delta','Theta','Alpha','Beta','Hi/Low','BB'});
    else
        sbplot(row,col,pl)
        %ph = plot([0 length(clsmems)+1],[100 100],'r-'); hold on;
        for s = 1:length(clsmems)            
            ph = bar([s-.2 s+.2],clsmems(:,s)/maxmems(s)); hold on;
            set(ph,'facecolor',cols(s,:));
            set(ph,'edgecolor',cols(s,:));
        end;
        set(gca,'xlim',[0 length(clsmems)+1]);%set(gca,'ylim',[-maxy maxy]);
        set(gca,'xtick',[1:length(clsmems)]);
        set(gca,'xticklabel',{'Delta1','Delta2','Theta1','Theta2','Alpha1','Alpha2','Alpha3','Beta1','Beta2','Beta3','HiLow','BB','Peaked'});
   end;        
