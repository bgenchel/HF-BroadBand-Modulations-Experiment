% takes the IMs found in Emo Space and quantifies for each spectral cluster:
%
%
% [clsmems,nims] = SortbyClusters(subjfacs,row,col,pl)
%
%
% subjfacs -- [cell array] of [subj im] indicators for a given emotion
% row -- number of rows in the current figure
% col -- number of columns in the current figure
% pl -- sub-plot number to plot into
% if row is [], will pop a new figure
% OUTPUTS:
% clsmems -- [vector] giving the number of IMs from 'subjfacs' in each 
%                     cluster (Delta,Theta,Alpha,Beta,Gamma) 
% nims -- [vector] giving the number of total IMs in each cluster 
%                 (Delta,Theta,Alpha,Beta,Gamma)

function [clsmems,nims,allidx,allidx2,maxy] = SortbyClusters(subjfacs,row,col,pl,maxy)
    
    
    for f = 1:length(subjfacs)
        maxnsubj(1,f) = max(subjfacs{f}(:,1));
    end;
    maxnsubj = max(maxnsubj);
    
    strs = {
    'load /data/common2/emotion/DeltaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/ThetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/AlphaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/BetaClust.mat finaltempls finalidx finalmean freqs',
    'load /data/common2/emotion/GammaClust.mat finaltempls finalidx finalmean freqs'};
    
    nims = zeros(1,length(strs)); allidx2 = cell(length(strs),4);
    for s = 1:length(strs) % spectral clusters
        eval(strs{s}); 
        facvec = cell(1,length(finalidx));
        for cls = 1:length(finalidx)
            facvec{cls} = cell(1,maxnsubj);
            for nx = 1:max(finalidx{cls}(:,1))
                facvec{cls}{nx} = unique(finalidx{cls}(find(finalidx{cls}(:,1) == nx),2));
                nims(1,s) = nims(1,s)+ length(unique(finalidx{cls}(find(finalidx{cls}(:,1) == nx),2)));
            end;
        end;
        
        clear emovec
        for f = 1:length(subjfacs) %  for each emo space dim, subj idxes for one emo
            for nx = 1:max(subjfacs{f}(:,1))
                emovec{f}{nx} = subjfacs{f}(find(subjfacs{f}(:,1) == nx),:); 
            end;
        end;
        tmpidx = []; 
        numpos = 0;numneg = 0;
        for f = 1:length(emovec) % 'emo space' dimensions involved 
            for cls = 1:length(facvec)% spectral sub-category
                if s == length(strs)
                    numpos = 0;numneg = 0;
                end;
                tmpidx2 = [];
                for nx = 1:length(emovec{f})
                    % these will find only how many IMs were involved (regardless of how many ICs)
                    numpos = numpos+length(find(ismember(emovec{f}{nx}(:,2),facvec{cls}{nx})));
                    numneg = numneg+length(find(ismember(-1*emovec{f}{nx}(:,2),facvec{cls}{nx})));
                    clustidxs = finalidx{cls}(find(finalidx{cls}(:,1)==nx),:);
                    for im = 1:size(emovec{f}{nx},1)
                        % this will find all instances of the im (all ICs)
                        subidxs = clustidxs(find(clustidxs(:,2)==abs(emovec{f}{nx}(im,2))),:);
                        if emovec{f}{nx}(im,2) < 0
                            subidxs(:,2) = subidxs(:,2)*-1;
                        end;
                        if size(emovec{f}{nx},2) == 3 % for the decile approach
                            subidxs(:,4) = emovec{f}{nx}(im,3);
                        end;
                        tmpidx = [tmpidx;subidxs];
                        tmpidx2 = [tmpidx2;subidxs];

                    end;
                end;
                allidx2{s,cls} = tmpidx2;
                if s == length(strs)
                    clsmems(1,s+(cls-1)) = numpos;
                    clsmems(2,s+(cls-1)) = numneg*-1;
                end;
            end;
        end;
        if s ~= length(strs)
            clsmems(1,s) = numpos;
            clsmems(2,s) = numneg*-1;
        end;         
        allidx{s} = tmpidx;
    end;  
    if isempty(maxy)
        maxy = max(clsmems)+4;
    end;
    cols = jet(length(clsmems));
    if isempty(row)
        figure; 
        %ph = plot([0 length(clsmems)+1],[100 100],'r-'); hold on;
        for s = 1:length(clsmems)
            ph = bar([s-.2 s+.2],clsmems(:,s)); hold on;
            set(ph,'facecolor',cols(s,:));
            set(ph,'edgecolor',cols(s,:));
        end;
        set(gca,'xlim',[0 length(clsmems)+1]); %set(gca,'ylim',[0 maxy]);
        set(gca,'xtick',[1:length(clsmems)]);
        set(gca,'xticklabel',{'Delta','Theta','Alpha','Beta','Hi/Low','Gamma'});
    else
        sbplot(row,col,pl)
        %ph = plot([0 length(clsmems)+1],[100 100],'r-'); hold on;
        for s = 1:length(clsmems)            
            ph = bar([s-.2 s+.2],clsmems(:,s)); hold on;
            set(ph,'facecolor',cols(s,:));
            set(ph,'edgecolor',cols(s,:));
        end;
        set(gca,'xlim',[0 length(clsmems)+1]);%set(gca,'ylim',[0 maxy]);
        set(gca,'xtick',[1:length(clsmems)]);
        set(gca,'xticklabel',{'Delta','Theta','Alpha','Beta','Hi/Low','Gamma'});
    end;        
