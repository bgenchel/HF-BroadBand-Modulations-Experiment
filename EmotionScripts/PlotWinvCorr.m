% takes two "EEG.icawinv" matrices, runs matcorr to find best fitting ICs and plots IC pairs
% this function compares the output of topoplot, not the winv itself, so # of chans and 
% # of components can differ between datasets.
% 
% [corr] = PlotWinvCorr(wv1,wv2,chanlocs1,chanlocs2,plotto);
%
% INPUTS:
% wv1 -- EEG.icawinv from first dataset
% wv2 -- EEG.icawinv from second dataset
% chanlocs1 -- EEG.chanlocs from first dataset
% chanlocs2 -- EEG.chanlocs from second dataset
% plotto -- [integer] how many components from each set to plot (max of 32 pairs are plotted per page); [] plots all components
% If 'plotto' is specified, it must be less than or equal to the # of components for the smallest wv.
%
% OUTPUT: 
% corr -- [1 X #chan] vector containing correlation coefficient for each pair


function [corr] = PlotWinvCorr(wv1,wv2,chanlocs1,chanlocs2,plotto);
    
  
    if isempty(plotto)
        if size(wv1,2) <= size(wv2,2)
            calcto = size(wv1,2);
        elseif size(wv1,2) >= size(wv2,2)
            calcto = size(wv2,2);
        end;        
    else
      calcto = plotto;
    end;    
    for ic = 1:calcto
        [h grid1(:,:,ic) plotrad, xmesh, ymesh]= topoplot(wv1(:,ic),chanlocs1,'noplot','on');
        [h grid2(:,:,ic) plotrad, xmesh, ymesh]= topoplot(wv2(:,ic),chanlocs2,'noplot','on');
    end;
       
    for rw = 1:size(grid1,1)
        nonan{rw} = find(grid1(rw,:,1) > 0 | grid1(rw,:,1) < 0);
    end;
    for ic = 1:size(grid1,3)
        ng1 = zeros(1,0);
        ng2 = zeros(1,0);    
        for rw = 2:size(grid1,1)-1
            ng1(1,end+1:end+length(nonan{rw})-4) = grid1(rw,[nonan{rw}(1) + 2:nonan{rw}(end) - 2],ic);
            ng2(1,end+1:end+length(nonan{rw})-4) = grid2(rw,[nonan{rw}(1) + 2:nonan{rw}(end) - 2],ic);
        end;
        newgrid1(ic,:) = ng1;
        newgrid2(ic,:) = ng2;
    end;   
    
    [corr,indx,indy,corrs] = matcorr(newgrid1,newgrid2);
    
    if isempty(plotto)
        col = round(sqrt(length(indx)*2));
        row = ceil(sqrt(length(indx)*2));
    else
        col = round(sqrt(plotto*2));
        row = ceil(sqrt(plotto*2));
    end;        
    if ~iseven(col)
        col = col+1;
        row = row-1;
    end;    
    if col > 8
        col = 8;
        row = 8;
    end;
    
    figure; pl = 1;
    for pr = 1:length(indx)
        if pl > (row*col)-1
          sbplot(1,1,1)
          if col == 8
              plot([.5 .5],[0 1],'k-','linewidth',2);hold on;
              set(gca,'xlim',[-1 1]);          
              plot([0 0],[0 1],'k-','linewidth',2);
              plot([-.5 -.5],[0 1],'k-','linewidth',2);
          elseif col == 6
              plot([.35 .35],[0 1],'k-','linewidth',2);hold on;
              set(gca,'xlim',[-1 1]);          
              plot([-.35 -.35],[0 1],'k-','linewidth',2);
          end;          
          axis('off')          
          set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            figure; pl = 1;
        end;
        sbplot(row,col,pl)
        topoplot(wv1(:,indx(pr)),chanlocs1,'electrodes','off','plotrad',.65);
        if pl <= col
          title(['Set 1; IC ',int2str(indx(pr))]);
        else
          title(['IC ',int2str(indx(pr))]);
        end;        pl = pl+1;
        sbplot(row,col,pl)
        topoplot(wv2(:,indy(pr)),chanlocs2,'electrodes','off','plotrad',.65);
        if pl <= col
          title(['Set 2; IC ',int2str(indy(pr))]);
        else
          title(['IC ',int2str(indy(pr))]);
        end;   pl = pl+1;     
    end;
    sbplot(1,1,1)
    if col == 8
        plot([.5 .5],[0 1],'k-','linewidth',2);hold on;
        set(gca,'xlim',[-1 1]);          
        plot([0 0],[0 1],'k-','linewidth',2);
        plot([-.5 -.5],[0 1],'k-','linewidth',2);
    elseif col == 6
        plot([.35 .35],[0 1],'k-','linewidth',2);hold on;
        set(gca,'xlim',[-1 1]);          
        plot([-.35 -.35],[0 1],'k-','linewidth',2);
    end;          
    axis('off')          
    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    
    corr = corr';
