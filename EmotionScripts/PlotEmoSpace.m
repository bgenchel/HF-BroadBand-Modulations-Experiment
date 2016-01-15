% plots emo spaces for single subject
%
%
%
%




function PlotEmoSpace(emomeans,maxfacs,nx,plot2d,cutline);

    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
    cols = jet(15);cols(10,:) = [.9 .9 0];
    dims = zeros(15,0);
    for dim = 1:size(emomeans{nx},1)
        dims(:,end+1) = emomeans{nx}(dim,:);
        dims(:,end) = dims(:,end)/std(dims(:,end));
    end;
    
    if isempty(plot2d)
        figure;  % for one point per emotion (mean)
        for wv = 1:size(dims,2)
            sbplot(round(sqrt(size(dims,2))),ceil(sqrt(size(dims,2))),wv)
            ph = plot([0 size(emomeans{nx},2)+3],[0 0],'k:');
            set(ph,'color',[.7 .7 .7]); hold on;
            ph = plot([0 size(emomeans{nx},2)+3],[-cutline -cutline],'k--');
            set(ph,'color',[.5 .5 .5]); hold on;
            ph = plot([0 size(emomeans{nx},2)+3],[cutline cutline],'k--');
            set(ph,'color',[.5 .5 .5]); hold on;
            for e = 1:size(dims,1)
                ph=plot(e,dims(e,wv),'.');hold on;
                set(ph,'markersize',20);set(ph,'color',cols(e,:));
                ph = text(e,dims(e,wv),emo2{e});
                set(ph,'color',cols(e,:)); 
            end;            
            set(gca,'ylim',[-4 4]);
            set(gca,'xlim',[0 size(emomeans{nx},2)+3]);
            title(['Sbj ',int2str(nx),' Dim ',int2str(wv)]);
        end;axcopy
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        textsc(['Subject ',int2str(nx)],'title');
    elseif length(plot2d) < 3
        %%  Plot 2 Dims vs each other:
        c1 = plot2d(1); c2 = plot2d(2);
        figure;  % for one point per emotion (mean)
        for e = 1:size(dims,1)
            ph=plot(dims(e,c1),dims(e,c2),'.');hold on;
            set(ph,'markersize',20);set(ph,'color',cols(e,:));
            ph = text(dims(e,c1),dims(e,c2),emo2{e});
            set(ph,'color',cols(e,:)); 
        end;
        xlabel(['dim ',int2str(c1)]);ylabel(['dim ',int2str(c2)]);
        
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        textsc(['Subject ',int2str(nx)],'title');
    elseif length(plot2d) ==3
        %%  Plot 3 Dims vs each other:
        figure; % just 3  dims vs each other
        c1 = plot2d(1); c2 = plot2d(2);c3 = plot2d(3);
        for e = 1:size(dims,1)
            ph=plot3(dims(e,c1),dims(e,c2),dims(e,c3),'.');hold on;
            set(ph,'markersize',25);                set(ph,'color',cols(e,:));
            ph = text(dims(e,c1),dims(e,c2),dims(e,c3),emo2{e});
            set(ph,'color',cols(e,:)); set(ph,'fontsize',14); 
        end;
        zl = get(gca,'zlim');
        for e = 1:size(dims,1)
            pl =plot3([dims(e,c1) dims(e,c1)],[dims(e,c2) dims(e,c2)],[zl(1)  dims(e,c3)]);
            set(pl,'color',cols(e,:)); set(pl,'linewidth',2)             
        end;
        set(gca,'xgrid','on');  set(gca,'ygrid','on');set(gca,'zgrid','on');
        xlabel(['dim ',int2str(c1)]);ylabel(['dim ',int2str(c2)]);zlabel(['dim ',int2str(c3)]);
        textsc(['Subject ',int2str(nx)],'title');
        
    end; 
    
