% for talk figures, plots all window weights for specified subject
%
%
% ims -- [vector] of independent modulators to plot time course for
% smoothon -- [0|1] if 1, then will smooth the weights over trials to see trends
% row/col/pl -- for plotting into an existing figure with rows/cols
% pl tells it where to start plotting
% labels -- [0|1] if 1, will plot emotion labels, none is 0

function [newpl] = PlotSubjBackProjs(fullpath,savedat,emos,ims,smoothon,emoorder,row,col,pl,labels)
    
     
    
    s = load([fullpath,savedat,'.mat']);  
    sph=floatread([fullpath,savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpath,savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    ws = wts*sph;     icawinv = pinv(ws);   %icawinv = icawinv';
    clear wts sph ws
    cols = jet(15);cols(10,:) = [.9 .9 0];
    %figure; 
    
    %row = round(sqrt(size(icawinv,2)));col = ceil(sqrt(size(icawinv,2)));
    for im = 1:length(ims)
        dim = ims(im);
        sbplot(row,col,pl); pl = pl+1;
        if smoothon == 1
			 xpl= 0;
            for e = 1:length(s.dstrials)
                em = find(strcmp(emoorder{e},emos));
                xwidth = round(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/40);
                [outdata,outx] = movav(icawinv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),dim),0,xwidth,0);
                %newdata(end+1:end+size(outdata,2),dim) = outdata;
                %intvls(emm,:) = [size(outdata,2) length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em))) - size(outdata,2)];
                datlen = length(outdata);
                ph = plot([xpl+1:xpl+datlen],outdata,'k-','linewidth',1); hold on;
                set(ph,'color',cols(em,:));
                set(gca,'ylim',[-2 4]);yl = get(gca,'ylim');
                ph = plot([xpl+1],[yl],'k-');
                if labels == 1
                    ph = text(xpl+1,2.5,emoorder{e});
                    set(ph,'rotation',45);set(ph,'color',cols(em,:)); 
                end; 
               ph = plot([xpl+1 xpl+1],[yl],'k-');
                xpl = xpl+datlen;
            end;
            set(gca,'xlim',[0 xpl+1]);  
            %for e = 1:size(intvls,1)
            %    emm =  find(strcmp(emoorder{e},emos));
            %    datlen = intvls(em,1);
            %    datdiff = intvls(em,2);
                %dattot = intvls(em,3);               
                
            %    ph = plot([sum(intvls(1:(e-1),1))+1:sum(intvls(1:(e-1),1))+datlen],newdata(sum(intvls(1:(e-1),1))+1:sum(intvls(1:(e-1),1))+datlen,dim),'k-','linewidth',1); hold on;
            %    set(ph,'color',cols(emm,:));
            %    set(gca,'ylim',[-2 4]);            yl = get(gca,'ylim');
            %    ph = plot([sum(intvls(1:(e-1),1)) sum(intvls(1:(e-1),1))],[yl],'k-');
            %end;
            %set(gca,'xlim',[0 sum(intvls(:,1))]);  
            %if labels == 1
            %    for e = 1:length(emos)
            %        emm =  find(strcmp(emoorder{e},emos));
            %        ph = text(sum(intvls(1:(e-1),1))+1,2.5,emoorder{e});
            %        set(ph,'rotation',45);set(ph,'color',cols(emm,:)); 
            %    end;                
            %end; 
        else   
         xpl= 0;
            for e = 1:length(emos)                
                em =  find(strcmp(emoorder{e},emos));
                datlen = length([sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em))]);
				ph = plot([xpl+1:xpl+datlen],icawinv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),dim));hold on;
                set(gca,'ylim',[-2 4]);        yl = get(gca,'ylim');
                if labels == 1
                    ph = text(xpl+1,2.5,emoorder{e});
                    set(ph,'rotation',45); set(ph,'color',cols(em,:)); 
                end;
                ph = plot([xpl+1 xpl+1],[yl],'k-');
                xpl = xpl+datlen;
            end; 
            set(gca,'xlim',[0 sum(s.dstrials)+1]);
        end;
        ph = plot([get(gca,'xlim')],[0 0],'k--','linewidth',2);
        set(ph,'color','y');
        title([fullpath(end-4:end-1),' IM ',int2str(dim)]);
    end;
    %axcopy
    ph=textsc(['Subject ',fullpath(end-4:end-1)],'title');set(ph,'fontsize',16);
    newpl = pl;
