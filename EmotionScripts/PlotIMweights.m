% for talk figures, plots all window weights for specified subject
%
%
% ims -- [vector] of independent modulators to plot time course for
% smoothon -- ['on' or 'off'] if 'on', then will smooth the weights over trials to see trends
% row/col/pl -- for plotting into an existing figure with rows/cols
% pl tells it where to start plotting
% labels -- [cell array] label for each dataset entered into analysis (must
% be the same length as s.dstrials
% condcols -- [cell array] should be a collection of strings or ? x 3 matrix 
% specifying for each corresponding 'label' what color the wts and labels
% should be. [] will create 'jet' color coding across conditions.

function [newpl] = PlotIMweights(fullpath,savedat,labels,condcols,ims,smoothon,row,col,pl)

s = load([fullpath,savedat,'.mat']);
sph=floatread([fullpath,[savedat],'.sph'],[s.pcs s.pcs],[],0);
wts=floatread([fullpath,[savedat],'.wts'],[s.pcs s.pcs],[],0);
icamatall = floatread([fullpath,[savedat],'.fdt'],[s.pcs s.numframes],[],0);
ws = wts*sph;    acts = ws*icamatall;    winv = pinv(ws);
clear wts sph ws icamatall
speceig = floatread([fullpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);
specwts = speceig*winv;
icawinv = specwts;

xwidth = round(size(icawinv,1)/150);
            
if isempty(condcols)
    condcols = jet(length(s.dstrials));%cols(10,:) = [.9 .9 0];
end;
if isempty(ims)
    ims = [1:s.pcs];
end;

for im = 1:length(ims)
    dim = ims(im);
    if pl > row*col
        ph=textsc(['Subject ',fullpath(end-4:end-1),': IM weights across conditions'],'title');set(ph,'fontsize',16);
        figure; pl = 1;
    end;
    sbplot(row,col,pl); pl = pl+1;
    if strcmp(smoothon,'on')
        xup = 0;
        for em = 1:length(unique(s.keeptrack(:,1)))
            [outdata,outx] = movav(icawinv([s.keeptrack(em,1):s.keeptrack(em,2)],dim),0,xwidth,0);
            if iscell(condcols)
                ph = plot([xup+1:xup+length(outx)],outdata,'-','linewidth',1,'color',condcols{em}); hold on;
            else
                ph = plot([xup+1:xup+length(outx)],outdata,'-','linewidth',1,'color',condcols(em,:)); hold on;
            end;
            set(gca,'ylim',[-2 4]);yl = get(gca,'ylim');
            if ~isempty(labels)
                ph = text(xup+1,2.5,labels{em});set(ph,'rotation',45);
                if iscell(condcols)
                    set(ph,'color',condcols{em});
                else
                    set(ph,'color',condcols(em,:));
                end;
            end;
            ph = plot([xup+1 xup+1],[yl],'k-');
            xup = xup + length(outx);
        end;
        set(gca,'xlim',[1 xup]);
    else
        for em = 1:length(unique(s.keeptrack(:,1)))
            if iscell(condcols)
                ph = plot([s.keeptrack(em,1):s.keeptrack(em,2)],icawinv(s.keeptrack(em,1):s.keeptrack(em,2),dim),'linewidth',1,'color',condcols{em});
            else
                ph = plot([s.keeptrack(em,1):s.keeptrack(em,2)],icawinv(s.keeptrack(em,1):s.keeptrack(em,2),dim),'linewidth',1,'color',condcols(em,:));
            end;
            hold on;
            set(gca,'ylim',[-2 4]);        yl = get(gca,'ylim');
            if ~isempty(labels)
                ph = text(s.keeptrack(em,1),2.5,labels{em});set(ph,'rotation',45);
                if iscell(condcols)
                    set(ph,'color',condcols{em});
                else
                    set(ph,'color',condcols(em,:));
                end;
            end;
            ph = plot([s.keeptrack(em,1) s.keeptrack(em,1)],[yl],'k-');
        end;
        set(gca,'xlim',[1 s.keeptrack(em,2)]);
    end;
    ph = plot([get(gca,'xlim')],[0 0],'k--','linewidth',2);
    set(ph,'color','y'); set(gca,'xticklabel',[]);
    title([fullpath(end-4:end-1),' IM ',int2str(dim)]);
    if pl == 2
        xlabel('Time points'); ylabel('IM weights');
    end;
    title(['IM ',int2str(dim)]);
end;
%axcopy
ph=textsc(['Subject ',fullpath(end-4:end-1),': IM weights across conditions'],'title');set(ph,'fontsize',16);
newpl = pl;
