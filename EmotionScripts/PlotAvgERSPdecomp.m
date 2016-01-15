% plots results from AvgERSPdecomp()
%
% [acts,row,col,pl] = PlotAvgERSPdecomp(wtsphname,writepath,plottype,groupidx);
%
%
% INPUTS:
% wtsphname -- [string] same filename title given to AvgERSPdecomp()
% writepath -- [string] same full directory path as given to AvgERSPdecomp()
% plottype -- ['winv' or 'wtdmean'] to plot factor templates of weighted mean, respectively
% plotfacs -- [vector] if not [], will plot only specified factors. 
% groupidx -- [vector of integers] each integer value indicates a specific a priori group
%             assignment. Start at 1 and increment by 1 with every additional group.
%
% OUTPUTS:
% acts -- [matrix] of weights of size (ndims x nICs)
% grpdims -- [vector] one value for each dimension. 1 indicates a high loading on the 
%                     group indicator input. -1 means that this dimension is not 
%                     associated with the group assignment, but has a significant ERSP
%                     pattern beyond chance. 0 means that neither of these cases was
%                     true. [] is returned if 'plottype' is 'winv'
% row, col and pl are for plotting purposes.


function [acts,grpdims,grpP,grpvals,plotidx,row,col,pl] = PlotAvgERSPdecomp(wtsphname,writepath,plottype,plotfacs,groupidx ,row,col,pl);

    
    ninputdims = 5; % number of independent factors from each subj
    plotscatter = 'off';
    grpdims = [];
    groupcut = .5;
    grpvals = [];
    plotidx = [];
    
    alpha = .01;
    if min(groupidx) < 0
        groupidx = groupidx + 2;
    end;
    maxsubj1 = length(find(groupidx==1))/ninputdims;
    maxsubj2 = length(find(groupidx==2))/ninputdims;
    
    shuffnum = 500; alpha = .01;
    cols = {'r','b','g','m','c'};
    s = load ([writepath,wtsphname,'.mat']);
    sph=floatread([writepath,wtsphname,'.sph'],[s.channo s.channo],[],0); 
    wts=floatread([writepath,wtsphname,'.wts'],[s.channo s.channo],[],0); 
    icamat=floatread([writepath,wtsphname,'.fdt'],[s.channo s.numframes],[],0);    
    ws = wts*sph;winv = pinv(ws);    acts = ws*icamat;
    
    %dat=floatread([writepath,wtsphname,'DAT.fdt'],[length(s.freqs)*length(s.freqs)*length(s.datafiles)+length(s.freqs) inf],[],0);   
 
    if isempty(plotfacs)
        plotfacs = [1:size(winv,2)];
    end;

    if ~isempty(plotfacs)
        winv = winv(:,plotfacs);
        acts = acts(plotfacs,:);
    end;
    
    erspeig = floatread([writepath,wtsphname,'ERSPEIG.fdt'],[length(s.freqs)*length(s.times)*length(s.datafiles) s.pcs],[],0);
    erspdat = erspeig*winv(1:s.pcs,:);
    erspback = erspeig*icamat(1:s.pcs,:);
    
    if ~isempty(s.pcb)
        speceig = floatread([writepath,wtsphname,'SPECEIG.fdt'],[length(s.freqs) s.pcb],[],0);
        specdat = speceig*winv(s.pcs+1:s.pcs+s.pcb,:); 
        specback = speceig*icamat(s.pcs+1:s.pcs+s.pcb,:); 
    end;
    
    if strcmp(plottype,'wtdmean')% calculated weighted means and bootstrapping distributions:
        fprintf('\nAccumulating bootstrap statistics...\n')
        clear limersp limmat 
        %%%%%%%  calculate weighted ERSPs for each dim:
        wtdersp = (erspback*acts')/size(erspback,2);
        % multiply each trial element,sum then divide by ntrials 
        bootwts = zeros(size(wtdersp,1),size(wtdersp,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (erspback*shuffle(acts,2)')/size(erspback,2);
        end;
        bootwts = sort(bootwts,3);
        limersp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limersp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        clear bootwts
        if ~isempty(s.pcb)
            %%%%%%%  calculate weighted mean spectra:
            wtdspec = (specback*acts')/size(erspback,2); 
            bootwts = zeros(size(wtdspec,1),size(wtdspec,2),shuffnum);
            for b= 1:shuffnum
                bootwts(:,:,b) = (specback*shuffle(acts,2)')/size(erspback,2);
            end;
            bootwts = sort(bootwts,3);
            limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
            limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        end;
    end;
    
    
    if strcmp(plottype,'winv')
        fprintf('\nPlotting factor templates...\n')
        if strcmp(plotscatter,'on')
            if isempty(row)
               figure; rw = size(winv,2)-1; col = size(winv,2)-1; 
            end;
            for fac1 = 1:size(winv,2)-1
                pl = (fac1-1)*rw+fac1;
                mnstd1 = mean(acts(fac1,:)) + std(acts(fac1,:));
                for fac2 = fac1+1:size(winv,2)
                    mnstd2 = mean(acts(fac2,:)) + std(acts(fac2,:));
                   sbplot(rw,col,pl); pl = pl+1;            
                    for gp = 1:max(groupidx)                
                        ph = plot(acts(fac1,find(groupidx == gp)),acts(fac2,find(groupidx == gp)),'k.'); hold on;
                        set(ph,'color',cols{gp});set(ph,'markersize',5);
                        set(gca,'box','off');
                    end; 
                    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
                    plot([get(gca,'xlim')],[0 0],'k-'); plot([0 0],[get(gca,'ylim')],'k-');  
                    title([int2str(fac1),'-',int2str(fac2)]);
                end;
            end; 
        end;
        if isempty(row)
            figure; pl = 1;
            if ~isempty(s.pcb) 
                row = ceil(size(winv,2)/3)+1;  col = (length(s.datafiles)+1)*3;
                if col > 6
                    col = 6; row = row+2;
                end;
            else
                row = ceil(size(winv,2)/4)+1;  col = (length(s.datafiles))*4;
                if col > 6
                    col = 6; row = row+1;
                end;            
            end;
        else
            close; % close the scatter plot in this case
        end;
        
        for dim = 1:size(winv,2)
            for sb = 1:size(acts,2)/ninputdims % for each IC in cluster
                [x y] = max(abs(acts(dim,(sb-1)*ninputdims+1:sb*ninputdims)));
                if groupidx((sb-1)*ninputdims+1) == 1
                    grpvals(dim,sb) = acts(dim,(sb-1)*ninputdims+y);
                    plotidx(dim,sb) = 1;
                elseif  groupidx((sb-1)*ninputdims+1) == 2
                    grpvals(dim,sb) = acts(dim,(sb-1)*ninputdims+y);
                    plotidx(dim,sb) = 2;
                end;
            end;
            
            [grpdims(1,dim),grpP(1,dim)] = ttest2(grpvals(find(plotidx == 1)),grpvals(find(plotidx == 2)),.01);
                
            hival = find(acts(dim,:) > mean(acts(dim,:))+1.5*std(acts(dim,:)));
            loval = find(acts(dim,:) < mean(acts(dim,:))-1.5*std(acts(dim,:)));
            %grpskew(1,1,dim) = length(find(hival == 1))/maxsubj1;
            %grpskew(1,2,dim) = length(find(hival == 2))/maxsubj2;
            %grpskew(2,1,dim) = length(find(loval == 1))/maxsubj1;
            %grpskew(2,2,dim) = length(find(loval == 2))/maxsubj2;
             
            
            for df = 1:length(s.datafiles)
                sbplot(row,col,pl); pl = pl+1;
                tmpcomp = erspdat((df-1)*length(s.freqs)*length(s.times)+1:df*length(s.freqs)*length(s.times),dim);
                tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
                lim = max(abs(tmpcomp(:)));
                if strcmp(s.freqscale,'quad')
                    quadimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                elseif  strcmp(s.freqscale,'log')                
                    mylogimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                else
                    imagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                    set(gca,'ytick',[10:10:s.freqs(end)]);
                end;
                set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                set(gca,'xticklabel',[]);
                title(['Winv Fac ',int2str(dim)]);
            end;
            if ~isempty(s.pcb)
                sbplot(row,col,pl); pl = pl+1;
                ph = plot(s.freqs,specdat(:,dim),'b-','linewidth',2);hold on;
                ph = plot([10 10],[get(gca,'ylim')],'r-');
                set(gca,'xlim',[s.freqs(1) s.freqs(end)]);
                set(gca,'xgrid','on');  
            end;
        end;
        if ~isempty(s.groupidx)
            sbplot(row,col,[pl pl+1])
            ph = plot(winv(end,:),'m.-','markersize',20); 
        end;
    else
        fprintf('\n	Plotting factor weighted means...\n')
        if isempty(row)
            figure; pl =	1;
            if ~isempty(s.pcb) 	
                row = ceil(size(winv,2)/3)+1;  col = (length(s.datafiles)+1)*3;
                if col > 6	
                    col = 6;	row = row+1;
                end;	
            else	
                row = ceil(size(winv,2)/4)+1;  col = (length(s.datafiles))*4;
                if col > 6
                    col = 6;	row = row+1;
                end;            	
            end;	
        end;	
        
        for dim = 1:size(winv,2)
            hival = find(acts(dim,:) > mean(acts(dim,:))+1.5*std(acts(dim,:)));
            loval = find(acts(dim,:) < mean(acts(dim,:))-1.5*std(acts(dim,:)));
            hivalid = groupidx(hival);
            lovalid = groupidx(loval);
            %grpskew(1,1,dim) = length(find(hival == 1))/maxsubj1;
            %grpskew(1,2,dim) = length(find(hival == 2))/maxsubj2;
            %grpskew(2,1,dim) = length(find(loval == 1))/maxsubj1;
            %grpskew(2,2,dim) = length(find(loval == 2))/maxsubj2;
            for sb = 1:size(acts,2)/ninputdims % for each IC in cluster
                [x y] = max(abs(acts(dim,(sb-1)*ninputdims+1:sb*ninputdims)));
                if groupidx((sb-1)*ninputdims+1) == 1
                    grpvals(dim,sb) = acts(dim,(sb-1)*ninputdims+y);
                    plotidx(dim,sb) = 1;
                elseif  groupidx((sb-1)*ninputdims+1) == 2
                    grpvals(dim,sb) = acts(dim,(sb-1)*ninputdims+y);
                    plotidx(dim,sb) = 2;
                end;
            end;
            
           [grpdims(1,dim),grpP(1,dim)] = ttest2(grpvals(find(plotidx == 1)),grpvals(find(plotidx == 2)),.01);

            for df = 1:length(s.datafiles)
                sbplot(row,col,pl); pl = pl+1;
                tmpcomp = wtdersp((df-1)*length(s.freqs)*length(s.times)+1:df*length(s.freqs)*length(s.times),dim);
                % mask the wtd mean by the bootstrap values:
                tmpcomp(find(tmpcomp>limersp((df-1)*length(s.freqs)*length(s.times)+1:df*length(s.freqs)*length(s.times),dim,1)&tmpcomp<limersp((df-1)*length(s.freqs)*length(s.times)+1:df*length(s.freqs)*length(s.times),dim,2))) = 0;
                tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
                
                lim = max(abs(tmpcomp(:)));
                if lim == 0
                    lim = 1;
                end;
                if strcmp(s.freqscale,'quad')
                    quadimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                elseif  strcmp(s.freqscale,'log')                
                    mylogimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                else
                    imagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
                    set(gca,'ytick',[10:10:s.freqs(end)]);
                end;
                set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                set(gca,'xticklabel',[]);
                title(['Factor ',int2str(dim)]);
                totvox1 = size(tmpcomp,1)*size(tmpcomp,2);
            end;
            if ~isempty(s.pcb)
                sbplot(row,col,pl); pl = pl+1;
                for ee = 1:size(wtdspec,1)
                    ph = plot([ee ee],[limmat(ee,dim,1) limmat(ee,dim,2)],'k-','linewidth',3);
                    set(ph,'color',[.8 .8 .8]);hold on;
                end;
                tmpcomp2 = wtdspec(:,dim);
                % mask the wtd mean by the bootstrap values:
                tmpcomp2(find(tmpcomp2>limmat(:,dim,1)&tmpcomp2<limmat(:,dim,2))) = 0;           
                ph = plot(s.freqs,tmpcomp2,'b-','linewidth',2);hold on;
                set(gca,'xlim',[s.freqs(1) s.freqs(end)]);
                set(gca,'xgrid','on');  
                ph = plot([10 10],[get(gca,'ylim')],'r-');  
                totvox2 = size(tmpcomp2,1)*size(tmpcomp2,2);
            end;
% $$$             if winv(end,dim) > groupcut & length(find(tmpcomp(:)))>totvox1/40 
% $$$                 grpdims(1,dim) = 1;
% $$$             elseif winv(end,dim) <= groupcut & length(find(tmpcomp(:)))>totvox1/20
% $$$                 grpdims(1,dim) = -1;
% $$$             else
% $$$                 grpdims(1,dim) = 0;
% $$$             end;                
        end;
        if ~isempty(s.groupidx)
        sbplot(row,col,[pl pl+1])
            ph = plot([0 size(winv,2)+1],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);hold on;
            ph = plot(winv(end,:),'m.-','markersize',20); 
            set(gca,'xlim',[0 size(winv,2)+1]);
            xlabel('Dimensions')
            title('Group indicator value over dimensions');
        end;
    end;        
    
