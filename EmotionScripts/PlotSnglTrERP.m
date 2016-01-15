% Plot results from SnglTrialICAClustering.m using the 'wt' option and NO context
%
% PlotWTsnglTrialERSP(pathname,savedat,savename,complist,toplot,whichfacs);
%
% savedat: [string] title stem used in SnglTrialERSPDecomp.m
% pathname: [string] full path where wts, sph can be found and savename will be saved
% wts: [string] Name of weights float file that can be found in pathname 
% sph: [string] Name of sphere float file that can be found in pathname 
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% pcs: [integer] Number of PCA dimensions retained before ICA on full, combined matrix
% toplot: [integer] Max number of factors to plot on one page.
% plotwts: if not [], will plot weights over the course of the session
%          0 will plot straight weights, any other number will plot weights smoothed by this factor
% whichfacs: vector of integers corresponding to factors to plot. [] for all

function [pl,cxtout] = PlotSnglTrERP(datset,pathname,stem,ic,plottype,labels,whichtmps,row,col,place,savettl,pltopt);

    
    
    pvafcalc = 0; % plot pvaf of context instead
    cxtout = [];
    erpim = 0; % 1 plots erpimage instead of erp
    distr = 0; % only if plottype == 'winv'
    tickl = .05;
    shuffnum = 500;
    alpha = .01;
    cols = hsv(length(labels));
     cols(2,:)=[.75 .75 0]; % make yellow darker
     fs = 10; % fontsize
             %%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_loadset(datset, pathname);
    savedat = [stem,int2str(ic)];
    s = load([pathname,savedat,'.mat']);     
    wts = floatread([pathname,savedat,'.wts'],[s.pcerp+s.pcctx s.numrows],[],0);
    sph = floatread([pathname,savedat,'.sph'],[s.numrows s.numrows],[],0);     
    dat = floatread([pathname,savedat,'.fdt'],[s.numrows s.numframes],[],0);     
    %realdat = floatread([pathname,savedat,'DAT.fdt'],[length(s.times) s.numframes],[],0);     
    ws = wts*sph;winv = pinv(ws); acts = ws*dat;
    
    if isempty(whichtmps)
        whichtmps = [1:size(winv,2)];
    end;

    if isfield(s,'pcerp')
        if s.pcerp > 0
            erpeig = floatread([pathname,savedat,'EIGVEC.fdt'],[length(s.times) s.pcerp],[],0);
            erpdat = erpeig*dat(1:s.pcerp,:);
            erptempl = erpeig*winv(1:s.pcerp,:);
            %erpproj = winv*acts;erpproj = erpeig*erpproj(1:s.pcerp,:);
        end;
    end;
    if isfield(s,'pcctx')
        if s.pcctx > 0
            cxteig = floatread([pathname,savedat,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
            cxtdat = cxteig*dat(s.pcerp+1:end,:);  
            cxttmpl = cxteig*winv(s.pcerp+1:end,:);  
        else
            cxtdat = [];
        end;
    end;    
    if strcmp(plottype,'wtdmean')
        clear limerp limmat 
        %%%%%%%  calculate weighted ERSPs for each dim:
        wtderp = (erpdat*acts')/size(erpdat,2);
        % multiply each trial element,sum then divide by ntrials 
        bootwts = zeros(size(wtderp,1),size(wtderp,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (erpdat*shuffle(acts,2)')/size(erpdat,2);
        end;
        bootwts = sort(bootwts,3);
        limerp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limerp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        clear bootwts
        if ~isempty(cxtdat) % only if context was included
            %%%%%%%  calculate context vectors for each dim:
            wtdctx = (cxtdat*acts')/size(erpdat,2); 
            bootwts = zeros(size(wtdctx,1),size(wtdctx,2),shuffnum);
            for b= 1:shuffnum
                bootwts(:,:,b) = (cxtdat*shuffle(acts,2)')/size(erpdat,2);            
            end;
            bootwts = sort(bootwts,3);
            limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
            limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
%mnctx = mean(bootwts,3);
%stdctx = std(bootwts,0,3);
        end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%% Plot Reg ERP factors onto one page %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(row)
        row = ceil(sqrt(length(whichtmps)));
    end;
    if isempty(col)
        col = ceil(sqrt(length(whichtmps)));
    end;
    if isempty(place)
        figure;   pl = 1; pg = 1;
    else 
        if place == 1
            if strcmp(pltopt,'on')
                figure;
            else
                clf; pl = 1;
            end;
        end;
        pl = place; 
    end;
    
    sbplot(row,col,pl);pl = pl+1;
    topoplot(EEG.icawinv(:,ic),EEG.chanlocs(EEG.icachansind),'electrodes','off');            
    title(['IC ',int2str(ic)]);
% $$$     % plot erp:-----------------------------------------------------
% $$$     sbplot(row,col,pl);pl = pl+1;
% $$$     ph = plot(s.times,s.icmean,'r-','linewidth',3);
% $$$     set(gca,'ticklength',[tickl tickl]);
% $$$     hold on; plot([0 0],[get(gca,'ylim')],'k-');
% $$$     hold on; plot([get(gca,'xlim')],[0 0],'k-');
% $$$     set(gca,'fontsize',fs);  
% $$$     if pl < col+1
% $$$         title(['Mean ERP']); 
% $$$     end;
% $$$     set(gca,'xlim',[s.times(1) s.times(end)]);
% $$$     if pl < (row-1)*col
% $$$         set(gca,'xticklabel',[]);
% $$$     end;
    for tp = 1:length(whichtmps)
        template = whichtmps(tp);       
        if pl > row * col
            if ~isempty(savettl)
                makesave = [pathname,savettl];    eval(makesave);   
            end;
            if strcmp(pltopt,'on')
                figure; pl=1;
            else
                clf; pl = 1;
            end;
        end;            
        if ~isempty(cxtdat)% only if context was included
            sbplot(row,col,pl);pl = pl+1;    
            if strcmp(plottype,'winv')
                allctxvals = cxttmpl(:,template);
                for ee = 1:length(allctxvals)
                    ph=plot(ee,allctxvals(ee),'.'); hold on;
                    set(ph,'markersize',18);                         
                    set(ph,'color',cols(ee,:));
                    ph=text(ee,allctxvals(ee),labels{ee});;
                    set(ph,'color',cols(ee,:));set(ph,'fontsize',6);
                    set(ph,'rotation',90);
                end;
                if isempty(cxtout)
                    cxtout = cxttmpl(:,template)';
                else
                    cxtout(end+1,:) = cxttmpl(:,template)';
                end;
                set(gca,'xlim',[0 length(allctxvals)+1]);
            elseif strcmp(plottype,'wtdmean')
                
                if pvafcalc == 1% calculate pvaf instead:---                
                                % first make back proj of single template/factor and only cxt pcs:
                    bp = winv(s.pcerp+1:end,template) * acts(template,:);
                    bpq = cxteig * bp;
                    % calculate pvaf of PCA-reduced data
                    for q = 1:size(bpq,1)
                        pv(1,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
                    end;    
                    pv = pv*100;
                    for b = 1:shuffnum
                        bp = winv(s.pcersp+1:end,template) * shuffle(acts(template,:),2);
                        bpq = cxteig * bp;
                        %bpq = cxteig * shuffle(bp,2);
                        for q = 1:size(bpq,1)
                            pvboot(b,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
                        end;   
                        pvboot(b,:) = pvboot(b,:).* (cxttmpl(:,template)'./abs(cxttmpl(:,template)'));
                    end;                
                    pvboot = pvboot*100;
                    pvboot = sort(abs(pvboot),1);
                    limmat(2,:) = pvboot(end-ceil(shuffnum*alpha),:); % max boot
                    limmat(1,:) = pvboot(ceil(shuffnum*alpha),:);  % min boot
                    tmpcomp = pv .* (cxttmpl(:,template)'./abs(cxttmpl(:,template)'));
                    for ee = 1:length(labels)
                        ph = plot([ee ee],[limmat(1,ee) limmat(2,ee)],'k-','linewidth',3);
                        set(ph,'color',[.8 .8 .8]);hold on;
                    end;
                else
                    tmpcomp = wtdctx(:,template);
                    for ee = 1:length(labels)
                        ph = plot([ee ee],[limmat(ee,template,1) limmat(ee,template,2)],'k-','linewidth',3);
                        set(ph,'color',[.8 .8 .8]);hold on;
                    end;
                end;
                % mask the wtd mean by the bootstrap values:
                tmpcomp(find(tmpcomp>limmat(:,template,1)&tmpcomp<limmat(:,template,2))) = 0;
                if isempty(cxtout)
                    cxtout = cxttmpl(:,template)';
                else
                    cxtout(end+1,:) = cxttmpl(:,template)';
                end;
                cxttmpl(:,template) = cxttmpl(:,template)*(sqrt(mean(tmpcomp.^2)))/(sqrt(mean(cxttmpl(:,template).^2))); % normalize template to rms of wtdmean
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.75 .75 .75]);   hold on;         
                for ee = 1:size(wtdctx,1)
                    if tmpcomp(ee) == 0
                        ph=plot(ee,tmpcomp(ee),'.');hold on; 
                        set(ph,'markersize',8);
                    else
                        ph=plot(ee,tmpcomp(ee),'*');hold on; 
                        set(ph,'markersize',5);
                    end;
                    set(ph,'color',cols(ee,:));
                    ph=text(ee,tmpcomp(ee),labels{ee});;
                    set(ph,'color',cols(ee,:));set(ph,'fontsize',6);
                    set(ph,'rotation',90);
                end;            
            end;
            set(gca,'xlim',[0 length(labels)+1]);            
            plot([1:size(cxttmpl,1)], cxttmpl(:,template)','k-','linewidth',1);hold on;             
            plot([get(gca,'xlim')],[0 0],'k-');
            set(gca,'ticklength',[tickl tickl]);
            set(gca,'fontsize',fs);
            if pl <= col+1
                title(['Dim ',int2str(template),' Context']);
            end;        
            if pl < (row-1)*col
                set(gca,'xticklabel',[]);
            end;
        end;
        
        sbplot(row,col,pl); pl = pl+1;
        if strcmp(plottype,'winv')
            if distr == 1
                % plot distribution:------------------------------------
                hist(acts(template,:),100); hold on;plot([0 0],[get(gca,'ylim')],'r-');
                title('Wts Distr.');
            else
                erptempl(:,template) = erptempl(:,template)*(sqrt(mean(s.mnrowpwr.^2)))/(sqrt(mean(erptempl(:,template).^2))); % normalize template to rms of mean erp
                ph = plot(s.times,s.mnrowpwr,'g-','linewidth',1.5); hold on;
                
                plot(s.times,s.mnrowpwr + erptempl(:,template),'r-','linewidth',2);hold on;               
                set(gca,'xlim',[s.times(1) s.times(end)]);
                set(gca,'ticklength',[tickl tickl]);
                hold on; plot([0 0],[get(gca,'ylim')],'k-');
                hold on; plot([get(gca,'xlim')],[0 0],'k-');
                set(gca,'fontsize',fs);  
                if pl <= col+1
                    title(['ERP']); 
                end;
            end;
        elseif strcmp(plottype,'wtdmean') 
            backproj = winv(:,template)*acts(template,:);
            backproj = erpeig*backproj(1:s.pcerp,:);
            mnerp = repmat(s.mnrowpwr,[1 size(backproj,2)]);
            %erpproj = mnerp + erpproj; % full back projection
            backproj = mnerp + backproj; % single dim back proj
            [x y] = sort(acts(template,:)); backproj = backproj(:,y);
            erptempl(:,template) = erptempl(:,template)*(sqrt(mean(wtderp(:,template).^2)))/(sqrt(mean(erptempl(:,template).^2))); % normalize template to rms of wtdmean
            if erpim == 1
                imlim = max(abs(backproj(:)));
                imagesc(s.times,[1:size(backproj,2)],backproj',[-imlim imlim]); 
                if pl < (row-1)*col
                    set(gca,'xticklabel',[]);
                end;
% $$$                 sbplot(row,col,pl); pl = pl+1;
% $$$                 for pt = 1:size(limerp,1)
% $$$                     ph = plot([s.times(pt) s.times(pt)],[limerp(pt,template,1) limerp(pt,template,2)],'k-','linewidth',3);
% $$$                     set(ph,'color',[.8 .8 .8]);hold on;
% $$$                 end;                
% $$$                 ph = plot(s.times,s.icmean,'g-','linewidth',1.5); hold on;% mean ERP
% $$$                 plot(s.times,  s.icmean+backproj(:,end)','linewidth',1.5);hold on; % max template proj            
% $$$                 plot(s.times,  erptempl(:,template)','r-','linewidth',1);hold on;             
% $$$                 set(gca,'ylim',[-max(abs(s.icmean+backproj(:,end)')) max(abs(s.icmean+backproj(:,end)'))]);
% $$$                 set(gca,'xlim',[s.times(1) s.times(end)]);
            else
% $$$                 for pt = 1:size(limerp,1)
% $$$                     ph = plot([s.times(pt) s.times(pt)],[limerp(pt,template,1) limerp(pt,template,2)],'k-','linewidth',3);
% $$$                     set(ph,'color',[.8 .8 .8]);hold on;
% $$$                 end;                
                
                z = env(backproj'); 
                t1 = [s.times,[s.times(end:-1:1)]];
                t2 = [z(1,:),[z(2,end:-1:1)]];
                ph = fill(t1,t2,[.75 .75 .75]);hold on;
                set(ph,'edgecolor',[.75 .75 .75]);set(ph,'edgealpha',0)
                ph = plot(s.times,s.mnrowpwr,'g-','linewidth',1.5); hold on;% mean ERP
                ph = plot(s.times,backproj(:,1)','b-','linewidth',1.5); hold on;% max template proj  
                ph = plot(s.times,backproj(:,end)','r-','linewidth',1.5); hold on;% max template proj  
%plot(s.times,  wtderp(:,template)','linewidth',2);hold on;             
%plot(s.times,  erptempl(:,template)','r-','linewidth',1);hold on; 
                set(gca,'ylim',[-max(abs(t2)) max(abs(t2))]);
                set(gca,'xlim',[s.times(1) s.times(end)]);
            end;
            set(gca,'ticklength',[tickl tickl]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            hold on; plot([get(gca,'xlim')],[0 0],'k-');
            set(gca,'fontsize',fs);  
            if pl <= col+1
                title(['ERP']); 
            end;
        end;        
        if pl <= (row-1)*col
            set(gca,'xticklabel',[]);
        end;
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 


