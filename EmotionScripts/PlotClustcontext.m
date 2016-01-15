% takes context info from TW weightings for individual ICs and plots scalp map,template and context
%
%
%
%
% plottype -- ['winv' or 'wtdmean'] for 'comb' decomp only
% alpha -- [pvalue] for bootstrap masking of weighted mean data (ERSP and context)
% OUTPUT:
% stdevs -- [3d matrix] gives standard deviations from the mean for each question (rows) and dim (col)
%                       quests/dims that are within stdcut are zero. third dimension is ICs

function [cxtout] = PlotClustcontext(filename,fullpath,plottype,labels,alpha,savettl)
    
    
    pvafcalc = 0; % plot pvaf of context instead
    
    cxtout = [];
    shuffnum = 500; % number of bootstrap iterations for wtd mean masking

    col = 5;   row = 6;
    
    s = load([fullpath,filename,'.mat']); s.pcs = s.channo;    
    if s.pcs==0
        s.pcs= s.channo;    
    end;
    wts = floatread([fullpath,filename,'.wts'],[s.pcs s.channo],[],0);
    sph = floatread([fullpath,filename,'.sph'],[s.channo s.channo],[],0);     
    dat = floatread([fullpath,filename,'.fdt'],[s.channo s.numframes],[],0);   
    ws = wts*sph;winv = pinv(ws);acts = ws*dat;
    if isfield(s,'pcersp')
        if s.pcersp > 0               
            erspeig = floatread([fullpath,filename,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcersp],[],0);
            erspdat = erspeig*dat(1:s.pcersp,:);
            ersptmpl = erspeig*winv(1:s.pcersp,:);            
        end;
    end;
    if isfield(s,'pcctx')
        if s.pcctx > 0
            cxteig = floatread([fullpath,filename,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
            cxtdat = cxteig*dat(s.pcersp+1:end,:);  
            cxttmpl = cxteig*winv(s.pcersp+1:end,:);  
        end;
    end;
    
    pl = 1; pg = 1;
    cols = hsv(length(labels)); cols(3,:) = [.8 .8 0];
    figure; 
    if strcmp(plottype,'wtdmean')
        fprintf('\nAccumulating bootstrap distribution ...\n');
        clear limersp limmat 
        %%%%%%%  calculate weighted ERSPs for each dim:
        wtdersp = (erspdat*acts')/size(erspdat,2);% multiply each trial element,sum then divide by ntrials 
        bootwts = zeros(size(wtdersp,1),size(wtdersp,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (erspdat*shuffle(acts,2)')/size(erspdat,2);
        end;
        bootwts = sort(bootwts,3);
        limersp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limersp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
        clear bootwts
        %%%%%%%  calculate context vectors for each dim:
        wtdctx = (cxtdat*acts')/size(erspdat,2); 
        %wtdctx = (addmat*acts')/size(erspdat,2); 
        bootwts = zeros(size(wtdctx,1),size(wtdctx,2),shuffnum);
        for b= 1:shuffnum
            bootwts(:,:,b) = (cxtdat*shuffle(acts,2)')/size(erspdat,2);
        end;
        bootwts = sort(bootwts,3);
        limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
        limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
                                                            %mnctx = mean(bootwts,3);
                                                            %stdctx = std(bootwts,0,3);
    end;
    clear sph wts ws dat bootwts 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for d = 1:size(winv,2)
        if strcmp(plottype,'wtdmean')
            goahead = 0;
            tmpcomp = wtdctx(:,d);
            tmpcomp(find(tmpcomp>limmat(:,d,1)&tmpcomp<limmat(:,d,2))) = 0;
            sigchk = (limmat(:,d,1)-tmpcomp).*(limmat(:,d,2)-tmpcomp);
            if ~isempty(find(sigchk>0))
                goahead = 1;
            else
                goahead = 0;
            end;
        else
            goahead = 1;
        end;
        if goahead == 1
            if pl > row * col
                textsc([filename,'; ',fullpath(end-8:end-1)],'title');
                set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                if ~isempty(savettl)
                    str = ['print ',fullpath,savettl,int2str(pg),'.jpg -djpeg']; eval(str)
                    clf;
                else
                    figure;
                end;
                pl = 1; pg = pg+1;
            end;
            
            sbplot(row,col,pl); pl = pl+1;
            if strcmp(plottype,'winv')
                tmpcomp =  ersptmpl(:,d);
            elseif strcmp(plottype,'wtdmean')
                tmpcomp = wtdersp(:,d);
                % mask the wtd mean by the bootstrap values:
                tmpcomp(find(tmpcomp>limersp(:,d,1)&tmpcomp<limersp(:,d,2))) = 0;
            end; 
            lim = max(abs(tmpcomp));
            if lim == 0
                lim = 1;% if no sig points, make lim a valid value
            end;
            tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  
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
            if isfield(s,'medwarpevs')
                wcols = lines(length(s.medwarpevs));
                for wp = 1:length(s.medwarpevs)
                    ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                    set(ph,'color',wcols(wp,:));
                end;
            end;  
            if pl < (row-1)*col
                set(gca,'xticklabel',[]);
            end;
            title(['Dim ',int2str(d)]);
            
            sbplot(row,col,[pl pl+3]); pl = pl+4;
            if strcmp(plottype,'winv')
                allctxvals = cxttmpl(:,d);
                ph=plot(allctxvals,'k-'); hold on;
                for ee = 1:length(allctxvals)
                    ph=plot(ee,allctxvals(ee),'.'); hold on;
                    set(ph,'markersize',15);                         
                    set(ph,'color',cols(ee,:));
                    ph=text(ee,allctxvals(ee),labels{ee});;
                    set(ph,'color',cols(ee,:));set(ph,'fontsize',14);
                    set(ph,'rotation',45);
                end;
                set(gca,'xlim',[0 length(allctxvals)+1]);
            elseif strcmp(plottype,'wtdmean')
                if pvafcalc == 1% calculate pvaf instead:---
                                % first make back proj of single template/factor and only cxt pcs:
                    bp = winv(s.pcersp+1:end,d) * acts(d,:);
                    bpq = cxteig * bp;
                    % calculate pvaf of PCA-reduced data
                    for q = 1:size(bpq,1)
                        pv(1,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
                    end;    
                    pv = pv*100;
                    tmpcomp = pv .* (cxttmpl(:,d)'./abs(cxttmpl(:,d)'));
                    for b = 1:shuffnum
                        bp = winv(s.pcersp+1:end,d) * shuffle(acts(d,:),2);
                        bpq = cxteig * bp;
                        %bpq = cxteig * shuffle(bp,2);
                        for q = 1:size(bpq,1)
                            pvboot(b,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
                        end;   
                        pvboot(b,:) = pvboot(b,:).* (cxttmpl(:,d)'./abs(cxttmpl(:,d)'));
                    end;    clear limmat             
                    pvboot = pvboot*100;
                    pvboot = sort(pvboot,1);
                    limmat(2,:) = pvboot(end-ceil(shuffnum*alpha),:); % max boot
                    limmat(1,:) = pvboot(ceil(shuffnum*alpha),:);  % min boot
                    for ee = 1:size(wtdctx,1)
                        ph = plot([ee ee],[limmat(1,ee) limmat(2,ee)],'k-','linewidth',3);
                        set(ph,'color',[.8 .8 .8]);hold on;
                    end;
                    sigchk = (limmat(1,:)-tmpcomp).*(limmat(2,:)-tmpcomp);
                    title(['Pvaf of Context']);
                else
                    tmpcomp = wtdctx(:,d);
                    tmpcomp(find(tmpcomp>limmat(:,d,1)&tmpcomp<limmat(:,d,2))) = 0;
                    for ee = 1:size(wtdctx,1)
                        ph = plot([ee ee],[limmat(ee,d,1) limmat(ee,d,2)],'k-','linewidth',3);
                        set(ph,'color',[.8 .8 .8]);hold on;
                    end;
                    sigchk = (limmat(:,d,1)-tmpcomp).*(limmat(:,d,2)-tmpcomp);
                    title(['Wtd Mean Context']);
                end;
                ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.75 .75 .75]);                
                cxttmpl(:,d) = cxttmpl(:,d)*(sqrt(mean(tmpcomp.^2)))/(sqrt(mean(cxttmpl(:,d).^2))); % normalize template to rms of wtdmean
                plot([1:size(cxttmpl,1)], cxttmpl(:,d)','k-','linewidth',1);hold on;             
                for ee = 1:size(wtdctx,1)
                    ph=plot(ee,tmpcomp(ee),'.');
                    set(ph,'markersize',18); 
                    set(ph,'color',cols(ee,:));
                    ph=text(ee,tmpcomp(ee),labels{ee});
                    set(ph,'color',cols(ee,:));set(ph,'fontsize',9);
                    set(ph,'rotation',45);
                end;                
                set(gca,'xlim',[0 size(wtdctx,1)+1]);
                
                if ~isempty(find(sigchk>0))
                    cxtout(:,d) = tmpcomp';
                else
                    cxtout(:,d) = zeros(length(tmpcomp),1); % if no context out of pvaf range
                end;                
            end;            
            plot([get(gca,'xlim')],[0 0],'k-');
            if pl < (row-1)*col
                set(gca,'xticklabel',[]);
            end;
        end;                
    end;
    textsc([filename,'; ',fullpath(end-8:end-1)],'title');
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    if ~isempty(savettl)
        str = ['print ',fullpath,savettl,int2str(pg),'.jpg -djpeg']; eval(str)
    end;
        
