% Plot results from SnglTrialERSPDecomp.m with standard context conditions (Reward Twoback)
%
% PlotNewCtxWtdTrials(origersp,origctx,finalwts,times,freqs)
%
% origersp -- [matrix] full data matrix with ERSP information (t/f points X trials)
% origctx -- [matrix] full data matrix with context information (questions X trials)
% finalwts -- [matrix] resulting matrix from ICA decomposition. (64(dims) X trials)

function PlotNewCtxWtdTrials(origersp,origctx,finalwts,times,freqs)
    
    
    fs = 11;  % fontsize
    ms = 3;   % markersize 
    nlett = 6; nfb = 21; nresp = 6;   % numbers of questions in the 3 categories
    ntot = 33; % number of questions total
    row = 12; col = 6;
    fs = 7;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% PCA on  ERSP AND Context %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure; pl = 1; ms = 4;
    % wt each trial of ersp matrix by finalwts
    erspmat = origersp *finalwts' ; % t/f points X 64 matrix
    ctxmat = origctx * finalwts' ; % 33 questions X 64 matrix
   
    %%%%%%%%%%%%%%%%%
    % make signif masks for ersp trial plotting
    shuffnum = 500;
    randvec = [1:size(finalwts,2)]; 
    for rep = 1:shuffnum
        randwts = shuffle(randvec); 
        shufferspmat = origersp *finalwts(:,randwts)' ; % t/f points X 64 matrix
        shuffctxmat = origctx * finalwts(:,randwts)' ; % 33 questions X 64 matrix
        randersps(:,:,rep) = shufferspmat;
        randctxs(:,:,rep) = shuffctxmat;
    end;
    
    for template = 1:size(fullwts,2)
        for tf = 1:size(randersps,1)        
            currvec = squeeze(randersps(tf,template,:));currvec = sort(currvec);
            emask(1,:) = currvec(round(shuffnum*.01)); % min 1% of data
            emask(2,:) = currvec(length(currvec) - round(shuffnum*.01));% max 1% of data
            currvec = squeeze(randctxs(tf,template,:));currvec = sort(currvec);
            cmask(1,:) = currvec(round(shuffnum*.01)); % min 1% of data
            cmask(2,:) = currvec(length(currvec) - round(shuffnum*.01));% max 1% of data
        end;
        
        if pl > row * col
            ph = textsc(['ERSP/Context decomposition: ERSP and context jointly decomposed before ICA'],'title'); set(ph,'fontsize',14);
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            axcopy
            figure; pl = 1;            
        end;        
        tmpcomp = erspmat(:,template); 
        tmpcomp(find(tmpcomp>emask(1,:)&tmpcomp<emask(2,:))) = 0;

        tmpcomp = reshape(erspmat(:,template),length(freqs),length(times));  %Good
        lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
        sbplot(row,col,[pl pl+1]);
        imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+2; %colormap('gray');
        set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
        hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
        set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
        set(gca,'fontsize',fs);  title(['Dim: ',int2str(template)]); colorbar;

        newwinv = ctxmat(:,template);
        newwinv(find(newwinv>emask(1,:)&newwinv<emask(2,:))) = 0;
        sbplot(row,col,[pl pl+3])
        ph = plot([1:nlett],newwinv(1:nlett),'ro-','linewidth',2);hold on;
        set(ph,'markersize',ms);
        ph = plot([nlett+1:nlett+nfb],newwinv(nlett+1:nlett+nfb),'bs-','linewidth',2); hold on;
        set(ph,'markersize',ms);
        ph = plot([nlett+nfb+1:nlett+nfb+nresp],newwinv(nlett+nfb+1:nlett+nfb+nresp),'g^-','linewidth',2);
        set(ph,'markersize',ms+ms/2); pl = pl+4;hold on;
        set(gca,'xlim',[0 ntot+1]);%set(gca,'ylim',[-ctlim ctlim]);
        ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.5 .5 .5]);
        yl = get(gca,'ylim');

        if pl > (row-1)*col 
            CnxtLabelRewFB(max(yl),.3,fs); 
        end;

    end;
    ph = textsc(['ERSP/Context decomposition: ERSP and context jointly decomposed before ICA'],'title'); set(ph,'fontsize',14);
    set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    axcopy
