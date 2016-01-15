% Plot results from SnglTrialERSPDecomp.m with standard context conditions (Reward Twoback)
%
% PlotNewCtxResults(fullwts,erspPCAwts,erspwts,ctxwts,times,freqs)
%
% fullwts -- [matrix] rank*2 X rank*2?
% erspPCAwts -- [matrix] the time/freq points X PCA dim matrix from original dimension
%                        reduction of ERSP matrix
% erspwts -- [matrix] ERSP time/freq points X rank number
% ctxwts -- [matrix] context questions X rank number


function PlotNewCtxResults(fullwts,erspPCAwts,erspwts,ctxwts,times,freqs)
    
    
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
    ctxmat = ctxwts * fullwts(size(fullwts,1)/2+1:end,:); % 33 X 64 matrix ready to plot (64 = dimensions)
    
    erspmat = erspwts * fullwts(1:size(fullwts,1)/2,:); % take top half of matrix
    erspmat = erspPCAwts * erspmat; % 3024 (t-f points) X 64 (dimensions); ready to plot
    
    for template = 1:size(fullwts,2)
        if pl > row * col
            ph = textsc(['ERSP/Context decomposition: ERSP and context jointly decomposed before ICA'],'title'); set(ph,'fontsize',14);
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            axcopy
            figure; pl = 1;            
        end;        

        tmpcomp = reshape(erspmat(:,template),length(freqs),length(times));  %Good
        lim = max(max(tmpcomp))-.1*max(max(tmpcomp));
        sbplot(row,col,[pl pl+1]);
        imagesc(times,freqs,tmpcomp,[-lim lim]); pl = pl+2; %colormap('gray');
        set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
        hold on; plot([0 0],[get(gca,'ylim')],'k-');set(gca,'yscale','log');
        set(gca,'ytick',[5:5:freqs(end)]); set(gca,'yticklabel',{5 10 [] 20 [] 30 [] [] [] []});
        set(gca,'fontsize',fs);  title(['Dim: ',int2str(template)]); colorbar;

        newwinv = ctxmat(:,template);
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
