% made to plot results of SpecDecomp() with no scalp maps
%
%  PlotSpecDecomp(datpath,complist,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone);
%
% INPUTS:
% datpath: full directory path where dataset can be found
% complist: full list of components/channels in orginal decomposition
% plotcomps: list of components/channels to plot (must be in original 'complist')(default: all)
% factors: list of factor numbers to plot (default: all)
% savedat: filename that float files were saved as (same as input to SpecCoModAnal.m)
% frqlim: [minfreq maxfreq] to plot
% superpos -- ['y' or 'n'] 'y' will plot a single axis with all comps superimposed.
%             (can only be used for plotting all comps (justone needs to be []) 
% meanvar -- [0 | 1] if 1, will plot mean spectra and trial backprojections, time intensive
% justone: (integer) IM number to plot alone. [] to plot all. Plots mean with spectral variations on top

function PlotSpecDecomp(datpath,complist,plotcomps,factors,savedat,frqlim,superpos,meanvar,justone);
    
    
    nlim = []; mlim = [];
    s = load([datpath,savedat,'.mat']);     
    for ff = 1:length(s.freqs)-1
        fspace(1,ff) = s.freqs(ff+1) - s.freqs(ff);
    end;
    if fspace(4) > fspace(1) & fspace(10) > fspace(1)
        logyes = 1;
    else
        logyes = 0;
    end;
    sph=floatread([datpath,savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    icamatall = floatread([datpath,savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws); clear wts sph ws
    if isempty(plotcomps)
        plotcomps = complist;
    end;
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    fr = find(s.freqs > frqlim(1) & s.freqs < frqlim(end));
    if isempty(justone)
        minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
        maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
        if superpos == 'y'
            figure;row = round(sqrt(length(factors))); col = ceil(sqrt(length(factors)));
            cols = jet(length(plotcomps)); pl = 1;
            for tpp = 1:length(factors)
                lnwdth = 2.5;
                lnfac = (lnwdth-2)/length(plotcomps);
                tp = factors(tpp);               
                sbplot(row,col,pl)
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == complist);
                    ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                    hold on; set(ph,'color',cols(cp,:));                    
                    lnwdth = lnwdth - lnfac;
                end;pl = pl+1;
                set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                set(gca,'ylim',[minl maxl]); title(['Fac ',int2str(factors(tpp))]);
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                set(gca,'box','off');
                    if logyes == 1
                        set(gca,'xscale','log');set(gca,'xtick',[4 6 10 20 40]);
                        set(gca,'xtick',[5,10,30:20:frqlim(2)]);  
                        set(gca,'xticklabel',{5 10 30 [] []}); 
                    else
                        set(gca,'xtick',[10:5:frqlim(2)]);
                        %set(gca,'xticklabel',{10 [] 30 [] []});
                    end;
                set(gca,'xgrid','on');
                if tpp == length(factors)
                    xlabel('Frequency (Hz)');
                    ylabel('Relative Power');
                end;                
            end; textsc(['Superimposed Factor Spectra: ',datpath],'title');
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            set(gcf,'color','w');
            axcopy
        else  
            figure;row = length(factors); 
            if row > 15
                row = 15;
            end;            
            col = length(plotcomps)+1;
            pl = 1;           
            for tpp = 1:length(factors)
                tp = factors(tpp);
                if pl > row*col
                    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                    ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                    set(ph,'fontsize',14);
                    figure; pl = 1;
                end;
                sbplot(row,col,pl)
                hist(winv(:,tp),75);pl = pl+1;hold on;
                set(gca,'fontsize',7);set(gca,'xlim',[-2 2]);
                plot([0 0],[get(gca,'ylim')],'r-');
                set(gca,'yticklabel',[]);                
                set(gca,'xticklabel',[]); 
                title(['IM ',int2str(tp)]);
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == complist);
                    sbplot(row,col,pl)
                    plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); pl = pl+1;hold on;
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'ylim',[minl maxl]); 
                    if logyes == 1
                        set(gca,'xscale','log');set(gca,'xtick',[4 6 10 20 40]);
                        set(gca,'xtick',[5,10,30:20:frqlim(2)]);  
                        set(gca,'xticklabel',{5 10 30 [] []}); 
                    else
                        set(gca,'xtick',[10:5:frqlim(2)]);
                        %set(gca,'xticklabel',{10 [] 30 [] []});
                    end;
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                    %set(gca,'xgrid','on');
                    set(gca,'ticklength',[.03 .03]);
                    if pl <= (col+1)
                        title(int2str(rcp));
                    end;                    
                    if tpp ~= length(factors) %| tpp ~= row-1
                        set(gca,'xticklabel',[]);
                        set(gca,'yticklabel',[]);
                    end;                    
                end;
            end;
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates for ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
        end;
    else
        minl = min(min(activations(justone,:)))-abs(min(min(activations(justone,:))))*.01;
        maxl = max(max(activations(justone,:)))+abs(max(max(activations(justone,:))))*.01;
       figure;row = round(sqrt(length(plotcomps)*2)); col = round(sqrt(length(plotcomps)*2))+1;pl = 1;
        if ~iseven(col)
            col = col+1; row = row-1;
        end;       
        pl = 2;
        if meanvar == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %z = load('/data/common4/emotion/KmeansClustCoMods.mat' );
            %s = load([datpath,'SpecCoModStuff.mat']);
            cols = jet(length(plotcomps));
            x=[1:size(activations,1)];x(justone) = [];                
            acts = activations; acts(x,:) = 0;
            backproj = winv*acts ; clear plotprj
            for cmp = 1:length(plotcomps)
                rcp =  find(plotcomps(cmp) == complist);
                sbplot(row,col,pl)
                onebkprj = backproj(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                for trls = 1:size(onebkprj,1)
                    plotprj(trls,:) = onebkprj(trls,:) + s.meanpwr(rcp,:);
                end;
                if isempty(nlim)
                    nlim = min(plotprj(:))+(min(plotprj(:))*.1); 
                    mlim = max(plotprj(:))+(max(plotprj(:))*.1); 
                end;
                ph = plot(freqs,plotprj,'k-'); hold on;
                set(ph,'color',cols(cmp,:));
                ph = plot(s.freqs,s.meanpwr(rcp,:),'k-'); 
                set(ph,'linewidth',2);set(gca,'xlim',[3 50]);
                set(gca,'ylim',[nlim mlim]);
                    if logyes == 1
                        set(gca,'xscale','log');set(gca,'xtick',[4 6 10 20 40]);
                        set(gca,'xtick',[5,10,30:20:frqlim(2)]);  
                        set(gca,'xticklabel',{5 10 30 [] []}); 
                    else
                        set(gca,'xtick',[10:10:frqlim(2)]);
                        set(gca,'xticklabel',{10 [] 30 [] []});
                    end;
                fprintf('.');pl = pl+2; 
                if cmp ~= length(plotcomps)
                    set(gca,'xticklabel',[]);
                end;            
                title(['IM ',int2str(justone),';CH ',int2str(plotcomps(cmp))]);
            end;
            if pl < (row * col)+1
                sbplot(row,col,pl)
                hist(winv(:,justone),75);hold on; pl = pl+1;
                set(gca,'fontsize',10);%set(gca,'xlim',[-1.5 1.5]);
                plot([0 0],[get(gca,'ylim')],'r-'); 
                %set(gca,'yticklabel',[]); set(gca,'xticklabel',[]);
                title('Factor Weights');
            end;
            if pl < (row * col)+1
                cols = jet(length(plotcomps));      
                sbplot(row,col,pl)
                for cp = 1:length(plotcomps)
                    rcp =  find(plotcomps(cp) == complist);
                    ph = plot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                    hold on; set(ph,'color',cols(cp,:));
                    set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                end;title(['IM ',int2str(justone)]); 
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);
                set(gca,'ylim',[min(activations(justone,:))-.1 max(activations(justone,:))+.1]);
            end;
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates for IM ',int2str(justone),'; ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % OLD style:      
        else
            
            tp = justone;
            ymin = min(min(activations(tp,:)))-.05*min(min(activations(tp,:)));
            ymax = max(max(activations(tp,:)))+.05*max(max(activations(tp,:)));
            for cp = 1:length(plotcomps)
                rcp =  find(plotcomps(cp) == complist);
                sbplot(row,col,pl)
                ph = plot(freqs,activations(tp,length(freqs)*(rcp-1)+1:length(freqs)*rcp)); pl = pl+2;hold on;
                set(ph,'linewidth',2);
                %plot([0 0],[get(gca,'ylim')],'g'); 
                set(gca,'fontsize',10);set(gca,'box','off');
                set(gca,'xlim',frqlim); set(gca,'ylim',[ymin ymax]);   
                    if logyes == 1
                        set(gca,'xscale','log');set(gca,'xtick',[4 6 10 20 40]);
                        set(gca,'xtick',[5,10,30:20:frqlim(2)]);  
                        set(gca,'xticklabel',{5 10 30 [] []}); 
                    else
                        set(gca,'xtick',[10:10:frqlim(2)]);
                        set(gca,'xticklabel',{10 [] 30 [] []});
                    end;
                set(gca,'xgrid','on');
            end;
            if pl < (row * col)+1
                sbplot(row,col,pl)
                hist(winv(:,tp),75);hold on; pl = pl+2;
                set(gca,'fontsize',10);%set(gca,'xlim',[-1.5 1.5]);
                plot([0 0],[get(gca,'ylim')],'r-'); 
                %set(gca,'yticklabel',[]); set(gca,'xticklabel',[]);
                title('Factor Weights');
            end;
            if pl < (row * col)+1
                cols = jet(length(plotcomps));      
                sbplot(row,col,pl:pl+1)
                for cp = 1:length(plotcomps)
                    rcp =  find(plotcomps(cp) == complist);
                    ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)); 
                    
                    hold on; set(ph,'color',cols(rcp,:));
                end;title(['Fac ',int2str(tp)]); 
                set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);
                set(gca,'ylim',[minl maxl]); 
            end;
            set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates Factor ',int2str(justone),' for ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end;
    end;
    
