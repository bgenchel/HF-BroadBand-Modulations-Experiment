% Plots output from SpecDecompWT.m
%
%PlotSpecDecompWT(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,justone,auxpath,maps);
%    
% auxpath -- [string] if .mat is saved in a different folder from .set
%             this string should be the one pointing to the .mat
% maps -- ['maps','off'] if 'maps', then will plot scalp maps. 
% superpos -- ['ims', 'ics' or 'off'] 'ics' will plot a single axis with all comps 
%             superimposed.Not available when 'justone' is used.
%             'ims' will show each IC scalp map with super imposed IM templates.
%             'off' will plot all ICs and IMs plotted separately. 


function PlotSpecDecompWT(datset,datpath,plotcomps,factors,savedat,frqlim,superpos,justone,auxpath,maps);
    
    
    plotbackproj = 'off'; % for single IM plotting, plots all trials vs activations
    maxrows = 11; % max # of rows to plot before starting new fig
    lnwdth = 2;
    if ~exist('maps')
        maps = 'off';
    end;
    
    nlim = []; mlim = [];
    fprintf('\nLoading subject spectral decomposition data...\n');
    if strcmp(maps,'maps')% plot maps if requested
        EEG = pop_loadset(datset, datpath);
        EEG = eeg_checkset(EEG);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('auxpath')
        if ~isempty(auxpath)
            datpath = auxpath;
        end;
    end;
    s = load([datpath,savedat,'.mat']);   
    sph=floatread([datpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([datpath,savedat,'.wts'],[s.pcs s.pcs],[],0);           
    icamatall = floatread([datpath,savedat,'.fdt'],[s.pcs s.numframes],[],0);      
    ws = wts*sph;  
    winv = pinv(ws);    
    activations = ws*icamatall;      
    speceig = floatread([datpath,s.eigfile],[length(s.rowmeans) s.pcs],[],0);    
    specwts = speceig*winv;      
    winv = specwts; % overwrite ICA winv with ICA/PCA winv
    clear speceig specwts   
    
    if s.complist == 0 % channels
      s.complist = [1:size(activations,2)/length(s.freqs)];
    end;
    clear wts sph ws icamatall
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(plotcomps)
        plotcomps = s.complist;
    end;
    if isempty(factors)
        factors = [1:size(activations,1)];
    end;
    fr = find(s.freqs >= frqlim(1) & s.freqs <= frqlim(end));
    
    if isempty(justone)
        minl = min(min(activations(factors,:)))-abs(min(min(activations(factors,:))))*.01;
        maxl = max(max(activations(factors,:)))+abs(max(max(activations(factors,:))))*.01;
        if strcmp(superpos,'ics') % superimpose ic templates
            figure;row = round(sqrt(length(factors))); col = ceil(sqrt(length(factors)));
            cols = lines(length(plotcomps)); pl = 1;
            for tpp = 1:length(factors)
                lnfac = (lnwdth-2)/length(plotcomps);
                tp = factors(tpp);               
                sbplot(row,col,pl)
                for cp = 1:length(plotcomps)
                    rcp = find(plotcomps(cp) == s.complist); 
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(cp,:)); hold on;
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'b'); hold on;
                        set(ph,'color',cols(cp,:)); 
                    else % otherwise linear
                        ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                        hold on; 
                        set(ph,'color',cols(cp,:)); 
                    end;
                    lnwdth = lnwdth - lnfac;
                end;
                set(gca,'ylim',[minl maxl]); title(['IM ',int2str(factors(tpp))]);
                set(gca,'box','off');
                set(gca,'xgrid','on');
                if pl == (row-1)*col+1
                    xlabel('Frequency (Hz)'); ylabel('Relative Power');
                elseif pl > (row-1)*col+1
                    xlabel('Frequency (Hz)');
                end;  
                if pl <= col*(row-1)
                    set(gca,'xticklabel',[]);
                end;pl = pl+1;
            end; textsc(['Superimposed IC Templates for single IMs: ',datpath],'title');
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            set(gcf,'color','w');
            axcopy
        elseif strcmp(superpos,'ims')  % superimpose IM templates for each IC
            figure;
            row = round(sqrt(length(plotcomps)*2)); 
            col = ceil(sqrt(length(plotcomps)*2));
            if ~iseven(col)
                col = col+1; row = row-1;
            end;
            cols = lines(length(factors)); 
            pl = 1;
            for cp = 1:length(plotcomps)
                rcp = find(plotcomps(cp) == s.complist); 
                sbplot(row,col,pl); pl = pl+1;
                topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
                set(gca,'fontsize',12);  title(['IC ',int2str(plotcomps(cp))]);
                sbplot(row,col,pl); 
                for tpp = 1:length(factors)
                    lnfac = (lnwdth-2)/length(plotcomps);
                    tp = factors(tpp);               
                    if strcmp(s.freqscale,'quad') % quadratic spaced freqs
                        ph = quadplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,cols(tpp,:)); hold on;
                    elseif strcmp(s.freqscale,'log') % log spaced
                        ph = logplot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'b'); hold on;
                        set(ph,'color',cols(tp,:)); 
                    else % otherwise linear
                        ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); 
                        hold on; 
                        set(ph,'color',cols(tpp,:)); 
                    end;
                    %lnwdth = lnwdth - lnfac;
                end;
                set(gca,'ylim',[minl maxl]); title(['IM templates']);
                set(gca,'xgrid','on');
                if pl == (row-1)*col+2
                    xlabel('Frequency (Hz)'); ylabel('Relative Power');
                elseif pl > (row-1)*col+1
                    xlabel('Frequency (Hz)');
                end;  
                if pl <= col*(row-1)
                    set(gca,'xticklabel',[]);
                end;
                pl = pl+1;
            end; 
            %textsc(['Superimposed IM templates for single ICs: ',datpath],'title');
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            set(gcf,'color','w');
            axcopy
        else  
            figure;row = length(factors)+1; 
            if row > maxrows
                row = round(row/2);
                if row > maxrows
                    row = maxrows;
                end;             
            end;            
            col = length(plotcomps)+1;
            if strcmp(maps,'maps')% plot maps if requested
                pl = 2;           
                for cp = 1:length(plotcomps)
                    sbplot(row,col,pl)
                    topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off','plotrad',.7); pl = pl+1;
                    set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
                end;
            else
                pl = 1;
            end;
            
            for tpp = 1:length(factors)
                tp = factors(tpp);
                if pl == row*col+1
                    set(gcf,'Position',[100 300 1400 900]);
                    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                    axcopy
                    ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
                    set(ph,'fontsize',14);
                    figure;
                    if strcmp(maps,'maps')% plot maps if requested
                        pl = 2;
                        for cp = 1:length(plotcomps)
                            sbplot(row,col,pl)
                            topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
                            set(gca,'fontsize',7);  title(int2str(plotcomps(cp)));
                        end;
                    else
                        pl = 1;
                    end;
                end;
                sbplot(row,col,pl)
                hist(winv(:,tp),75);pl = pl+1;hold on;
                set(gca,'fontsize',7);   %set(gca,'xlim',[-2 2]);
                plot([0 0],[get(gca,'ylim')],'r-');
                set(gca,'yticklabel',[]);   set(gca,'xticklabel',[]); 
                title(['IM ',int2str(tp)]);
                for cp = 1:length(plotcomps)
                    rcp = find(ismember(s.complist,plotcomps(cp)));
                    sbplot(row,col,pl);                    
                    
                    if strcmp(s.freqscale,'quad') % quadratic spacing
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        quadplot(s.freqs(fr),tpact(:,fr),lnwdth,'b');
                        pl = pl+1;hold on;
                    elseif strcmp(s.freqscale,'log')  % log spacing
                        tpact = activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp);
                        logplot(s.freqs(fr),tpact(:,fr),lnwdth,'b');
                        pl = pl+1;hold on;
                    else % otherwise linear
                        plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',lnwdth); pl = pl+1;hold on;
                        %set(gca,'xtick',[10,40,100,200]);
                        %set(gca,'xtick',[10:10:frqlim(2)]);
                        %set(gca,'xticklabel',{10 [] 30 [] []});
                        set(gca,'xlim',[s.freqs(fr(1)) s.freqs(fr(end))]);   
                        % temp added:
                    end;                    
                    
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    set(gca,'ylim',[minl maxl]); 
                    set(gca,'xgrid','on');
                    set(gca,'fontsize',7);set(gca,'box','off');
                    set(gca,'ticklength',[.03 .03]);
                    plot([get(gca,'xlim')],[0 0],'r-');
                    if pl <= (row-1)*col+1
                        if tpp ~= length(factors)
                            set(gca,'xticklabel',[]);
                            set(gca,'yticklabel',[]);
                        end;
                    end;  
                    if ~strcmp(maps,'maps') & pl <= (col+1)
                        title(int2str(plotcomps(cp)));
                    end;
                end;
            end;
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            ph=textsc(['Spectral Templates for all good comps from ',datpath,': ',savedat],'title');
            set(ph,'fontsize',14);
            axcopy
        end;
    else % if plotting back-projection of just one IM:---
        
        figure;pl = 1;
        row = round(sqrt(length(plotcomps)*2))+1; col = round(sqrt(length(plotcomps)*2))+1;
        
        if ~iseven(col)
            col = col+1; row = row-1;
        end;  
        if strcmp(maps,'maps')% plot maps if requested
            
            for cp = 1:length(plotcomps)
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,plotcomps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+2;
                set(gca,'fontsize',10);
                title(int2str(plotcomps(cp)));
            end;
        end;
        pl = 2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cols = jet(length(plotcomps));
        
        if strcmp(plotbackproj,'on')
            x=[1:size(activations,1)];x(justone) = [];                
            activations(x,:) = 0;
            backproj = winv*activations ; 
            
            minl = min(backproj(:))-max(abs(s.meanpwr(:)));
            maxl = max(backproj(:))+max(abs(s.meanpwr(:)));
        else
            minl = floor(min(activations(justone,:))-abs(min(activations(justone,:)))*.01);
            maxl = ceil(max(activations(justone,:))+abs(max(activations(justone,:)))*.01);            
        end;
        
        for cmp = 1:length(plotcomps)
            rcp =  find(plotcomps(cmp) == s.complist);
            sbplot(row,col,pl)
            if strcmp(plotbackproj,'on')
                plotprj = backproj(:,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                plotprj = plotprj + repmat(s.meanpwr(rcp,:),[size(plotprj,1) 1]);
                if strcmp(s.freqscale,'quad') % quadratic spacing
                    ph = quadplot(s.freqs,plotprj,1,'k'); hold on;
                    set(ph,'color',cols(cmp,:));
                    ph = quadplot(s.freqs,s.meanpwr(rcp,:),lnwdth,'k'); 
                elseif strcmp(s.freqscale,'log') % quadratic spacing
                    ph = logplot(s.freqs,plotprj,1,'k'); hold on;
                    set(ph,'color',cols(cmp,:));
                    ph = logplot(s.freqs,s.meanpwr(rcp,:),lnwdth,'k'); 
                elseif strcmp(s.freqscale,'linear') % quadratic spacing
                    ph = plot(s.freqs,plotprj,'k-','linewidth',1); hold on;
                    set(ph,'color',cols(cmp,:));
                    ph = plot(s.freqs,s.meanpwr(rcp,:),'k-','linewidth',lnwdth); 
                end;
            else
                plotprj = activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp); 
                if strcmp(s.freqscale,'quad') % quadratic spacing
                    ph = quadplot(s.freqs,plotprj,lnwdth,'b'); hold on;
                elseif strcmp(s.freqscale,'log') % quadratic spacing
                    ph = logplot(s.freqs,plotprj,lnwdth,'b'); hold on;
                elseif strcmp(s.freqscale,'linear') % quadratic spacing
                    ph = plot(s.freqs,plotprj,'b-','linewidth',lnwdth); hold on;
                end;  
                set(gca,'xgrid','on')
            end;
            set(gca,'ylim',[minl maxl]);
% $$$             if logyes == 1
% $$$                 set(gca,'xscale','log');
% $$$             else
% $$$             end;
            fprintf('.');pl = pl+2; 
            if cmp ~= length(plotcomps)
                set(gca,'xticklabel',[]);
            end;            
            title(['IM ',int2str(justone),'; IC ',int2str(plotcomps(cmp))]);
        end;
        if pl < (row * col)+1
            sbplot(row,col,pl-1)
            hist(winv(:,justone),75);hold on; pl = pl+1;
            set(gca,'fontsize',10);
            plot([0 0],[get(gca,'ylim')],'r-'); 
            title('Factor Weights');
        end;
        if pl < (row * col)+1  % plot IMs superimposed at the end
            cols = jet(length(plotcomps));      
            sbplot(row,col,pl-1)

            for cp = 1:length(plotcomps)
                rcp =  find(plotcomps(cp) == s.complist);
                
                if strcmp(s.freqscale,'quad') % quadratic spacing
                    ph = quadplot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'k'); 
                elseif strcmp(s.freqscale,'log') % quadratic spacing
                    ph = logplot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),lnwdth,'k'); 
                elseif strcmp(s.freqscale,'linear') % quadratic spacing
                    ph = plot(s.freqs,activations(justone,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'k-','linewidth',lnwdth); 
                    
                end;                
                hold on; set(ph,'color',cols(cp,:));
            end;title(['IM ',int2str(justone)]); 
            set(gca,'ylim',[minl maxl]);
        end;
        set(gcf,'Position',[100 300 1400 900]);
        set(gcf,'PaperOrientation','landscape');   set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        ph=textsc(['Spectral Templates Factor ',int2str(justone),' for ',datpath,': ',savedat],'title');
        set(ph,'fontsize',14);
        %axcopy            
    end;
    
