% Plot results from SnglTrialERSPDecomp.m using NO context
%
% PlotTWsnglTrialERSP(pathname,savedat,savename,plotmax);
%
% INPUTS:
% pathname: [string] full path where wts, sph can be found and savename will be saved
% wtsphname: [string] Name of weights/sphere float file that can be found in pathname (same as name
%                     given to SnglTrialERSPDecomp.m)
% savename: [string] or [] If not empty, will save the generated plot as a jpg in the pathname directory
% plotfacs: [vector] Indices of dimensions to plot. [] for all.
% plotpcs: ['on' or 'off'] if 'on', will plot PCA ERSP templates for comparison {'off'}
%
%

function PlotTWsnglTrialERSP(pathname,savedat,savename,plotfacs,plotpcs);
    
    fs = 14; % fontsize
%%%  Do this for all  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = load([pathname,savedat,'.mat']);     
    wts = floatread([pathname,savedat,'.wts'],[s.pcmat s.pcmat],[],0);
    sph = floatread([pathname,savedat,'.sph'],[s.pcmat s.pcmat],[],0);     
    ws = wts*sph;winv = pinv(ws);
    if isfield(s,'pcmat')
        if s.pcmat > 0               
            erspeig = floatread([pathname,savedat,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcmat],[],0);
            winv = erspeig*winv(1:s.pcmat,:);  % just ERSP templates   
        end;
    end;
    %%%%%%  Check for log frequency spacing  %%%%%%%%%%%%%%%%
    
    if isfield(s,'freqscale')
        freqscale = s.freqscale;
    else 
        fprintf('\nError: .mat files must include a "freqscale" field\n'); return;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if ~exist('plotfacs')
        plotfacs = [1:size(winv,2)];
    elseif isempty(plotfacs)
        plotfacs = [1:size(winv,2)];
     end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if ~exist('plotpcs')
        plotpcs = 'off';
    elseif isempty(plotpcs)
        plotpcs = 'off';
     end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%% Plot Reg ERSP factors onto one page %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    row = ceil(sqrt(length(plotfacs)));
    if row > 5
        row = 6; 
    end;
    col = round(sqrt(length(plotfacs)));
    if col > 4
        col = 5;
    end;
    figure;pl = 1; ms = 4;
    for tts = 1:length(plotfacs)
        template = plotfacs(tts);
        if pl > row * col 
            ph = textsc([savedat,'; ICA ERSP Templates'],'title'); set(ph,'fontsize',fs);
            figure; pl=1;
        end;            
        cnt = 0;
        tmpcomp = winv(1:length(s.freqs)*length(s.times),template)';  %makes a freqs*times X 1
        lim = max(abs(tmpcomp))+.1*max(abs(tmpcomp));
        tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
        sbplot(row,col,pl); 
        if strcmp(freqscale,'quad')
            quadimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
        elseif strcmp(freqscale,'log')                
            mylogimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
        else
            imagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;               
            set(gca,'ytick',[10:10:s.freqs(end)]);
        end;
        set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
        hold on; plot([0 0],[get(gca,'ylim')],'k-');
        set(gca,'fontsize',fs);  title(['IC ',int2str(template)]); 
        if isfield(s,'medwarpevs')
            cols = lines(length(s.medwarpevs));
            for wp = 1:length(s.medwarpevs)
                ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                set(ph,'color',cols(wp,:));
            end;
        end;  
        if tts == row*col | tts == length(plotfacs)
            cbar;
        end;pl = pl+1;
    end;
    ph = textsc([savedat,'; ICA ERSP Templates'],'title'); set(ph,'fontsize',fs);
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
    if ~isempty(savename)
        makesave = [pathname,savename];    eval(makesave);   
    end;

    if strcmp(plotpcs,'on')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%%%% Plot PCA ERSP templates onto one page %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        row = ceil(sqrt(length(plotfacs)));
        if row > 5
            row = 6; 
        end;
        col = round(sqrt(length(plotfacs)));
        if col > 4
            col = 5;
        end;
        figure;pl = 1; ms = 4;
        for tts = 1:length(plotfacs)
            template = plotfacs(tts);
            if pl > row * col 
                ph = textsc([savedat,'; PCA ERSP Templates'],'title'); set(ph,'fontsize',fs);
                figure; pl=1;
            end;            
            cnt = 0;
            tmpcomp = erspeig(1:length(s.freqs)*length(s.times),template)';  %makes a freqs*times X 1
            lim = max(abs(tmpcomp))+.1*max(abs(tmpcomp));
            tmpcomp = reshape(tmpcomp,length(s.freqs),length(s.times));  %Good
            sbplot(row,col,pl); 
            if strcmp(freqscale,'quad')
                quadimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
            elseif strcmp(freqscale,'log')                
                mylogimagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;
            else
                imagesc(s.times, s.freqs, tmpcomp,[-lim lim]);hold on;               
                set(gca,'ytick',[10:10:s.freqs(end)]);
            end;
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            set(gca,'fontsize',fs);  title(['PC ',int2str(template)]); 
            if isfield(s,'medwarpevs')
                cols = lines(length(s.medwarpevs));
                for wp = 1:length(s.medwarpevs)
                    ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                    set(ph,'color',cols(wp,:));
                end;
            end;  
            if tts == row*col | tts == length(plotfacs)
                cbar;
            end;pl = pl+1;
        end;
        ph = textsc([savedat,'; PCA ERSP Templates'],'title'); set(ph,'fontsize',fs);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
        
    end;
