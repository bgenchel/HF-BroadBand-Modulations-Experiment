% calculates the power spectrum for weights/time for each subject included




function CoModPwrTrends(savedat,fullpaths,overlap,smoothfac);
    
    
    if overlap == 2
        srate = .5;
    elseif overlap == 1
        srate = 1;
    end;
    
% $$$     for nx = 1:length(fullpaths)
% $$$         if ~isempty(fullpaths{nx})
% $$$             s = load([fullpaths,savedat,'.mat']);  
% $$$             sph=floatread([fullpaths,savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
% $$$             wts=floatread([fullpaths,savedat,'.wts'],[s.pcs s.numtrials],[],0); 
% $$$             ws = wts*sph;     winv = pinv(ws);   
% $$$             clear wts sph ws
% $$$     
% $$$             win = 256; ovrlap = 128; fftwin = 256;
% $$$ 
% $$$             for im = 1:size(winv,2)
% $$$                 [pwr(:,im),freqs] = pwelch(winv(:,im),win,ovrlap,fftwin,srate);
% $$$                 pwr(:,im) = 10*log10(pwr(:,im));
% $$$             end;
% $$$             str = [fullpaths{nx},'CoModSpectra.mat pwr freqs']; eval(str);
% $$$         end;
% $$$     end;

    s = load([fullpaths,savedat,'.mat']);  
    sph=floatread([fullpaths,savedat,'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths,savedat,'.wts'],[s.pcs s.numtrials],[],0); 
    ws = wts*sph;     winv = pinv(ws);   
    cols = hsv(s.pcs);
    
    figure; row = round(sqrt(s.pcs)); col = ceil(sqrt(s.pcs));
    for im = 1:size(winv,2)
        if ~isempty(smoothfac)
            if smoothfac ~= 1
                xwidth = round(size(winv,1)/smoothfac);
            else 
                xwidth = 0;
            end;
            [outdata,outx] = movav(winv(:,im),0,xwidth,0);
            outdata = winv(:,im);
            outx = [1:size(winv,1)];
            sbplot(row,col,im)
            ph = plot([1:length(outdata)],outdata,'k-','linewidth',1); hold on;
            set(ph,'color',cols(im,:));
            set(gca,'xlim',[0 size(winv,1)]);
        else            
            sbplot(row,col,im); 
            ph = plot(winv(:,im));hold on;
            set(gca,'xlim',[0 size(winv,1)]);
            set(ph,'color',cols(im,:));            
        end;
        ph = plot([get(gca,'xlim')],[0 0],'k--','linewidth',2);
        set(ph,'color','y');
        title([fullpaths(end-4:end-1),' IM ',int2str(im)]);
    end;
    
