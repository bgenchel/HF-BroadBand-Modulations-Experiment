% takes CoMod decomposition data from specified clusters and plots mean
% weights/time spectral decompositions
%
%
%
%
%
%
% ems -- [0|1] if 1 will plot all emotions separately for each cluster
%

function PlotClustWtsSpectra(savedat,fullpaths,facvec,ems,row,col,place)
    
    
    emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
    %win = 256; ovrlap = 128; fftwin = 256;srate = .5;
    win = 128; ovrlap = 64; fftwin = 128;srate = .5;
    nx = 1; im = 1;
    s = load([fullpaths{nx},savedat,'.mat']);  
    sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
    wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
    ws = wts*sph;     winv = pinv(ws);   %icawinv = icawinv';
    clear wts sph ws
    [Px,frqs] = pwelch(winv(:,im),win,ovrlap,fftwin,1);
    Px = zeros(0,length(frqs));
    Pxemos = zeros(0,length(frqs),length(s.dstrials));
    for nx = 1:length(facvec)
        if ~isempty(facvec{nx})
            fprintf('.');
            s = load([fullpaths{nx},savedat,'.mat']);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            ws = wts*sph;     winv = pinv(ws);   %icawinv = icawinv';
            clear wts sph ws
            allims = unique(facvec{nx});
            for im = 1:length(allims)
                if ems == 1
                    for em = 1:length(s.dstrials)
                        if length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em))) > win
                        [Pxx(:,em),frqs] = pwelch(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),allims(im)),win,ovrlap,fftwin,srate);
                        Pxx(:,em) = 10*log10(Pxx(:,em));
                        else
                            Pxx(:,em) = zeros(length(frqs),1);
                        end;
                    end;
                    Pxemos(end+1,:,:) = Pxx;
                else                    
                    [Px(end+1,:),frqs] = pwelch(winv(:,allims(im)),win,ovrlap,fftwin,srate);
                    Px(end,:) = 10*log10(Px(end,:));
                end;
                
            end;
        end;
    end;
    if ems == 1
        cols = jet(length(s.dstrials));
        for em = 1:length(s.dstrials)
            sbplot(row,col,place)
            plot(frqs,Pxemos(:,:,em),'k-'); hold on;
            ph = plot(frqs,mean(Pxemos(:,:,em),1),'r-','linewidth',2.5);
            set(ph,'color',cols(em,:));
            set(gca,'xlim',[frqs(1) frqs(end)]); set(gca,'xscale','log');  
            if place < 16
                title(emos{em});
            end;
            place = place + 1;
        end;
    else        
        sbplot(row,col,place)
        plot(frqs,Px,'g-'); hold on;
        plot(frqs,mean(Px,1),'r-','linewidth',2);
        set(gca,'xlim',[frqs(1) frqs(end)]);
        set(gca,'xscale','log');
    end;
    
        
