% Calls in ITC info from each subject for selected components and plots avg
% includes optional masking for binomial probability across subjects
%
% [eallersps,times,freqs,colorlim] = PlotSigAvgERSP(filename,paths,origlist,clustcomps,freqlims,timelims,clim,binomalpha);
% 
% datafile: name of datafile in respective subject paths (.mat). (string)
% paths: cell array of strings with full paths to datafiles.
% clustcomps: cell array of component numbers to include in average from each subject.
% freqlims: [minfreq maxfreq]. frequencies in Hz to plot. (integers)
% timelims: [mintime maxtime]. times in ms to plot. If 0 is included, a black line will be drawn at 0. (integers)
% binomalpha: if not empty, will mask resulting avg ITC using binomial
%            probability at the specified alpha level (p value: decimal) 
%

function [eallersps,times,freqs,colorlim] = PlotSigAvgERSP(filename,paths,origlist,clustcomps,freqlims,timelims,clim,binomalpha);
    
cd (paths{1});   s = load (filename);
     fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
    tms=find(s.times > timelims(1) & s.times<timelims(2));

    eallersps = zeros(length(fr),length(tms),0);
    fprintf('\nsubjs done: '); 
    for nx = 1:length(paths)
        fprintf('..%s..',int2str(nx));            
        if~isempty(clustcomps{nx})
        s = load ([paths{nx},filename]);
           for k=1:length(clustcomps{nx})
                iminmask = s.itc_boot(fr,find(origlist{nx} == clustcomps{nx}(k)));
                iminmask = repmat(iminmask,[1 length(tms)]);
                iersp = s.comp_itc(fr,tms,find(origlist{nx} == clustcomps{nx}(k)));
                iersp = sqrt(real(iersp).^2 + imag(iersp).^2);
                iersp(find(iminmask >= iersp)) = 0;
                eallersps(:,:,end+1) = iersp; 
            end;            
        end;
    end;
    if ~isempty(binomalpha)
        % mask for binomial probability across components
        for k = 1:30
            W = factorial(size(eallersps,3))/(factorial(k)*factorial(size(eallersps,3)-k))*(.01^k)*(1-.01^(size(eallersps,3)-k) );
            if W< binomalpha
                break
            end;
        end;
        plotmat = eallersps;
        for r = 1:size(plotmat,1)
            for ss = 1:size(plotmat,2)
                nonsig = plotmat(find(plotmat(r,ss,:)==0));
                if length(nonsig) > size(plotmat,3) - k
                    plotmat(r,ss,:) = 0;
                end;
            end;        
        end;
    else
        plotmat = eallersps;
    end;    
    if isempty(clim)
        colorlim = max(max( abs(mean(plotmat,3))));
    else
        colorlim = clim;
    end;
    imagesc(s.times(tms),s.freqs(fr), mean(plotmat,3),[-colorlim colorlim]);
    set(gca,'ydir','norm'); hold on;
    plot([0 0],[get(gca,'ylim')],'k-');
    
    times = s.times(tms); freqs=s.freqs(fr);