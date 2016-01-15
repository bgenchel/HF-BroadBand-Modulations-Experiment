% Calls in ERSP info from each subject for selected components and plots avg
% includes optional masking for binomial probability across subjects
%
% [eallersps,times, freqs,colorlim] = PlotSigAvgERSP(filename,paths,origlist,clustcomps,freqlims,timelims,clim,binomalpha,sigdiff);
% 
% datafile: name of datafile in respective subject paths (.mat). (string)
% paths: cell array of strings with full paths to datafiles.
% origclust: cell array of vectors with original component lists for each subject (provides index into matrices)
% clustcomps: cell array of component numbers to include in average from each subject.
% freqlims: [minfreq maxfreq]. frequencies in Hz to plot. (integers)
% timelims: [mintime maxtime]. times in ms to plot. If 0 is included, a black line will be drawn at 0. (integers)
% clim -- [integer] max and abs(min) color scale value (only enter one number and will scale around 0.
% binomalpha: if not empty, will mask resulting avg ERSP using binomial
%            probability at the specified alpha level (p value: decimal) 
% sigdiff: 1 if dealing with two condition ERSP data, 0 otherwise

function [allersps,times, freqs,colorlim] = PlotSigAvgERSP(filename,paths,origlist,clustcomps,freqlims,timelims,clim,binomalpha,sigdiff,vertlines);
    
    
    
    if length(clim) > 1
        fprintf('Colorlim must be a single value (will scale around zero automatically)');
    end;
    
    
    if ~exist('sigdiff')
        sigdiff = 0; % by default, single condition ersp
    end;    
    if sigdiff == 1
         cd (paths{1});   s = load (filename);
         fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
         tms=find(s.times > timelims(1) & s.times<timelims(2));
       fprintf('\nsubjs done: '); 
        for cond = 1:3
            eallersps = zeros(length(fr),length(tms),0);
            for nx = 1:length(paths)
              fprintf('..%s..',int2str(nx))           
                if ~isempty(clustcomps{nx})
                  cd (paths{nx});   s = load (filename);
                  fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
                  tms=find(s.times > timelims(1) & s.times<timelims(2));

                for k=1:length(clustcomps{nx})
                    eformask = s.ersp_boot{find(origlist{nx} == clustcomps{nx}(k)),cond};% fr X 2
                    if cond ~= 3
                        eminmask = eformask(fr,1);
                        emaxmask = eformask(fr,2);
                        eminmask = repmat(eminmask,[1 length(tms)]);
                        emaxmask = repmat(emaxmask,[1 length(tms)]);
                    else
                        eminmask = eformask(fr,tms,1);
                        emaxmask = eformask(fr,tms,2);
                    end;
                    ersp = s.comp_ersp{find(origlist{nx} == clustcomps{nx}(k)),cond};
                    ersp = ersp(fr,tms);
                    ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;          
                    eallersps(:,:,end+1) = ersp; 
                end;
                condersps{cond} = eallersps;                
            end;                
                end;
        end;        
        % mask for binomial probability across components
        for cond = 1:3
            eallersps = condersps{cond};

            if ~isempty(binomalpha)
             [plotmat] = GroupSig(eallersps,.01,binomalpha);
               %for k = 1:size(eallersps,3)
                    % .01 is the individual ERSP alpha
                    % W is the binomial alpha
               %     W = factorial(size(eallersps,3))/(factorial(k)*factorial(size(eallersps,3)-k))*(.01^k)*(1-.01^(size(eallersps,3)-k) );
               %     if W< binomalpha
               %         break
               %     end;
               % end;
               % plotmat = eallersps;
               % for r = 1:size(plotmat,1)
               %     for ss = 1:size(plotmat,2)
               %         nonsig = plotmat(find(plotmat(r,ss,:)==0));
               %         if length(nonsig) > size(plotmat,3) - k
               %             plotmat(r,ss,:) = 0;
               %         end;
               %     end;        
               % end;                
            else
                plotmat = mean(eallersps,3);
            end;             
            if isempty(clim)
                colorlim = max(max( abs(plotmat,3)));
            else
                colorlim = clim;
            end;
            subplot(2,2,cond)
            if cond == 3
            imagesc(s.times(tms),s.freqs(fr), plotmat,[-colorlim colorlim]);
            else
            imagesc(s.times(tms),s.freqs(fr), plotmat,[-colorlim colorlim]);
            end;
            
            set(gca,'ydir','norm'); hold on;
            plot([0 0],[get(gca,'ylim')],'k-');
        end;        
                allersps = condersps;
    else   % single condition (no condition differences)     
        s = load ([paths{1},filename]);
        fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
        tms=find(s.times > timelims(1) & s.times<timelims(2));
        eallersps = zeros(length(fr),length(tms),0);
        fprintf('\nsubjs done: '); 
        for nx = 1:length(clustcomps)            
            if ~isempty(clustcomps{nx})                
                fprintf('\t%s',int2str(nx)); 
                s = load ([paths{nx},filename]);
                if isfield(s,'erspboot')
                    s.ersp_boot = s.erspboot;
                    s.comp_ersp = s.ersp;
                end;
                fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
                tms=find(s.times > timelims(1) & s.times<timelims(2));
                for k=1:length(clustcomps{nx})
                    eminmask = s.ersp_boot(fr,1,origlist{nx}== clustcomps{nx}(k));
                    emaxmask = s.ersp_boot(fr,2,origlist{nx}== clustcomps{nx}(k));
                    eminmask = repmat(eminmask,[1 length(tms)]);
                    emaxmask = repmat(emaxmask,[1 length(tms)]);
                    ersp = s.comp_ersp(fr,tms,origlist{nx}== clustcomps{nx}(k));
                    ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;          
                    eallersps(:,:,end+1) = ersp; 
                end;            
            end;            
        end;
        if ~isempty(binomalpha)
            [plotmat] = GroupSig(eallersps,.01,binomalpha,'binom',2);
            % mask for binomial probability across components
            %for k = 1:size(eallersps,3)
            %    W = factorial(size(eallersps,3))/(factorial(k)*factorial(size(eallersps,3)-k))*(.01^k)*(1-.01^(size(eallersps,3)-k) );
            %    if W< binomalpha
            %        break
            %    end;
            %end;
            %plotmat = eallersps;
            %for r = 1:size(plotmat,1)
            %    for ss = 1:size(plotmat,2)
            %        nonsig = plotmat(find(plotmat(r,ss,:)==0));
            %        if length(nonsig) > size(plotmat,3) - k
            %            plotmat(r,ss,:) = 0;
            %        end;
            %    end;        
            %end;
        else
            plotmat = mean(eallersps,3);
        end;    
        if isempty(clim)
            colorlim = max(max( abs(plotmat,3)));
        else
            colorlim = clim;
        end;
        quadimagesc(s.times(tms),s.freqs(fr), plotmat,[-colorlim colorlim]);
        set(gca,'ydir','norm'); hold on;
        plot([0 0],[get(gca,'ylim')],'k-','linewidth',3);  
        for v = 1:length(vertlines)
            plot([vertlines(v) vertlines(v)],[get(gca,'ylim')],'k-');
        end;
        allersps = eallersps;
    end;
    times = s.times(tms);  freqs = s.freqs(fr);
