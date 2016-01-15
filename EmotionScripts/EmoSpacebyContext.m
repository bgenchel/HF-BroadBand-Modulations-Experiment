% use context ICA to find emotion space for all subjects
%
%
%
%
%
%
%
%

function EmoSpacebyContext(savedat,fullpaths,subjlist);
    
    keyboard
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
    cols = jet(15);cols(10,:) = [.9 .9 0];
    allsubjmat = zeros(15,0);
    for nxs = 1:length(subjlist)
        nx = subjlist(nxs);
        s = load([fullpaths{nx},savedat,'.mat']);  
        sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
        ws = wts*sph;     icawinv = pinv(ws); icawinv = icawinv';
        % make the addmat matrix
        addmat = ones(15,sum(s.dstrials)); addmat = addmat*-1;
        for e = 1:15
            addmat(e,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))) = 1;
        end;     
        %%%%%%%%%%
        icawinv(end+1:end+size(addmat,1),:) = addmat;
        
        pcdims = 29;
        [weights,sphere,compvars,bias,signs,lrates,activations] = runica(icawinv,'extended',1,'stop',1e-7,'verbose','on','maxsteps',1000,'pca',pcdims);
        [weights,sphere,compvars,bias,signs,lrates,activations] = runica(icawinv,'extended',1,'stop',1e-7,'verbose','on','maxsteps',1000);
        
        ws = weights*sphere; newwinv = pinv(ws); clear ws weights sphere 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot back-projections of each emotion dimension
        clear alldat
        for dim = 1:size(activations,1)
            zerout = [1:size(activations,1)];
            zerout(dim) = [];   acts = activations;
            acts(zerout,:) = 0; % take out all dims but one
            
            newdat = newwinv * acts;
            
            %alldat(:,:,dim) = newdat(16:30,:);
            alldat(:,:,dim) = newdat;
        end;
        figure;  dimcols = hsv(size(alldat,3));
        for e = 1:size(alldat,1)
            for dim = 12:12%size(alldat,3)
                %ph = plot(sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),alldat(e,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e)),dim),'k-');
                %ph = plot(alldat(e,:,dim),'k-');
                ph = plot(icawinv(e,:),'k-');
                set(ph,'color',dimcols(dim,:)); hold on;
            end;
            ph = text(sum(s.dstrials(1:e-1))+1,1.75,emo2{e});
            set(ph,'color',cols(e,:)); 
            
        end;        
        set(gca,'xlim',[0 size(newdat,2)]);
        set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
          
                
        
        
        
        winv = newwinv(16:30,:);
        
        cols = jet(15);cols(10,:) = [.9 .9 0];
        figure;  % for one point per emotion (mean)
        for wv = 1:size(winv,2)
            subplot(round(sqrt(size(winv,2))),ceil(sqrt(size(winv,2))),wv)
            for e = 1:size(winv,1)
                ph=plot(e,winv(e,wv),'.');hold on;
                set(ph,'markersize',20);set(ph,'color',cols(e,:));
                ph = text(e,winv(e,wv),emo2{e});
                set(ph,'color',cols(e,:)); 
            end;
        end;axcopy
        textsc(fullpaths{nx},'title');
        
        winv2 = winv;
        clear realwts
        for dim = 1:size(activations,1)
            for e = 1:length(emo2)  % start with 2 for2 straight nums (not diffs)
                realwts(dim,e) = mean(activations(dim,sum(s.dstrials(1:e-1))+1:sum(s.dstrials(1:e))));
            end;
        end;
        winv=realwts';

        fprintf('\n One More SUBJECT Done: %i',nx);
        clear sph winv wts ws activations icamatall
    end;


