% plots results from EmoSpacebyICA()
%
%


function PlotEmoSpacebyICA(savedat,fullpaths,subjlist,whichfacs,emomeans);


    q = load('/data/common4/emotion/clustfacs.mat'); 
    faccols = {'r','b','g','m','c','k','y'};
    cols = jet(15);cols(10,:) = [.9 .9 0];
    dimcols = hsv(size(alldat,3));
    for nx = 1:length(whichfacs)
        if ~isempty(whichfacs{nx})
            EEG = pop_loadset('sources1.set', fullpaths{nx});
            s = load([fullpaths{nx},savedat,'.mat']);  
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            icamatall = floatread([fullpaths{nx},savedat,'.fdt'],[s.numtrials s.numframes],[],0);    
            ws = wts*sph;    activations = ws*icamatall;   clear wts sph ws
            pg=0; close; close; close; close;  close;  close;  close; 
            for dm= 1:length(whichfacs{nx})
                figure; pl = 1; 
                if length(whichfacs{nx}{dm}) == 1
                    row = 2; col=2; zoom = 1;
                elseif length(whichfacs{nx}{dm}) == 2
                    row = 3;  col = 2; zoom = 1.3;
                elseif length(whichfacs{nx}{dm}) == 3
                    row = 2;  col = 4;zoom = .75;
                elseif length(whichfacs{nx}{dm}) == 4|length(whichfacs{nx}{dm}) == 5
                    row = 3; col = 4;zoom = .85;
                elseif length(whichfacs{nx}{dm}) == 6|length(whichfacs{nx}{dm}) == 7
                    row = 4; col = 4;zoom = .9;
                end;                
                maxsz = 25;
                for fac = 1:length(whichfacs{nx}{dm})
                    sbplot(row,col,pl);            
                    cmpcols = lines(length(gdcomps{nx}));
                    tp = whichfacs{nx}{dm}(fac);               
                    %sbplot(row,col,pl)
                    for cp = 1:length(gdcomps{nx})
                        rcp =cp;
                        if tp < 0
                            ph = plot(s.freqs,activations(abs(tp),length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp)*-1,'linewidth',1.5); 
                        else
                            ph = plot(s.freqs,activations(tp,length(s.freqs)*(rcp-1)+1:length(s.freqs)*rcp),'linewidth',1.5); 
                        end;
                        
                        hold on; %set(ph,'color',cmpcols(cp,:));                    
                        set(ph,'color',faccols{fac});
                    end;pl = pl+1;
                    set(gca,'ytick',[-10:5:15]);set(gca,'yticklabel',[-10:5:15]);
                    title(['Fac ',int2str(whichfacs{nx}{dm}(fac))]);
                    set(gca,'xlim',[s.freqs(1) s.freqs(end)]); 
                    set(gca,'ylim',[-10 10]);
                    set(gca,'box','off');
                    set(gca,'xgrid','on');
                    set(gca,'ticklength',[.05 .05]);
                    if fac == length(whichfacs{nx}{dm})
                        xlabel('Frequency (Hz)');
                        ylabel('Relative Power');
                    end; 
                    facfac = maxsz - 5;
                    if ~isempty(q.allbigs{nx})
                        sbplot(row,col,pl);  
                        mydipplot(EEG.dipfit.model(cell2mat(q.allbigs{nx}(abs(whichfacs{nx}{dm}(fac))))),'image','mri','gui','off','dipolelength',0,'dipolesize',facfac,'normlen','on','spheres','on','color',{faccols{fac}},'projlines','on','projimg','on','coordformat','spherical');                hold on;
                        view(60,30); pl = pl+1;
                    end;
                    camzoom(zoom); 
                end;
                set(gcf,'color','w');
                if col < 3
                    pl = pl+1;
                end;                    
                sbplot(row,col,pl);
                ph = plot([0 length(emos)+1],[0 0],'k-');set(ph,'color',[.5 .5 .5]);hold on;
                for e = 1:length(emos)
                    dim = nxgdims{nx}(dm);
                    ph = plot(e,emomeans(e,dim,nx),'k.','markersize',15);
                    set(ph,'color',dimcols(dim,:)); hold on;
                end;
                yl = get(gca,'ylim');
                set(gca,'ylim',[yl]);
                for e = 1:length(emos)
                    ph = text(e,yl(1)+.01,emos{e});
                    set(ph,'rotation',90);
                    set(ph,'color',cols(e,:));             
                    set(gca,'xlim',[0 length(emos)+1]);
                end;
                pl = pl+1;title(['Subj ',int2str(nx),' Valence Order']);
                textsc(['Subject ',fullpaths{nx}(end-4:end-1)],'title'); 
                pg = pg+1;
                str = ['print ',fullpaths{nx},'EmoDimFacs',int2str(pg),'.jpg -djpeg']; eval(str);
            end;
        end;
    end;
    
