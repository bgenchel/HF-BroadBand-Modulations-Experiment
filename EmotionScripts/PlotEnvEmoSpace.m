% plots the low-pass back-projections of emo dimensions with background of full origwinv
%
%
%
%
%

function PlotEnvEmoSpace(savedat,fullpaths,origwinvs,spaceacts,savewinvs,subjlist);
    
    for nxx = 1:length(subjlist)
        nx = subjlist(nx);
        s = load([fullpaths{nx},savedat,'.mat']);  
        alldat = zeros(size(origwinvs{nx},1),size(spaceacts{nx},2),size(spaceacts{nx},1));
        for dim = 1:size(spaceacts{nx},1)
            zerout = [1:size(spaceacts{nx},1)];
            zerout(dim) = [];   acts = spaceacts{nx};
            acts(zerout,:) = 0; % take out all dims but one            
            alldat(:,:,dim) = savewinvs(:,:,nx) * acts; % orig Fac X windows X back-proj (new)dim
            
            newdata = zeros(size(origwinvs{nx},1),0,size(spaceacts{nx},1));
            for em = 1:length(s.dstrials)
                xwidth = round(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/40);
                [outdata,outx] = movav(alldat(:,sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),dim),0,xwidth,0);
                newdata(:,end+1:end+size(outdata,2),dim) = outdata;
                intvls(em,:) = [size(outdata,2) length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em))) - size(outdata,2)];
            end;
            
            envdata{dim} = env(newdata(:,:,dim)); 
            
        end;

        
        origsmooth = zeros(size(origwinvs{nx},1),0);
        for em = 1:length(s.dstrials)
            
            xwidth = round(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/40);
            [outdata,outx] = movav(origwinvs{nx}(:,sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em))),0,xwidth,0);
            origsmooth(:,end+1:end+size(outdata,2)) = outdata;
            
        end;
        totalenv = env(origsmooth); 
        
      
        figure; cols = jet(size(origwinvs{nx},1));
        row = round(sqrt(size(spaceacts{nx},1))); col = round(sqrt(size(spaceacts{nx},1)));
        f1 = [[1:size(totalenv,2)],[size(totalenv,2):-1:1]];
        f2 = [totalenv(1,:),[totalenv(2,size(totalenv,2):-1:1)]];
        
        for dim = 1:length(envdata)
            sbplot(row,col,dim)
            ph = fill(f1,f2,[.92 .92 .92]);hold on;
            set(ph,'edgecolor',[.92 .92 .92]);
            for em = 1:size(intvls,1)
                datlen = intvls(em,1);
                datdiff = intvls(em,2);
                dattot = intvls(em,3);
                
                
                %ph = plot([sum(intvls(1:(em-1),3))+round(datdiff/2)+1:sum(intvls(1:(em-1),3))+datlen+round(datdiff/2)],envdata{dim}(1,sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen),'k-','linewidth',1);
                ph = plot([sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen],envdata{dim}(1,sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen),'k-','linewidth',1);
                set(ph,'color',cols(em,:));
                %ph = plot([sum(intvls(1:(em-1),3))+round(datdiff/2)+1:sum(intvls(1:(em-1),3))+datlen+round(datdiff/2)],envdata{dim}(2,sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen),'k-','linewidth',1);
                ph = plot([sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen],envdata{dim}(2,sum(intvls(1:(em-1),1))+1:sum(intvls(1:(em-1),1))+datlen),'k-','linewidth',1);
                set(ph,'color',cols(em,:));
                
            end;
            
            set(gca,'xlim',[0 size(origwinvs{nx},2)]);
            set(gca,'ylim',[-1.5 1.5]);
            title(['Dim ',int2str(dim)]);
        end;
        axcopy
        textsc(['Subject ',int2str(nx),':  Envelopes of each dimension (smoothed by #pnts/40 (~10)'],'title');
     set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
       
