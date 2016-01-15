% finds percent trials on upper and lower edges of all-emotion weight distribution
% ploton -- [0 | 1] 1 to plot results, one page/subj. 0 otherwise

function [allposperc,allnegperc] = PercDistShift(subjlist,savedat,fullpaths,emos,ploton);
    
    
    cols = jet(length(emos));
    allposperc = cell(1,subjlist(end));
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        s = load([fullpaths{nx},savedat,'.mat']);         
        sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
        wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
        ws = wts*sph;    winv = pinv(ws); clear wts sph ws 
        if ploton == 1
            row = round(sqrt(size(winv,2)));col= round(sqrt(size(winv,2)));
            figure;
        end;        
        for fac = 1:size(winv,2)
            fulldist = winv(:,fac);
            fullpos = length(find(winv(:,fac) > 0))/sum(s.dstrials);
            fullneg = length(find(winv(:,fac) < 0))/sum(s.dstrials);
            if ploton == 1
                sbplot(row,col,fac)
            end;            
            for em = 1:length(emos)
                %emdist = winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),fac);
                emopos(em,fac) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),fac) > 0));  
                emoneg(em,fac) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),fac) < 0));
                %emopos(em,fac) = emopos(em,fac)/s.dstrials(1,em);
                %emoneg(em,fac) = emoneg(em,fac)/s.dstrials(1,em);
                emopos(em,fac) = emopos(em,fac)/s.dstrials(1,em) - fullpos;
                emoneg(em,fac) = emoneg(em,fac)/s.dstrials(1,em) - fullneg;
                if ploton == 1
                    ph=plot(em,emopos(em,fac),'k.'); hold on; 
                    set(ph,'markersize',15);
                    set(ph,'color',cols(em,:));
                    ph=plot(em,emoneg(em,fac),'k.'); hold on; 
                    set(ph,'markersize',15);
                    set(ph,'color',cols(em,:));
                    yl = get(gca,'ylim');
                    ph = text(em,yl(1),emos{em});hold on;
                    set(ph,'rotation',90);set(ph,'color',cols(em,:));
                end;
                set(gca,'xlim',[1 15]);
            end;
            if ploton == 1
                plot(emopos(:,fac),'k-'); hold on; 
                plot(emoneg(:,fac),'k-'); hold on; 
                %plot(emoneg(:,fac),'c-');
                %set(gca,'ylim',[-.5 .5]);
                plot([get(gca,'xlim')],[0 0],'g-');
                %plot([get(gca,'xlim')],[fullpos fullpos],'r-'); hold on; 
                %plot([get(gca,'xlim')],[fullneg fullneg],'b-');
            end;
        end;
        allposperc{nx} = emopos;  
        allnegperc{nx} = emoneg;  
        fprintf('.');
    end;
