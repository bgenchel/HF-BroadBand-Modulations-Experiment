% plots, for one subject, the pairwise scatter plots of specified IM weights
%
%
%
% INPUTS:
% savedat -- [string] file name of wts and sph files
% fullpath -- [string] full directory path where wts and sph files are saved
% ims -- [vector] of IM indices to include in the pairwise plot
% trialtype -- [0 or 1] if 1, will plot different trial types as different color
% plottype -- ['trial' or 'median'] former = scatter plot, latter= median weight
% btstrap -- [] or integer. If not [], will run shuffle permutation stats on 
%            weights interactions.

function PlotIMwtsInteractions(savedat,fullpath,ims,trialtype,plottype,btstrap);
    
    
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'}; % for plotting purposes

    s = load([fullpath,savedat,'.mat']);     
    wts = floatread([fullpath,savedat,'.wts'],[s.pcs s.numtrials],[],0);
    sph = floatread([fullpath,savedat,'.sph'],[s.numtrials s.numtrials],[],0);  
    ws = wts*sph;  winv = pinv(ws); % winv is what you need.
    
    cols = jet(length(s.dstrials));
    if iscell(ims) % paired between cells
        figure(2); pl = 1; nxt = 3;
        row = round(sqrt(length(ims{1}))); col = ceil(sqrt(length(ims{1}))); 
        for imm1 = 1:length(ims{1})
            figure(2);
            im1 = ims{1}(imm1);
            im2 = ims{2}(imm1);
            sbplot(row,col,pl);pl = pl+1;
            if ~isempty(btstrap)
                for b = 1:btstrap                        
                    perms(:,:,b) = [winv(randperm(size(winv,1)),im1),winv(randperm(size(winv,1)),im2)];
                    %ph = plot(squeeze(perms(:,1,b)),squeeze(perms(:,2,b)),'b.','markersize',8);hold on;
                    quadvals(b,1) = length(find(perms(:,1,b) < 0 & perms(:,2,b) > 0)')/(size(winv,1)/4);
                    quadvals(b,2) = length(find(perms(:,1,b) > 0 & perms(:,2,b) > 0)')/(size(winv,1)/4);
                    quadvals(b,3) = length(find(perms(:,1,b) < 0 & perms(:,2,b) < 0)')/(size(winv,1)/4); 
                    quadvals(b,4) = length(find(perms(:,1,b) > 0 & perms(:,2,b) < 0)')/(size(winv,1)/4); 
                end;
                quadprobs(1,:) = max(quadvals,[],1);
                quadprobs(2,:) = min(quadvals,[],1);
                quadwinv(1,1) = length(find(winv(:,im1) < 0 &  winv(:,im2) > 0))/(size(winv,1)/4);
                quadwinv(1,2) = length(find(winv(:,im1) > 0 &  winv(:,im2) > 0))/(size(winv,1)/4);
                quadwinv(1,3) = length(find(winv(:,im1) < 0 &  winv(:,im2) < 0))/(size(winv,1)/4);
                quadwinv(1,4) = length(find(winv(:,im1) > 0 &  winv(:,im2) < 0))/(size(winv,1)/4);            
            end;
            if trialtype == 1
                for em = 1:length(s.dstrials)
                    if strcmp(plottype,'trial')
                    ph = plot(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1),winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2),'b.','markersize',8);hold on;
                    quaddens(em,1) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) < 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) > 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
                    quaddens(em,2) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) > 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) > 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
                    quaddens(em,3) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) < 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) < 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
                    quaddens(em,4) = length(find(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1) > 0 & winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2) < 0 ))/(length(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)))/4);
                    elseif strcmp(plottype,'median')
                        ph = plot(median(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1)),median(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2)),'b.','markersize',15);hold on;
                    end;
                    set(ph,'color',cols(em,:));
                end;
                figure(nxt);
                for q = 1:4
                    sbplot(2,2,q)
                    plot([0 length(s.dstrials)+1],[quadwinv(1,q) quadwinv(1,q)],'k-');hold on; 
                    plot([0 length(s.dstrials)+1],[quadprobs(1,q) quadprobs(1,q)],'m-');hold on; 
                    plot([0 length(s.dstrials)+1],[quadprobs(2,q) quadprobs(2,q)],'m-');hold on; 
                    for e = 1:length(s.dstrials)
                        ph = bar(e,quaddens(e,q)); 
                        set(ph,'facecolor',cols(e,:));  
                        ph = text(e,0,emo2{e});
                        set(ph,'rotation',90);
                    end;
                    set(gca,'xlim',[0 length(s.dstrials)+1]);
                end;
                nxt = nxt+1;                    
            else                
                ph = plot(winv(:,im1),winv(:,im2),'b.','markersize',10);hold on;
            end;
            ph = plot([0 0],[get(gca,'ylim')],'k-');
            ph = plot([get(gca,'xlim')],[0 0],'k-');
            xl = get(gca,'xlim');
            set(gca,'xlim',xl);
            %set(gca,'xticklabel',[]);
            %set(gca,'yticklabel',[]);
            xlabel(['IM ',int2str(im1)]);
            ylabel(['IM ',int2str(im2)]);
        end;
    else      % pairwise between all   
        figure; pl = 1;
        row = length(ims)-1; col = length(ims)-1;
        for imm1 = 1:length(ims)-1
            im1 = ims(imm1); pl = pl+(imm1-1);
            for imm2 = imm1+1:length(ims)
                im2 = ims(imm2);
                sbplot(row,col,pl);pl = pl+1;
                if trialtype == 1
                    for em = 1:length(s.dstrials)
                        if strcmp(plottype,'trial')
                            ph = plot(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1),winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2),'b.','markersize',5);hold on;
                        elseif strcmp(plottype,'median')
                            ph = plot(median(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im1)),median(winv(sum(s.dstrials(1:em-1))+1:sum(s.dstrials(1:em)),im2)),'b.','markersize',15);hold on;

                        end;       set(ph,'color',cols(em,:));
                        end;
                else                
                    ph = plot(winv(:,im1),winv(:,im2),'b.','markersize',2);hold on;
                end;
                ph = plot([0 0],[get(gca,'ylim')],'r-');
                ph = plot([get(gca,'xlim')],[0 0],'r-');
                xl = get(gca,'xlim');
                set(gca,'xlim',xl);
                if imm2 == imm1+1
                    text((xl(1) - 5),0,['IM ',int2str(im1)]);
                end;
                if imm1 == 1
                    title(['IM ',int2str(im2)]);
                end;
            end;
        end;
    end;
    set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
