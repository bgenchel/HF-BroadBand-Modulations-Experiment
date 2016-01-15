% plots 'emotion space' weightings on relevant dipoles and spectral clusters
%
%
%
%
% fullwts -- [vector] of weights from the relevant spectral cluster 
%            and dimension of the 'emotion space'
% keeptrack -- [matrix] number of rows must equal the length(fullwts),
%              the 1st column is the subject index and the 2nd is the IM.
% specclust -- ['string'] full path and filename of spectral cluster to load.
% dens -- ['string'] 'neg' or 'pos' weights to plot as dipole density.
%         Input [] to plot dipoles instead of density.


function EmoDimDipoles(subjlist,fullpaths,fullwts,fullmd,keeptrack,pics,cut,specclust,subclust,clustinfo,dens);
    

    
    dipargs = {'image','mri','gui','off','dipolelength',0,'normlen','on','spheres','on','projlines','off','projimg','off','coordformat','spherical'};
    densargs = {'mrislices',[50:-15:-30],'mriview','top','geom',[4 3]};%,'cmax',4.5
    densargs = {'mrislices',[63:-5:-27],'mriview','top'};%many slices
    densargs = {'mrislices',[33,13,-4,-15,-24],'mriview','top'};
    
    dipslices = densargs{2};
    
    emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sadness','  compassion','  love','  relief','  contentedness','  awe','  happiness','  joyfulness','  excitment'}; % for plotting purposes
    
    cols = jet(15);cols(10,:) = [1 .9 0];
    
    
    clear  mnspecs minmask maxmask kptk emomeans finaltempls finalmeans finalmean   
    if ~isempty(cut) % if only highest weighted to be plotted
        cutval = cut*std(fullwts); % set cut value here****
        ifacs = find(abs(fullwts) > cutval);
        selwts = fullwts(ifacs);
        subjfacs = keeptrack(ifacs,:);       
    else
        selwts = fullwts;
        subjfacs = keeptrack;       
    end;
    eval(specclust)% call in respective cluster info
    eval(clustinfo)% call in all IM info

    
    [facvec,comods,wtsmat1,justcomps,jcwts,denslist,denswts] = Var4DipPlot(finalidx,allbigs,bigwts,orivec);

    
    figure;  
    row = 2; col = 2; place = 1; czoom = 1.3;
    clear allics allwts
    newlist = cell(1,max(subjlist));
    allics = cell(1,max(subjlist));
    nallics = cell(1,max(subjlist));
    pallics = cell(1,max(subjlist));
    allwts = cell(1,max(subjlist));
    nallwts = cell(1,max(subjlist));
    pallwts = cell(1,max(subjlist));
    counts = 0;
    for nxx = 1:length(subjlist)
        nx = subjlist(nxx);
        if length(denslist{subclust}) >= nx
            newlist{nx} = denslist{subclust}{nx};
            subidxs = subjfacs(find(subjfacs(:,1)==nx),:);
            subwts = selwts(find(subjfacs(:,1)==nx));
            ics = []; negics = [];posics = [];  pt=[];swts = [];dswts = [];pwts=[];nwts=[];
            for im = 1:size(subidxs,1)
                if ~isempty(find(ismember(abs(finalidx{subclust}(find(finalidx{subclust}(:,1)==nx),2)),subidxs(im,2))))
                    assocics = pics{nx}(find(abs(pics{nx}(:,2))==abs(subidxs(im,2))),:);
                    
                    clsics = finalidx{subclust}(find(finalidx{subclust}(:,1) == nx&abs(finalidx{subclust}(:,2))==subidxs(im,2)),3);
                    keepidx = find(ismember(assocics(:,3),clsics));
                    assocics = assocics(keepidx,:);
                    for ic = 1:size(assocics,1)
                        if assocics(ic,2) < 0
                            ss = subwts(im)*-1;
                            swts = [swts subwts(im)*-1];
                        else     
                            ss = subwts(im);
                            swts = [swts subwts(im)];
                        end;
                        if swts(end) < 0 % negative weight
                            negics = [negics assocics(ic,3)];
                            nwts = [nwts swts(end)];
                        else
                            posics = [posics assocics(ic,3)];                                    
                            pwts = [pwts swts(end)];
                        end;
                        ics = [ics assocics(ic,3)]; % undifferentiated for pos/neg wts
                                                    %swts = [swts ss];
                    end;                    
                end;
            end;
            [x y] =  unique(ics);
            allics{nx} =  ics(y);
            allwts{nx} =  swts(y);
            % if you wanna separate pos and neg wts:----
            [x y] =  unique(negics);
            nallics{nx} = negics(y);  
            nallwts{nx} =  abs(nwts(y)); % make abs and then do pos - neg
            [x y] =  unique(posics);
            pallics{nx} = posics(y);  
            pallwts{nx} =  pwts(y);
            %-----------------------
            counts = counts + length(allics{nx});
        end;   
    end;
    
    if ~isempty(dens) % if 'on'
%%-----  this is for dipole density differences from full cluster:-----------
        if counts > 0 % if there are any dipoles to plot
            clear  mnspecs minmask maxmask kptk emomeans finaltempls finalmeans finalmean templcell onebig orivec modcorr clustfacs allbigs bigwts
        %[minmask maxmask] = SigDensityDiff('sources.set','sources.set',fullpaths,fullpaths,pallics,nallics,pallwts,nallwts,300,.01,[]);

        %PlotSigDensDiff('sources.set','sources.set',fullpaths,fullpaths,pallics,nallics,pallwts,nallwts,minmask,maxmask,densargs,[]);

        
        %figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[],densargs,'yred',[],'on'); % no mask
            %figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,[],[],[],densargs,'yred',[],'on'); % no mask, un-wted
            %figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,newlist,[],densargs,'yred',[],'on');% compare to all gdcomps AND shuffle wts
            figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,allics,[],densargs,'yred',[],'on'); % compare to shuffle weights
            set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
            
            density(find(density>minmask & density<maxmask)) = 0; 
            str = ['save /home/julie/EmoCorrDens.mat density']; eval(str)
        end;
        %%--------------------
    else % otherwise plot color-coded dipoles:---------------
        if counts > 0 % only if it found dipoles
            %iclist{subclust} =  allics;
            %cwts{subclust} = allwts;         
            %sbplot(row,col,place);place = place +1;
            ph = plot([0 20],[0 0],'k-'); hold on;
            for e = 1:length(fullmd)
                ph=plot(e,fullmd(e),'.');hold on;
                set(ph,'markersize',18);set(ph,'color',cols(e,:));
                ph = text(e,fullmd(e),emo2{e});
                set(ph,'color',cols(e,:)); set(ph,'fontsize',12); 
            end; title(['Dimension Template']);
            figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[dipargs,'view',[0 0 1]],[],[],[],dipslices); 
            figure; [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[dipargs,'view',[0 0 1]],[],[],[],[]); 
            
            
% $$$             sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[dipargs,'view',[0 0 1]],[],[],[],[]);place = place +1;camzoom(czoom)
% $$$             
% $$$             sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[dipargs,'view',[1 0 0]],[],[],[],[]);place = place +1;camzoom(czoom)
% $$$             
% $$$             
% $$$             sbplot(row,col,place); [density,minmask,maxmask] = PlotDipoles('sources.set', fullpaths, allics,allwts,[],[dipargs,'view',[0 -1 0]],[],[],[],[]);place = place +1;camzoom(czoom)
            
        else
            place = place + 4;
        end;
        %%--------------------
    end;

    set(gcf,'color','w');set(gcf,'PaperOrientation','landscape');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]); 
