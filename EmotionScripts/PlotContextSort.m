% plots context vectors sorted by dimension weight in a 'TW' context decomposition
% 
% [outimages] = PlotContextSort(wtsphname,writepath,idxmat,clabels,plotfacs,imgtype,imgval,alpha)
%
% 
% filename -- name of .mat file with context decomp variables 
% datset -- name of a .set file that contains ICA weights for scalp map plotting
% fullpath -- full data directory path
% wtsphname -- name of .mat file with context decomp variables 
% writepath -- full data directory path
% clabels -- [cell array of strings] label names for each context question.
% delqs -- [numbers] indices of context questions to delete from the saved 
%          context matrix. [] if no questions to delete.
% plotfacs -- [vector] TW context decomp dimensions to plot
% imgtype -- [string] 'oldview'(hatch marks) or 'imgview'(color-coded moving
%             average) or 'wtmodes'(color-coded moving average of HISTOGRAM)
% imgval -- ['percent' or 'raw'] percent will plot 'imgview' or 'wtmodes' as
%           percents of randomly weighted context values by bootstrap method.
% alpha -- [decimal] p value desired for masking of 'imgview' or 'wtmodes' in 
%          'percent' imgval mode.
%

function [outimages,hicontigs,lowcontigs,sigout] = PlotContextSort(wtsphname,writepath,clabels,delqs,plotfacs,imgtype,imgval,alpha)

stem2 = 'TWwarpSubQsIC';
truecolor = 'off';
if ~exist('alpha')
  alpha = [];
end;
sigout = [];
tiles = [0:.1:100]; % for percentiles
shuffnum = 50; % don't need a lot cuz it pools trials as well
outimages = []; % 
movingavg = 'off';
% color-code the contexts (twoback):
% $$$     cols(1,:) = [1 0 0]; % Letter
% $$$     cols(2:9,:) = repmat([0 .7 0],[8 1]); % correct
% $$$     cols(10:11,:) = repmat([0 .8 .8],[2 1]); % bonus
% $$$     cols(12:13,:) = repmat([.9 .5 0],[2 1]); % neutral
% $$$     cols(14:15,:) = repmat([1 0 1],[2 1]); % penalty
% $$$     cols(16:21,:) = repmat([0 0 1],[6 1]); % response
cols =  hsv(length(clabels));
s = load ([writepath,wtsphname,'.mat']);
sph=floatread([writepath,wtsphname,'.sph'],[s.numrows s.numrows],[],0); 
wts=floatread([writepath,wtsphname,'.wts'],[s.numrows s.numrows],[],0); 
icamat=floatread([writepath,wtsphname,'.fdt'],[s.numrows s.numframes],[],0);    
ws = wts*sph;winv = pinv(ws);    acts = ws*icamat;  
idxmat = s.idxmat;
if ~isempty(s.pcmat) & ~isempty(s.cxtmean) % if not context or ERSP only
  erspeig = floatread([writepath,s.eigfile],[length(s.freqs)*length(s.times) s.pcmat],[],0);
  ersptempl = erspeig*winv(1:s.pcmat,:);% just makes ERSP templates
  erspdat = erspeig*winv(1:s.pcmat,:);% just makes ERSP templates
  if ~isempty(s.cxteigfile)
    idxeig = floatread([writepath,s.cxteigfile],[length(s.cxtmean) s.pcctx],[],0);
    idxtempl = idxeig*winv(s.pcmat+1:end,:);  % templates
  elseif isempty(s.cxteigfile) & ~isempty(s.anaeigfile)
    anaeig = floatread([writepath,s.anaeigfile],[length(s.anarows) inf],[],0);
    bineig = floatread([writepath,s.bineigfile],[length(s.binrows) inf],[],0);
    tmpwinv = winv(s.pcmat+1:end,:); % separate to use ana/bin rows
    anawinv = tmpwinv(s.anarows,:);
    binwinv = tmpwinv(s.binrows,:);
    anatmpl = anaeig*anawinv;
    bintmpl = bineig*binwinv;
    for a = 1:length(s.anarows)
      idxtempl(s.anarows(a),:) = anatmpl(a,:);
    end;
    for b = 1:length(s.binrows)
      idxtempl(s.binrows(b),:) = bintmpl(b,:);
    end;
  end;
end;
if isempty(s.cxteigfile) & isempty(s.anaeigfile) & s.pcmat ~= s.pctot % ERSP or context only
  ss = load([writepath,stem2,wtsphname(end),'.mat']); % hijacking for context only decomp
  dat2 = floatread([writepath,stem2,wtsphname(end),'.fdt'],[ss.numrows ss.numframes],[],0);
  if s.pcctx == s.pctot % context only
    idxeig = floatread([writepath,s.cxteigfile],[length(s.cxtmean) s.pcctx],[],0);          
    idxtempl = idxeig*winv;  % context only decomp
    erspeig = floatread([writepath,stem2,wtsphname(end),'EIGVEC.fdt'],[length(ss.freqs)*length(ss.times) ss.pcmat],[],0);
    erspdat = erspeig*dat2(1:ss.pcmat,:);% back-proj to orig data
    ersptmpl = erspeig*dat2(1:ss.pcmat,:);  % same as back proj cuz no template
  elseif s.pcmat == s.pctot % ERSP only          
    erspeig = floatread([writepath,s.eigfile],[length(s.freqs)*length(s.times) s.pcmat],[],0);
    ersptempl = erspeig*winv(1:s.pcmat,:);% just makes ERSP templates
    erspdat = erspeig*winv(1:s.pcmat,:);% just makes ERSP templates
    cxtdat = idxmat; % straight context vectors
    idxtempl = cxtdat*acts'; % weighted mean
  end;
end;

if isempty(plotfacs)
  plotfacs = [1:size(acts,1)];
end;

col = 1;   row = size(idxmat,1)+3;
hicontigs = cell(1,max(plotfacs)); % cell for each question
lowcontigs = cell(1,max(plotfacs));

for d = 1:length(plotfacs)
  dim = plotfacs(d);
  figure; pl = 1;  pg = 1;    
  
  % first plot the context vector:-----------------
  if strcmp(imgtype,'oldview')
    sbplot(row,col,pl);pl = pl+1;
  else
    sbplot(row,col,[pl pl+1]);pl = pl+2;
  end;
  ph = plot([1:size(idxtempl,1)],idxtempl(:,dim),'k-','linewidth',1);
  set(ph,'color',[.5 .5 .5]);
  set(gca,'xlim',[0 size(idxtempl,1)+1]); hold on;
  ph = plot([get(gca,'xlim')],[0 0],'r-');
  for q = 1:size(idxtempl,1)
    ph = plot(q,idxtempl(q,dim),'b.');
    set(ph,'color',cols(q,:));set(ph,'markersize',10)
  end;
  clim = max(abs(idxtempl(:,dim)))+max(abs(idxtempl(:,dim)))*.05;
  set(gca,'ylim',[-clim clim]);
  for q = 1:size(idxtempl,1)
    %ph = text(q,-clim+abs(-clim*.05),clabels{q}); 
    ph = text(q,idxtempl(q,dim)+idxtempl(q,dim)*.01,clabels{q}); 
    set(ph,'color',cols(q,:));
    set(ph,'rotation',90); set(ph,'fontsize',10);
  end;  set(gca,'xticklabel',[]);
  title(['Dim ',int2str(dim),'-Context Template ']); axis('off')
  
  % sort the weights:-----------------
  [vals idx] = sort(acts(dim,:));
  sortcon = idxmat(:,idx);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(imgtype,'oldview') |  strcmp(imgtype,'imgview')
    % variables:-------------------
    imgxwidth = round(size(sortcon,2)/200);% larger = less quantized
    imgxadv = .5;
    gausswin = 1000; gaussstd = .5;
    if strcmp(movingavg,'on')
      % plot the sorted weights-----------------
      [histout,imgoutx] = movav(acts(dim,idx),[],imgxwidth,imgxadv,[],[]);
      %[histout,imgoutx] = movav(acts(dim,idx),[],imgxwidth,imgxadv,[],[],gauss(gausswin,gaussstd));
    else
      histout = acts(dim,idx); imgoutx = [1:length(histout)];
    end;
    if strcmp(imgtype,'oldview')
      sbplot(row,col,pl);pl = pl+1;
    else
      sbplot(row,col,[pl pl+1]);pl = pl+2;
    end;
    
    %plot(acts(dim,idx),'linewidth',2);  hold on; 
    %set(gca,'xlim',[1 size(acts,2)]);
    plot(imgoutx,histout,'linewidth',2);  hold on; 
    set(gca,'xlim',[min(imgoutx) max(imgoutx)]);
    set(gca,'ylim',[min(acts(dim,:)) max(acts(dim,:))]);
    set(gca,'xticklabel',[]);set(gca,'ticklength',[0 0]);
    ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.3 .3 .3]);       
  elseif strcmp(imgtype,'wtmodes')
    % variables:-------------------
    winsize = .005;  %.015
    xadv = .0004; %.004
    mnval = min(vals); mxval = max(vals); 
    gausswin = 500; gaussstd = .5;
    % Plot the weights histogram-----------------
    [hh bins] = hist(vals,size(acts,2)); 
    if strcmp(movingavg,'on')
      %[histout,imgoutx] = movav(hh,bins,winsize,xadv,mnval,mxval,gauss(gausswin,gaussstd));
      [histout,imgoutx] = movav(hh,bins,winsize,xadv,mnval,mxval);
    else                
      histout = hh; imgoutx = bins;
    end;
    sbplot(row,col,[pl pl+1]);pl = pl+2;
    ph = bar(imgoutx,histout); hold on;
    %ph = plot(imgoutx,histout,'linewidth',2);hold on; 
    ph = plot([0 0],[get(gca,'ylim')],'r-');set(gca,'ticklength',[0 0]);
    set(gca,'xlim',[min(vals) max(vals)]);axis('off');
    %set(gca,'xlim',[min(vals) mxval]);
  end;
  
  %--------------------------------------------------------------------------------
  % Begin plotting question (ie, answers)
  if strcmp(imgtype,'oldview')
    if strcmp(imgval,'percent')
      if strcmp(movingavg,'on')
        [outdata,outx] = movav(sortcon,[],imgxwidth,imgxadv,[],[],[]);
      else
        outdata=sortcon;outx=[1:size(outdata,2)];
      end;
      
      clear bootout 
      % collect random limits to plot as percentiles:---------------
      fprintf('\nCollecting bootstrap limits for percentiles...\n');
      for b = 1:shuffnum
        randidx = ceil(rand(1,size(sortcon,2))*size(sortcon,2));% with substitution
        if strcmp(movingavg,'on')
          [bootout(:,:,b),bootx] = movav(sortcon(:,randidx),[],imgxwidth,imgxadv,[],[]);
        else
          bootout = sortcon(:,randidx); bootx = [1:size(sortcon,2)];
        end;
      end;
      poolboot = reshape(bootout,[size(bootout,1) size(bootout,2)*size(bootout,3)]);
      y=prctile(poolboot,tiles,2); %clear poolboot
    else % plot moving avg w/mean removed
      if strcmp(truecolor,'on')
        scon = sortcon;
      else
        scon = sortcon - repmat(median(sortcon,2),[1 size(sortcon,2)]);
      end;
      if strcmp(movingavg,'on')
        [outdata,outx] = movav(scon,[],imgxwidth,imgxadv,[],[],[]);
      else
        outdata=scon;outx=[1:size(outdata,2)];
      end;
      outdata = outdata/max(abs(outdata(:)));               
    end;
    for qq = 1:size(idxmat,1)
      if pl > row*col % start a new page:
        set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 8 10.5]);
        textsc([wtsphname,'; Subject ',writepath(end-4:end-1)],'title');  
        figure; pl = 1; pg = pg + 1;
        % first plot the context vector:
        if strcmp(imgtype,'oldview')
          sbplot(row,col,pl);pl = pl+1;
        else
          sbplot(row,col,[pl pl+1]);pl = pl+2;
        end;
        ph = plot([1:size(idxtempl,1)],idxtempl(:,dim),'k-','linewidth',1);
        set(ph,'color',[.5 .5 .5]);
        set(gca,'xlim',[0 size(idxtempl,1)+1]); hold on;
        ph = plot([get(gca,'xlim')],[0 0],'r-');
        for q = 1:size(idxtempl,1)
          ph = plot(q,idxtempl(q,dim),'b.');
          set(ph,'color',cols(q,:));set(ph,'markersize',10)
        end;
        clim = max(abs(idxtempl(:,dim)))+max(abs(idxtempl(:,dim)))*.05;
        set(gca,'ylim',[-clim clim]);
        for q = 1:size(idxtempl,1)
          ph = text(q,idxtempl(q,dim)+(idxtempl(q,dim)*.01),clabels{q}); 
          set(ph,'color',cols(q,:));
          set(ph,'rotation',90); set(ph,'fontsize',10);
        end;  set(gca,'xticklabel',[]);
        title(['Dim ',int2str(dim),'-Context Template ']);axis('off')
      end;
      sbplot(row,col,pl);pl = pl+1;
      if strcmp(imgval,'percent')
        z=squeeze(y(qq,:));
        % take difference to find the closest bin (tile)
        [v id] = min(abs(repmat(outdata(qq,:)',[1 size(z,2)]) - repmat(z,[size(outdata,2) 1])),[],2);
        percdata = tiles(id); % change into percentiles
                              % mask out insig percentiles
        if ~isempty(alpha)
          %percdata(q,find(percdata>100*alpha&percdata<100-(100*alpha))) = 50; 
        end;
        out{dim} = percdata;
        ph = plot(outx,percdata,'-k');hold on;
        set(ph,'color',[.7 .7 .7]);
        ph = plot([get(gca,'xlim')],[50 50],'k-');%set(ph,'color',[.3 .3 .3]);
        ints = linspace(log(.1),log(50),500);
        ints = exp(ints);
        highints = 100-ints; highints = highints(end:-1:1);
        qintervals = [ints,highints];
        %qintervals = [0:.1:100];
        avgcols = jet(length(qintervals));
        for w = 1:size(outdata,2)
          ph = plot(outx(w),percdata(1,w),'.k','markersize',4);hold on;
          tmpdat = percdata(1,w); 
          [v wint] = min(abs(qintervals - tmpdat));                       
          %wint = find(ismember(qintervals,tmpdat));
          set(ph,'color',avgcols(wint,:));
        end;                    
        set(gca,'xlim',[outx(1) outx(end)]);
        ph = text(-1,50,clabels{qq});
        %ph = text(-(.12*length(outx)),50,clabels{qq});
        set(ph,'color',cols(qq,:));hold on;
        axis('off')
        
        ttl = ['Trial view; Percent of expected random answers; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
      elseif strcmp(imgval,'raw')
        if strcmp(truecolor,'off')
          outdata = outdata - repmat(mean(outdata,2),[1 size(outdata,2)]);
        end;
        ph = plot(outx,outdata(qq,:),'-k');hold on;
        set(ph,'color',[.7 .7 .7]);
        ph = plot([get(gca,'xlim')],[0 0],'k-');%set(ph,'color',[.3 .3 .3]);
        qintervals = [-10:1:10];
        avgcols = jet(length(qintervals));
        for w = 1:size(outdata,2)
          ph = plot(outx(w),outdata(qq,w),'.k','markersize',4);hold on;
          tmpdat = outdata(qq,w); tmpdat = round(tmpdat*10);
          wint = find(ismember(qintervals,tmpdat));
          set(ph,'color',avgcols(wint,:));
        end;
        
        %ph = plot(outx,outdata(qq,:),'k-','linewidth',2); 
        %ph = plot(outx,outdata(qq,:) - mean(outdata(qq,:)),'k-','linewidth',2); 
        %set(ph,'color',[.7 .7 .7]);
        set(gca,'xlim',[outx(1) outx(end)]);
        ph = text(-(.12*size(sortcon,2)),.5,clabels{qq});set(ph,'color',cols(qq,:));hold on;
        axis('off')
        ttl = ['Trial view; Raw,smooth,mean removed;',wtsphname,'; Subject ',writepath(end-4:end-1)];
      end;
    end;% to qq loop
    ylabel('Context Questions');
    xlabel('Trials');
    
  elseif strcmp(imgtype,'imgview')
    sbplot(row,col,[pl row*col]);
    textcols =  gray(length(clabels)+15);textcols(end-14:end,:) =[];
    if strcmp(imgval,'percent')
      sigout = zeros(size(sortcon));
      % find probability for each trial/window
      wsize = 9; 
      for q = 1:size(sortcon,1)                    
        nyestot = length(find(sortcon(q,:)==1));
        qperc = nyestot/size(sortcon,2);     clear sigboot               
        randidx = ceil(rand(1,size(sortcon,2))*size(sortcon,2));% with substitution
        bb=1;clear sigboot
        sigboot = zeros(1,shuffnum*(size(sortcon,2)-wsize));
        for b=1:shuffnum
          bootimg = sortcon(q,randidx);
          for w = 1:size(bootimg,2)-wsize
            nyes = length(find(bootimg(1,w:w+(wsize-1))==1));
            sigboot(1,bb) = nyes/wsize; bb=bb+1;
          end;
        end;
        qlim = max(sigboot);
        [hh bin] = hist(sigboot,wsize);
        for w = 1:size(sortcon,2)-wsize
          nyes = length(find(sortcon(q,w:w+wsize)==1));
          if nyes/wsize > max(sigboot)
            altmat(q,w) = 100;
          elseif nyes/wsize < min(sigboot)
            altmat(q,w) = 0;
          else
            [v altmat(q,w)] = min(abs(bin-nyes/wsize));% percentile
            altmat(q,w) = altmat(q,w)*10;
          end;
        end;
        % find sig stretches of trials 
        hiprob = find(altmat(q,:) == 100);contig = [];
        tmpcontigs=cell(1,0); g=1;
        for h = 2:length(hiprob)
          if hiprob(h) == hiprob(h-1)+1
            contig = [contig idx(hiprob(h))]; % save trial #
            sigout(q,hiprob(h)) = 1;
          elseif length(contig) > (nyestot/4)/wsize
            tmpcontigs{g} = contig; contig = [];g=g+1;
          end;
        end;  
        if length(contig) > (nyestot/4)/wsize
          tmpcontigs{g} = contig; 
          sigout(q,hiprob(h)) = 1;
        end;
        hicontigs{dim}{q} = tmpcontigs;
        hiprob = find(altmat(q,:) == 0);contig = [];
        tmpcontigs=cell(1,0); g=1;
        for h = 2:length(hiprob)
          if hiprob(h) == hiprob(h-1)+1
            contig = [contig idx(hiprob(h))]; % save trial #
            sigout(q,hiprob(h)) = -1;
          elseif length(contig) > (nyestot/4)/wsize
            tmpcontigs{g} = contig; contig = [];g=g+1;
          end;
        end;  
        if length(contig) > (nyestot/4)/wsize
          tmpcontigs{g} = contig;
          sigout(q,hiprob(h)) = -1;
        end;
        lowcontigs{dim}{q} = tmpcontigs;
      end;
      percimg = altmat; 
      %                 clear bootout 
      %                 % collect random limits to plot as percentiles:---------------
      %                 fprintf('\nCollecting bootstrap limits for percentiles...\n');
      %                 shuffnum = 100;
      %                 for b = 1:shuffnum
      %                     randidx = ceil(rand(1,size(sortcon,2))*size(sortcon,2));% with substitution
      %                     %randidx = randperm(size(sortcon,2));% without substitution
      %                     if strcmp(movingavg,'on')
      %                         [bootout(:,:,b),bootx] = movav(sortcon(:,randidx),[],imgxwidth,imgxadv,[],[]);
      %                     else
      %                         bootout = sortcon(:,randidx); bootx = [1:size(sortcon,2)];
      %                     end;
      %                 end;
      %                 poolboot = reshape(bootout,[size(bootout,1) size(bootout,2)*size(bootout,3)]);
      %                 y=prctile(poolboot,tiles,2);       
      %                 clear percimg
      %                 for q = 1:size(imgout,1)
      %                     z=squeeze(y(q,:));
      %                     % if repeated percentiles:----
      %                     if length(find(z == z(1))) > 1
      %                         [vv ii] = find(z == z(1));
      %                         z(ii(2:end)) = 100; % really big
      %                     end;
      %                     if length(find(z == z(end))) > 1
      %                         [vv ii] = find(z == z(end));
      %                         z(ii(1:end-1)) = 100; % really big
      %                     end;
      %                     
      %                     [v id] = min(abs(repmat(imgout(q,:)',[1 size(z,2)]) - repmat(z,[size(imgout,2) 1])),[],2);                    
      %                     percimg(:,q) = tiles(id); % change into percentiles
      %                 end;
      %                 percimg = percimg';
      % mask out insig percentiles
      %                 if ~isempty(alpha)
      %                     for q = 1:size(percimg,1)
      %                         percimg(q,find(percimg(q,:)>100*alpha&percimg(q,:)<100-(100*alpha))) = 50; %mask 
      %                     end;
      %                 end;
      outimages{dim} = percimg;
      if strcmp(movingavg,'on')
        [imgout,imgoutx] = movav(percimg,[],imgxwidth,imgxadv,[],[]);
      else
        imgout=percimg;imgoutx=[1:size(percimg,2)];
      end;
      %imagesc([100/size(imgout,2):100/size(imgout,2):100],[1:size(imgout,1)],imgout); hold on; % for group plotting
      imagesc(imgoutx,[1:size(imgout,1)],imgout); hold on;     
      ttl = ['Trial view; Percents; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
    else % raw, not percents
      ttl = ['Trial view; Raw data; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
      if strcmp(truecolor,'off')
        sortcon = sortcon - repmat(mean(sortcon,2),[1 size(sortcon,2)]);
      ttl = ['Trial view; Raw data; mean removed; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
      end;
      if strcmp(movingavg,'on')
        [imgout,imgoutx] = movav(sortcon,[],imgxwidth,imgxadv,[],[]);
      else
        imgout=sortcon;  imgoutx=[1:size(sortcon,2)];
      end;
      outimages{dim} = imgout;
      imagesc(imgoutx,[1:size(imgout,1)],imgout); hold on;
    end            
    for q = 1:size(idxtempl,1)
      %ph = text(imgoutx(end)+imgoutx(end)*.01,q,clabels{q}); 
      ph = text(imgoutx(round(length(imgoutx)/2)),q,clabels{q}); set(ph,'color','w');
      %set(ph,'color',textcols(q,:));
      set(ph,'fontsize',10);
    end;              
    set(gca,'xlim',[min(imgoutx) max(imgoutx)]);
    set(gca,'yticklabel',[]);
    ylabel('Context Questions');
    xlabel('Trials');
    
  elseif strcmp(imgtype,'wtmodes')%----------------------------               
    newvals = [vals(1):.001:vals(end)];
    clear imgout imgoutx avgimgout
    for qs = 1:size(idxmat,1)                
      for c = 1:length(newvals)-1
        if ~isempty(idxmat(qs,find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))))
          imgout(qs,c) = (mean(idxmat(qs,find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))))*length(find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))));
        else
          imgout(qs,c) = 0;
        end;
        imgoutx(1,c) = mean(newvals(c:c+1));
      end;
      imgout(qs,c+1) =(mean(idxmat(qs,find(acts(dim,:)>=newvals(c))))*length(find(acts(dim,:)>=newvals(c))));
      imgoutx(1,c+1) = newvals(c+1);
      if strcmp(movingavg,'on')
        imgxwidth = round(size(imgout,2)/200);% larger = less quantized
        [avgimgout(qs,:),avgimgoutx] = movav(imgout(qs,:),[],imgxwidth,1,1,size(imgout,2));
        if qs == size(idxmat,1)
          imgout = avgimgout; 
          xdiff = length(imgoutx) - length(avgimgoutx);
          imgoutx(1:xdiff/2) = []; imgoutx(end-(xdiff/2-1):end) = [];                        
        end;
      end;
    end;
    textcols =  gray(length(clabels)+15);textcols(end-14:end,:) =[];
    sbplot(row,col,[pl row*col]);
    if strcmp(imgval,'percent')
      clear bootout 
      % collect random limits to plot as percentiles:---------------
      fprintf('\nCollecting bootstrap limits for percentiles...\n');
      for b = 1:shuffnum
        randidx = ceil(rand(1,size(sortcon,2))*size(sortcon,2));
        clear imgout imgoutx
        for qs = 1:size(idxmat,1)
          for c = 1:length(newvals)-1
            if ~isempty(idxmat(qs,find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))))
              bootout(qs,c) = (mean(idxmat(qs,find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))))*length(find(acts(dim,:)>=newvals(c)&acts(dim,:)<newvals(c+1))));
            else
              bootout(qs,c) = 0;
            end;
            bootx(1,c) = mean(newvals(c:c+1));
          end;
          imgout(qs,c+1) =(mean(idxmat(qs,find(acts(dim,:)>=newvals(c))))*length(find(acts(dim,:)>=newvals(c))));
          imgoutx(1,c+1) = newvals(c+1);
          if strcmp(movingavg,'on')
            imgxwidth = round(size(bootout,2)/300);% larger = less quantized
            [avgbootout(qs,:),avgbootoutx] = movav(imgout(qs,:),[],imgxwidth,1,1,size(imgout,2));
            if qs == size(idxmat,1)
              bootout = avgbootout; 
              xdiff = length(bootx) - length(avgbootoutx);
              bootx(1:xdiff/2) = []; bootx(end-(xdiff/2-1):end) = [];
            end;
          end;
        end;
        %                     if strcmp(movingavg,'on')
        %                         [bootout(:,:,b),bootx] = movav(sortcon(:,randidx),vals,winsize,xadv,mnval,mxval);
        %                     else
        %                         bootout=sortcon(:,randidx);bootx=vals;
        %                     end;
      end;
      poolboot = reshape(bootout,[size(bootout,1) size(bootout,2)*size(bootout,3)]);
      y=prctile(poolboot,tiles,2);       
      clear percimg
      for q = 1:size(imgout,1)
        z=squeeze(y(q,:));
        % if repeated percentiles:----
        if length(find(z == z(1))) > 1
          [vv ii] = find(z == z(1));
          z(ii(2:end)) = 100; % really big
        end;
        if length(find(z == z(end))) > 1
          [vv ii] = find(z == z(end));
          z(ii(1:end-1)) = 100; % really big
        end;
        
        [v id] = min(abs(repmat(imgout(q,:)',[1 size(z,2)]) - repmat(z,[size(imgout,2) 1])),[],2);
        %z=squeeze(y(q,:,:));
        %[v id] = min(abs(repmat(imgout(q,:)',[1 size(z,2)]) - z),[],2);
        
        percimg(:,q) = tiles(id); % change into percentiles
      end;
      percimg = percimg';
      % mask out insig percentiles
      if ~isempty(alpha)
        for q = 1:size(percimg,1)
          percimg(q,find(percimg(q,:)>100*alpha&percimg(q,:)<100-(100*alpha))) = 50; %mask 
        end;
      end;
      imagesc(imgoutx,[1:size(imgout,1)],percimg); hold on;            
      outimages{dim} = percimg;
      ttl = ['Histogram view; Percents; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
    else                
      outimages{dim} = imgout;
      imagesc(imgoutx,[1:size(imgout,1)],imgout); hold on;            
      ttl = ['Histogram view; Raw data; ',wtsphname,'; Subject ',writepath(end-4:end-1)];
    end;
    set(gca,'yticklabel',[]);
    for q = 1:size(idxtempl,1)
      %ph = text(imgoutx(end)+imgoutx(end)*.01,q,clabels{q}); 
      %set(ph,'color',textcols(q,:));
      ph = text(0,q,clabels{q}); set(ph,'color','k');
      set(ph,'fontsize',10);
    end;            
    set(gca,'xlim',[min(imgoutx) max(imgoutx)]);
    xlabel('Weights');
    ylabel('Context Questions');
  end;
  set(gcf,'PaperOrientation','portrait');  set(gcf,'PaperPosition',[0.25 0.25 8 10.5]);
  textsc(ttl,'title');
end;

