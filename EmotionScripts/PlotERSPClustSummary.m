% plots newtimef ERSP results with a dipole fig of clustered dipoles
%
% [pl,allbase,freqs] = PlotERSPClustSummary(filename,dipset,paths,origlist,clustcomps,freqlims,timelims,clim,binomalpha,row,col,dipcolor,startpl);
% 
% INPUTS:
% filename -- [string] name of data file with ERSP information
% dipset -- [string] name of dataset to load with dipole information
% paths -- [cell array of strings] full data directory path for each subject
% clustcomps -- [cell array] vector of IC indices for each subject to plot
% freqlims -- [minfreq maxfreq] in Hz
% timelims -- [mintime maxtime] in ms
% clim -- [integer] color limit (dB)
% binomalpha -- [decimal] level of significance for binomial probability masking (p < .00001)
%               if [], will plot ERSPs with no masking for significance (individual means
%               as well as the mean across subjects.
% plotall -- [0 or 1] if 1, will plot all individual ERSPs for all ICs in the cluster
% have to specify row and columns of figure that YOU open before running function.
% plotall -- [0 | 1] input 1 to plot individual subject ERSPs
%            in addition to cluter means
% plotitc -- ['itc' or 'noitc'] if 'itc', when will plot additional panels with masked 
%            ITC values. 'noitc' will plot ERSP only. Default: 'noitc'

function [pl,allbase,freqs] = PlotERSPClustSummary(filename,dipset,paths,clustcomps,freqlims,timelims,clim,binomalpha,row,col,startpl,plotall,plotitc);
     
ilim = []; % itc colorlim
if ~exist('plotitc')
  plotitc = 'noitc';
elseif isempty(plotitc)
  plotitc = 'noitc';
end;
if ~exist('plotall')
  plotall = 0;
elseif isempty(plotall)
  plotall = 0;
end;
fg = gcf;
pl = startpl;
if length(clim) > 1
  fprintf('Colorlim must be a single value (will scale around zero automatically)');
end;
if iscell(filename)
  s = load ([paths{1},filename{1}]);
else
  s = load ([paths{1},filename]);
end;
fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
tms=find(s.times > timelims(1) & s.times<timelims(2));
%%%%%%  Check for log frequency spacing  %%%%%%%%%%%%%%%%
%%%  need to fix this
if isfield(s,'freqscale') % then this is easy
  freqscale = s.freqscale;
else
  tst1 = sqrt(s.freqs);
  tst2 = log(s.freqs);
  lfr1 = linspace(sqrt(s.freqs(1)),sqrt(s.freqs(end)),length(s.freqs));
  lfr2 = linspace(log(s.freqs(1)),log(s.freqs(end)),length(s.freqs));
  
  if tst1 == lfr1
    freqscale = 'quad'; % quadratic spacing
  elseif tst2 == lfr2 
    freqscale = 'log'; % log spacing
  else
    freqscale = 'linear'; % linear spacing
  end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if iscell(s.comp_ersp) % newtimef of two expt'l conditions
  allbase = zeros(0,length(fr));
  fprintf('\nsubjs done: '); 
  %cord = [2,1,3];
  for cond = 1:3
    %cond=cord(cc);
    eallersps = zeros(length(fr),length(tms),0);
    iallersps = zeros(length(fr),length(tms),0);
    for nx = 1:length(paths)
      if ~isempty(clustcomps{nx})
        s = load ([paths{nx},filename]);
        fprintf('..%s..',int2str(nx));  
        %r = load ([paths{nx},'taskbase.mat']);
        fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
        tms=find(s.times > timelims(1) & s.times<timelims(2));
        gdcomps{nx} = s.complist;
        for k=1:length(clustcomps{nx})
          eformask = s.ersp_boot{find(s.complist == clustcomps{nx}(k)),cond};% fr X 2
          if cond ~= 3
            eminmask = eformask(fr,1);
            emaxmask = eformask(fr,2);
            eminmask = repmat(eminmask,[1 length(tms)]);
            emaxmask = repmat(emaxmask,[1 length(tms)]);
          else
            eminmask = eformask(fr,tms,1);
            emaxmask = eformask(fr,tms,2);
          end;
          ersp = s.comp_ersp{find(s.complist == clustcomps{nx}(k)),cond};
          ersp = ersp(fr,tms);
          if ~isempty(binomalpha)
            ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;   
          end;
          eallersps(:,:,end+1) = ersp; 
          allbase(end+1,:) = s.baseline(find(s.complist == clustcomps{nx}(k)),fr);
          if strcmp(plotitc,'itc') % do plot ITC as well
            iformask = s.itc_boot{find(s.complist == clustcomps{nx}(k)),cond};% fr X 2
            if cond ~= 3
              imaxmask = iformask(fr,1);
              imaxmask = repmat(imaxmask,[1 length(tms)]);
            else
              iminmask = iformask(fr,tms,1);
              imaxmask = iformask(fr,tms,2);
            end;
            itc = s.comp_itc{find(s.complist == clustcomps{nx}(k)),cond};
            itc = itc(fr,tms);
            if ~isempty(binomalpha)
              if cond ~= 3
                itc(find(itc <= imaxmask)) = 0; 
              else % for difference, mask is min and mask
                itc(find(iminmask <= itc& itc <= imaxmask)) = 0;   
              end;
            end;
            iallersps(:,:,end+1) = itc; 
          end;
        end;     % for k           
      end;% for ~isempty
    end;  % for nx
    conditcs{cond} = iallersps;                
    condersps{cond} = eallersps;                
  end;% for cond
  for cond = 1:3
    if ~isempty(binomalpha)
      [plotmat] = GroupSig(condersps{cond},s.alpha,binomalpha,'binom');    
    else
      plotmat = condersps{cond};
    end;           
    if isempty(clim)
      clim = max(max( abs(mean(plotmat,3))));
    end;
    if plotall == 1
      figure; rw = round(sqrt(size(condersps{cond},3)));
      cl = ceil(sqrt(size(condersps{cond},3)));
      for m = 1:size(condersps{cond},3)
        sbplot(rw,cl,m);
        if strcmp(freqscale,'quad')
          quadimagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); 
          hold on;
        elseif strcmp(freqscale,'log')
          mylogimagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); 
          hold on;
        else
          imagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); hold on;
        end;            
        plot([0 0],[get(gca,'ylim')],'k-');pl = pl+1;
      end; 
    end;
    figure(fg);     
    sbplot(row,col,pl);pl = pl+1;
    if strcmp(freqscale,'quad')
      quadimagesc(s.times(tms),s.freqs(fr), mean(plotmat,3),[-clim clim]); 
      hold on;
    elseif strcmp(freqscale,'log')
      mylogimagesc(s.times(tms),s.freqs(fr),mean(plotmat,3),[-clim clim]); 
      hold on;
    else
      imagesc(s.times(tms),s.freqs(fr),mean(plotmat,3),[-clim clim]); hold on;
    end;
   set(gca,'fontsize',12);
    %set(gca,'ydir','norm'); hold on;
    plot([0 0],[get(gca,'ylim')],'k-');
    if cond == 3
      title('Difference');
    end;      
  end;   
  if strcmp(plotitc,'itc')
    for cond = 1:3
      if ~isempty(binomalpha)
        [itcmat] = GroupSig(conditcs{cond},s.alpha,binomalpha,'binom');  
      else
        itcmat = condersps{cond};
      end;           
      if isempty(ilim)
        ilim = max(max( abs(mean(itcmat,3))));
      end;
      if plotall == 1
        figure; rw = round(sqrt(size(conditcs{cond},3)));
        cl = ceil(sqrt(size(conditcs{cond},3)));
        for m = 1:size(conditcs{cond},3)
          sbplot(rw,cl,m);
          if strcmp(freqscale,'quad')
            quadimagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); 
            hold on;
          elseif strcmp(freqscale,'log')
            mylogimagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); 
            hold on;
          else
            imagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); hold on;
          end;            
          plot([0 0],[get(gca,'ylim')],'k-');pl = pl+1;
        end; 
      end;
      figure(fg);
      sbplot(row,col,pl);pl = pl+1;
      if strcmp(freqscale,'quad')
        quadimagesc(s.times(tms),s.freqs(fr), mean(itcmat,3),[-ilim ilim]); 
        hold on;
      elseif strcmp(freqscale,'log')
        mylogimagesc(s.times(tms),s.freqs(fr),mean(itcmat,3),[-ilim ilim]); 
        hold on;
      else
        imagesc(s.times(tms),s.freqs(fr),mean(itcmat,3),[-ilim ilim]); hold on;
      end;
      set(gca,'fontsize',12);
      %set(gca,'ydir','norm'); hold on;
      plot([0 0],[get(gca,'ylim')],'k-');
      if cond == 3
        title('Difference');
      end; 
    end;
  end;   
else% single condition newtimef 
  
  allbase = zeros(0,length(fr));
  plotallttls = [];ploticttls = [];
  if iscell(filename) % plot several conditions from diff files
    for fn = 1:length(filename)
      eallersps = zeros(length(fr),length(tms),0);
      iallersps = zeros(length(fr),length(tms),0);
      for nx = 1:length(paths)
        if ~isempty(clustcomps{nx})
          s = load ([paths{nx},filename{fn}]);
          %r = load ([paths{nx},'taskbase.mat']);
          fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
          tms=find(s.times > timelims(1) & s.times<timelims(2));
          gdcomps{nx} = s.complist;
          for k=1:length(clustcomps{nx})
            eformask = s.ersp_boot(fr,:,find(s.complist == clustcomps{nx}(k)));% fr X 2
            eminmask = eformask(:,1);
            emaxmask = eformask(:,2);
            eminmask = repmat(eminmask,[1 length(tms)]);
            emaxmask = repmat(emaxmask,[1 length(tms)]);
            ersp = s.comp_ersp(fr,tms,find(s.complist == clustcomps{nx}(k)));
            if ~isempty(binomalpha) % don't mask indiv ersps
              ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;          
            end;
            eallersps(:,:,end+1) = ersp; 
            if strcmp(plotitc,'itc') % plot ITC too
              iformask = s.itc_boot(fr,:,find(s.complist == clustcomps{nx}(k)));% fr X 2
              imaxmask = iformask(:,1);
              imaxmask = repmat(imaxmask,[1 length(tms)]);
              itc = s.comp_ersp(fr,tms,find(s.complist == clustcomps{nx}(k)));
              if ~isempty(binomalpha) % don't mask indiv ersps
                itc(find(itc <= imaxmask)) = 0; 
              end;
              iallersps(:,:,end+1) = itc; 
            end;
            allbase(end+1,:) = s.baseline(find(s.complist == clustcomps{nx}(k)),fr);
            if fn == 1
              plotallttls = [plotallttls nx];
              ploticttls = [ploticttls s.complist(k)];
            end;
          end;% for k
        end;  % for ~isempty       
      end;% for nx
      condersps{fn} = eallersps; 
    end;  % for filename      
  else % filename is not a cell (single filename)
    eallersps = zeros(length(fr),length(tms),0);
    iallersps = zeros(length(fr),length(tms),0);
    for nx = 1:length(paths)
      if ~isempty(clustcomps{nx})
        s = load([paths{nx},filename]);    
        %r = load ([paths{nx},'taskbase.mat']);
        fr=find(s.freqs >= freqlims(1) & s.freqs<=freqlims(2));  
        tms=find(s.times > timelims(1) & s.times<timelims(2));
        gdcomps{nx} = s.complist;
        for k=1:length(clustcomps{nx})
          eformask = s.ersp_boot(fr,:,find(s.complist == clustcomps{nx}(k)));% fr X 2
          eminmask = eformask(:,1);
          emaxmask = eformask(:,2);
          eminmask = repmat(eminmask,[1 length(tms)]);
          emaxmask = repmat(emaxmask,[1 length(tms)]);
          ersp = s.comp_ersp(fr,tms,find(s.complist == clustcomps{nx}(k)));
          ersp(find(eminmask <= ersp& ersp <= emaxmask)) = 0;          
          eallersps(:,:,end+1) = ersp; 
          if strcmp(plotitc,'itc') % plot ITC too
            iformask = s.itc_boot(fr,:,find(s.complist == clustcomps{nx}(k)));% fr X 2
            imaxmask = iformask(:,1);
            imaxmask = repmat(imaxmask,[1 length(tms)]);
            itc = s.comp_ersp(fr,tms,find(s.complist == clustcomps{nx}(k)));
            itc(find(itc <= imaxmask)) = 0;          
            iallersps(:,:,end+1) = itc; 
          end;
          allbase(end+1,:) = s.baseline(find(s.complist == clustcomps{nx}(k)),fr);
          plotallttls = [plotallttls nx];
        end;
        fprintf('..%s..',int2str(nx));
      end;                
    end
    condersps{1} = eallersps;            
   
    fnm= filename;clear filename
    filename{1} = fnm;
  end;        
  for cond = 1:length(condersps)
    if ~isempty(binomalpha)
      [plotmat] = GroupSig(condersps{cond},s.alpha,binomalpha,'binom');                  
    else
      plotmat = condersps{cond};
    end;    
    if isempty(clim)
      clim = max(max(abs(mean(plotmat,3))));
    end;
    if plotall == 1 
      figure; rw = round(sqrt(size(condersps{cond},3)));
      cl = ceil(sqrt(size(condersps{cond},3)));
      for m = 1:size(condersps{cond},3)
        sbplot(rw,cl,m);
        if strcmp(freqscale,'quad')
          quadimagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); 
          hold on;
        elseif strcmp(freqscale,'log')
          mylogimagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); 
          hold on;
        else
          imagesc(s.times(tms),s.freqs(fr),condersps{cond}(:,:,m),[-clim clim]); hold on;
        end;  
        plot([0 0],[get(gca,'ylim')],'k-');
        title(['Subj ',int2str(plotallttls(m)),'; IC ',int2str(ploticttls(m))]);
      end; cbar;
      ph = textsc(['Condition ',int2str(cond)],'title');
      set(ph,'fontsize',14);
    end;
    figure(fg);
    sbplot(row,col,pl);pl = pl+1;
    if strcmp(freqscale,'quad')
      quadimagesc(s.times(tms),s.freqs(fr), mean(plotmat,3),[-clim clim]); 
      hold on;
    elseif strcmp(freqscale,'log')
      mylogimagesc(s.times(tms),s.freqs(fr),mean(plotmat,3),[-clim clim]); 
      hold on;
    else
      imagesc(s.times(tms),s.freqs(fr),mean(plotmat,3),[-clim clim]); hold on;
    end;            
    set(gca,'fontsize',12);
    %set(gca,'ydir','norm'); hold on;
    plot([0 0],[get(gca,'ylim')],'k-');
    title(filename{cond}(1:end-4));
    if cond == length(condersps)
      cbar;
    end;
  end;  
  if strcmp(plotitc,'itc')
    for cond = 1:3
      if ~isempty(binomalpha)
        [itcmat] = GroupSig(conditcs{cond},s.alpha,binomalpha,'binom');  
      else
        itcmat = condersps{cond};
      end;           
      if isempty(ilim)
        ilim = max(max( abs(mean(itcmat,3))));
      end;
      if plotall == 1
        figure; rw = round(sqrt(size(conditcs{cond},3)));
        cl = ceil(sqrt(size(conditcs{cond},3)));
        for m = 1:size(conditcs{cond},3)
          sbplot(rw,cl,m);
          if strcmp(freqscale,'quad')
            quadimagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); 
            hold on;
          elseif strcmp(freqscale,'log')
            mylogimagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); 
            hold on;
          else
            imagesc(s.times(tms),s.freqs(fr),conditcs{cond}(:,:,m),[-clim clim]); hold on;
          end;            
          plot([0 0],[get(gca,'ylim')],'k-');pl = pl+1;
        end; 
      end;
      figure(fg);
      sbplot(row,col,pl);pl = pl+1;
      if strcmp(freqscale,'quad')
        quadimagesc(s.times(tms),s.freqs(fr), mean(itcmat,3),[-ilim ilim]); 
        hold on;
      elseif strcmp(freqscale,'log')
        mylogimagesc(s.times(tms),s.freqs(fr),mean(itcmat,3),[-ilim ilim]); 
        hold on;
      else
        imagesc(s.times(tms),s.freqs(fr),mean(itcmat,3),[-ilim ilim]); hold on;
      end;
      set(gca,'fontsize',12);
      %set(gca,'ydir','norm'); hold on;
      plot([0 0],[get(gca,'ylim')],'k-');
      if cond == 3
        title('Difference');
      end; 
    end;
  end;   
  end;  
  
    clustcps{1} = clustcomps;
    ttl = '';
    viewnum=[1,2,3]; % top,sag,rear
    PlotDipoleClusters(dipset,paths,gdcomps,clustcps,1,row,col,pl,ttl,viewnum,[1 0 0],'off','off');
    pl = pl+length(viewnum);
    freqs = s.freqs(fr);
    set(gcf,'color','w');
