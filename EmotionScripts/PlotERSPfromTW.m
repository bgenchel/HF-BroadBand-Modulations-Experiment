% takes data from ERSP TW decomp and plots
% will check for whether one or more ICs were included in decomp
% correl -- [cell array of strings]
% sortfac -- [string] single sorting factor. Field of EEG.event

function [corrs,boots,sigs] = PlotERSPfromTW(savedat,fullpath,dims,backproj,correl,sortfac)

    shuffnum = 500; % for bootstrap for correlation
    corrs=[]; boots = []; sigs = []; % only if correlation requested
    if ~exist('correl')
      correl = [];
    elseif isempty(correl)
      correl = [];
    end;
    if ~exist('backproj')
      backproj = 'off';
    elseif isempty(backproj)
      backproj = 'off';
    end;
      
    s = load ([fullpath,savedat,'.mat']);
    if ~isfield(s,'pcs')
      s.pcs = s.numrows;
    end;
    dat=floatread([fullpath,savedat,'.fdt'],[s.numrows s.numcols],[],0);        
    sph=floatread([fullpath,savedat,'.sph'],[s.pcs s.pcs],[],0); 
    wts=floatread([fullpath,savedat,'.wts'],[s.pcs s.pcs],[],0);        
    ws = wts*sph;   winv = pinv(ws);   acts = ws*dat;    
    
    speceig = floatread([fullpath,s.eigfile],[s.origdims inf],[],0); % PCA eigvecs
    % do the following below in case there are multiple ICs
    %specwts = speceig*winv;  % mix with ICA winv
    %winv = specwts;    % t*f*nics x npcs, new winv (templates) 
    
    if ~exist('dims')
      dims = [];
    elseif isempty(dims)
      dims = [1:size(winv,2)];
    end;
    acts = ws*dat;  % npcs x nepochs  (weights)

    if ~isempty(correl) % mutually exclusive with ERSP plotting
      EEG = pop_loadset('filename',s.datset,'filepath',fullpath);
      cols = hsv(length(correl));
      allz=[]; % use only time-lockin events
      for ep = 1:length(EEG.epoch)
        zlat = cell2mat(EEG.epoch(ep).eventlatency);
        zlat = find(zlat==0);
        allz = [allz,EEG.epoch(ep).event(zlat)];
      end;
      %figure; 
      row = round(sqrt(length(dims)));col = ceil(sqrt(length(dims)));
      for d = 1:length(dims)
        dim = dims(d);
        sbplot(row,col,d);
        pres = cell2mat({EEG.event(allz).OrderPresented});
        for e = 1:length(correl)        
          str =['fac = cell2mat({EEG.event(allz).',correl{e},'});']; eval(str)
          dd = dims(d);
          cc = corrcoef([fac',acts(dd,:)']);
          %cc = corrcoef([fac(find(pres==1))',acts(dd,find(pres==1))']);
          corrs(e,d) = cc(1,2);
          for b = 1:shuffnum
            f = shuffle(fac); a = shuffle(acts(dd,:));
            %f = shuffle(fac(find(pres==1))); a = shuffle(acts(dd,find(pres==1)));
            bc = corrcoef([f',a']);bcorr(1,b) = bc(1,2);
          end;
          bcorr = abs(bcorr);
          boots(e,1) = min(bcorr); boots(e,2) = max(bcorr);
          if corrs(e,d) > boots(1,2)
            sigs(e,d) = 1;
          else
            sigs(e,d) = 0;
          end;
          plot(acts(dd,:),fac,'.','color',cols(e,:)); hold on;
        end;
        title(['D ',int2str(dd),'; r = ',num2str(corrs(e,d)),' [',num2str(boots(1,1)),' ',num2str(boots(1,2)),']']);
      end;
    
    else    
      
      % determine rows/columns
    if s.nics > 3
      row = s.nics;
      col = 6;
    elseif s.nics == 1
      col = ceil(sqrt(s.pcs)); row = round(sqrt(s.pcs));
    else
      row = s.nics*2;
      col = 6;
    end;
    if length(dims) > 6 & s.nics > 1
      col = 6;
    elseif  length(dims) > 6 & s.nics == 1
      col = ceil(sqrt(s.pcs)); row = round(sqrt(s.pcs));
    else
      col = length(dims);
    end;
   
    if strcmp(backproj,'on')
      col = ceil(sqrt(s.numcols)); row = round(sqrt(s.numcols));
      bp = winv*dat;
      c=0;
      for ic = 1:s.nics
        for dm = 1:length(dims)
          dim = dims(dm);
          onewinv = winv(c+1:c+s.pcmat,dim);
          bp = onewinv*acts(dim,:);
          ersp = speceig(:,c+1:c+s.pcmat)*bp; % ERSP x ICs
          ersp = reshape(ersp,[length(s.freqs) length(s.times) size(ersp,2)]); % f x t x trials
          if ~isempty(sortfac)
            EEG = pop_loadset('filename',s.datset,'filepath',fullpath);
            allz=[]; % use only time-lockin events
            for ep = 1:length(EEG.epoch)
              zlat = cell2mat(EEG.epoch(ep).eventlatency);
              zlat = find(zlat==0);
              allz = [allz,EEG.epoch(ep).event(zlat)];
            end;
            str =['fac = cell2mat({EEG.event(allz).',sortfac,'});']; eval(str)
           [v i] = sort(fac);
           ersp = ersp(:,:,i); % sort based on event variable
          end;
          clim = max(max(abs(ersp(:))));clear bp
          figure; pl = 1; 
          for trl = 1:size(ersp,3)
            sbplot(row,col,pl);pl = pl+1;
            if strcmp(s.freqscale,'quad')
              quadimagesc(s.times,s.freqs,ersp(:,:,trl),[-clim clim]);
            elseif strcmp(s.freqscale,'log')
              mylogimagesc(s.times,s.freqs,ersp(:,:,trl),[-clim clim]);
            else
              imagesc(s.times,s.freqs,ersp(:,:,trl),[-clim clim]);
            end;
            hold on;set(gca,'ydir','norm'); plot([0 0],[get(gca,'ylim')],'k-');
            title(['trial ',int2str(trl)]); 
          end;
          textsc(['IC ',int2str(s.complist(ic)),'; Dim ',int2str(dim)],'title');
        end;      
        c = c + s.pcmat;
      end;
        
    else
      
    clim = max(max(abs(speceig*winv)));
    %clim = max(max(abs(speceig)));
    d = 0; last = 0;
    for pg = 1:ceil((length(dims)*s.nics)/(row*col))
      if s.nics > 1
        if d+col > s.pcs & last == 0
          enddim = s.pcs; col = enddim - d;last = 1;
        elseif d+col > s.pcs & last == 1
          break;
        else
          enddim = d+col;
        end;
      else
        enddim = length(dims);
      end;
      figure; pl = 1; c=0; 
     
      for ic = 1:s.nics
        for dm = d+1:enddim
          dim = dims(dm);
          onewinv = winv(c+1:c+s.pcmat,dim);
          %ersp = speceig(:,dim); % this will show the PCA dims and magnitude
          ersp = speceig(:,c+1:c+s.pcmat)*onewinv;
          ersp = reshape(ersp,[length(s.freqs) length(s.times)]);
          
          sbplot(row,col,pl);pl = pl+1;
          if strcmp(s.freqscale,'quad')
            quadimagesc(s.times,s.freqs,ersp,[-clim clim]);
          elseif strcmp(s.freqscale,'log')
            mylogimagesc(s.times,s.freqs,ersp,[-clim clim]);
          else
            imagesc(s.times,s.freqs,ersp,[-clim clim]);
          end;
          hold on;set(gca,'ydir','norm'); plot([0 0],[get(gca,'ylim')],'k-');
          title(['IC ',int2str(s.complist(ic)),'; Dim ',int2str(dim)]); 
        end;      
        c = c + s.pcmat;
      end;
      d = d+col;
    end;
    end;
    end;