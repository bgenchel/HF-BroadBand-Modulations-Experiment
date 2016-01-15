% takes context info from TW weightings for individual ICs and plots scalp map,template and context
%
% [cxtout] = PlotTWcontext(filename,datset,fullpath,gdcomps,plottype,labels,alpha,savettl,plton);
%
%
% filename -- name of .mat file with context decomp variables
% datset -- name of a .set file that contains ICA weights for scalp map plotting
% fullpath -- full data directory path
% gdcomps -- vector of ICs to plot (must be included in 'complist' in filename.mat
% plottype -- ['winv' or 'wtdmean'] 'winv' plots the ICA/PCA templates, 'wtdmean' plots
%             the weighted means of each trial's ERSP/context according to ICA weights
%             and masks for significance.
% labels -- [cell array of strings] label names for each context question.
% alpha -- [pvalue] for bootstrap masking of weighted mean data (ERSP and context)
% savettl -- [string] If not [], will save all figures with this filename and page number as
%            jpeg images.
% plton -- ['on' or 'off'] if 'off', will not create figures, but return the 'cxtout'
% addmat - [matrix] straight context answers for only epochs decomposed.
% Only need if ERSP only or no context PCA
% OUTPUT:
% cxtout -- [matrix] matrix of context templates or weighted means for all requested ICs
%           matrix sizes: (# context questions) x (# Modes) x (# ICs)
%
%
% Author Julie Onton, Swartz Center for Computational Neuroscience.
%

function [cxtout,erspout] = PlotTWcontext(filename,datset,fullpath,gdcomps,plottype,labels,alpha,savettl,plton,plotdims)

checksig = 'on'; % 'off': outputs unmasked ERSP/context and does not skip dims with non-sigs
lfontsz = 14; % font size for context labels
if ~exist('plotdims')
   plotdims = [];
end;
stem2 = 'TWwarpSubQsIC';
cntxtonly = 'off'; % 'on' turns off ERSP bootstrap for speed
percersp = .0001; % percent of ERSP values that need to be sig for output

pvafcalc = 0; % plot pvaf of context instead

shuffnum = 500; % number of bootstrap iterations for wtd mean masking
EEG = pop_loadset(datset, fullpath);

col = 6;   row = 5;

name = [filename,int2str(gdcomps(1))];
s = load([fullpath,name,'.mat']);
cxtout = zeros(length(labels),s.numrows,length(gdcomps));
erspout = zeros(length(s.times)*length(s.freqs),s.numrows,length(gdcomps));
addmat = s.idxmat;

pl = 1; pg = 1;
cols = hsv(length(labels)); cols(3,:) = [.8 .8 0];
figure;
for ic = 1:length(gdcomps)
   %%%  load decomp data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   name = [filename,int2str(gdcomps(ic))];
   s = load([fullpath,name,'.mat']);
   wts = floatread([fullpath,name,'.wts'],[s.numrows s.numrows]);
   sph = floatread([fullpath,name,'.sph'],[s.numrows s.numrows]);
   dat = floatread([fullpath,name,'.fdt'],[s.numrows s.numframes]);
   ws = wts*sph;winv = pinv(ws);acts = ws*dat; acts = acts*-1;
   if isempty(plotdims)
      plotdims = [1:s.numrows];
   end;
   if ~isempty(s.eigfile)
      if s.pcmat > 0
         if s.pcctx == s.pctot % context only
            ss = load([fullpath,stem2,int2str(gdcomps(ic)),'.mat']); % hijacking for context only decomp
            dat2 = floatread([fullpath,stem2,int2str(gdcomps(ic)),'.fdt'],[ss.numrows ss.numframes]);
            erspeig = floatread([fullpath,stem2,int2str(gdcomps(ic)),'EIGVEC.fdt'],[length(ss.freqs)*length(ss.times) ss.pcmat]);
            erspdat = erspeig*dat2(1:ss.pcmat,:);% back-proj to orig PCA-reduced data
            ersptmpl = erspeig*winv(1:ss.pcmat,:);  % (garbage) templates
         else
            erspeig = floatread([fullpath,name,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcmat]);
            erspdat = erspeig*dat(1:s.pcmat,:);% back-proj to orig PCA-reduced data
            ersptmpl = erspeig*winv(1:s.pcmat,:);  % templates
         end;
      end;
   end;
   if ~isempty(s.cxteigfile) % if PCA performed on context
      cxteig = floatread([fullpath,name,'ADDEIG.fdt'],[length(labels) s.pcctx],[],0);
      if s.pcctx == s.pctot % context only
         cxttmpl = cxteig*winv;  % context only decomp
         cxtdat = cxteig*dat;  % back-proj to orig data
         cxtdat = cxtdat + repmat(s.cxtmean,[1 size(cxtdat,2)]); % add back mean
      else
         cxttmpl = cxteig*winv(s.pcmat+1:end,:);  % templates
         cxtdat = cxteig*dat(s.pcmat+1:end,:);  % back-proj to orig data
         cxtdat = cxtdat + repmat(s.cxtmean,[1 size(cxtdat,2)]); % add back mean
      end;
   elseif isempty(s.cxteigfile) & ~isempty(s.anaeigfile) % split matrices between analog and binary
      anaeig = floatread([fullpath,name,'ANAEIG.fdt'],[length(s.anarows) inf],[],0);
      bineig = floatread([fullpath,name,'BINEIG.fdt'],[length(s.binrows) inf],[],0);
      tmpwinv = winv(s.pcmat+1:end,:); % separate to use ana/bin rows
      anawinv = tmpwinv(s.anarows,:);
      binwinv = tmpwinv(s.binrows,:);
      anatmpl = anaeig*anawinv;
      bintmpl = bineig*binwinv;
      tmpdat = dat(s.pcmat+1:end,:);
      anadat = tmpdat(s.anarows,:);anadat = anaeig*anadat;
      anadat = anadat + repmat(s.cxtmean(s.anarows,1),[1 size(anadat,2)]);
      bindat = tmpdat(s.binrows,:);bindat = bineig*bindat;
      bindat = bindat + repmat(s.cxtmean(s.binrows,1),[1 size(bindat,2)]);
      
      for a = 1:length(s.anarows)
        cxttmpl(s.anarows(a),:) = anatmpl(a,:);
        cxtdat(s.anarows(a),:) = anadat(a,:);
      end;
      for b = 1:length(s.binrows)
        cxttmpl(s.binrows(b),:) = bintmpl(b,:);
        cxtdat(s.binrows(b),:) = bindat(b,:);
      end;
      
   elseif isempty(s.cxteigfile) & s.pcmat == s.pctot % ERSP only
      cxtdat = s.idxmat; % straight context vectors
      %ss = load([fullpath,stem2,int2str(gdcomps(ic)),'.mat']); % hijacking for context only decomp
      %dat2 = floatread([fullpath,stem2,int2str(gdcomps(ic)),'.fdt'],[ss.numrows ss.numframes],[],0);
      %cxteig = floatread([fullpath,stem2,int2str(gdcomps(ic)),'ADDEIG.fdt'],[length(labels) ss.pcctx],[],0);
      %cxtdat = cxteig*dat2(ss.pcmat+1:end,:);  % back-proj to orig data
      %cxtdat = cxtdat + repmat(ss.cxtmean,[1 size(cxtdat,2)]); % add back mean
   else % no PCA on context, or no context
      if isempty(s.pcctx) % no context
         cxtdat = s.idxmat; % straight context vectors
      else
         cxtdat = s.idxmat; % raw context
         cxttmpl = winv(s.pcmat+1:end,:); % just templates
      end;
   end;

   if strcmp(plottype,'wtdmean')
      fprintf('\nAccumulating bootstrap distribution ...\n');
      wtdersp = (erspdat*acts')/size(erspdat,2);% multiply each trial element,sum then divide by ntrials
      if strcmp(cntxtonly,'off')
         %%%%%%%  calculate weighted ERSPs for each dim:
         clear limersp limmat
         bootwts = zeros(size(wtdersp,1),size(wtdersp,2),shuffnum);
         for b= 1:shuffnum
            bootwts(:,:,b) = (erspdat*shuffle(acts,2)')/size(erspdat,2);
         end;
         bootwts = sort(bootwts,3);
         limersp(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
         limersp(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
         clear bootwts
      end;
      %%%%%%%  calculate context vectors for each dim:
      wtdctx = (cxtdat*acts')/size(erspdat,2);
      bootwts = zeros(size(wtdctx,1),size(wtdctx,2),shuffnum);
      for b= 1:shuffnum
         bootwts(:,:,b) = (cxtdat*shuffle(acts,2)')/size(erspdat,2);
      end;
      bootwts = sort(bootwts,3);
      limmat(:,:,2) = bootwts(:,:,end-ceil(shuffnum*alpha)); % max boot
      limmat(:,:,1) = bootwts(:,:,ceil(shuffnum*alpha));  % min boot
      %mnctx = mean(bootwts,3);
      %stdctx = std(bootwts,0,3);
   end;
   clear sph wts ws dat bootwts
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for dd = 1:length(plotdims)
      d = plotdims(dd);
      if strcmp(plottype,'wtdmean')% only plot sig contexts
         goahead = 0;
         tmpcomp = wtdctx(:,d);
         tmpcomp(find(tmpcomp>limmat(:,d,1)&tmpcomp<limmat(:,d,2))) = 0;
         sigchk = (limmat(:,d,1)-tmpcomp).*(limmat(:,d,2)-tmpcomp);
         if ~isempty(find(sigchk>0))
            goahead = 1;
         else
            if strcmp(checksig,'off')
               goahead = 1; % 0
            else
               goahead = 0; % 0
            end;
         end;
      else
         goahead = 1; % if winv, then plot all
      end;
      if goahead == 1
         if strcmp(plton,'on')
            if pl > row * col
               textsc([filename,'; Subject ',fullpath(end-8:end-1)],'title');
               set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
               if ~isempty(savettl)
                  str = ['print ',fullpath,savettl,int2str(pg),'.eps -depsc2 -adobe -painters']; eval(str)
                  clf;
               else
                  figure;
               end;
               pl = 1; pg = pg+1;
            end;
         else
            clf; % clear figure and continue plotting on same fig
         end;

         if strcmp(plottype,'winv')
            tmpersp =  ersptmpl(:,d);
            erspout(:,d,ic) = tmpersp;
            erspsig = 'sig'; % for later, allows output
         elseif strcmp(plottype,'wtdmean')
            tmpersp = wtdersp(:,d);
            erspout(:,d,ic) = tmpersp;
            if strcmp(cntxtonly,'off')
               % mask the wtd mean by the bootstrap values:
               tmpersp(find(tmpersp>limersp(:,d,1)&tmpersp<limersp(:,d,2))) = 0;
               erspout(:,d,ic) = tmpersp;
               if length(find(tmpersp~=0)) > round(percersp*(length(s.freqs)*length(s.times)))
                  erspsig = 'sig';
               else
                  if strcmp(checksig,'on')
                     erspsig = 'nonsig'; % 'nonsig'
                  else
                     erspsig = 'sig'; % 'nonsig'
                  end;
               end;
            else
               erspsig = 'sig'; % for later, allows output
            end;
         end;
         if strcmp(erspsig,'sig')
            sbplot(row,col,pl);pl = pl+1; % plot histogram of weights
            hist(acts(d,:),100); hold on; plot([0 0],[get(gca,'ylim')],'r-');

            sbplot(row,col,pl);pl = pl+1;
            topoplot(EEG.icawinv(:,gdcomps(ic)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
            title(['IC ',int2str(gdcomps(ic))]);
            sbplot(row,col,pl); pl = pl+1;
            lim = max(abs(tmpersp));
            if lim == 0
               lim = 1;% if no sig points, make lim a valid value
            end;
            tmpersp = reshape(tmpersp,length(s.freqs),length(s.times));
            if strcmp(s.freqscale,'quad')
               quadimagesc(s.times, s.freqs, tmpersp,[-lim lim]);hold on;
            elseif  strcmp(s.freqscale,'log')
               mylogimagesc(s.times, s.freqs, tmpersp,[-lim lim]);hold on;
            else
               imagesc(s.times, s.freqs, tmpersp,[-lim lim]);hold on;
               set(gca,'ytick',[10:10:s.freqs(end)]);
            end;
            set(gca,'ydir','norm'); set(gca,'ticklength',[.02 .02]);
            hold on; plot([0 0],[get(gca,'ylim')],'k-');
            if isfield(s,'medwarpevs')
               wcols = lines(length(s.medwarpevs));
               for wp = 1:length(s.medwarpevs)
                  ph = plot([s.medwarpevs(wp) s.medwarpevs(wp)],[get(gca,'ylim')],'k-');
                  set(ph,'color',wcols(wp,:));
               end;
            end;
            if pl < (row-1)*col
               set(gca,'xticklabel',[]);
            end;
            title(['Dim ',int2str(d)]);
         end;

         if strcmp(plottype,'winv')
            allctxvals = cxttmpl(:,d);
            if strcmp(erspsig,'sig')
               sbplot(row,col,[pl pl+2]); pl = pl+3; hold on;
               plot([1:size(cxttmpl,1)], cxttmpl(:,d)','k-','linewidth',1);hold on;
               ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.75 .75 .75]);
               for ee = 1:length(allctxvals)
                  ph=plot(ee,allctxvals(ee),'.'); hold on;
                  set(ph,'markersize',15);
                  set(ph,'color',cols(ee,:));
                  ph=text(ee,allctxvals(ee),labels{ee});
                  set(ph,'color',cols(ee,:));set(ph,'fontsize',lfontsz);
                  set(ph,'rotation',90);
               end;
               set(gca,'xlim',[0 length(allctxvals)+1]);
            end;
            cxtout(:,d,ic) = allctxvals';%  template
         elseif strcmp(plottype,'wtdmean')
            if pvafcalc == 1% calculate pvaf instead:---
               %keyboard
               % first make back proj of single template/factor and only cxt pcs:
               bp = winv(s.pcersp+1:end,d) * acts(d,:);
               bpq = cxteig * bp;
               % calculate pvaf of PCA-reduced data
               for q = 1:size(bpq,1)
                  pv(1,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
               end;
               pv = pv*100;
               tmpcomp = pv .* (cxttmpl(:,d)'./abs(cxttmpl(:,d)'));
               for b = 1:shuffnum
                  bp = winv(s.pcersp+1:end,d) * shuffle(acts(d,:),2);
                  bpq = cxteig * bp;
                  %bpq = cxteig * shuffle(bp,2);
                  for q = 1:size(bpq,1)
                     pvboot(b,q) = 1 - (var(cxtdat(q,:)' - bpq(q,:)')/var(cxtdat(q,:)'));
                  end;
                  pvboot(b,:) = pvboot(b,:).* (cxttmpl(:,d)'./abs(cxttmpl(:,d)'));
               end;    clear limmat
               pvboot = pvboot*100;
               pvboot = sort(pvboot,1);
               limmat(2,:) = pvboot(end-ceil(shuffnum*alpha),:); % max boot
               limmat(1,:) = pvboot(ceil(shuffnum*alpha),:);  % min boot
               if strcmp(erspsig,'sig')
                  sbplot(row,col,[pl pl+2]); pl = pl+3;
                  for ee = 1:size(wtdctx,1)
                     ph = plot([ee ee],[limmat(1,ee) limmat(2,ee)],'k-','linewidth',3);
                     set(ph,'color',[.8 .8 .8]);hold on;
                  end;
                  title(['Pvaf of Context']);
                  plot([get(gca,'xlim')],[0 0],'k-');
                  if pl < (row-1)*col
                     set(gca,'xticklabel',[]);
                  end;
               end;
               sigchk = (limmat(1,:)-tmpcomp).*(limmat(2,:)-tmpcomp);
            else
               tmpcomp = wtdctx(:,d);
               if strcmp(checksig,'on')
                  tmpcomp(find(tmpcomp>limmat(:,d,1)&tmpcomp<limmat(:,d,2))) = 0;
               end;
               if strcmp(erspsig,'sig')
                  sbplot(row,col,[pl pl+2]); pl = pl+3;
                  for ee = 1:size(wtdctx,1)
                     ph = plot([ee ee],[limmat(ee,d,1) limmat(ee,d,2)],'k-','linewidth',3);
                     set(ph,'color',[.8 .8 .8]);hold on;
                  end;
                  title(['Wtd Mean Context']);
                  plot([get(gca,'xlim')],[0 0],'k-');
                  set(gca,'box','off');
                  if pl < (row-1)*col
                     set(gca,'xticklabel',[]);
                  end;
               end;
               if strcmp(checksig,'on')
                  sigchk = (limmat(:,d,1)-tmpcomp).*(limmat(:,d,2)-tmpcomp);
               else
                  sigchk = 1;
               end;
            end;
            if strcmp(erspsig,'sig')
               ph = plot([get(gca,'xlim')],[0 0],'k-'); set(ph,'color',[.75 .75 .75]);
               plot([1:size(wtdctx,1)], wtdctx(:,d)','k-','linewidth',1);hold on;

               for ee = 1:size(wtdctx,1)
                  ph=plot(ee,tmpcomp(ee),'.');
                  set(ph,'markersize',18);
                  set(ph,'color',cols(ee,:));
                  ph=text(ee,tmpcomp(ee),labels{ee});
                  set(ph,'color',cols(ee,:));set(ph,'fontsize',lfontsz);
                  set(ph,'rotation',45);
               end;
               set(gca,'xlim',[0 size(wtdctx,1)+1]);
               plot([get(gca,'xlim')],[0 0],'k-');
               if pl < (row-1)*col
                  set(gca,'xticklabel',[]);
               end;
            end;
            if ~isempty(find(sigchk>0))
               if strcmp(erspsig,'sig')% only if ersp also sig
                  cxtout(:,d,ic) = tmpcomp';% masked template
                  %cxtout(:,d,ic) = cxttmpl(:,d)';% output template
               else
                  cxtout(:,d,ic) = zeros(1,size(cxtout,1));% nonsig ersp, so zeros
               end;
            else
               cxtout(:,d,ic) = zeros(length(tmpcomp),1); % if no context out of pvaf range
            end;
         end;
      end;
   end;
end;
textsc([filename,'; Subject ',fullpath(end-8:end-1)],'title');
set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
if ~isempty(savettl)
   str = ['print ',fullpath,savettl,int2str(pg),'.eps -depsc2 -adobe -painters']; eval(str)
end;

