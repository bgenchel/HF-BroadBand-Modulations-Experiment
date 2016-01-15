% Plot trial by trial profiles of context/ERSP based on context decomp
%
%
%
%
%

function PlotSnglTrialContext(name,fullpath);


s = load([fullpath,name,'.mat']);  
wts = floatread([fullpath,name,'.wts'],[s.numrows s.numrows],[],0);
sph = floatread([fullpath,name,'.sph'],[s.numrows s.numrows],[],0);     
dat = floatread([fullpath,name,'.fdt'],[s.numrows inf],[],0);   
erspeig = floatread([fullpath,s.eigfile],[length(s.freqs)*length(s.times) s.pcmat],[],0);
cxteig = floatread([fullpath,s.cxteigfile],[length(s.cxtmean) s.pcctx],[],0);
ws = wts*sph;winv = pinv(ws);acts = ws*dat;
ersptmpl = erspeig*winv(1:s.pcmat,:);  % templates (context winv)
cxttmpl = cxteig*winv(s.pcmat+1:end,:);  % templates (ersp winv)
bpcontext = zeros(size(cxttmpl,1),size(acts,2),size(acts,1));
bpersp = zeros(size(ersptmpl,1),size(acts,2),size(acts,1));
bp = [ersptmpl;cxttmpl]*acts;
for dim = 1:size(acts,1) % back proj each dim individually
   bpcontext(:,:,dim) = cxttmpl(:,dim)*acts(dim,:);
   bpersp(:,:,dim) = ersptmpl(:,dim)*acts(dim,:);
end;

lim = max(max(abs(bpcontext(:,:,dim))));
figure;pl=1;
row = 4; col = 4;
for t = 1:size(acts,2) % for each trial
   if pl > row*col
      input('Plot more?')
      figure; pl=1;
   end;
   sbplot(row,col,pl); pl=pl+1; hold on;
   plot([1:size(bpcontext,1)],bp(length(s.times)*length(s.freqs)+1:end,t),'k-','linewidth',2);
   plot([1:size(bpcontext,1)],squeeze(bpcontext(:,t,:)),'linewidth',1);
   set(gca,'xlim',[1 size(bpcontext,1)]);
   set(gca,'xticklabel',[]); set(gca,'ytick',[-lim lim]);
   set(gca,'yticklabel',[]);
   title(int2str(t));
 end;
row = 4; col = 4;
for t = 1:size(acts,2) % for each trial
   lim = max(max(abs(bpersp(:,t,:))));
   figure;pl=1;
   for dim = 1:size(acts,1)
      if pl > row*col
         figure; pl=1;
      end;
      sbplot(row,col,pl); pl=pl+1;
      tmpersp = reshape(bpersp(:,t,dim),length(s.freqs),length(s.times));
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
      title(['Dim ',int2str(dim)]);
   end;
   sbplot(row,col,pl); pl=pl+1;
   tmpersp = reshape(bp(1:length(s.times)*length(s.freqs),t),length(s.freqs),length(s.times));
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
   title(['Full Back Proj']);
   input('Plot more?')
end;





