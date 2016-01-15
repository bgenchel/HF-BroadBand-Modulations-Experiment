% plots specified scalp maps from multiple subjects
%
% PlotScalpMaps(datset,paths,complist,ttl,savettl,rowcolpl,specs)
%
%
%
% savettl -- [string] fullpath plus savename and extension.
%            ex: '/data/common2/emotion/MyScalpMaps'
%            function will add page number and .jpg extension
% specs -- ['on' or 'off'] 'on' will plot spectra along with scalp map. 
%          Default is 'off'

function PlotScalpMaps(datset,paths,complist,ttl,savettl,rowcolpl,specs)
    
    if ~exist('specs')
        specs = 'off';
    end;

    if ~exist('rowcolpl')
        rowcolpl = [];
    end;
    if ~exist('savettl')
        savettl = [];
    end;
    howmany= 0;
    for nx = 1:length(complist)
        if ~isempty(complist{nx})
            hm = length(complist{nx});
            howmany = howmany+hm;
        end;
    end;
    if ~isempty(savettl)
        pg = 1; % if saving figures, start a page count
    end;
    if isempty(rowcolpl)
       figure; pl = 1;
       if strcmp(specs,'on')
          frqlim = [2 70]; % 70
          freqs = linspace(frqlim(1), frqlim(end), frqlim(2));
          col = 12; % 3 ICs per row
          row = 5; 
       else          
          row=round(sqrt(howmany)); col=ceil(sqrt(howmany));
          if row > 9
             row = 4; col = 4;
          end;
       end;
    else% plot into existing figure;
      if strcmp(specs,'on')
        frqlim = [2 70]; % 70
        freqs = linspace(frqlim(1), frqlim(end), frqlim(2));
        col = 12; % 3 ICs per row
        row = 5; 
      end          
      row = rowcolpl(1); col = rowcolpl(2); pl = rowcolpl(3);
    end;
    for nx = 1:length(complist)
        if ~isempty(complist{nx})
            EEG = pop_loadset(datset ,paths{nx}); 
            if strcmp(specs,'on')
               pwr = zeros(length(complist{nx}),length(freqs));
               wavelets = computepsifamilyQodd(freqs,1/EEG.srate,frqlim(1),frqlim(2),EEG.srate+1); %
            end;
            for c = 1:length(complist{nx})
                if pl > row*col
                    ph = textsc(ttl,'title');  set(ph,'fontsize',14);
                    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
                    if ~isempty(savettl)
                        str = ['print ',savettl,'-',int2str(pg),'.jpg -djpeg']; eval(str);pg = pg+1;
                    end;
                    figure; pl = 1;
                end;
                sbplot(row,col,pl)
                topoplot(EEG.icawinv(:,complist{nx}(c)),EEG.chanlocs,'electrodes','off');
                hold on;     set(gca,'fontsize',14);
                title([int2str(nx),'-',int2str(complist{nx}(c))]); 
                pl = pl+1; %cbar;
                if strcmp(specs,'on') % calculates from continuous data
                   alldat = squeeze(EEG.icaact(complist{nx}(c),:));
                   %alldat = squeeze(EEG.data(complist{nx}(c),:));
                   xx =floor(length(alldat)/(EEG.srate+1));
                   alldat = alldat(1,1:xx*(EEG.srate+1));
                   alldat = reshape(alldat,[(EEG.srate+1) xx]);
                   pwr1 = wavelets * alldat;
                   pwr1 = pwr1.*conj(pwr1);% pwr1 is now non-complex
                   pwr1 = 10*log10(pwr1);% convert to log
                   pwr = mean(pwr1,2)'; clear pwr1
                   sbplot(row,col,[pl pl+2]); pl=pl+3;
                   plot(freqs,pwr,'b-','linewidth',2);title(['IC ',int2str(complist{nx}(c))]);
                   set(gca,'xlim',[freqs(1) freqs(end)]); set(gca,'xgrid','on');
                end;
            end;
        end;
    end;
    ph = textsc(ttl,'title');  set(ph,'fontsize',14);
    set(gcf,'PaperOrientation','portrait');set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    if ~isempty(savettl)
        str = ['print ',savettl,'-',int2str(pg),'.jpg -djpeg']; eval(str);
    end;
    

