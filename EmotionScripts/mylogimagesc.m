% plots a color images that has log-spaced frequencies
%
% Usage:  >> [realy] = mylogimagesc(times,freqs,data,lim,ticks);
%
% Input:
%   times = vector of x-axis values
%   freqs = vector of y-axis values
%   data  = matrix of size (freqs,times)
%   lim   = [min max] for color axis. Default is abs(max(data(:)))
%   ticks = [vector] of y-axis values to create ticks and ticklabels for
%
% Author Julie Onton, Mar 1, 2007


function [realy,labely] = mylogimagesc(times,freqs,data,lim,ticks)
    
    if ~exist('ticks')
      yt = [4,10,20,40,100,200];
    elseif isempty(ticks)
      if freqs(end) < 40
        yt = [2,4,10,20,30]; % ytick labels
      elseif freqs(end) > 200
        yt = [10,70,200,400]; % ytick labels
      else 
        yt = [4,10,20,40,100,200]; % ytick labels
      end;
    else
      yt = ticks;      
    end;
    yt = yt(find(yt > freqs(1) & yt < freqs(end)));

    if ~exist('lim')
        if max(abs(data(:))) == 0
            lim = [-1 1];
        else
            lim = [-max(abs(data(:))) max(abs(data(:)))];
        end;
    elseif isempty(lim)
        if max(abs(data(:))) == 0
            lim = [-1 1];
        else
            lim =  [-max(abs(data(:))) max(abs(data(:)))];
        end;
    end;
    if length(lim) < 2
        lim = [-lim lim];
    end;
    
    
    lfreqs = log(freqs); % find linear freqs which is what is plotted by default
    
    imagesc(times,lfreqs,data,lim); % plot the image (linear-spaced freqs)
    set(gca,'ydir','normal')
    
    
    for f = 1:length(yt)
        [val mch(f)] = min(abs(freqs-yt(f)));
    end;
    realy = lfreqs(mch); % find the closest freq in the linearly spaced freqs      
    
    set(gca,'ytick',realy); % appropriately label
    set(gca,'yticklabel',yt);
      labely = yt;
 
