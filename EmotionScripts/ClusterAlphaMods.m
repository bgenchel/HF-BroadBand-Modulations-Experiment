% only to be used for modulator templates that have their highest
% point in the alpha range (8-12 Hz)
%
%
% [finaltempls, finalmeans, finalidx] = ClusterAlphaMods(alphaclust,freqs,frqlims,stdcut,minval);
%
% Categorizes alpha templates as 1) below, 2) at or 3) above mean peak
% alpha frequency. This algorithm assumes that the templates have already
% been stretched so that the mean peak alpha is 10 Hz. 
%


function [finaltempls, finalmeans, finalidx] = ClusterAlphaMods(alphaclust,freqs,frqlims,freqscale,stdcut,minval);
    
    
    lowfrq = frqlims(1); hifrq = frqlims(2);
    lowalph = find(freqs < lowfrq);lowalph = lowalph(end);
    hialph = find(freqs > hifrq); hialph = hialph(1);
    
    if ~isempty(minval)% delete templates with low max val
        fr60 = find(freqs > 57 & freqs < 63);
        delmem = [];
        for mem = 1:size(alphaclust{1},1)
            [maxval maxidx] = max(abs(alphaclust{1}(mem,:)));
            if maxval < minval
                delmem = [delmem mem];
            end;
        end;
        alphaclust{1}(delmem,:) = [];
        alphaclust{2}(delmem,:) = [];
        alphaclust{3}(delmem,:) = [];
    end;  

    lowclust = [];medclust=[];   hiclust = [];
    lowmeans = [];medmeans = []; himeans = [];
    lowidx = [];  medidx=[];     hiidx=[];
    
    for t = 1:size(alphaclust{3},1)
        [val idx] = max(alphaclust{1}(t,:));
        if idx <= lowalph
            lowclust = [lowclust; alphaclust{1}(t,:)];
            lowmeans = [lowmeans; alphaclust{2}(t,:)];
            lowidx = [lowidx; alphaclust{3}(t,:)];
        elseif idx >= hialph
            hiclust = [hiclust; alphaclust{1}(t,:)];
            himeans = [himeans; alphaclust{2}(t,:)];
            hiidx = [hiidx; alphaclust{3}(t,:)];
        else            
            medclust = [medclust; alphaclust{1}(t,:)];
            medmeans = [medmeans; alphaclust{2}(t,:)];
            medidx = [medidx; alphaclust{3}(t,:)];
        end;
    end;
    finalidx{1} = lowidx;
    finalidx{2} = medidx;
    finalidx{3} = hiidx;
    
    finalmeans{1} = lowmeans;
    finalmeans{2} = medmeans;
    finalmeans{3} = himeans;
    
    finaltempls{1} = lowclust;
    finaltempls{2} = medclust;
    finaltempls{3} = hiclust;
    
    if ~isempty(stdcut)
        for cls = 1:length(finaltempls)
            zs = zscore(finaltempls{cls});
            delmem = [];
            for mem = 1:size(finaltempls{cls},1)
                if find(mean(abs(zs(mem,:))) > stdcut)
                    delmem = [delmem mem];
                end;
            end;
            if ~isempty(delmem)
                finaltempls{cls}(delmem,:) = [];
                finalmeans{cls}(delmem,:) = [];
                finalidx{cls}(delmem,:) = [];
            end;
        end;        
    end;
    fr = find(freqs < 25);
    figure; 
    for cls = 1:length(finalidx)
        sbplot(2,2,cls);
        if ~isempty(finalidx{cls}) % can be empty depending on freq limits
        if strcmp(freqscale,'quad')
            [handle realx labelx] = quadplot(freqs(fr),finaltempls{cls}(:,fr)',2,'m'); hold on;
            [handle realx labelx] = quadplot(freqs(fr),mean(finaltempls{cls}(:,fr),1)',2,'b');
            plot([realx(find(labelx==10)) realx(find(labelx==10))],[get(gca,'ylim')],'g-');
        elseif strcmp(freqscale,'log')
            [handle realx labelx] = logplot(freqs(fr),finaltempls{cls}(:,fr)',2,'m'); hold on;
            [handle realx labelx] = logplot(freqs(fr),mean(finaltempls{cls}(:,fr),1)',2,'b');
            plot([realx(find(labelx==10)) realx(find(labelx==10))],[get(gca,'ylim')],'g-');
            set(gca,'xlim',[freqs(fr(1)) freqs(fr(end))]);
        else
            plot( freqs,finaltempls{cls}','m','linewidth',2); hold on;
            plot( freqs,mean(finaltempls{cls},1)','b','linewidth',2);
            plot([10 10],[get(gca,'ylim')],'g-');
            set(gca,'xlim',[freqs(fr(1)) freqs(fr(end))]);
        end;        
        plot([get(gca,'xlim')],[0 0],'k-'); hold on;
        title([int2str(cls),'-',int2str(size(finaltempls{cls},1))]);
        end;
    end;
