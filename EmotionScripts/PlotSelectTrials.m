% selects and plots single 'trials' from emotion clustering listed below
% results 
wts = floatread('/data/common1/emotion/tfXcmpSngltrPC30.wts',[30 5500]);
sph = floatread('/data/common1/emotion/tfXcmpSngltrPC30.sph',[5500 5500]);
ws = wts*sph;   clear wts sph
winv = pinv(ws);
alltrials = floatread('/home/julie/Emos/alltrials.fdt',[ 5500 24481]); % rand epochs tf
activations = ws*alltrials;
load /home/julie/Emos/alltrialStuff.mat

% load a dataset and run one iteration of timef for times/freqs (see beg of script)
subj1 =  [5,6,7,10,11,12,13,15,18,19,22,23,24,25,40,46]; % jo74  for ICA on all time/freq
gdcomps = {subj1};
clusts = size(ws,1);
icaresults = zeros (length(freqs),length(tm),clusts);
for n = 1:size(winv,2)
    tmpcomp = winv (:,n)';  %makes a totalength X 1
    tmpcomp = reshape(tmpcomp,length(freqs),length(tm));  %Good
    icaresults(:,:,n) = tmpcomp;
end; 
figure;
for n = 1:12%clusts
    subplot(3,4,n)
    imagesc(times,freqs,icaresults(:,:,n),[-5 5]);hold on;
    set(gca,'ydir','norm');
    set(gca,'ticklength',[.02 .02]);    
    if n < 9
    set(gca,'xticklabel',[]);
    end;   
    set(gca,'fontsize',14);
    title(int2str(n));
end; colorbar;

% find where each emotion starts and stops for correlated trials only run
load /data/common1/emotion/Corridx.mat 
emos = {'awe', 'frustration','joy','anger','sad','happy','fear','love' ,'jealousy','compassion','emabarrass','content','grief','relief'};
for nx = 1:length(gdcomps)
    for em = 1:length(emos)
        oneemocorrs = allcorridx{em};
        for cp = 1:length(gdcomps{nx})
        ntcomp(1,cp) = size(oneemocorrs{cp},2);
        end;
        numtrials{em} = ntcomp;
    end;
end;

% Find > zero activations for clusters 1:12
for n = 1:10%size(activations,1)
    hiact{n} = find(activations(n,:) > std(activations(n,:)));
end;
% find sum of positive activations for all components
for clust = 1:10
    numtot = 0; clear sums
    sc = hiact{clust};
    for nx = 1:length(gdcomps)
        if ~isempty(gdcomps{nx})
            if nx>1
                skipact = allsubjix{nx-1};
            end;
            for em = 1:length(emos)
                for cp = 1:length(gdcomps{nx})
                    ntrials = numtrials{em}(cp); % will have to change to cell array for more subj
                    ft = find(sc>=numtot+1 & sc<=numtot+ntrials);
                    sumcomp = sum(activations(clust,sc(ft)));
                    sums(cp) = sumcomp/ntrials*100;
                    numtot = ntrials+numtot;
                end;             
                emosumsonesubj{em} = sums;
            end;            
        end;
        allsubjsums{nx} = emosumsonesubj;
    end;
    allclustsums{clust} = allsubjsums;
end;


% Select trials that are positive or negative for a cluster and plot
% Find > zero activations for clusters 1:12
figure;    % ,8,9
clustset = [1,4,5,6];
for clust = 1:4
    p = clust;
    clust = clustset(clust);
    for cond = 1:4
        if cond == 1
            hiact{clust} = find(activations(clust,:) > 1);
            lim = 15;
        elseif cond == 2
            hiact{clust} = find(activations(clust,:) < -1);
            lim = 10;
            
        elseif cond == 3
            hiact{clust} = find(activations(clust,:) < 1 & activations(n,:) > -1);            
            lim = 5;
        elseif cond == 4
            hiact{clust} = find(activations(clust,:) < .2 & activations(n,:) > -.2);            
            lim = 5;    
        end;
        % find trials selected by +- activation above
        clear allsubjtris allclusttris
        sc = hiact{clust};
        numtot = 0;
        for nx = 1:length(gdcomps)
            if ~isempty(gdcomps{nx})
                clear alltris tris
                    for cp = 1:length(gdcomps{nx})
                        ntrials = numtrials{em}(cp); % will have to change to cell array for more subj
                        ft = find(sc>=numtot+1 & sc<=numtot+ntrials);
                        tris = alltrials(:,sc(ft));
                        alltris{cp} = tris;  % putting in a 2800 X ? matrix
                        numtot = ntrials+numtot;
                    end;             
                allsubjtris{nx} = alltris;
            end;
        end;
        % plot results:

        thisclust = zeros(size(alltrials,1),0);
        for nx = 1:length(gdcomps)
            if ~isempty(allsubjtris{nx})
                for cp = 1:length(gdcomps{nx})
                    thisclust(:,end+1:end+size(allsubjtris{nx}{cp},2)) = allsubjtris{nx}{cp};
                end;
            end;
        end;
        thisclust = reshape(thisclust,length(freqs),length(tm),size(thisclust,2));
        subplot(4,4,p)
        imagesc(times,freqs,mean(thisclust,3),[-lim lim]); p = p+4;
        set(gca,'ydir','norm'); 
        set(gca,'fontsize',14);
        set(gca,'xticklabel',[]);
        if cond == 1
            title(['Cluster: ',int2str(clust),'; ',int2str(size(thisclust,3)),' Trials; ',int2str(size(thisclust,3)/size(alltrials,2)*100),'%']);
        else
            title([int2str(size(thisclust,3)),' Trials; ',int2str(size(thisclust,3)/size(alltrials,2)*100),'%']);
        end;
        
        if p == 4
            colorbar
        end;
    end;
end;
textsc(['Avg Trials Selected by criteria at left. ICA from Emotion data (jo74); ',int2str(size(alltrials,2)),' Tot Trials'],'title');
