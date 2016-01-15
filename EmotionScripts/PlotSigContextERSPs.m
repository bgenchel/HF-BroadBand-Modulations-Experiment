% plots specified trial ERSPs (from PCA'd back-proj) from PlotContexSort
% output
% 
%PlotSigContextERSPs(name,name2,fullpath,hicontigs,lowcontigs,questions);
%
% name -- [string] stem+ic number of decomposition to load
% name2 -- [string] only used if loading a 'context-only' decomposition;
% This file should have ERSP data from a decomp that is otherwise identical
% to 'name'
% fullpath - [string] data  directory where name and name2 are
% hicontigs,lowcontigs: output from PlotContextSort (cell array with a cell
% for each question. Within each question are cell array(s) with contiguous
% trials that were identified by this decomposition. 
% questions -- [vector] index of context questions for which to plot mean
% trial ERSPs that were isolated by the specified decomposition. If [],
% will go through all questions and plot any selected trials. 
%
% Currently this function is only for plotting a single decomp dimension.
%


function [allersps] = PlotSigContextERSPs(name,name2,fullpath,hicontigs,lowcontigs,questions);

s = load([fullpath,name,'.mat']);
wts = floatread([fullpath,name,'.wts'],[s.numrows s.numrows],[],0);
sph = floatread([fullpath,name,'.sph'],[s.numrows s.numrows],[],0);
dat = floatread([fullpath,name,'.fdt'],[s.numrows s.numframes],[],0);
ws = wts*sph;winv = pinv(ws);acts = ws*dat;
if ~isempty(s.eigfile)
    if s.pcmat > 0
        if s.pcctx == s.pctot % context only
            ss = load([fullpath,name2,'.mat']); % hijacking for context only decomp
            erspeig = floatread([fullpath,name2,'EIGVEC.fdt'],[length(ss.freqs)*length(ss.times) ss.pcmat],[],0);
            dat = floatread([fullpath,name2,'.fdt'],[ss.numrows ss.numframes],[],0);
            erspdat = erspeig*dat(1:ss.pcmat,:);% back-proj to orig data
        else
            erspeig = floatread([fullpath,name,'EIGVEC.fdt'],[length(s.freqs)*length(s.times) s.pcmat],[],0);
            erspdat = erspeig*dat(1:s.pcmat,:);% back-proj to orig data
        end;
    end;
end;
if isempty(questions)
    questions = [1:length(s.cxtmean)];
end;

for dim = 1:length(hicontigs)
    if ~isempty(hicontigs{dim})
        hiersps = zeros(length(s.freqs),length(s.times),0);
        figure; pl=1;
        row = 4; col = 4;
        for q = 1:length(questions)
            if ~isempty(hicontigs{dim}{questions(q)})
                for h=1:length(hicontigs{dim}{questions(q)})
                    if pl>row*col
                        textsc([name,', Dim ',int2str(dim),'; Mean ERSPs for ''yes'' trials'],'title');
                        figure; pl=1;
                    end;
                    oneersp = erspdat(:,hicontigs{dim}{questions(q)}{h});
                    ntrials = size(oneersp,2);
                    oneersp = reshape(oneersp,length(s.freqs),length(s.times),size(oneersp,2));
                    oneersp = mean(oneersp,3);lim = max(abs(oneersp(:)));
                    hiersps(:,:,end+1) = oneersp;
                    sbplot(row,col,pl);pl = pl+1;
                    quadimagesc(s.times,s.freqs,oneersp,lim); hold on;
                    title(['Q ',int2str(questions(q)),'; Run ',int2str(h),'; ',int2str(ntrials),' trials']);
                    plot([0 0],[get(gca,'ylim')],'k-');cbar;
                end;
            end;
        end;
        allersps{1}{dim} = hiersps;
        textsc([name,', Dim ',int2str(dim),'; Mean ERSPs for ''yes'' trials'],'title');
    end;
end;



for dim = 1:length(lowcontigs)
    if ~isempty(lowcontigs{dim})
        lowersps = zeros(length(s.freqs),length(s.times),0);
        figure; pl=1;
        row = 4; col = 4;
        for q = 1:length(questions)
            if ~isempty(lowcontigs{dim}{questions(q)})
                for h=1:length(lowcontigs{dim}{questions(q)})
                    if pl>row*col
                        textsc([name,', Dim ',int2str(dim),'; Mean ERSPs for ''no'' trials'],'title');
                        figure; pl=1;
                    end;
                    oneersp = erspdat(:,lowcontigs{dim}{questions(q)}{h});
                    ntrials = size(oneersp,2);
                    oneersp = reshape(oneersp,length(s.freqs),length(s.times),size(oneersp,2));
                    oneersp = mean(oneersp,3);lim = max(abs(oneersp(:)));
                    lowersps(:,:,end+1) = oneersp;
                    sbplot(row,col,pl);pl = pl+1;
                    quadimagesc(s.times,s.freqs,oneersp,lim); hold on;
                    title(['Q ',int2str(questions(q)),'; Run ',int2str(h),'; ',int2str(ntrials),' trials']);
                    plot([0 0],[get(gca,'ylim')],'k-');cbar;
                end;
            end;
        end;
        allersps{2}{dim} = lowersps;
        textsc([name,', Dim ',int2str(dim),'; Mean ERSPs for ''no'' trials'],'title');
    end;
end;


