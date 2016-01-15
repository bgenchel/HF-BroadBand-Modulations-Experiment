% make emotion IM movies
eeglab
datpath = '/data/common1/emotion/mi83/';
datset = 'sources.set';
savedat = 'SpecCoModMoreFreqs';

% Plot one IM, sorted, for demonstration:
im = [1,4];
comp = 5; % mid occ
frqlim = [3 125];
moviename = 'GammaSorted.avi';
moviename = 'AGdiffEmos.avi';
trialset = [3,4,10,12]; % jealousy
weights = []; % no normalization
sorttrls = 'on';
setcols = []; % make a rainbow
ttls = {'jealousy','fear','relief','awe'}; % no titles
ttls = []; % no titles
setcols = [];
[M] = MakeIMmovie(datset,savedat,datpath,im,comp,frqlim,moviename,trialset,weights,sorttrls,setcols,ttls);
s = load([datpath,savedat,'.mat']); 
SpecCoModPlot('sources.set',datpath,[],[1:15],savedat,[3 125],'n',0,[]);

% subj mi83, IMs 1 and 4 were in movies for talks
% Scatter plots:

