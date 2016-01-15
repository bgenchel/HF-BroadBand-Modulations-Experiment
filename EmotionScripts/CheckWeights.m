% go through specified datasets and return which datasets are different
% takes first dataset listed as the 'correct' weights
%
% wtcorrs: row 1: weight correlation, row 2: sphere correlation
%          NaN is if EEG.ica* is empty, 0 is when size is different
% problem is 1 if correlation is ever < 1

function [wtcorrs,filenames,problem] = CheckWeights(filenames,fullpath);

if isempty(filenames)
  % then take all EEGLAB datasets in the directory
  [stat,ff] = system(['cd ',fullpath,'; ls *.set']);
  filenames = textscan(ff,'%s','collectoutput',1);
  filenames = filenames{1};
end;
% load first datraset to get weights
EEG = pop_loadset('filename',[filenames{1}(1:end-4),'.set'],'filepath',fullpath);
wts = EEG.icaweights;
sph = EEG.icasphere;
wtcorrs(1:2,1) = 1; % by definition
for f = 2:length(filenames)
  ALLEEG=[];EEG=[]; 
  EEG = pop_loadset('filename',[filenames{f}(1:end-4),'.set'],'filepath',fullpath);
  if isempty(EEG.icaweights)
    wtcorrs(1:2,f) = NaN;
  elseif size(EEG.icaweights,1) ~= size(wts,1)
    wtcorrs(1:2,f) = 0; % worse than NaN, implies different decomp
  else
    [wtcorrs(1,f)] = corr2(wts,EEG.icaweights);
    [wtcorrs(2,f)] = corr2(sph,EEG.icasphere);
  end;
end;
if ~isempty(find(wtcorrs<1))
  problem = 1;
else
  problem = 0;
end;