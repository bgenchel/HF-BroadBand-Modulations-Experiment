clear,close all

InputPath=[pwd '\Individual Results\'];

%% Load Add-one-feature-in results of IM-, IC- or Channel-based data
InputFileName='YuanPinEmoWeights1sec';                                   % IM variable
% InputFileName='YuanPinICSpecPwr1Sec';                                  % IC variable
% InputFileName='YuanPinChanSpecPwr1Sec';                             % Channel variable

%% Load Feature selection results  (Fscore and Random feature selection for add-one-feature-in scheme)
load([InputPath 'FeatSelRawData_' InputFileName '.mat']);      
% FeatSel{NumOfSubjects,NumOfEmoCate}.Strategy{NumOfStrategy}.SortedIndex(NumOfFeature,1)
% NumOfSubjects: 35
% NumOfEmoCate: 3 (1: All, 2: Valence, 3: Arousal)
% NumOfStrategy: 2 (1: Fscore, 2: Random)
% NumOfFeature: dependents on subject

%% 
SubjectStartIndex=1;
SubjectEndIndex=8;
EmCateIndex=1;                                                                                       % 1: All, 2: Valence, 3: Arousal
StrategyIndex=1;                                                                                     % 1: Fscore, 2: Random

NumOfSubject=SubjectEndIndex-SubjectStartIndex+1;
NameOfStrategy{1}='Fscore';
NameOfStrategy{2}='Random';
NameOfEmoCate{1}='All';
NameOfEmoCate{2}='Val';
NameOfEmoCate{3}='Aro';

for SubjectIndex=SubjectStartIndex:SubjectEndIndex
 try
    % load individual results. Ex: YuanPinEmoWeights1sec_SingleSub_S1_All_Fscore
    load([InputPath InputFileName '_SingleSub_S' num2str(SubjectIndex) '_' NameOfEmoCate{EmCateIndex} '_' NameOfStrategy{StrategyIndex} '.mat']);
    % SingleSub.Accuracy.Mean: (NumOfFeat,1)
    % SingleSub.Accuracy.PerIteration: (NumOfFeat,NumOfAvg,NumOfFold)
    % SingleSub.Accuracy.PerIterationStr='(NumOfFeature, NumOfAvg,NumOfFold)';
    % SingleSub.ConfusionTablePerIteration: {NumOfFeat,1}
    
end
end