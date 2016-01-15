

eeglab
%addpath('/home/julie/MatlabScripts/')
addpath('E:/')

% WT spectral IM clusters:
strs = {
    'load /data/common1/emotion/ThetaClust.mat',
    'load /data/common1/emotion/AlphaClust.mat',
    'load /data/common1/emotion/BetaClust.mat',
    'load /data/common1/emotion/GammaClust.mat'
};

% fullpaths{1} =  '/data/cta/emotion/tl81/';
% fullpaths{2} =  '/data/cta/emotion/mi83/';
% fullpaths{3} =  '/data/cta/emotion/ms82/';
% fullpaths{4} =  '/data/cta/emotion/js75/';
% fullpaths{5} =  '/data/cta/emotion/kw78/';
% fullpaths{6} =  '/data/cta/emotion/jo82/';
% fullpaths{7} =  '/data/cta/emotion/kl80/';
% fullpaths{8} =  '/data/cta/emotion/ar81/';
% fullpaths{9} =  '/data/cta/emotion/eb79/';
% fullpaths{10} = '/data/cta/emotion/dg75/';
% fullpaths{11} = '/data/cta/emotion/an82/';
% fullpaths{12} = '/data/cta/emotion/jw84/';
% fullpaths{13} = '/data/cta/emotion/tv81/';
% fullpaths{14} = '/data/cta/emotion/sr81/';
% fullpaths{15} = '/data/cta/emotion/an70/';
% fullpaths{16} = '/data/cta/emotion/sg75/';
% fullpaths{17} = '/data/cta/emotion/mr72/';
% fullpaths{18} = '/data/cta/emotion/dk74/';
% fullpaths{19} = '/data/cta/emotion/dn86/';
% fullpaths{20} = '/data/cta/emotion/mr71/';
% fullpaths{21} = '/data/cta/emotion/md85/';
% fullpaths{22} = '/data/cta/emotion/mr72-2/';
% fullpaths{23} = '/data/cta/emotion/cj82/';
% fullpaths{24} = '/data/cta/emotion/kc66/';
% fullpaths{25} = '/data/cta/emotion/ts79/';
% fullpaths{26} = '/data/cta/emotion/es76/';
% fullpaths{27} = '/data/cta/emotion/mm78/';
% fullpaths{28} = '/data/cta/emotion/ab75/';
% fullpaths{29} = '/data/cta/emotion/hs83/';
% fullpaths{30} = '/data/cta/emotion/ps82/';
% fullpaths{31} = '/data/cta/emotion/as82/';
% fullpaths{32} = '/data/cta/emotion/ef76/';
% fullpaths{33} = '/data/cta/emotion/jl83/';
% fullpaths{34} = '/data/cta/emotion/rr83/';
% fullpaths{35} = '/data/cta/emotion/jw77/';
% paths{36} = 'hf45/'; % no EXG5 and 6 on mastoid, used old linked mastoids
% paths{37} = 'js78/';
% paths{38} = 'jc66/';
% fullpaths{36} = '/data/common4/emotion/hf45/'; % no EXG5 and 6 on mastoid, used old linked mastoids
% fullpaths{37} = '/data/common4/emotion/js78/';
% fullpaths{38} = '/data/common4/emotion/jc66/';
% oldpaths = fullpaths;

savename = 'Emo-HP';

emos = {'anger','frustration','jealousy','fear' ,'disgust','grief','sad','compassion','love','relief','content','awe','happy','joy','excite'}; % for all new ones
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sadness','  compassion','  love','  relief','  contentment','  awe','  happiness','  joy','  excitement'}; % for plotting purposes
 cols = jet(15);cols(10,:) = [1 .9 0];% for plotting purposes also
 
datset = {'anger.set', 'frustration.set', 'jealousy.set', 'fear.set', 'disgust.set','grief.set','sad.set','compassion.set','love.set','relief.set','content.set','awe.set','happy.set','joy.set','excite.set'}; 

fullpaths{1} =  '/data/common1/emotion/tl81/';
fullpaths{2} =  '/data/common1/emotion/mi83/';
fullpaths{3} =  '/data/common1/emotion/ms82/';
fullpaths{4} =  '/data/common1/emotion/js75/';
fullpaths{5} =  '/data/common1/emotion/kw78/';
fullpaths{6} =  '/data/common1/emotion/jo82/';
fullpaths{7} =  '/data/common1/emotion/kl80/';
fullpaths{8} =  '/data/common1/emotion/ar81/';
fullpaths{9} =  '/data/common1/emotion/eb79/';
fullpaths{10} = '/data/common1/emotion/dg75/';
fullpaths{11} = '/data/common1/emotion/an82/';
fullpaths{12} = '/data/common1/emotion/jw84/';
fullpaths{13} = '/data/common1/emotion/tv81/';
fullpaths{14} = '/data/common1/emotion/sr81/';
fullpaths{15} = '/data/common1/emotion/an70/';
fullpaths{16} = '/data/common1/emotion/sg75/';
fullpaths{17} = '/data/common1/emotion/mr72/';
fullpaths{18} = '/data/common1/emotion/dk74/';
fullpaths{19} = '/data/common1/emotion/dn86/';
fullpaths{20} = '/data/common1/emotion/mr71/';
fullpaths{21} = '/data/common1/emotion/md85/';
fullpaths{22} = '/data/common1/emotion/mr72-2/';
fullpaths{23} = '/data/common1/emotion/cj82/';
fullpaths{24} = '/data/common1/emotion/kc66/';
fullpaths{25} = '/data/common1/emotion/ts79/';
fullpaths{26} = '/data/common1/emotion/es76/';
fullpaths{27} = '/data/common1/emotion/mm78/';
fullpaths{28} = '/data/common1/emotion/ab75/';
fullpaths{29} = '/data/common1/emotion/hs83/';
fullpaths{30} = '/data/common1/emotion/ps82/';
fullpaths{31} = '/data/common1/emotion/as82/';
fullpaths{32} = '/data/common1/emotion/ef76/';
fullpaths{33} = '/data/common1/emotion/jl83/';
fullpaths{34} = '/data/common1/emotion/rr83/';
fullpaths{35} = '/data/common1/emotion/jw77/';
%fullpaths{36} = '/data/common21/emotion/hf45/'; % no EXG5 and 6 on mastoid, used old linked mastoids
%fullpaths{37} = '/data/common1/emotion/js78/';
% newpaths = fullpaths;
% fullpaths = oldpaths;

gdcomps=[];
gdcomps{1} = [1   2   3   4   5   6   7   8  10  11  12  13  14  15  16  18  19  21  22  23  24  27  28  29  30  31  32  33  37];
gdcomps{2} = [ 2     3     4     5   7    10    11    12    13    15    17    18 19    20    24    25];
gdcomps{3} = [3   4   5   6   7   8  10  12  14  15  22  23];
gdcomps{4} = [1   2   3   5   7  15  16  17  18  19  20  21  29];
gdcomps{5} = [1   2   3   5   6   8   9  10  11  13  14  15  18  19  21  23  24  25  26  27  28];
gdcomps{6} = [3    4    6    7    9   10   12   14   16   24   25   27   28];
gdcomps{7} = [2   3   5   6   7   9  15  18  55];
gdcomps{8} = [2   5   7   8  10  14  17  21  24  30];
gdcomps{9} = [2    3    4    5    9   10   12   13   14   15   16   17   18   19   22   24   28   31   33   34   35   36   39   41   48   55];

gdcomps{10} = [3   5   6   7   8   9  11  16  18  19  20  23  25];

gdcomps{11} = [1   3   4   5   6   8   9  10  11  14  15  16  20  22  23  25  29  32  33  48];

gdcomps{12} = [2   4   6   8  10  12  13  16  17  21  23  24  25  26];
gdcomps{13} = [1   2   3   4   5   6   9  10  11  13  14  17  19  23];
gdcomps{14} = [1    2    3    5    6    7    8    9   10   11   12   13   14   15   17   19   20   21   22   26   33   42];
gdcomps{15} = [1   2   3   4   5   6   8  22  28  29  37];

gdcomps{16} = [3   4   6   7   8   9  10  11  12  13  15  16];

gdcomps{17} = [1   3   8   9  10  11  12  13  14  16  17  18  22];
gdcomps{18} = [2   4   5   6   7   8  10  11  13  14  15  17  18  19  20 ];
gdcomps{19} = [5   6   7   8   9  10  11  12  13  14  15  16  18  20  22  23  24  25  26  27  28  29  30  32  33  34  35  38  41  42  43];

gdcomps{20} = [1   2   4   6   8   9  11  12  14  15  16];

gdcomps{21} = [1   5   7   9  11  12  13  14  15  17  18  20  22  24  25  30];
gdcomps{22} = [1    2    5    7    8    9   10   11   13   14   16   17   19   20   21   22   24   26   27   28   29   37   55];
gdcomps{23} = [2    3    4    5    9   10   11   12   13    16   19   21   23   27   28   36   38   39   43   49   51];

gdcomps{24} = [4   5   6   7   8   9  10  11  12  13  15  16  18  19  20  22  23  26];
gdcomps{25} = [9  10  12  14  15  16  18  19  24];

gdcomps{26} = [4   5   7   9  10  12  13  15  16  19  20  31];
gdcomps{27} = [2    4    6    7    8   10   15   20   23];
gdcomps{28} = [2    4    6    8    9   10   11   13   14   15   16   17   18   20   21  26];
gdcomps{29} = [2    3    4    5    6    7    8   10   11   12   13   16   18   19   20   22   24   25   27   29   30   33   35   36   37   39   42   43   44];

gdcomps{30} = [2    3    7    8    9   10   11   12   13   14   15   17   19   21  26];
gdcomps{31} = [3   5   7   8  10  11  14  15  16  17  18];
gdcomps{32} = [5   6   10  12  15  16  22  24  28  30  38];
gdcomps{33} = [4,6,7,9,11,16,17,20,25];
gdcomps{34} = [1   2   3   4   5   6   7   9  10  11  12  14  16  17  18  22  25  26  27  33];
gdcomps{35} = [3   5   6   8   9  12  13  14  16  18  19  20  22  25  27  30]; 



%gdcomps{36} = [5,6,8,9,10,11,12,14,15,16,20,32];% hf45
%gdcomps{36} = [];% take this out for now
%gdcomps{37} = [];% js78

%gdcomps{38} = [4    13    17    20]; % jc66


%blocks{36} = {[0 2100],[2101 4123]}; %hf45
%blocks{37} = {[0 2100],[2101 4123]}; %hf45

% Noisy channels identified by hand:
% $$$ badchan{1}={'D13','EXG4'}; %
% $$$ badchan{2}={'C18','F28'}; %
% $$$ badchan{3} = { 'A1'  'A11'  'A15'  'A17'  'A2'  'A22'  'A23'  'A24'  'A25'  'A29'  'A3'  'A32'  'A8'  'B29'  'C16'  'C17'  'C20'  'C21'  'C23'  'C24'  'C25'  'C26'  'C28'  'C3'  'D32'  'E12'  'E16'  'E17'  'E25'  'E28'  'E29'  'E3'  'E4'  'E5'  'EXG1'  'EXG2'  'EXG4'  'EXG5'  'F2'  'F24'  'F5'  'G12'  'G13'  'G14'  'G15'  'G17'  'G18'  'G21'  'G23'  'G26'  'G31'  'G4'  'G5'  'G6'  'G7'  'G9'  'H18'  'H19'  'H21'  'H23'  'H24','C14','C29','F25','C12','C30','D27','E1' };
% $$$ badchan{5}={'C23'}; %
% $$$ badchan{6}={'D18'};
% $$$ badchan{8}={'A30','D21','E28','H9'};
% $$$ badchan{9}={'C2','D18','EXG2'};
% $$$ badchan{11}={'D28'};                    % 
% $$$ badchan{12}={'A5','A31'}; % 
% $$$ badchan{13}={'C23','G21','G32'}; % 
% $$$ badchan{14}={'B25','C2'}; % 
% $$$ badchan{15}={'B25'}; % 
% $$$ badchan{16}={'A6','C2','D6','D18','D25','G31'}; % 
% $$$ badchan{18}={'C11'}; % 
% $$$ badchan{19}={'D18'}; % 
% $$$ badchan{20}={'A32','C10','E11','E13'}; % 
% $$$ badchan{21}={'A26','B13','B20','B27','D7','D21','G19','G25','H1'}; % 
% $$$ badchan{22}={'C5'}; % 
% $$$ badchan{23} = {'A32' 'D20' 'A14'  'A16'  'A23'  'A26'  'A27'  'A29'  'A31'  'B13'  'B16'  'C21'  'C3'  'C5'  'D14' 'D25'  'D28'  'D29'  'E12'  'E3'  'E5'  'F1'  'F24'  'F25'  'G10'  'G12'  'G13'  'G18'  'G21'  'G23'  'G25'  'G26'  'G27'  'H19'  'H21' };
% $$$ badchan{24} = {'B22','F3' 'A12'  'A19'  'A24'  'A25'  'A30'  'B13'  'B16'  'C1'  'C16'  'C18'  'C7'  'D16'  'D21'  'D23'  'D26'  'D28'  'D8'  'E3'  'E7'  'EXG4'  'EXG5'  'F20'  'F25'  'G10'  'G15'  'G19'  'G23'  'G28'  'G31' 'A31' 'A20' 'F21'};
% $$$ badchan{25}={'A6','A19','A31','G19','G32'}; % 
% $$$ badchan{26}={'A29','A31','F13'}; % 
% $$$ badchan{27}={'C28','G25'}; % 
% $$$ badchan{28}={'C10','D1'}; % 
% $$$ badchan{29}={''}; % 
% $$$ badchan{30} = { 'A12'  'A13'  'A14'  'A20'  'A21'  'A22'  'A27'  'A32'  'A6'  'A7'  'B13'  'B14'  'B17'  'B20'  'B26'  'C16'  'C19'  'C26'  'C7'  'D1'  'D11'  'D23'  'D25'  'D4'  'D5'  'E1'  'E2'  'E3'  'E32'  'EXG1'  'EXG4'  'EXG5'  'F1'  'F11'  'F14'  'F18'  'F31'  'G11'  'G23'  'G5'  'G9'  'H10'  'H24'  'H7'  'H9' 'C5'};
% $$$ badchan{31}={'A31','G21'}; % 
% $$$ badchan{32}={'B9','B27'}; % 
% $$$ badchan{33} = { 'A1'  'A18'  'A31'  'A6'  'A7'  'C18'  'C8'  'D1'  'D13'  'D17'  'D18'  'D19'  'D24'  'D25'  'E19'  'E20'  'E23'  'E24'  'E26'  'E29'  'E3'  'E30'  'E5'  'E7'  'EXG1'  'EXG4'  'EXG5'  'F23'  'F24'  'F29'  'F30'  'F32'  'G10'  'G11'  'G13'  'G23'  'G27'  'H11'  'H16'  'H18'  'H19'  'H20'  'H21'  'H23'  'H24' 'F4' 'F22' 'F28' 'C16' 'F15' 'H9' 'H7' 'B19'};% uh, E3 is bad, meaning don't use it as a reference! (reimported)
% $$$ 
% $$$ badchan{34}={ 'A11'    'A21'    'A30'    'A32'    'A8'    'B14'    'C11'    'C16'  'D13'    'D16'    'D18'    'D19'    'D2'    'D21'    'D25'    'D28'    'D3'    'D31'    'D4'    'D9'    'E1'    'E13'    'E15'    'E17'    'E21'    'E28'    'E29'    'E4'    'E5'    'E6'    'F19'    'F20'    'F22'    'F27'    'G1'    'G13'    'G16'    'G20'    'G27'    'G32'    'G5'    'H11'    'H13'    'H18'    'H20'    'H21'    'H24'    'H5'    'H6'    'H8'};
% $$$ 
% $$$ badchan{35}={'A32','C10','D20','G29'}; % 
% $$$ badchan{36}={'A1','A13','A19','A30','A31','A32','C10','C11','C16','C21','D3','D6','D14','D22','E3','E4','E5','E16','E32','F4','F9','F27','G11','G12','G24','G27'}; % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combination of by hand and automated rejection:
badchan{15} = { 'A26'  'A27'  'A31'  'B25'  'B27'  'B30'  'C10'  'C27'  'C30'  'C31'  'C32'  'C5'  'D19'  'E12'  'E2'  'E27'  'E3'  'E9'  'G12'  'G15'  'G2'  'G21'  'G23'  'G27'  'G29'  'G30'  'G4'  'H12'  'H19'  'H20'  'H22' };
badchan{11} = { 'B21'  'B25'  'D14'  'D20'  'D27'  'D28'  'D29'  'E15'  'E2'  'E20'  'E3'  'G1'  'G10'  'G23'  'H12'  'H16'  'H18'  'H19'  'H20'  'H21'  'H22'  'H23'  'H24' };

badchan{6} = { 'A11'  'A18'  'B27'  'B32'  'D1'  'D14'  'D15'  'D9'  'F9'  'G25' };
badchan{22} = {};
badchan{23} = {};
badchan{7} = { 'A11'  'A18'  'B27'  'B32'  'D1'  'D14'  'D15'  'D9'  'F9'  'G25' };
badchan{24} = {};
badchan{2} = { 'A24'  'A28'  'A9'  'B15'  'D11'  'D20'  'E3'  'E8'  'E9'  'F20'  'G10'  'G15'  'G17'  'G19'  'G23'  'G26'  'G28'  'G29'  'G30'  'G32'  'G4'  'G7'  'H20'  'H24'  'H6' };
badchan{4} = { 'A21'  'A31'  'B1'  'B14'  'B25'  'C13'  'C19'  'C26'  'D20'  'D4'  'D8'  'E16'  'E19'  'E23'  'E24'  'E3'  'E9'  'EXG1'  'EXG4'  'F24'  'G1'  'G10'  'G2'  'G23'  'G3'  'G4'  'G5'  'G6'  'G8'  'G9'  'H2'  'H24' };
badchan{7} = { 'A14'  'A18'  'A19'  'A21'  'A25'  'A26'  'A30'  'C22'  'C5'  'D14'  'D22'  'E3'  'E30'  'F25'  'G14'  'G23'  'G24'  'G26'  'G27' };
badchan{21} = { 'A1'  'A21'  'A29'  'B15'  'B5'  'C18'  'C23'  'C24'  'C28'  'C30'  'D13'  'D16'  'D20'  'D29'  'D8'  'E17'  'E22'  'E3'  'E30'  'E5'  'EXG4'  'F11'  'F19'  'F3'  'G13'  'G14'  'G18'  'G22'  'G23'  'H23'  'H24' };
badchan{14} = { 'A2'  'A24'  'A3'  'A30'  'A8'  'B17'  'C1'  'C22'  'C27'  'C5'  'D1'  'D13'  'E3'  'E9'  'F12'  'F28'  'G1'  'G13'  'G14'  'G15'  'G16'  'G19'  'G2'  'G23'  'G5'  'H13'  'H14'  'H17'  'H19'  'H2'  'H21'  'H4'  'H6' };
badchan{10} = { 'A1'  'A16'  'A17'  'A20'  'A24'  'A30'  'A8'  'A9'  'C11'  'C12'  'C25'  'C28'  'D17'  'D18'  'D23'  'D24'  'D25'  'D26'  'D28'  'D29'  'E10'  'E12'  'E16'  'E2'  'E20'  'E23'  'E24'  'E3'  'E5'  'E9'  'EXG2'  'EXG3'  'G11'  'G18'  'G2'  'G22'  'G23'  'G24'  'G25'  'G28'  'G32'  'G5'  'G9'  'H21'  'H23' };
badchan{32} = { 'A1'  'A13'  'A14'  'A2'  'A26'  'A5'  'A9'  'B26'  'C10'  'D11'  'D12'  'D18'  'D25'  'E11'  'E12'  'E13'  'E3'  'E7'  'E8'  'E9'  'G15'  'G23'  'G32'  'G8'  'G9'  'H13'  'H18'  'H19'  'H20' };
badchan{1} = { 'A24'  'A28'  'A30'  'A9'  'B15'  'D11'  'D20'  'E3'  'E8'  'E9'  'EXG5'  'F20'  'G10'  'G15'  'G16'  'G17'  'G19'  'G23'  'G26'  'G28'  'G29'  'G30'  'G32'  'G4'  'G7'  'H20'  'H24'  'H6' };
badchan{26} = { 'A11'  'A16'  'A17'  'A18'  'A20'  'A21'  'A23'  'A25'  'A9'  'B10'  'B26'  'B28'  'C1'  'C17'  'C23'  'C27'  'C29'  'C5'  'D13'  'D14'  'D19'  'D23'  'D24'  'D25'  'E15'  'E3'  'E31'  'F27'  'F5'  'G18'  'G19'  'G20'  'G23'  'G26'  'G27'  'G28'  'G29'  'G31'  'G32'  'H2' };
badchan{29} = { 'A31'  'B26'  'C26'  'C7'  'D16'  'D18'  'D2'  'D7'  'E26'  'E27'  'E29'  'E3'  'F12'  'F24'  'G10'  'G11'  'G23'  'G30'  'H2'  'H23'  'H24'  'H8' };
badchan{22} = { 'A20'  'A27'  'B13'  'C21'  'C27'  'C29'  'D15'  'D19'  'D22'  'D26'  'D27'  'D28'  'D6'  'E12'  'E27'  'E3'  'E7'  'EXG4'  'F14'  'F22'  'F23'  'F31'  'G17'  'G21'  'G22'  'G23'  'G24'  'G25'  'G26'  'G30'  'G4'  'G6'  'G7'  'G9'  'H9' };
badchan{12} = { 'A24'  'A29'  'A30'  'A32'  'C1'  'C12'  'C2'  'C26'  'C27'  'D25'  'D28'  'E3'  'EXG2'  'EXG3'  'EXG6'  'G11'  'G13'  'G15'  'G19'  'G23'  'G28'  'G29'  'G30'  'G31'  'G32'  'H16'  'H17'  'H18'  'H20'  'H22'  'H23'  'H24' };
badchan{2} = { 'A20'  'A21'  'A28'  'A30'  'B1'  'C17'  'C2'  'C22'  'C23'  'C26'  'C28'  'C30'  'C5'  'D13'  'D15'  'D19'  'D21'  'D5'  'D6'  'D8'  'E24'  'E27'  'E3'  'E30'  'E5'  'E9'  'F13'  'G11'  'G13'  'G14'  'G16'  'G17'  'G22'  'G23'  'G28'  'G30'  'G32'  'G8' };
badchan{8} = { 'A13'  'A21'  'A24'  'A31'  'C15'  'C17'  'C21'  'C3'  'D14'  'D17'  'D18'  'D26'  'E2'  'E24'  'E29'  'E3'  'E30'  'E5'  'EXG4'  'F12'  'F28'  'G11'  'G13'  'G18'  'G21'  'G23'  'G24'  'G5'  'G6' };
badchan{20} = { 'A1'  'A11'  'A19'  'A2'  'B10'  'B16'  'B4'  'C27'  'C29'  'C31'  'D11'  'D15'  'D22'  'D24'  'D28'  'D30'  'E10'  'E12'  'E24'  'E25'  'E3'  'E30'  'EXG3'  'G12'  'G13'  'G20'  'G23'  'G29'  'G3'  'G7'  'H14'  'H15' };
badchan{13} = { 'A1'  'A13'  'A21'  'A25'  'A30'  'C11'  'C12'  'C17'  'C4'  'C5'  'D25'  'E10'  'E15'  'E20'  'E25'  'E28'  'E29'  'E3'  'E30'  'E8'  'F12'  'G1'  'G11'  'G13'  'G16'  'G17'  'G18'  'G23'  'G26'  'G27'  'G30'  'G9'  'H1'  'H18'  'H21'  'H24' };
badchan{5} = { 'A25'  'A30'  'C12'  'C17'  'C24'  'C25'  'C6'  'C7'  'D15'  'D20'  'D21'  'D26'  'D6'  'D7'  'D9'  'E24'  'E27'  'E3'  'E5'  'F27'  'F28'  'F5'  'G15'  'G21'  'G23'  'G25'  'G26'  'G5'  'H12'  'H18'  'H20'  'H21'  'H23'  'H24' };
badchan{3} = { 'A1'  'A11'  'A15'  'A17'  'A2'  'A22'  'A23'  'A24'  'A25'  'A29'  'A3'  'A32'  'A8'  'B29'  'C16'  'C17'  'C20'  'C21'  'C23'  'C24'  'C25'  'C26'  'C28'  'C3'  'D32'  'E12'  'E16'  'E17'  'E25'  'E28'  'E29'  'E3'  'E4'  'E5'  'EXG1'  'EXG2'  'EXG4'  'EXG5'  'F2'  'F24'  'F5'  'G12'  'G13'  'G14'  'G15'  'G17'  'G18'  'G21'  'G23'  'G26'  'G31'  'G4'  'G5'  'G6'  'G7'  'G9'  'H18'  'H19'  'H21'  'H23'  'H24' };
badchan{9} = { 'A1'  'A16'  'A17'  'A32'  'A8'  'C11'  'C27'  'D10'  'D11'  'D13'  'D16'  'D20'  'D22'  'D24'  'D7'  'D9'  'E3'  'E4'  'E9'  'EXG1'  'G13'  'G23'  'G32'  'G5' };
badchan{16} = { 'A12'  'A13'  'A20'  'A21'  'A25'  'A28'  'A30'  'B25'  'B28'  'C10'  'C11'  'C17'  'C23'  'C24'  'C27'  'C28'  'D14'  'D27'  'D9'  'E13'  'E3'  'E30'  'EXG1'  'EXG4'  'EXG6'  'F2'  'F25'  'G15'  'G20'  'G21'  'G23'  'G26'  'G27'  'G4'  'G9'  'H24' };
badchan{25} = { 'A16'  'A17'  'A32'  'A8'  'B2'  'C10'  'C11'  'C17'  'C18'  'C2'  'C22'  'C26'  'C27'  'C30'  'C31'  'D13'  'D18'  'D25'  'D28'  'E10'  'E11'  'E12'  'E13'  'E14'  'E28'  'E3'  'E30'  'E7'  'E8'  'E9'  'G1'  'G12'  'G14'  'G17'  'G20'  'G23'  'G27'  'H2'  'H20'  'H22'  'H23'  'H24' };
badchan{31} = { 'A1'  'A16'  'A17'  'A2'  'A20'  'A22'  'A23'  'A8'  'A9'  'C13'  'C15'  'D16'  'D18'  'D19'  'D21'  'D24'  'D25'  'D27'  'D29'  'D9'  'E15'  'E20'  'E21'  'E22'  'E26'  'E29'  'E3'  'E5'  'E7'  'E8'  'F31'  'G17'  'G23'  'G7'  'H19'  'H20'  'H21'  'H23'  'H24' };
badchan{6} = { 'A11'  'A18'  'B27'  'B32'  'D1'  'D14'  'D15'  'D9'  'F9'  'G25' };
badchan{17} = { 'A20'  'A21'  'A24'  'A25'  'A30'  'A31'  'B7'  'C1'  'C2'  'C22'  'C23'  'C27'  'D13'  'D14'  'D20'  'D26'  'D5'  'D7'  'E21'  'E22'  'E23'  'E24'  'E26'  'E27'  'E3'  'E30'  'EXG1'  'F13'  'F24'  'G10'  'G11'  'G12'  'G13'  'G14'  'G17'  'G21'  'G22'  'G23'  'G24'  'G28'  'G5'  'G9' };
badchan{27} = { 'A12'  'A13'  'A14'  'A24'  'A28'  'A30'  'A7'  'B19'  'D11'  'D30'  'D6'  'D7'  'E25'  'E26'  'E27'  'E29'  'E3'  'E30'  'E4'  'G13'  'G21'  'G22'  'G23'  'G31'  'G5'  'G9' };
badchan{18} = { 'A13'  'A14'  'A22'  'A27'  'B23'  'C1'  'C10'  'C2'  'C26'  'C30'  'D19'  'D24'  'D28'  'D29'  'E10'  'E3'  'E6'  'E7'  'EXG6'  'G1'  'G10'  'G13'  'G14'  'G23'  'G28'  'G3'  'G30'  'G32'  'G9' };
badchan{19} = { 'A1'  'A20'  'C1'  'C11'  'C12'  'C26'  'C27'  'D16'  'D26'  'D28'  'D32'  'D8'  'E1'  'E13'  'E14'  'E15'  'E3'  'E5'  'EXG1'  'F28'  'F4'  'G18'  'G20'  'G23'  'G24'  'G25'  'G26'  'G28'  'G9' };
badchan{28} = { 'A18'  'A32'  'B10'  'B13'  'B14'  'B19'  'B25'  'C12'  'C29'  'C30'  'C5'  'D18'  'D25'  'D30'  'E16'  'E23'  'E24'  'E27'  'E3'  'E30'  'E5'  'E7'  'E8'  'E9'  'F14'  'F31'  'G13'  'G17'  'G23'  'G9'  'H23'  'H24'  'H9' };
badchan{30} = { 'A19'  'A24'  'A26'  'A30'  'C16'  'C19'  'C7'  'D11'  'D23'  'D4'  'D5'  'E1'  'E2'  'E23'  'E24'  'E26'  'E27'  'F14'  'F18'  'F23'  'F24'  'G6'  'H10'  'H7' };
badchan{35} = { 'A1'  'A15'  'A31'  'A32'  'A8'  'B20'  'B31'  'C1'  'C10'  'C5'  'D18'  'D20'  'D25'  'D28'  'D29'  'D7'  'E10'  'E12'  'E18'  'E29'  'E3'  'E32'  'EXG1'  'EXG3'  'EXG4'  'EXG5'  'F1'  'F15'  'F24'  'F28'  'F4'  'G1'  'G14'  'G16'  'G20'  'G23'  'G29'  'G31'  'G5'  'G6'  'G7'  'G8'  'H19'  'H21'  'H23'  'H24' };

badchan{23} = { 'A13'  'A21'  'D18'  'D24'  'D25'  'E4'  'E7'  'F23'  'G15'  'G19'  'G22'  'G24'  'G31'  'H24' };

badchan{24} = { 'A1'  'A11'  'A13'  'A14'  'A20'  'A21'  'A26'  'A31'  'C10'  'C13'  'C30'  'D12'  'D14'  'D22'  'D27'  'D9'  'E15'  'E16'  'EXG6'  'F1'  'F18'  'F21'  'F23'  'F24'  'G13'  'G21'  'H21' };

badchan{3} = { 'A10'  'C1'  'C10'  'C11'  'C2'  'D23'  'D24'  'D28'  'D29'  'E11'  'E13'  'E14'  'E15'  'E21'  'E22'  'E23'  'E24'  'E26'  'E27'  'E30'  'E6'  'E9'  'EXG3'  'G10'  'G11'  'G16'  'G19'  'G29'  'G32'  'G8' };

badchan{34} = { 'A11'  'A12'  'A21'  'A26'  'A27'  'A30'  'A31'  'A32'  'A8'  'B14'  'C10'  'C11'  'C16'  'C30'  'C9'  'D13'  'D16'  'D18'  'D19'  'D2'  'D20'  'D21'  'D25'  'D28'  'D29'  'D3'  'D30'  'D31'  'D4'  'D9'  'E1'  'E13'  'E15'  'E17'  'E20'  'E21'  'E22'  'E28'  'E29'  'E3'  'E30'  'E32'  'E4'  'E5'  'E6'  'EXG3'  'EXG4'  'EXG5'  'F19'  'F20'  'F22'  'F24'  'F27'  'G1'  'G13'  'G14'  'G16'  'G20'  'G23'  'G27'  'G32'  'G5'  'G6'  'H11'  'H13'  'H18'  'H2'  'H20'  'H21'  'H24'  'H4'  'H5'  'H6'  'H8' };
badchan{30} = { 'A12'  'A13'  'A14'  'A19'  'A20'  'A21'  'A22'  'A24'  'A26'  'A27'  'A30'  'A32'  'A6'  'A7'  'B13'  'B14'  'B17'  'B20'  'B26'  'C16'  'C19'  'C26'  'C5'  'C7'  'D1'  'D11'  'D23'  'D25'  'D4'  'D5'  'E1'  'E14'  'E2'  'E3'  'E32'  'EXG1'  'EXG4'  'EXG5'  'F1'  'F11'  'F14'  'F18'  'F23'  'F24'  'F31'  'G11'  'G23'  'G5'  'G9'  'H10'  'H24'  'H7'  'H9' };
badchan{33} = { 'A1'  'A18'  'A2'  'A31'  'A6'  'A7'  'B19'  'C16'  'C18'  'C8'  'D1'  'D13'  'D17'  'D18'  'D19'  'D20'  'D24'  'D25'  'E19'  'E20'  'E23'  'E24'  'E25'  'E26'  'E27'  'E29'  'E3'  'E30'  'E4'  'E5'  'E7'  'EXG1'  'EXG4'  'EXG5'  'F1'  'F15'  'F19'  'F2'  'F22'  'F23'  'F24'  'F28'  'F29'  'F30'  'F32'  'F4'  'G1'  'G10'  'G11'  'G13'  'G23'  'G27'  'G5'  'G9'  'H11'  'H16'  'H17'  'H18'  'H19'  'H20'  'H21'  'H23'  'H24'  'H7'  'H9' };
badchan{37} = { 'A1'  'A11'  'A12'  'A13'  'A8'  'A9'  'D11'  'D12'  'D18'  'D24'  'D25'  'D29'  'D32'  'E3'  'E7'  'EXG7'  'EXG8'  'G20'  'G22'  'G23'  'H18'  'H21'  'H24'  'H7' };
badchan{38} = { 'A1'  'A12'  'A15'  'A16'  'A20'  'A31'  'A32'  'A6'  'B1'  'B2'  'B25'  'B3'  'C1'  'C13'  'C19'  'C2'  'C20'  'C31'  'C6'  'C8'  'D10'  'D12'  'D15'  'D17'  'D2'  'D23'  'D27'  'D3'  'D32'  'D4'  'D6'  'D7'  'E16'  'E22'  'E23'  'E26'  'E27'  'E3'  'E30'  'E5'  'EXG5'  'EXG6'  'F32'  'F4'  'G1'  'G10'  'G14'  'G17'  'G20'  'G21'  'G23'  'G3'  'G32'  'G5'  'G6'  'G7'  'G9'  'H12'  'H17'  'H21' };

% Number of channels in datasets

nchans(15) = 229;
nchans(11) = 231;
nchans(4) = 222;
nchans(7) = 235;
nchans(21) = 214;
nchans(14) = 219;
nchans(10) = 209;
nchans(32) = 223;
nchans(1) = 224;
nchans(26) = 211;
nchans(29) = 232;
nchans(22) = 218;
nchans(12) = 220;
nchans(2) = 214;
nchans(8) = 221;
nchans(20) = 218;
nchans(13) = 215;
nchans(5) = 219;
nchans(24) = 223;
nchans(3) = 183;
nchans(9) = 227;
nchans(16) = 212;
nchans(25) = 207;
nchans(31) = 213;
nchans(6) = 221;
nchans(17) = 212;
nchans(27) = 226;
nchans(18) = 224;
nchans(19) = 224;
nchans(28) = 219;
nchans(35) = 208;
nchans(23) = 206;
nchans(24) = 196;
nchans(3) = 134;
nchans(34) = 180;
nchans(30) = 201;
nchans(33) = 189;

nchans(36) = 144;
nchans(37) = 215;
nchans(38) = 194;

% Number of datasets with continuous data
nsets(1) = 5;
nsets(2) = 6;
nsets(3) = 5;% was 6??
nsets(4) = 6;
nsets(5) = 5;
nsets(6) = 5;
nsets(7) = 4;
nsets(8) = 5;
nsets(9) = 9;
nsets(10) = 7;
nsets(11) = 6;
nsets(12) = 5;
nsets(13) = 5;
nsets(14) = 5;
nsets(15) = 5;
nsets(16) = 5;
nsets(17) = 5;
nsets(18) = 7;
nsets(19) = 4;
nsets(20) = 5;
nsets(21) = 5;
nsets(22) = 5;
nsets(23) = 9;
nsets(24) = 6;
nsets(25) = 6;
nsets(26) = 6;
nsets(27) = 5;
nsets(28) = 5;
nsets(29) = 6;
nsets(30) = 7;
nsets(31) = 7;
nsets(32) = 5;
nsets(33) = 7;
nsets(34) = 5;
nsets(35) = 6;
nsets(36) = 5;
nsets(37) = 3;
nsets(38) = 5;
% made by (manually input subjs 3(6),22(5),23(9)
% $$$ for nx = 1:35
% $$$     x=openbdf([fullpaths{nx},fullpaths{nx}(end-4:end-1),'.bdf']);
% $$$     numframes = x.Head.NRec;
% $$$     nsets(nx) = ceil(numframes/1000);
% $$$     fid = fopen([fullpaths{1}(1:end-5),'Nsets'],'a');
% $$$     fprintf(fid, '\nnsets(%s) = %s;', int2str(nx),int2str(nsets(nx)));
% $$$     fclose(fid);
% $$$ end;


