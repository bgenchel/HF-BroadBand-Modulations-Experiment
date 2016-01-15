% find accuracy for each emotion by comparison with known labels

function [corrperc] = FindEmoAccuracy(AllClassLabels,AllTestLabels);


corrperc = zeros(size(AllClassLabels,2),size(AllClassLabels,3));
labeltypes = unique(AllTestLabels(:,1,1))';

for t = 1:size(AllClassLabels,3) % for each iteration
   for f = 1:size(AllClassLabels,2) % for each # features (IMs)
      for lb = 1:length(labeltypes) % for all test types
         onetype = find(AllTestLabels(:,f,t) == labeltypes(lb));
         corrmatch = length(find(ismember(AllClassLabels(onetype,f,t),labeltypes(lb))));
         corrperc(f,t,lb) = 100*(corrmatch/length(onetype)); % percent correct
      end;
   end;
end;
 