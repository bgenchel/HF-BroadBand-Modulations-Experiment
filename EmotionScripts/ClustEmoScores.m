% Spin off from PCAforEmos.m to cluster Fact*subj X Trials matrix

button = [1:12,21:26]; % all button presses, early and 'only when you feel it' subjects
button = [13:20]; % no button press (apart from the first one)
button = [1,2,4:6,8:12,14,17:21,23,25,26,27]; % all 'good' subjects (ones that said they got into it)
button = [1:8,10:21,23:27];  % not mr72-2

pcamatall = zeros(15,0);clear mnemodiff forstats
alldat = zeros(45,0);moredat = zeros(0,300);
for nx = 1:length(button)
    sph=floatread(['/data/common2/emotion/clusters/',sphfile{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(1)]); 
    wts=floatread(['/data/common2/emotion/clusters/',wtsfile{button(nx)}],[pcadims{button(nx)} subjdims{button(nx)}(1)]); 
    if length(gdcomps{1}) < 20
    icamatall = floatread(['/data/common2/emotion/clusters/',Frontsubjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    else
    icamatall = floatread(['/data/common2/emotion/clusters/',subjspecs{button(nx)}],[subjdims{button(nx)}(1) subjdims{button(nx)}(2)]);
    end;
    ws = wts*sph;    activations = ws*icamatall;    winv = pinv(ws);    emomap = ones(1,1);
    for e = 2:length(numtrials{button(nx)})+1
        emomap(1,e) = emomap(e-1) + numtrials{button(nx)}(e-1); % marks where each emotion STARTS
    end; 
    %alltps = zeros(45,0);
    percscores = zeros(300,0);
    for tp = 1:size(winv,2)
        tpwts =  winv(:,tp)'; clear allemos
        allemos = zeros(0,1);pscores = zeros(0,1);
        for e = 1:length(numtrials{button(nx)})  % start with 2 for straight nums (not diffs)
            clear newmat newmat2
            tempmat = tpwts(emomap(e):emomap(e+1)-1); tempmat = sort(tempmat);
            %pl=1;      % for inputing more than just mean           
            % for pk = .25:.25:.75
            %    newmat(1,pl) = tempmat(round(size(tempmat,2)*pk));pl = pl+1;
            %end;
            %allemos(end+1:end+3,1) = newmat';
            pl=1;      % for inputing more than just mean           
             for pk = .01:.05:.99
                newmat2(1,pl) = tempmat(ceil(size(tempmat,2)*pk));pl = pl+1;
            end;
            pscores(end+1:end+size(newmat2,2),1) = newmat2';
        %mnemo(e,tp) = mean(tempmat);%instead of diff
        end; 
        %alltps(:,end+1) = allemos;
         percscores(:,end+1) = pscores;
    end;
    moredat(end+1:end+size(percscores,2),:) = percscores';
    %alldat(:,end+1:end+size(alltps,2)) = alltps;
    %pcamatall(:,end+1:end+size(mnemo,2)) = mnemo; 
    fprintf('\n One More SUBJECT Done: %i',button(nx));
end;
 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(moredat,'pca',6,'extended',1,'stop',1e-7,'maxsteps',2000);
ws = weights*sphere; winv = pinv(ws);

% each bin of 20 is one emotion from one subject 
for dim = 1:size(activations,1)
for e = 1:length(emo2)
    ewts(e,:,dim) = activations(dim,(e-1)*20+1:e*20);
end;
end;
 figure; 
 for x=  1:size(activations,1)
     subplot(2,3,x)
     imagesc(ewts(:,:,x));
 end;
 
 mn2 = mean(ewts(:,:,2),2); mn3 = mean(ewts(:,:,3),2);
cols = jet(15);
emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
figure; 
for e = 1:length(emos)
    ph = plot(mn2(e),mn3(e),'k.');hold on;
    set(ph,'markersize',20); set(ph,'color',cols(e,:));
    ph = text(mn2(e),mn3(e),emo2{e});
    set(ph,'color',cols(e,:));     
end;
