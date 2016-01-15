






function PairEmos(emomeans,subjlist,useims,selemos,allemos);
    

    emoidx = find(ismember(allemos,selemos));

    

    %indiv = 0;
    %cols = jet(15);cols(10,:) = [.9 .9 0];
    %emo2 = {'  anger','  frustration','  jealousy','  fear' ,'  disgust','  grief','  sad','  compassion','  love','  relief','  content','  awe','  happy','  joy','  excited'};
    if length(size(emomeans{1})) > 2
        keeptrack = zeros(0,3);
    else
        keeptrack = zeros(0,2);
    end;
    % take instead IMs that were specified 
    collmeans = []; 
    for nx = 1:length(emomeans)
        if ~isempty(useims{nx}) & ~isempty(find(ismember(nx,subjlist)))
            for dm = 1:length(useims{nx})
                if length(size(emomeans{nx})) > 2
                    tmpmeans =  squeeze(emomeans{nx}(abs(useims{nx}(dm)),emoidx,:)); %  emos x deciles 
                    collmeans = [collmeans, tmpmeans]; % emos x IMs*ndeciles
                    keeptrack(end+1:end+size(tmpmeans,2),:) = [repmat(nx,[size(tmpmeans,2) 1]), repmat(useims{nx}(dm),[size(tmpmeans,2) 1]), [1:9]'];
                else              
                    %collmeans = [collmeans, emomeans{nx}(abs(useims{nx}(dm)),emoidx)'];% single subject emo space  
                    collmeans = [collmeans, emomeans{nx}(abs(useims{nx}(dm)),emoidx)'/std(emomeans{nx}(abs(useims{nx}(dm)),:))];% emos x IMs,/std     
                                                                                                                                %collmeans = [collmeans, zscore(emomeans{nx}(abs(useims{nx}(dm)),emoidx)')];% emos x IMs , zscore         
                    keeptrack(end+1,:) = [nx, useims{nx}(dm)];
                end;
            end;
        end;
    end;

    figure; 
    plot(collmeans(1,:),collmeans(2,:),'r.'); hold on;
    plot([0 0],[get(gca,'ylim')],'k-');
    plot([get(gca,'xlim')],[0 0],'k-');
    xlabel(selemos{1});
    ylabel(selemos{2});
    
    title([selemos{1},' vs ',selemos{2},' weights']);
