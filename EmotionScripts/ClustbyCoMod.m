% sorts a collection of templates by the number of comodulates and groups whole IMs together



function [modtempls,comodtempls,modmeans,comodmeans,modidx,comodidx] = ClustbyCoMod(gdidx,gdtempls,gdmeans,modcorr,corrcut);
    
    
    
    comodidx = zeros(0,size(gdidx,2));
    comodtempls = zeros(0,size(gdtempls,2));
    comodmeans = zeros(0,size(gdmeans,2));
    modidx = zeros(0,size(gdidx,2));
    modtempls = zeros(0,size(gdtempls,2));
    modmeans = zeros(0,size(gdmeans,2));
            
    for m = 1:size(gdidx,1)
        nx = gdidx(m,1);
        allcorrs = modcorr{nx}{gdidx(m,2)};
        ncomod = length(find(allcorrs > corrcut));
        if ncomod > 1            
            comodidx(end+1,:) = gdidx(m,:);
            comodtempls(end+1,:) = gdtempls(m,:);
            comodmeans(end+1,:) = gdmeans(m,:);
            
        else
            modidx(end+1,:) = gdidx(m,:);
            modtempls(end+1,:) = gdtempls(m,:);
            modmeans(end+1,:) = gdmeans(m,:);
        end;
    end;
    
    
