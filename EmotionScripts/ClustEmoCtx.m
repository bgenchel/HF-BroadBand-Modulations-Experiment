% clusters emotion indicator vectors from context single-trial decomposition across subjects
%
%
% subjlist -- [vector] list of subject indexes to cluster among
% fullpaths -- [cell array of strings] full paths to directories of ALL subjects (indexes point
%              to subjects of interest
% savedat -- [string] prefix under which data is saved for all subjects in resp. directories
% clustby -- [string] 'e'-emotion, 's'-spectra, 'b'-both
% dbl -- [1 | 0] 1 to double the matrix to match clustering matrix.
%

function [clustmat,keeptrack] = ClustEmoCtx(subjlist,fullpaths,savedat,clustby,dbl);

    
    keeptrack = zeros(0,2);
    if clustby == 'e'
        clustmat = zeros(0,15); keeptrack = zeros(0,2);
        for nxx = 1:length(subjlist)
            nx = subjlist(nxx);
            s = load([fullpaths{nx},savedat,'Stuff.mat']);         
            data = floatread([fullpaths{nx},savedat,'.fdt'],[15+s.pcnum inf],[],0);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            pcaeig = floatread([fullpaths{nx},savedat,'EIGVEC.fdt'],[length(s.freqs) s.pcnum],[],0);   
            ws = wts*sph;    
            if size(ws,1) == size(ws,2)
                winv = inv(ws); 
            else
                winv = pinv(ws);
            end;
            activations = ws*data;
            clear wts sph ws emowts
            cpc = data(size(pcaeig,2)+1:end,:);
            fprintf('Subj %s\n',int2str(nx))
            emowts = zeros(0,15);
            for fac = 1:size(activations,1)
            emowts(end+1,:) = winv(size(pcaeig,2)+1:end,fac)';
            keeptrack(end+1,:) = [nx,fac];
            if dbl == 1
                emowts(end+1,:) = winv(size(pcaeig,2)+1:end,fac)'*-1;
                keeptrack(end+1,:) = [nx,fac];
            end;
            %    for p = 1:size(cpc,1)
            %        cpc(p,:) = cpc(p,:).*activations(fac,:);
            %    end;
            %    emowts(end+1,:) = mean(cpc,2)';
            %    keeptrack(end+1,:) = [nx fac];
            %    if dbl == 1
            %    emowts(end+1,:) = mean(cpc,2)'*-1;
            %    keeptrack(end+1,:) = [nx fac];
            %    end;
            %backproj = winv(:,fac) * activations(fac,:);
            %emowts(fac,:) = mean(backproj(s.pcnum+1:end,find(activations(fac,:)> mean(activations(fac,:) + 3*std(activations(fac,:))))),2); %pos weights
            end;
            if dbl == 1
                clustmat(end+1:end+s.pcs*2,:) = emowts;
            else
                clustmat(end+1:end+s.pcs,:) = emowts;
            end;
        end;
        
    elseif clustby == 's'
        clustmat = zeros(0,99);
        for nxx = 1:length(subjlist)
            nx = subjlist(nxx);
            s = load([fullpaths{nx},savedat,'Stuff.mat']);         
            data = floatread([fullpaths{nx},savedat,'.fdt'],[15+s.pcnum inf],[],0);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            pcaeig = floatread([fullpaths{nx},savedat,'EIGVEC.fdt'],[length(s.freqs) s.pcnum],[],0);    
            ws = wts*sph;    
            if size(ws,1) == size(ws,2)
                winv = inv(ws); 
            else
                winv = pinv(ws);
            end;
            activations = ws*data;
            clear wts sph ws
            fprintf('Subj %s\n',int2str(nx))
            epc = data(1:size(pcaeig,2),:);
            for fac = 1:size(activations,1)
                fprintf('Fac %s\t',int2str(fac))
                for p = 1:size(epc,1)
                    epc(p,:) = epc(p,:).*activations(fac,:);
                end;
                onespec = pcaeig*epc;
                onespec = mean(onespec,2)';
                spectemps(fac,:) = onespec;
                
                %backproj = winv(:,fac) * activations(fac,:);
                %backproj = pcaeig*backproj(1:s.pcnum,:);
                %cutval = mean(activations(fac,:)) + 2*std(activations(fac,:));
                %spectemps(fac,:) = mean(backproj(:,find(activations(fac,:)> cutval)),2)';
            end;
            %spectemps = winv(1:size(pcaeig,2),:)'*pcaeig';
            clustmat(end+1:end+s.pcs,:) = spectemps;     
            keeptrack(end+1:end+s.pcs,:) = [repmat(nx,[s.pcs 1]) [1:s.pcs]'];
        end;
        
    elseif clustby == 'b'
        clustmat = zeros(0,114);
        for nxx = 1:length(subjlist)
            nx = subjlist(nxx);
            s = load([fullpaths{nx},savedat,'Stuff.mat']);         
            data = floatread([fullpaths{nx},savedat,'.fdt'],[15+s.pcnum inf],[],0);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            pcaeig = floatread([fullpaths{nx},savedat,'EIGVEC.fdt'],[length(s.freqs) s.pcnum],[],0);    
            ws = wts*sph;    
            if size(ws,1) == size(ws,2)
                winv = inv(ws); 
            else
                winv = pinv(ws);
            end;
            fprintf('Subj %s\n',int2str(nx))
            activations = ws*data;
            clear wts sph ws spectemps
            epc = data(1:size(pcaeig,2),:);
            cpc = data(size(pcaeig,2)+1:end,:);
            for fac = 1:size(activations,1)
                for p = 1:size(epc,1)
                    epc(p,:) = epc(p,:).*activations(fac,:);
                    cpc(p,:) = cpc(p,:).*activations(fac,:);
                end;
                onespec = pcaeig*epc;
                onespec = mean(onespec,2)';
                spectemps(fac,:) = onespec;
                emowts(fac,:) = cpc';
                
                %backproj = winv(:,fac) * activations(fac,:);
                %specs = pcaeig*backproj(1:s.pcnum,:);
                %backproj = [specs;backproj(s.pcnum+1:end,:)];
                %spectemps(fac,:) = mean(backproj(:,find(activations(fac,:)> mean(activations(fac,:) + 2*std(activations(fac,:))))),2)';
            end;            
            clustmat(end+1:end+s.pcs,:) = [spectemps emowts]; 
            %clustmat(end+1:end+s.pcs,:) = spectemps;     
            %spectemps = winv(1:size(pcaeig,2),:)'*pcaeig';
            %emowts = winv(size(pcaeig,2)+1:end,:)';
            %clustmat(end+1:end+s.pcs,:) = [spectemps emowts];     
            keeptrack(end+1:end+s.pcs,:) = [repmat(nx,[s.pcs 1]) [1:s.pcs]'];
        end;        
    elseif clustby == 't' % template (doubled to give both orientations)
        clustmat = zeros(0,99); keeptrack = zeros(0,2);
        for nxx = 1:length(subjlist)
            nx = subjlist(nxx);
            s = load([fullpaths{nx},savedat,'Stuff.mat']);         
            data = floatread([fullpaths{nx},savedat,'.fdt'],[15+s.pcnum inf],[],0);
            sph=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.sph'],[s.numtrials s.numtrials],[],0); 
            wts=floatread([fullpaths{nx},savedat,'PC',int2str(s.pcs),'.wts'],[s.pcs s.numtrials],[],0); 
            pcaeig = floatread([fullpaths{nx},savedat,'EIGVEC.fdt'],[length(s.freqs) s.pcnum],[],0);    
            ws = wts*sph;    
            if size(ws,1) == size(ws,2)
                winv = inv(ws); 
            else
                winv = pinv(ws);
            end;
            fprintf('Subj %s\n',int2str(nx))
            activations = ws*data;
            clear wts sph ws 
            spectemps = zeros(0,99);
            for fac = 1:size(activations,1)
                spectemps(end+1,:) = pcaeig*winv(1:s.pcnum,fac);
                keeptrack(end+1,:) = [nx fac];
                spectemps(end+1,:) = (pcaeig*winv(1:s.pcnum,fac))*-1;
                keeptrack(end+1,:) = [nx fac];
            end;            
            clustmat(end+1:end+(s.pcs*2),:) = spectemps; 
        end;        
    end;
    

