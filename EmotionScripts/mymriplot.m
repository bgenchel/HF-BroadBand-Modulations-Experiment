% mri3dplot() - plot 3-D array on top of the average MRI image
%                 from the dipplot() function. Plot from (default) 'top' (axial), 
%                 'side' (sagittal), or 'rear' (coronal).
% Usage:
%      >> [smoothed_3ddens, mriplanes] = mri3dplot(array3d, mri, 'key', 'val');
%
% Input: 
%   array3d  - 3-D array to plot translucently on top of MRI image planes
%   mri      - MRI structure (returned by dipoledensity for instance)
%
% Optional inputs:
%   'mriview'   - ['top'|'side'|rear'] MRI image projection to plot. 'Axial',
%                 'coronal', and 'saggital' are also recognized keywords.
%                 {default|[]: 'top'}
%   'mrislices' - [real] MNI coordinates of slices  {default: every 10 mm}
%   'kernel'    - 3-D smoothing. {default: 0 voxels}.
%   'geom'      - [rows, cols] of subplot for output figure.
%   'rotate'    - [0|90|180|270] rotate 2-D graphic. Default is 90.
%   'cmap'      - [float array] colormap for 3-D array. default is hot.
%   'cmax'      - [float] manual scaling of last color. Default: from data.
%
% Outputs:
%  smoothed_3ddens - smoothed and plotted 3-d density matrix
%  mriplanes - depth (in mri_view direction) of the plotted mri image slices
%
% Author: Arnaud Delorme & Scott Makeig, SCCN, 10 June 2003
%
% See also: plotmri()

 % Copyright (C) Arnaud Delorme, sccn, INC, UCSD, 2003-
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [mriplot, mriplanes] = mymriplot(prob3d, mri, varargin)
% REDUCEPATCH  Reduce number of patch faces.

    
    if nargin < 1
        help mri3dplot;
        return;
    end;

    translucency = 0.5;    % default alpha value
    mri_lim      = 85;     % axis limits of MNI image
    
    g = finputcheck( varargin, { 'mriview'   'string'   { 'sagital' 'axial' 'coronal' ...
                                                          'top' 'side' 'rear' }   'top';
                        'mrislices' 'float'    []                        [];
                        'view'      'float'    []                        [];
                        'geom'      'float'    []                        [];
                        'cmap'      'float'    []                        jet;
                        'cmax'      'float'    []                        [];
                        'rotate'    'integer'  { 0 90 180 270 }          90;
                        'kernel'    'float'    []                        0 });
    
    if isstr(g), error(g); end;
        
    if strcmpi(g.mriview,'sagittal'),    g.mriview = 'side'; 
    elseif strcmpi(g.mriview,'axial'),   g.mriview = 'top'; 
    elseif strcmpi(g.mriview,'coronal'), g.mriview = 'rear';
    end;
    
    
    % normalize prob3d for 1 to ncolors and create 3-D dim
    % -----------------------------------------------
    if g.kernel ~= 0
        disp('Smoothing...');
        smoothprob3d    = smooth3(prob3d, 'gaussian', g.kernel);
        prob3d          = smoothprob3d;
    end;
    maxdens = max(prob3d(:));
    mindens = min(prob3d(:));
    
     g.cmap=jet(128); % don't know why I had to add this
        
    if mindens < 0

        %g.cmap([58:70],:) = repmat([1 1 1],[13 1]); % make zero white
        g.cmap([58:70],:) = repmat([.75 .66 .54],[13 1]); % make zero brownish
    else    
        g.cmap([1:71],:) = repmat([.75 .66 .54],[71 1]); % make zero brownish
        g.cmap([1:60],:) = []; % lower to mask more green, increase to mask less
    end;
    ncolors = size(g.cmap,1);

    lim = max(abs(prob3d(:)));
    if lim == 0 % so function doesn't crash with no sig points
        lim = 1;
    end;
    
    if isempty(g.cmax), 
        if mindens < 0
            if abs(maxdens) > abs(mindens)
                prob3d = prob3d + maxdens;
            else                
                prob3d = prob3d + abs(mindens);
            end;
            g.cmax = 2*lim;
        else
            g.cmax = lim;
        end;
    else
        if mindens < 0
            prob3d = prob3d + g.cmax;
            g.cmax = 2*g.cmax;
            prob3d(find(prob3d < 0)) = 0;
            prob3d(find(prob3d > g.cmax*2)) = g.cmax*2;
        end;        
    end;        

    prob3d = round(prob3d/g.cmax*(ncolors-1))+1; 
    prob3d( find(prob3d > ncolors) ) = ncolors;
    newprob3d = zeros(size(prob3d,1), size(prob3d,2), size(prob3d,3), 3);
    tmp = g.cmap(prob3d,1); newprob3d(:,:,:,1) = reshape(tmp, size(prob3d));
    tmp = g.cmap(prob3d,2); newprob3d(:,:,:,2) = reshape(tmp, size(prob3d));
    tmp = g.cmap(prob3d,3); newprob3d(:,:,:,3) = reshape(tmp, size(prob3d));

    % plot MRI slices
    % ---------------
    if isempty(g.mrislices), g.mrislices = linspace(-50, 50, 11); end;
    if isempty(g.geom), 
        g.geom = ceil(sqrt(length(g.mrislices)+1)); 
        g.geom(2) = ceil((length(g.mrislices)+1)/g.geom);
    end;

    fig = figure;
    pos = get(gcf, 'position');
    set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/4*g.geom(1) pos(4)/3*g.geom(2) ]);% new
    %set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3) pos(4)/g.geom(2)*g.geom(1) ]);
    disp('Plotting...');
    for index = 1:length( g.mrislices )
        mysubplot(g.geom(1), g.geom(2), index);
        switch g.mriview
         case 'side', coord = [  g.mrislices(index) 0 0 1 ]; 
         case 'top' , coord = [  0 0 g.mrislices(index) 1 ]; 
         case 'rear', coord = [  0 g.mrislices(index) 0 1 ]; 
        end;
        
        coord = round( pinv(mri.transform)*coord' )';
        
        switch g.mriview
         case 'side', mriplot  = squeeze( mri.anatomy(coord(1), :, :) );
         case 'top' , mriplot  = squeeze( mri.anatomy(:, :, coord(3)) );
         case 'rear', mriplot  = squeeze( mri.anatomy(:, coord(2), :) );
        end;
        mriplot(:,:,2) = mriplot(:,:,1);
        mriplot(:,:,3) = mriplot(:,:,1);
        mriplot = rotatemat( mriplot, g.rotate );
        
        switch g.mriview
         case 'side', densplot = squeeze( newprob3d  (coord(1), :, :, :) );
         case 'top' , densplot = squeeze( newprob3d  (:, :, coord(3), :) );
         case 'rear', densplot = squeeze( newprob3d  (:, coord(2), :, :) );
        end;
        densplot = rotatemat( densplot, g.rotate );
        
        mriplot  = mriplot/2+densplot/2;
        mriplanes(:,:,:,index) = mriplot;
        
        imagesc(mriplot);
        axis off;
        xl = xlim;
        yl = ylim;
        zl = zlim;
        hold on;
        
        %options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
        %            'scaled','facelighting', 'none', 'facealpha', translucency};
        %h = surface( [xl(1) xl(2); xl(1) xl(2)], [yl(1) yl(1); yl(2) yl(2)], ...
        %             [1 1; 1 1], densplot, options{:});
        axis equal;
        if ~isempty(g.view), view(g.view); end;
        %title( [ int2str(g.mrislices(index)) ' mm' ], 'color', 'w');
    end;
    
    % plot colorbar
    % -------------
    if 1
        h = mysubplot(g.geom(1), g.geom(2), length(g.mrislices)+1);
        pos = get(h, 'position');
        pos(1) = pos(1)+pos(3)/3;
        pos(3) = pos(3)/6;
        pos(2) = pos(2)+pos(4)/5;
        pos(4) = pos(4)*3/5;
        axis off;
        h = axes('unit', 'normalized', 'position', pos);
        tmpmap = g.cmap/2 + ones(size(g.cmap))/4;
        colormap(tmpmap);
        if mindens < 0
            cbar(h, [1:length(g.cmap)], [-g.cmax/2 g.cmax/2]);
        else       
            cbar(h, [1:length(g.cmap)], [0 g.cmax]);
        end;
        box off;
        set(h, 'ycolor', [0.7 0.7 0.7]);
    end;
    
    fprintf('\n');
    %set(gcf,'color', g.cmap(1,:)/2);
    %set(gcf,'color', [.52 .52 .52]); % match the grey (= mri plus white)
    set(gcf,'color', [1 1 1]); % set background to white
return;

function mat = rotatemat(mat, angle);
    
    if angle == 0, return; end;
    if angle >= 90,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;
    if angle >= 180,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;
    if angle >= 270,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;

function h = mysubplot(geom1, geom2, coord);
    
    coord = coord-1;
    horiz_border = 0;
    vert_border  = 0.1;
    
    coordy = floor(coord/geom1);
    coordx = coord - coordy*geom1;
    
    posx   = coordx/geom1+horiz_border*1/geom1/2;
    posy   = 1-(coordy/geom2+vert_border*1/geom2/2)-1/geom2;
    width  = 1/geom1*(1-horiz_border);
    height = 1/geom2*(1- vert_border);
    
    h = axes('unit', 'normalized', 'position', [ posx posy width height ]);
    %h = axes('unit', 'normalized', 'position', [ coordx/geom1 1-coordy/geom2-1/geom2 1/geom1 1/geom2 ]);
    
