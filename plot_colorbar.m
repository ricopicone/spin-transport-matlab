function h = plot_colorbar(dims, orientation, title_string, cmap)
% PLOT_COLORBAR plot a standalone colorbar for inclusion in a publication
%       H = PLOT_COLORBAR(DIMS, ORIENTATION TITLE_STRING) Plot a colorbar for
%       inclusion in a publication. DIMS sets the length and width of the
%       colorbar (in vertical mode). DIMS(1) will be the size of the colormap
%       used and DIMS(2) will be the number of times it is repeated (thickness
%       of image). ORIENTATION sets the orientation of the bar -- 'h', or 'v'.
%       TITLE_STRING sets the title of the axis used.
%
%       H = PLOT_COLORBAR(DIMS, ORIENTATION TITLE_STRING, CMAP) Works as above,
%       except that CMAP is a handle to a function to generate the colormap.
%
%       Examples:
%               h1 = plot_colorbar([100, 5], 'h', 'Test Colormap')
%               h2 = plot_colorbar([150, 10], 'v', 'Test Colormap', @hsv)
%
%       Bugs:
%               May not work well with wide images. 
%               Feel free to send in patches etc for any problems you find.
%
% Matt Foster <ee1mpf@bath.ac.uk>

error(nargchk(1, 4, nargin, 'struct'))

% Extract the width froms dims, if there is one.
if length(dims) < 2
  width = 5;
else
  width = dims(2);
end

if nargin < 4
  cmap = @jet;
end

map = colormap(cmap(dims(1)));

switch lower(orientation)
case {'v', 'vert', 'vertical'}
  h = image(repmat(cat(3, map(:,1), map(:,2), map(:,3)), 1, width));
  
  % Remove ticks we dont want.
  set(gca, 'xtick', 0);
  
  ticks = get(gca, 'ytick');
  set(gca, 'ytick', [0.5, ticks]);
  set(gca, 'yticklabel', [0, ticks]);

case {'h', 'horiz', 'horizontal'}
  h = image(repmat(cat(3, map(:,1)', map(:,2)', map(:,3)'), width, 1));
  
  % Remove ticks we dont want.
  set(gca, 'ytick', 0);
  
  ticks = get(gca, 'xtick');
  set(gca, 'xtick', [0.5, ticks]);
  set(gca, 'xticklabel', [0, ticks]);

otherwise
  error('unknown colorbar orientation');
end

% Set up the axis
title(title_string)
axis equal
axis tight
axis xy
