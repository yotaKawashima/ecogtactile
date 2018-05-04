function [p, FV] = draw_biprref(values, ch_label, spatial_cfg, clim)

% draw_biprref: Calls the patch function to plot bipolar re-referenced data
%   draw_biprref(values, ch_label, spatial_cfg) draws the image.
%
%   draw_biprref(values, ch_label, spatial_cfg, clim) uses clim as the
%   limits of the colormap.
%
%   p = draw_biprref(...) returns the handle of the patch.
%
%   [p, FV] = draw_biprref(...) returns the handle and the FV data of the
%   patch.
%
%   INPUT FORMAT 
%   values: a 1D matrix holding the values of data for ONLY the 
%           bipolar-rereferenced channels.
%   ch_label: numel(values)*3 matrix. columns 1 and 2 contain the numbers
%           of the channels the re-referenced channel is made of. column 3
%           is the direction (0 for vertical, 1 for horizontal)
%   spatial_cfg: [rows, columns] of the spatial arrangement of the channels

% Written by Rannee Li, Jan 2015
% Change log:
%   Jul 2015:   fixed some typos, vectorised make_vertices
%               subfunction add_face replaced with subfunction all_faces, 
%               which is vectorised
%               all for loops have now been replaced with array operations
%               updated clim in case of NaN (assigned min value) and +/-Inf
%               (assigned max/min value)
%   Jul 2017:   added NaN-ignoring versions of min/max for CLim property

% call subfunctions
vertices = make_vertices(spatial_cfg(2), spatial_cfg(1));
faces = all_faces(ch_label, vertices);

% place into structure as requested by patch
FV.vertices = vertices;
FV.faces = faces;

clear faces vertices

% draw
% figure
p = patch(FV);

% color
clear cdata
if nargin>3 && ~isempty(clim)
    set(gca,'CLim',clim)
else % use limits of data if no limits provided (no arg or empty arg)
    % first remove infs and nans
    if ~isempty(find(isnan(values), 1))
        warning('NaNs detected in input values. Assigning NaN = min(values).')
        values(isnan(values)) = min(values);
    end
    if ~isempty(find(isinf(values), 1))
        warning('Infs detected in input values. Assigning +/-Inf = max/min(values).')
        noinfval = values;
        noinfval(isinf(values)) = NaN;
        values(values==Inf) = max(noinfval);
        values(values==-Inf) = min(noinfval);
    end
    
    set(gca,'CLim',[nanmin(values) nanmax(values)])
end
cdata = values;
set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled')
axis tight
axis equal
axis ij

% clear if no outputs
if nargout < 2
    clear FV
end
if nargout < 1
    clear p
end

%% END MAIN FUNCTION CALL. SUBFUNCTIONS FOLLOW

%% MAKE_VERTICES
% creates all vertices needed for the patch.
function vertices = make_vertices(spatial_width, spatial_height)

X_gv = 0.5:0.5:(spatial_width + 0.5);
Y_gv = 0.5:0.5:(spatial_height + 0.5);

Y_use = Y_gv(1:2:end);
nummidpts = numel(Y_use);

% make the points that represent the channels
tmp = {1:spatial_height, 1:spatial_width};
[o1, o2] = ndgrid(tmp{:});
vert_1 = [o2(:) o1(:)];

% make vertices in column 1
vert_2 = [X_gv(1)* ones(nummidpts-2, 1), Y_gv(3:2:end-2)'];

% make vertices in column 2 - end-1
tmp = {Y_use, X_gv(3:2:end-2)};
[o1, o2] = ndgrid(tmp{:});
vert_3 = [o2(:) o1(:)];

% make vertices in end column
vert_4 = [X_gv(end)* ones(nummidpts-2, 1), Y_gv(3:2:end-2)'];

vertices = [vert_1; vert_2; vert_3; vert_4];


%% ALL_FACES
% maps all of the bipolar channels onto the vertices
function face = all_faces(ch_label, vertices)
% find the midpoints of all the channels
midps = (vertices(ch_label(:, 1), :) + vertices(ch_label(:, 2), :)) /2;

% find all the horizontally re-referenced channels
horzinds = find(ch_label(:, 3));
vertinds = find(~ch_label(:, 3));

% calculate the location of the other two points
% split into horx and vert as this affects the result
point2 = zeros(size(ch_label, 2), 2); % preallocate
point4 = point2;

point2(horzinds, :) = [midps(horzinds, 1) midps(horzinds, 2)+0.5];
point4(horzinds, :) = [midps(horzinds, 1) midps(horzinds, 2)-0.5];

point2(vertinds, :) = [midps(vertinds, 1)+0.5 midps(vertinds, 2)];
point4(vertinds, :) = [midps(vertinds, 1)-0.5 midps(vertinds, 2)];

% locate the actual vertex numbers
[~, f_2] = ismember(point2, vertices, 'rows');
[~, f_4] = ismember(point4, vertices, 'rows');

face = [ch_label(:, 1) f_2 ch_label(:, 2) f_4];
% the first two points are the same as the two channels
%%

