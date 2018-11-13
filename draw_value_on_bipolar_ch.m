function figure = draw_value_on_bipolar_ch(values, ch_label, spatial_cfg, start_loc, direction, nColor, cLim)

% draw_value_on_bipolar_ch: Plot bipolar re-referenced data on bipolar
%                           channels
%   INPUT
%       values      : 1D matrix 
%                       the values of data for ONLY the
%                       bipolar-rereferenced channels
%       ch_label    : 2D marix
%                       columns 1 and 2 contain the numbers of the
%                       channels the re-referenced channel is made of
%       spatial_cfg : 2D matrix
%                       [rows, columns] of the spatial arrangement of the
%                       channels
%       start_loc   : String
%                       the location where electrode id starts 
%                       (UpLeft, DownLeft, UpRight, DownRight)
%       direction   : String
%                       the direction to wihch electrode id increases 
%                       (Up, Down, Left, Right)
%       --optional--
%       nColo       : Integer
%                       the number of colors by which values are represented  
%       cLim        : 1D matrix
%                       [min, max] min and max of colorbar

%   OUTPUT
%       figure      : Fig
%                       The figure
%   

% Written by Yota Kawashima, Nov 2018
% Change log:
%   Nov 20??:   Comment here
%
%

% Create electrode information and plot electrodes
dimRow = spatial_cfg(2);
dimColumn = spatial_cfg(1);
nElectrode = dimRow * dimColumn;
electrode_coordinate = ones(dimRow, dimColumn,2);
electrode_ids = ones(dimRow, dimColumn);
x_fake = dimRow:-1:1;
y_fake = dimColumn:-1:1;

ax1 = axes;
hold on;
for x = 1:dimRow
    for y = 1:dimColumn
        % Set electrode id
        switch start_loc
            case 'UpLeft'
                switch direction
                    case 'Down' 
                        electrode_ids(x,y) = (x-1)*dimColumn + y_fake(y); 
                    case 'Right'
                        electrode_ids(x,y) = x + (y_fake(y)-1)*dimRow; 
                    otherwise
                        disp('Set proper direction');
                        return;
                end
                    
            case 'DownLeft'
                switch direction
                    case 'Up'
                        electrode_ids(x,y) = (x-1)*dimColumn + y;
                    case 'Right'
                        electrode_ids(x,y) = x + (y-1)*dimRow;
                    otherwise
                        disp('Set proper direction');
                        return;
                end
                
            case 'UpRight'
                switch direction
                    case 'Down'
                        electrode_ids(x,y) = (x_fake(x)-1)*dimColumn + y_fake(y);
                    case 'Left'
                        electrode_ids(x,y) = x_fake(x) + (y_fake(y)-1)*dimRow;
                    otherwise
                        disp('Set proper direction');
                        return;
                end
                
            case 'DownRight'
                switch direction
                    case 'Up'
                        electrode_ids(x,y) = (x_fake(x)-1)*dimColumn + y;
                    case 'Left'
                        electrode_ids(x,y) = x_fake(x) + (y-1)*dimRow;
                    otherwise
                        disp('Set proper direction');
                        return;
                end
            otherwise 
                dips('Set proper start_loc');
                return;
                
        end
        
        electrode_coordinate(x, y, 1) = x;
        electrode_coordinate(x, y, 2) = y;
        scatter(x, y, 800, 'black');
        current_electrode_id = electrode_ids(x,y);
        text(x-0.05, y, num2str(current_electrode_id), 'Fontsize', 15);
        
    end
end

% Plot biploar values using electrode information
nBipolarch = size(ch_label,1);

% Create colormap

% Normalise data into [1,64] range. (Default colormap has [1,64] range.)
% First normalise data into [0,1] range, then renormalise it into [1,64]
% range. (which means the number of colors is 64).

% You can change the number of color.  
if nargin>5 && ~isempty(nColor) % When nColor is not identified.
    nColor = 64;
end

% Change cLim
if nargin>6 && ~isempty(cLim)
    % Get logical indices of values which are below min or beyond max
    less_than_min = values <= cLim(1);
    more_than_max = values >= cLim(2);
    
    % Update values 
    values(less_than_min) = cLim(1);
    values(more_than_max) = cLim(2);
    min_value = cLim(1);
    max_value = cLim(2);
else
    min_value = min(values);
    max_value = max(values);
end

color_ids = floor(((values-min_value)/(max_value-min_value))*(nColor-1)) + 1;
old_cmap = hot(nColor);
cmap=flipud(old_cmap);    

% Plot values on bipolar channel
for bipolar_ch_id = 1:nBipolarch
    
    % Get electrodes which this bipolar channel composes of.
    electrode1 = ch_label(bipolar_ch_id,1);
    electrode2 = ch_label(bipolar_ch_id,2);
    
    % Extract indices of electrodes
    [m1, n1] = find(electrode_ids == electrode1);
    [m2, n2] = find(electrode_ids == electrode2);
    
    % Get coordinates of electrodes
    x1 = electrode_coordinate(m1, n1, 1);
    y1 = electrode_coordinate(m1, n1, 2);
    x2 = electrode_coordinate(m2, n2, 1);
    y2 = electrode_coordinate(m2, n2, 2);
    
    % Get color id for the value of current bipolar channel
    color_id = color_ids(bipolar_ch_id);
    color_now = cmap(color_id,:);
    % Plot the value as line    
    plot(ax1,[x1,x2],[y1,y2], 'Color', color_now, 'LineWidth', 20);
end

hold off;

% Change object hierarchy 
hierarchy_ax1 = get(ax1, 'Children');
bipolar_plot_obj_ind_fake = zeros(nBipolarch,1);
electrode_plot_obj_ind_fake = ones(nElectrode*2,1);
electrode_plot_obj_ind = cat(1, bipolar_plot_obj_ind_fake, electrode_plot_obj_ind_fake);
electrode_plot_obj_ind = logical(electrode_plot_obj_ind);
new_hierarchy = [hierarchy_ax1(electrode_plot_obj_ind); hierarchy_ax1(~electrode_plot_obj_ind)];
set(ax1, 'Children', new_hierarchy);

% Make invisible plot just for creating colorbar
ax2 = axes;
colormap(cmap);
inv_p = scatter(ax2, ones(length(values),1), ones(length(values),1), 1, values);
set(inv_p, 'visible', 'off');

% Link ax1 and ax2 together
linkaxes([ax1, ax2]);
ax1.Visible ='off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

set([ax1, ax2],'Position',[.03 .03 .685 .815]);
colorbar(ax2,'Position',[.88 .11 .0675 .815]);
caxis(ax2, [min_value, max_value]);
figure = gcf;
