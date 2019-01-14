function correlation_fig = plot_grid_map(fig_name, map_data, bipolar_reref, frequency_type)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot correlation between log data and frequency array for plot on bipolar grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    fig_name                        : figure name
    %    correlation_data                : Correlation data
    %                                           dim = (numbe of bipolar channel id x 2)
    %    bipolar_reref                   : bipolar channel info
    %    data_type                       : 0/1  0-> power, 1-> SNR
    %    finger1/2                       : Finger1/2
    %    sweep_type                      : 'A' Ascending or 'D' Descending
    %    frequency_type                  : 'Fundamental', 'IM', 'Fundamental + IM' 
    % Output 
    %    correlation_fig                 : Correlation figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    title_text = '10Log10SNR';
    
    figure;
    set(gcf, 'Name', fig_name);
    cmap = colormap('hot');
    cmap = flipud(cmap);
    [p, FV] = draw_biprref(map_data, bipolar_reref, [5 4]);
    title({strcat('Correlation between', " ", title_text, ' spectrogram and stimulus', " ", frequency_type, " ", 'frequency array'), " "});
    colormap(cmap);
    cbar = colorbar;
    cbar.Label.String = 'R^{2} [-]';
    format_fig;
    set(gcf, 'Position', [1000, 390, 1080, 940]);
    correlation_fig = gcf;
end
