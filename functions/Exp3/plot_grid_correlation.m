function correlation_fig = plot_grid_correlation(fig_name, correlation_data, bipolar_reref, data_type, frequency_type)
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
    
    switch data_type
        case 0
            title_text = 'LogPower';
            %title_text = 'Log(Power+1)';
        case 1
            title_text = 'LogSNR';
            %title_text = 'LogSNR (computed from Log(Power+1)';
    end
    
    figure;
    set(gcf, 'Name', fig_name);
    cmap = colormap('hot');
    cmap = flipud(cmap);
    [p, FV] = draw_biprref(correlation_data, bipolar_reref, [5 4], [0 1]);
    title({strcat('Correlation between', " ", title_text, ' spectrogram and stimulus', " ", frequency_type, " ", 'frequency array'), " "});
    colormap(cmap);
    cbar = colorbar;
    cbar.Label.String = 'R^{2} [-]';
    format_fig;
    set(gcf, 'Position', [1000, 390, 1080, 940]);
    correlation_fig = gcf;
end
