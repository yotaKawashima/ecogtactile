
function spectrogram = plot_spectrogram_for_one_finger(fig_name, logdata_this_ch, faxis, trial, frequency_array_for_plot, data_type, finger, sweep_type, frequency_id)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot spectrogram for each bipolar channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    fig_name                        : figure name
    %    logdata                         : Log scaled data (Power/SNR)
    %                                           dim = (faxis x epochs)
    %    faxis                           : Frequency range
    %    trial                           : Trial id 
    %    frequency_array_for_plot        : Frequency array for plot   
    %                                           dim = (faxis x epochs x RBG) 
    %                                           value--- 1=frequency on, 0=frequency off
    %    data_type                       : 0/1  0-> power, 1-> SNR
    %    finger1/2                       : Finger1/2
    %    sweep_type                      : 'A' Ascending or 'D' Descending
    % Output 
    %    fig                             : Spectrogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    set(gcf,'Name',fig_name)
    
    % Plot frequency applied to each finger
    ax1 = subplot(1,2,1);
    imagesc(trial, faxis, frequency_array_for_plot);
    set(gca,'Ydir','Normal');
    
    % Set parameter 
    xlabel('Trial ID [-]');
    ylabel('Frequency [Hz]');
    title_freq_line1 = strcat('Frequency stimulated to', " ", finger);
    
    switch sweep_type
        case 'A'
            switch frequency_id 
                case 1
                    title_freq_line2 = '8:1:28 [Hz]';
                case 2
                    title_freq_line2 = '10:2:50 [Hz]';
            end
        case 'D'
            switch frequency_id
                case 1
                    title_freq_line2 = '50:-1:30 [Hz]';
                case 2
                    title_freq_line2 = '48:-2:8 [Hz]';
            end
    end
    
    title({title_freq_line1, title_freq_line2});
    format_fig;
    
    % Plot spectrogram
    ax2 = subplot(1,2,2);
    % Colormap
    cmap_2 = flipud(colormap('hot'));
    % Plot
    imagesc(trial, faxis, logdata_this_ch);
    set(gca,'Ydir','Normal');
    % Set parameter 
    xlabel('Trial ID [-]');
    ylabel('Frequency [Hz]');
    % Get colormap lim 
    %clim = [min(min(logdata_this_ch)), max(max(logdata_this_ch))];
    clim  = [0, max(max(logdata_this_ch))];
    
    % Set colormap limit
    caxis(clim); 
    % Get original figure size
    originalSize = get(gca, 'Position');
    % Set colorbar
    cbar = colorbar;
    switch data_type
        case 0 
            cbar_text = 'LogPower [dB]';
            %cbar_text = 'Log(Power+1) [dB]';
        case 1
            cbar_text = 'LogSNR [dB]';
            %cbar_text = 'LogSNR [dB] (computed from Log(Power+1))';
    end
    cbar.Label.String = cbar_text;
    title('Spectrogarm');
    set(gca, 'Position', originalSize);
    format_fig;
    
    % Colormap
    colormap(ax2, cmap_2);
    set(gcf, 'Position', get(0, 'Screensize'));
    %set(gcf, 'Position',[3000, 310, 940, 1000] )
    spectrogram = gcf;
    %pause;
    %close all;
    
end
