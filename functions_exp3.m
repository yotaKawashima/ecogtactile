%% Functions
% For experiment 3

% Set experiment conditions
function [ascending_frequency, descending_frequency] = make_frequency_array()
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    None
    % Output 
    %    ascending_frequency    : Frequency array for ascending condition 
    %    descending_frequency   : Frequency array for descending condition
    %    ----------
    %    dim = (21 x 2)
    %    [[frequency1, frequency2],
    %     [frequency1, frequency2],
    %              ...
    %     [frequency1, frequency2]]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ascending
    ascending_frequency = zeros(21,2);
    frequency1 = 8 ;
    frequency2 = 10 ;
    for f_id = 1:21
        % Store frequency
        ascending_frequency(f_id,:) = [frequency1, frequency2];

        % Update frequency
        frequency1 = frequency1 + 1;
        frequency2 = frequency2 + 2;
    end % f_id

    % Descending 
    descending_frequency = zeros(21,2);
    frequency1 = 50 ;
    frequency2 = 48 ;
    for f_id = 1:21
        % Store frequency
        descending_frequency(f_id,:) = [frequency1, frequency2];

        % Update frequency
        frequency1 = frequency1 - 1;
        frequency2 = frequency2 - 2;
    end % f_id
end

function frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    faxis                  : Frequency array from FFT
    %    frequency_array        : Frequency array for trial condition
    % Output 
    %    frequency_array_for_correlation   : Frequency array for plot
    %    ----------
    %    dim = (faxis x condition)
    %    if frequency on, then assign 1 to the element of the array
    %    else (frequency off), then assign 0 to the element of the array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation 
    frequency_array_for_correlation = zeros(length(faxis), 21);
    
    for condition_id = 1:21
        for frequency_id = 1:size(frequency_array, 2)
            % Get frequency
            freq_for_one_finger = frequency_array(condition_id, frequency_id);
            
            % Get index
            [~,idxF]=findclosest(faxis,freq_for_one_finger);
            
            % Assign 1 to the element
            frequency_array_for_correlation(idxF, condition_id) = 1;
            
        end
    end
    
end

function fundamental_IM_array_for_image = get_frequency_array_for_image(fundamental_array_for_correlation, IM_array_for_correlation)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    fundamental_array_for_correlation : Fundamental frequency array for correlation 
    %    IM_array_for_correlation          : IM frequency for correlation
    % Output 
    %    fundamental_IM_frequency_array_for_image   : Frequency array for image
    %    ----------
    %    dim = (faxis x condition x RGB)
    %    fundamental-> black, IM-> gray, background-> white
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation
    fundamental_IM_array_for_imageR = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageG = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageB = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    
    % Set black for fundamental frequency
    fundamental_IM_array_for_imageR(logical(fundamental_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageG(logical(fundamental_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageB(logical(fundamental_array_for_correlation)) = 0;
    
    % Set pink for IM frequency
    fundamental_IM_array_for_imageR(logical(IM_array_for_correlation)) = 240/256;
    fundamental_IM_array_for_imageG(logical(IM_array_for_correlation)) = 0/256;
    fundamental_IM_array_for_imageB(logical(IM_array_for_correlation)) = 120/256;
    
    fundamental_IM_array_for_image = cat(3, fundamental_IM_array_for_imageR, fundamental_IM_array_for_imageG, fundamental_IM_array_for_imageB); 
end

function fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger(target_array_for_correlation, fundamental_array_for_correlation, IM_array_for_correlation)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    target_array_for_correlation      : Fundamental frequency array, whose frequency is allocated to the target finger    
    %    fundamental_array_for_correlation : Fundamental frequency array for correlation 
    %    IM_array_for_correlation          : IM frequency for correlation
    % Output 
    %    fundamental_IM_frequency_array_for_image   : Frequency array for image
    %    ----------
    %    dim = (faxis x condition x RGB)
    %    fundamental-> black, IM-> gray, background-> white
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation
    fundamental_IM_array_for_imageR = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageG = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageB = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    
    % Set gray for fundamental frequency
    fundamental_IM_array_for_imageR(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageG(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageB(logical(fundamental_array_for_correlation)) = 180/256;
    
    % Set black for target fundamental frequency
    fundamental_IM_array_for_imageR(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageG(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageB(logical(target_array_for_correlation)) = 0;
    
    % Set pink for IM frequency
    if nargin > 3
        fundamental_IM_array_for_imageR(logical(IM_array_for_correlation)) = 240/256;
        fundamental_IM_array_for_imageG(logical(IM_array_for_correlation)) = 0/256;
        fundamental_IM_array_for_imageB(logical(IM_array_for_correlation)) = 120/256;
    end
    
    fundamental_IM_array_for_image = cat(3, fundamental_IM_array_for_imageR, fundamental_IM_array_for_imageG, fundamental_IM_array_for_imageB); 
end

function fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger_4IM(target_array_for_correlation, fundamental_array_for_correlation, IM1_array_for_correlation, IM2_array_for_correlation, IM3_array_for_correlation, IM4_array_for_correlation)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make frequency array for correspondence of group id to frequency
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    target_array_for_correlation      : Fundamental frequency array, whose frequency is allocated to the target finger    
    %    fundamental_array_for_correlation : Fundamental frequency array for correlation 
    %    IM_array_for_correlation          : IM frequency for correlation
    % Output 
    %    fundamental_IM_frequency_array_for_image   : Frequency array for image
    %    ----------
    %    dim = (faxis x condition x RGB)
    %    fundamental-> black, IM-> gray, background-> white
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialisation
    fundamental_IM_array_for_imageR = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageG = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    fundamental_IM_array_for_imageB = ones(size(fundamental_array_for_correlation,1), size(fundamental_array_for_correlation,2));
    
    % Set gray for fundamental frequency
    fundamental_IM_array_for_imageR(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageG(logical(fundamental_array_for_correlation)) = 180/256;
    fundamental_IM_array_for_imageB(logical(fundamental_array_for_correlation)) = 180/256;
    
    % Set pink for IM frequency
    fundamental_IM_array_for_imageR(logical(IM1_array_for_correlation)) = 137/256;
    fundamental_IM_array_for_imageG(logical(IM1_array_for_correlation)) = 15/256;
    fundamental_IM_array_for_imageB(logical(IM1_array_for_correlation)) = 80/256;
    
    fundamental_IM_array_for_imageR(logical(IM2_array_for_correlation)) = 256/256;
    fundamental_IM_array_for_imageG(logical(IM2_array_for_correlation)) = 24/256;
    fundamental_IM_array_for_imageB(logical(IM2_array_for_correlation)) = 69/256;
    
    fundamental_IM_array_for_imageR(logical(IM3_array_for_correlation)) = 245/256;
    fundamental_IM_array_for_imageG(logical(IM3_array_for_correlation)) = 144/256;
    fundamental_IM_array_for_imageB(logical(IM3_array_for_correlation)) = 178/256;
    
    fundamental_IM_array_for_imageR(logical(IM4_array_for_correlation)) = 253/256;
    fundamental_IM_array_for_imageG(logical(IM4_array_for_correlation)) = 229/256;
    fundamental_IM_array_for_imageB(logical(IM4_array_for_correlation)) = 237/256;
        
    % Set black for target fundamental frequency
    fundamental_IM_array_for_imageR(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageG(logical(target_array_for_correlation)) = 0;
    fundamental_IM_array_for_imageB(logical(target_array_for_correlation)) = 0;
    
    fundamental_IM_array_for_image = cat(3, fundamental_IM_array_for_imageR, fundamental_IM_array_for_imageG, fundamental_IM_array_for_imageB); 
end

function [finger1, finger2] = get_finger_pair(finger_pair_id)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get finger pair out of Thumb, Index, and Middle, using finger_pair_id  
    % Just for experiment 3 (Sweep experiment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    finger_pair_id         : a -> (finger1, finger2) = (Thumb, Index)
    %                             b -> (finger1, finger2) = (Index, Middle)
    %                             c -> (finger1, finger2) = (Middle, Thumb)
    % Output 
    %    finger1                : Being applied frequency1  
    %    finger2                : Being applied frequency2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch finger_pair_id
        case 'a'
            finger1 = 'Thumb';
            finger2 = 'Index';
        case 'b'
            finger1 = 'Index';
            finger2 = 'Middle';
        case 'c' 
            finger1 = 'Middle';
            finger2 = 'Thumb';
        otherwise 
            disp('Set a correct id for getting the finger pair!');
    end

end

function [file_name_prefix1, file_name_prefix2] = get_two_file_name_prefix(file_name_prefixs, sweep_id, frequency_id, finger)
    % Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get two file name prefix for the same finger, the same frequency, and the same sweep type   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input 
    %    file_name_prefixs          : all file name prefixs (cell array) 
    %    sweep_id                   : Sweep_id (Ascending = 1, Descending = 2)
    %    frequency_id               : Frequency_id (frequency1 = 1, frequency2 = 2)
    %    finger                     : Finger(Thumb/Index/Middle)
    % Output 
    %    file_name_prefix1/2        : File name prefix 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch frequency_id 
        case 1
            switch finger
                case 'Thumb'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 1),1}; % (Thumb, Index)  = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 6),1}; % (Thumb, Middle) = (frequency1, frequency2)
                case 'Index'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 2),1}; % (Index, Middle) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 4),1}; % (Index, Thumb)  = (frequency1, frequency2)
                case 'Middle'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 3),1}; % (Thumb, Index)  = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 5),1}; % (Thumb, Middle) = (frequency1, frequency2)
            end

        case 2
            switch finger
                case 'Thumb'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 3),1}; % (Middle, Thumb) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 4),1}; % (Index,  Thumb) = (frequency1, frequency2)
                case 'Index'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 1),1}; % (Thumb,  Index) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 5),1}; % (Middle, Index) = (frequency1, frequency2)
                case 'Middle'
                    file_name_prefix1 = file_name_prefixs{(6*(sweep_id -1) + 2),1}; % (Index, Middle) = (frequency1, frequency2)
                    file_name_prefix2 = file_name_prefixs{(6*(sweep_id -1) + 6),1}; % (Thumb, Middle) = (frequency1, frequency2)
            end
    end

end

% Plot 
function spectrogram = plot_spectrogram(fig_name, logdata_this_ch, faxis, trial, frequency_array_for_image, data_type, finger1, finger2, sweep_type)
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
    ax1 = subplot(2,1,1);
    image(trial, faxis, frequency_array_for_image);
    set(gca,'Ydir','Normal');
    
    % Set parameter 
    xlabel('Trial ID [-]');
    ylabel('Frequency [Hz]');
    title_freq_line1 = strcat('Frequency stimulated to', " ", finger1, " ", 'and',  " ", finger2);
    
    switch sweep_type
        case 'A'
        title_freq_line2 = strcat(finger1, '=8:1:28 [Hz],',  "  ", finger2, '=10:2:50 [Hz]');
        case 'D'
        title_freq_line2 = strcat(finger1, '=50:-1:30 [Hz],',  "  ", finger2, '=48:-2:8 [Hz]');
    end
    title({title_freq_line1, title_freq_line2});
    format_fig;
    
    % Plot spectrogram
    ax2 = subplot(2,1,2);
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
            cbar_text = 'LogPower [-]';
            %cbar_text = 'Log(Power+1) [-]';
        case 1
            cbar_text = 'LogSNR [-]';
            %cbar_text = 'LogSNR [-] (computed from Log(Power+1))';
    end
    cbar.Label.String = cbar_text;
    title('Spectrogarm');
    set(gca, 'Position', originalSize);
    format_fig;
   
    % Colormap
    cmap_2 = flipud(colormap('hot'));
    colormap(ax2, cmap_2);
    set(gcf, 'Position',[3000, 310, 940, 1000] )
    spectrogram = gcf;
    %pause;
    %close all;
    
end

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

function correlation_fig = plot_grid_correlation_YV(fig_name, correlation_data, bipolar_reref, data_type, frequency_type)
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
    %[p, FV] = draw_biprref(correlation_data, bipolar_reref, [5 4], [0 1]);
    draw_value_on_bipolar_ch(correlation_data, bipolar_reref, [5 4], 'UpLeft', 'Down', 64, [0 1]);
    set(gca, 'color', [220/256, 220/256, 220/256]);
    title({strcat('Correlation between', " ", title_text, ' spectrogram and stimulus', " ", frequency_type, " ", 'frequency array'), " "});
    %colormap(cmap);
    %cbar = colorbar;
    cbar.Label.String = 'R^{2} [-]';
    format_fig;
    set(gcf, 'Position', [1000, 390, 1080, 940]);
    correlation_fig = gcf;
end

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
            cbar_text = 'LogPower [-]';
            %cbar_text = 'Log(Power+1) [-]';
        case 1
            cbar_text = 'LogSNR [-]';
            %cbar_text = 'LogSNR [-] (computed from Log(Power+1))';
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


%% For experiment 1
% Plot 
function plot_spectrum(fig_name, logSNR, logPow, trials, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot spectrum (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    subplot(1,2,1); 
    plot(faxis,mean(logPow(bipolar_reref(:,3)==0,:,trials),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logPow(bipolar_reref(:,3)==1,:,trials),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(Power) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;
    
    subplot(1,2,2); 
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==0,:,trials),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==1,:,trials),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(SNR) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;

    % mark fundamentals
    subplot(1,2,2); hold on;
    lgd=[];
    for nF=1:length(LeftF0)
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        lgd(1)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','o','MarkerFaceColor','r','SizeData',100);
    end
    for nF=1:length(RightF0)
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        lgd(2)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','o','MarkerFaceColor','b','SizeData',100);
    end
    format_fig;

    % mark IM
    subplot(1,2,2); hold on;
    for nF=1:length(LeftIM)
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        lgd(3)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','*','MarkerEdgeColor','r','SizeData',100,'LineWidth',1.5);
    end
    for nF=1:length(RightIM)
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        lgd(4)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','*','MarkerEdgeColor','b','SizeData',100,'LineWidth',1.5);
    end
    legend(lgd,{'contraF','ipsiF','contraIM','ipsiIM'})
    format_fig;
    set(gcf, 'Position', get(0, 'Screensize'));
end

function plot_grid_logSNR_values_vs0(fig_name, logSNR, trials, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot grid longSNR values vs0 (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    cmap=colormap('hot');
    cmap=flipud(cmap);
    for nF=1:length(LeftF0)
        subplot(2,4,nF);
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraF (%g)',LeftF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(LeftIM)
        subplot(2,4,nF+2);
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraIM (%g)',LeftIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightF0)
        subplot(2,4,nF+4);
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiF (%g)',RightF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightIM)
        subplot(2,4,nF+6);
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiIM (%g)',RightIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

function plot_grid_logSNR_values_vs_trials(fig_name, logSNR, trials_1, trials_2, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot grid logSNR values vs trials (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    cmap=colormap('hot');
    cmap=flipud(cmap);
    for nF=1:length(LeftF0)
        subplot(2,4,nF);
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials_1)),squeeze(logSNR(:,idxF,trials_2)),'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraF (%g)',LeftF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(LeftIM)
        subplot(2,4,nF+2);
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials_1)),squeeze(logSNR(:,idxF,trials_2)),'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraIM (%g)',LeftIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightF0)
        subplot(2,4,nF+4);
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials_1)),squeeze(logSNR(:,idxF,trials_2)),'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiF (%g)',RightF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightIM)
        subplot(2,4,nF+6);
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials_1)),squeeze(logSNR(:,idxF,trials_2)),'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiIM (%g)',RightIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

function plot_grid_diff_t_values_of_logSNR_vs_trials(fig_name, logSNR, trials_1, trials_2, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot grid logSNR values vs trials (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    cmap=colormap('parula');
    cmap=flipud(cmap);
    for nF=1:length(LeftF0)
        subplot(2,4,nF);
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        [h_1, pV_1, ~, stats_1]=ttest(squeeze(logSNR(:,idxF,trials_1)),0,'dim',2);
        [h_2, pV_2, ~, stats_2]=ttest(squeeze(logSNR(:,idxF,trials_2)),0,'dim',2);
        diff_t_values = mean(stats_1.tstat,3) - mean(stats_2.tstat,3);
        [p, FV] = draw_biprref(diff_t_values, bipolar_reref, [5 4], [-9 0]);
        title(sprintf('contraF (%g)',LeftF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(LeftIM)
        subplot(2,4,nF+2);
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        [h_1, pV_1, ~, stats_1]=ttest(squeeze(logSNR(:,idxF,trials_1)),0,'dim',2);
        [h_2, pV_2, ~, stats_2]=ttest(squeeze(logSNR(:,idxF,trials_2)),0,'dim',2);
        diff_t_values = mean(stats_1.tstat,3) - mean(stats_2.tstat,3);
        [p, FV] = draw_biprref(diff_t_values, bipolar_reref, [5 4], [-9 0]);
        title(sprintf('contraIM (%g)',LeftIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightF0)
        subplot(2,4,nF+4);
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        [h_1, pV_1, ~, stats_1]=ttest(squeeze(logSNR(:,idxF,trials_1)),0,'dim',2);
        [h_2, pV_2, ~, stats_2]=ttest(squeeze(logSNR(:,idxF,trials_2)),0,'dim',2);
        diff_t_values = mean(stats_1.tstat,3) - mean(stats_2.tstat,3);
        [p, FV] = draw_biprref(diff_t_values, bipolar_reref, [5 4], [-9 0]);
        title(sprintf('ipsiF (%g)',RightF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightIM)
        subplot(2,4,nF+6);
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        [h_1, pV_1, ~, stats_1]=ttest(squeeze(logSNR(:,idxF,trials_1)),0,'dim',2);
        [h_2, pV_2, ~, stats_2]=ttest(squeeze(logSNR(:,idxF,trials_2)),0,'dim',2);
        diff_t_values = mean(stats_1.tstat,3) - mean(stats_2.tstat,3);
        [p, FV] = draw_biprref(diff_t_values, bipolar_reref, [5 4], [-9 0]);
        title(sprintf('ipsiIM (%g)',RightIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

