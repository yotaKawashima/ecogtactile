%% Init
clear all;
close all;

% EXPERIMENT 13 on SUBJECT 3
% Stimulation on both hands without attentional focus on hand
% 12 versions (same experiment but counterbalanced frequency mapping)
% Only the 20 fist channels (1-20) are on the somatosensory cortex

% Details for recordings (fron info.docx)
% % Sampling Rate: 1200 Hz
% % CH1: Time
% % CH2-101: ECoG 1-100 (bipolar vs. REF)
% % CH102: REF
% % CH103: GND (bipolar vs. REF)
% % CH104: DI
% % CH105: TaskId
% % CH106: GroupId
% % 
% % Variant 1 Frequency to  Finger Mapping (file name = Exp3_1*_*)
% % Exp3_1*_a:	[f1,f2] = [Thumb,Index]
% % Exp3_1*_b:	[f1,f2] = [Index,Middle]
% % Exp3_1*_c:	[f1,f2] = [Middle,Thumb]
% % 
% % Variant 2 Frequency to  Finger Mapping (file name = Exp3_2*_*)
% % Exp3_2*_a:	[f1,f2] = [Index,Thumb]
% % Exp3_2*_b:	[f1,f2] = [Middle,Index]
% % Exp3_2*_c:	[f1,f2] = [Thumb,Middle]
% %
% % Group Ids ASCENDING (file name = Exp3_*A_*):
% % 1	...	f1=8Hz		f2=10Hz
% % 2 	...	f1=9Hz		f2=12Hz
% % 3   ...	f1=10Hz     f2=14Hz
% %     ...
% % 21	…	f1=28Hz     f2=50Hz
% %
% % Group Ids DESCENDING (file name = Exp3_*D_*):
% % 1	...	f1=50Hz     f2=48Hz
% % 2 	...	f1=49Hz     f2=46Hz
% % 3	...	f1=48Hz     f2=44Hz
% % 	...
% % 21	...	f1=30Hz     f2=8Hz

% prepare path
localdef;
path_exp3data=[path_data filesep 'S3' filesep 'Exp3', filesep];
addpath(genpath(path_toolbox));
path_save_dir_exp3=[path_save_dir, filesep, 'Exp3', filesep, 'mean_power_relation', filesep];
addpath(genpath(path_functions));
%path_save_dir_exp3=[path_save_dir, filesep, 'Exp3', filesep, 'mean_power', filesep];
%path_save_dir=[path_save_dir, filesep, 'Exp3' filesep, 'normal', filesep];
if exist(path_save_dir_exp3)
    addpath(genpath(path_save_dir_exp3));
else
    mkdir(path_save_dir_exp3)
    addpath(genpath(path_save_dir_exp3));
end

% channels of interest
S1_channels = (2:21);   % CH1 is for time. CH2~Ch21 is for electrodes on the somatosensory cortex.
tot_channels = 100;     % CH2~CH101 is for electrodes.
cond_channels = 106;    % CH106 is for groupid.
load(path_bipolar_reref);

duration_stim = 4; % in seconds

% maximum frequency 
max_frequency = 400;

% preprocessing steps
%  Flag about re-referencing and filter
% make sure you pick one of the two! (bipolar or average)
preproc=[];
preproc.bipolar=1; % Whether bipolar-referencing or not
preproc.avref=0; % Whether average-referencing or not
preproc.notch=0; % Whether notch-filtering or not (Need chronux tool box)0
%% Create file name for each condition
% Freqyecny sweep type
% A -> Ascending
% D -> Descending
sweep_types = {'A', 'D'};            

% Finger pair
% a -> (finger1, finger2) = (Thumb, Index)
% b -> (finger1, finger2) = (Inde, Middle)
% c -> (finger1, finger2) = (Middle, Thumb)
finger_pair_ids = {'a', 'b', 'c'};   

% % Mapping frequency to fingers
% % 1 -> (finger1, finger2) = (Frequency1, frequency2)
% % 2 -> (finger1, finger2) = (Frequency2, Frequencty1)
variants = [1, 2];                   

% Initialisation 
sweep_type_array = {};
finger_pair_id_array = {};
variant_array = [];
file_name_prefixs = {};
% Counter for file name
counter_file_name = 1;

for sweep_id = 1:length(sweep_types)
    for finger_id = 1:length(finger_pair_ids)
        for variant_id = 1:length(variants)
            sweep_type_array{counter_file_name, 1} = sweep_types{1, sweep_id};
            finger_pair_id_array{counter_file_name, 1} = finger_pair_ids{1, finger_id};
            variant_array = [variant_array; variants(variant_id)];
            file_name_prefixs{counter_file_name, 1} = join(['Exp3_', string(variant_array(counter_file_name, 1)), sweep_type_array{counter_file_name, 1}, '_', finger_pair_id_array{counter_file_name, 1}],'');
            
            % Increment of counter 
            counter_file_name = counter_file_name + 1;
            
        end % variant_id
    end % finger_mapping_id
end % sweep_id

%% Correspondence of group id to frequency1,2 
[ascending_frequency, descending_frequency] = make_frequency_array();

%% Looping on the different versions of the experiment
%   
%   Plot data according to experiments
%


for nversion=1:length(file_name_prefixs)
    % Get condition
    sweep_type = sweep_type_array{nversion, 1};
    finger_pair_id = finger_pair_id_array{nversion, 1};
    variant = variant_array(nversion, 1);
    
    % Get data file name
    file_name_prefix = file_name_prefixs{nversion, 1};
    file_name_temp = join([path_exp3data filesep file_name_prefix, '*.mat'],'');
    
    % Get matched file name
    file_name_info = dir(file_name_temp{1}); % Here add {1}
    data_file = file_name_info.name;
    
    % Get frequency array for this version
    switch sweep_type 
        case 'A'
            frequency_array = ascending_frequency;
        case 'D'
            frequency_array = descending_frequency;
    end
    
    % Get finger pair for this version
    [finger1, finger2] = get_finger_pair(finger_pair_id);
    
    % Swap fingers if variant == 2
    if variant == 2
        % Swap fingers
        finger1_temp = finger1;
        finger2_temp = finger2;
        finger1 = finger2_temp;
        finger2 = finger1_temp;
    end
    
    % Get Intermodulation frequency array
    % Intermodulation array column=1 -> subtraction,  column=2 -> addition
    IM_frequency_array = zeros(size(frequency_array));
    IM_frequency_array(:,1) = abs(frequency_array(:,1) - frequency_array(:,2));
    IM_frequency_array(:,2) = abs(frequency_array(:,1) + frequency_array(:,2));
    
    % Load data 
    fprintf('... loading %s\n',data_file)
    load([path_exp3data filesep data_file])
    
    % Get groupid (frequency condition on each trial)
    fprintf('... finding epochs\n')
    groupid = squeeze(y(cond_channels,1,:));
    epochs = [];
    
    % Cut data into each trial
    for ngroup=1:21 
        temp_start=(groupid==ngroup);
        temp_end=(groupid==ngroup+1);
        start=find(diff(temp_start)==1)+1;
        ending=find(diff(temp_end)==1)+1;
        duration=ending-start;
        durationsec = (duration/SR);
        fprintf('... ... cond %g ... found %g epochs of length %2.2fs on average\n',ngroup,length(start), mean(durationsec))
        % matrix of epochs: condition; stimulus onset (index); epoch end; epoch start (onset -1s)
        %epochs=[epochs ; [ngroup*ones(length(start),1) start start+duration_stim*SR start-SR]];
        epochs=[epochs ; [ngroup*ones(length(start),1) start start+duration_stim*SR start]]; % Extract data during stimulus on 
        
    end
    
    % Re-reference (bi-polar referencing of the 60 electrodes of the 1st grid)
    if preproc.bipolar
        % data dim (bipolar-reref channels x 1 x sample points)
        fprintf('... bi-polar referencing on somatosensory channels\n')
        y=y(S1_channels,:,:);
        refy=nan(size(bipolar_reref,1),1,size(y,3)); 
        for nch=1:size(bipolar_reref,1)
            refy(nch,1,:)=y(bipolar_reref(nch,1),1,:)-y(bipolar_reref(nch,2),1,:); % Computation 
        end 
    elseif preproc.avref % reference to the average
        % data dim (channels in 1st grid x 1 x sample points)
        fprintf('... average referencing on somatosensory channels\n')
        y=y(S1_channels,:,:);  
        refy=y-repmat(mean(y,1),[size(y,1) 1 1]); % Computation
    else % no re-referencing
        % data dim (channels in all grids x 1 x sample points)
        fprintf('... no referencing and keep all channels\n')
        refy=y(2:tot_channels + 1,:,:);
    end
    
    % filtering
    if preproc.notch
        fprintf('... 50Hz notch-filter using chronux\n')
        data=squeeze(refy)'; % time*channel
        movingwin=[14 7]; % window length and step
        tau=10; %
        params=[];  %params has the following fields: tapers, Fs, fpass, pad
        params.tapers=[3 4];
        params.Fs=SR;
        f0= 50; % take out the 50Hz noise
        [dataf,~,~,~]=rmlinesmovingwinc(data,movingwin,tau,params,[],[],f0);
        %         f0= 100; % 
        %         [datac,~,~,~]=rmlinesmovingwinc(datac,movingwin,tau,params,[],[],f0);
        dataf=dataf'; % channel*time
    else
        dataf=squeeze(refy);
    end
    
    % Epoch Data
    fprintf('... epoching data (-1:11s)\n')
    times= 0:1/SR:duration_stim;
    sorted_data=nan(size(dataf,1),length(times),size(epochs,1)); % dim (channels x sample points x epochs)
    for ntr=1:size(epochs,1)
        sorted_data(:,:,ntr)=squeeze(dataf(:,epochs(ntr,4):epochs(ntr,3))); % epochs(:,4) - start index (sample points), epochs(:, 3) - end index (sample points)
    end
    
    % Remove erp part for experiment 3
    % erp=sorted_data-repmat(mean(sorted_data(:,times<0,:),2),[1 length(times) 1]); % Computation 
    
    % Compute log(SNR)
    fprintf('... extracting log(Pow) and log(SNR) by trial\n')
    snrparam=[];
    snrparam.method='fft'; % use FFT
    snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)
    [logSNR, faxis, logPow]=get_logSNR(sorted_data,SR,snrparam); % logSNR/logPow --dim (channels x faxis(frequency) x epochs)
    %[logSNR, faxis, logPow]=get_logSNR_from_power_plus_1(sorted_data,SR,snrparam); % logSNR/logPow --dim (channels x faxis(frequency) x epochs)
    
    %%% Plot Part
    
    %%% Get plot parameter
    % Get x-axis 
    group_id_for_plot = 1:21;
    
    % Get frequency lim for plot 
    [~,idxFmax] = findclosest(faxis,49);% max frequency = 50
    
    % Get frequency array for correlation
    fundamental_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array);
    IM_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, IM_frequency_array);
    fundamental_frequency_array_for_correlation = fundamental_frequency_array_for_correlation(2:idxFmax, :); % This avoids assigning 1 to the element corresponding more than 50Hz
    IM_frequency_array_for_correlation = IM_frequency_array_for_correlation(2:idxFmax, :);
    fundamental_IM_array_for_correlation = fundamental_frequency_array_for_correlation + IM_frequency_array_for_correlation;
    
    % Get frequency arrau for image
    % Background white, fundamental Black, IM 
    fundamental_IM_array_for_image = get_frequency_array_for_image(fundamental_frequency_array_for_correlation, IM_frequency_array_for_correlation);
    
    % Initialisation 
    correlation_array_fundamental = zeros(size(logSNR, 1), 2); % dim = number of bipolar channel id 
    correlation_array_IM = zeros(size(logSNR, 1), 2);
    correlation_array_fundamental_IM = zeros(size(logSNR, 1), 2);
    
    for bipolar_ch_id = 1:size(logSNR,1)
        for data_type=0:1
            
            if data_type == 0
                logdata_this_ch = squeeze(logPow(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                fig_file_prefix = 'spectrogram_logPower_';
                
            elseif data_type == 1
                logdata_this_ch = squeeze(logSNR(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                fig_file_prefix = 'spectrogram_logSNR_';
            end
            
            %%% Plot spectrogram
            fig_name = strcat(file_name_prefix, ' ', string(bipolar_ch_id), 'ch');
            spectrogram = plot_spectrogram(fig_name, logdata_this_ch(2:idxFmax, :), faxis(2:idxFmax), group_id_for_plot, fundamental_IM_array_for_image, data_type, finger1, finger2, sweep_type);
            saveas(spectrogram, strcat(path_save_dir_exp3, fig_file_prefix, file_name_prefix, '_', num2str(bipolar_ch_id) ,'ch.png'));
            %saveas(spectrogram, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, fig_file_prefix, file_name_prefix, '_', num2str(bipolar_ch_id) ,'ch.png'));
            close gcf;
            
            %%% Compute correlation between spectrogram and frequency_array_for_plot
            correlation_array_fundamental(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), fundamental_frequency_array_for_correlation); 
            correlation_array_IM(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), IM_frequency_array_for_correlation);
            correlation_array_fundamental_IM(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), fundamental_IM_array_for_correlation);
            
        end
    end
    
    % Plot correlation for bipolar channel map
    %%% Fundamental %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log Power
    grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental(:,1)),  bipolar_reref, 0, 'Fundamental');
    saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_fundamental_grid_logPower_', file_name_prefix, '.png'));
    %saveas(grid_logPow, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_grid_logPower_plus1_', file_name_prefix, '.png'));
    close gcf;
    % Log SNR
    grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental(:,2)),  bipolar_reref, 1, 'Fundamental');
    saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_fundamental_gird_logSNR_', file_name_prefix, '.png'));
    %saveas(grid_logSNR, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_gird_logSNR_from_plus1_power_', file_name_prefix, '.png'));
    close gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log Power
    grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM(:,1)),  bipolar_reref, 0, 'IM');
    saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM_grid_logPower_', file_name_prefix, '.png'));
    %saveas(grid_logPow, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_grid_logPower_plus1_', file_name_prefix, '.png'));
    close gcf;
    % Log SNR
    grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM(:,2)),  bipolar_reref, 1, 'IM');
    saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM_gird_logSNR_', file_name_prefix, '.png'));
    %saveas(grid_logSNR, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_gird_logSNR_from_plus1_power_', file_name_prefix, '.png'));
    close gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fundamental  +IM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log Power
    grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental_IM(:,1)),  bipolar_reref, 0, 'Fundamenta + IMl');
    saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_fundamental_IM_grid_logPower_', file_name_prefix, '.png'));
    %saveas(grid_logPow, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_grid_logPower_plus1_', file_name_prefix, '.png'));
    close gcf;
    % Log SNR
    grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental_IM(:,2)),  bipolar_reref, 1, 'Fundamental + IM');
    saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_fundamental_IM_gird_logSNR_', file_name_prefix, '.png'));
    %saveas(grid_logSNR, strcat(path_save_dir_exp3, filesep, 'power_plus_one', filesep, 'corr_gird_logSNR_from_plus1_power_', file_name_prefix, '.png'));
    close gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %{
    %%% Plot full spectrum
    % Use all trials
    trials_all = ones(length(epochs(:,1)),1);
    trials_all = logical(trials_all);
    plot_spectrum(file_name_prefix, logSNR, logPow, trials_all, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'spectrum_', file_name_prefix, '.png'])
    close gcf;

    % Plot grid logSNR values
    % Plot grid logSNR values
    plot_grid_logSNR_values_vs0(file_name_prefix, logSNR, trials_all, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'grid_logSNR_', file_name_prefix, '.png'])
    close gcf;
    
    % Separately use trials in which subject attend left land or right hand
    % Attend left hand (contralateral); groupid = 1
    trials_left = (epochs(:,1) ==  1);
    plot_spectrum([file_name_prefix, '_attend_left'], logSNR, logPow, trials_left, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, 'spectrum_', file_name_prefix, '_attend_left.png'])
    close gcf;
    
    % Plot grid logSNR values
    plot_grid_logSNR_values_vs0([file_name_prefix, '_attend_left'], logSNR, trials_left, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'grid_logSNR_', file_name_prefix, '_attend_left.png'])
    close gcf;
    
    % Attend right hand (ipsilateral); groupid = 3
    trials_right = (epochs(:,1) == 3);
    plot_spectrum([file_name_prefix, '_attend_right'], logSNR, logPow, trials_right, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'spectrum_', file_name_prefix, '_attend_right.png'])
    close gcf;

    % Plot grid logSNR values
    plot_grid_logSNR_values_vs0([file_name_prefix, '_attend_right'], logSNR, trials_right, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'grid_logSNR_', file_name_prefix, '_attend_right.png'])
    close gcf;
    
    % Plot grid logSNR values left hand trials vs right hand trials
    plot_grid_logSNR_values_vs_trials([file_name_prefix, '_attend_left_vs_attend_right'], logSNR, trials_left, trials_right, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis);
    saveas(gcf, [path_save_dir_exp1, filesep, 'grid_logSNR_', file_name_prefix, '_attend_left_vs_attend_right.png'])
    close gcf;
    
    % Plot grid differences of logSNR values left hand trials vs right hand
    % trials (Just subtract t-values)
    plot_grid_diff_t_values_of_logSNR_vs_trials([file_name_prefix, '_attend_left_vs_attend_right'], logSNR, trials_left, trials_right, bipolar_reref, LeftF0, IM, RightF0, RightIM, faxis)
    saveas(gcf, [path_save_dir_exp1, filesep, 'grid_diff_t_values_of_logSNR_', file_name_prefix, '_attend_left_vs_attend_right.png'])
    close gcf;
    %}
    
end

