%% 
% Plot relation between SNR vs Intermodulation of all bipolar-channels
% on the same fig for each experiment 
% Unit db --- 10log10(x)

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
addpath(genpath(path_functions));
path_save_dir_exp3=[path_save_dir, filesep, 'Exp3', filesep, 'IM_multiple_range_vs_SNR_condition', filesep];
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

% Frequency range 
frequency_range = [[0,100];[100,200];[200,300];[300,400]];  
trial_array = cat(2,1:21,1:21,1:21, 1:21, 1:21, 1:21);

% Create coefficience set for IM 
nOrder = 8; 
[coefficience_set, coefficience_set_start_id] = create_coefficience(nOrder);
nHarmonic = 8;

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
            sweep_name = 'Ascending';
        case 'D'
            frequency_array = descending_frequency;
            sweep_name = 'Descending';
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
    
    % Load data
    fprintf('... loading %s\n',data_file)
    load([path_exp3data filesep data_file])
        
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

    %%% Compute log(SNR) and plot part 
    fprintf('... extracting 10log10(Pow) and 10log10(SNR) by trial\n')
    snrparam=[];
    snrparam.method='fft'; % use FFT
    snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)  
    [SNRdB, faxis, PowdB]=get_SNRdB(sorted_data, SR, snrparam); % logSNR/logPow --dim (channels x faxis(frequency) x epochs)
    
    % Get target frequency 
    f1_array = frequency_array(:,1);
    f2_array = frequency_array(:,2);
    
    % Create IM array (1st~8th order) for a pair of (f1, f2) and store the
    % corresponding data(Pow or SNR)
    % Store the IM array using structure. 
    % struct(i).ims has ith-order ims
    % struct(i).power has ith-order ims' data (dim : num. of ims x num. of bipolar_ch)
    multiple_im_struct = struct();
    
    for trial_id = 1:size(frequency_array, 1)
        f1 = f1_array(trial_id);
        f2 = f2_array(trial_id);
        
        start_id_for_im_array = zeros(nOrder, 1);
        
        for i_order = 1:nOrder
            
            start_id_for_coeff = coefficience_set_start_id(i_order);
            
            if i_order ~= nOrder 
                end_id_for_coeff = coefficience_set_start_id(i_order + 1) -1;
            else
                end_id_for_coeff = size(coefficience_set, 1);
            end
            
            % Initialisation
            ims = zeros(size(start_id_for_coeff:end_id_for_coeff,1)*2, 1);
            
            counter_im = 1;
            
            if i_order == 1
                for i_set = start_id_for_coeff:end_id_for_coeff 
                    set_now = coefficience_set(i_set,:);
                    im = f1*set_now(1) + f2*set_now(2);
                    ims(counter_im) = im;
                    % Find closest point to current IM
                    [~, frequency_id_in_array] = findclosest(faxis, im);
                    % Get data
                    current_logpow = squeeze(PowdB(:, frequency_id_in_array, trial_id))'; 
                    current_logpow = (current_logpow / 10); % Change unit of Power
                    current_logSNR= squeeze(SNRdB(:, frequency_id_in_array, trial_id))';
                    
                    if trial_id == 1 && i_set == start_id_for_coeff
                        multiple_im_struct(i_order).log10power = current_logpow;
                        multiple_im_struct(i_order).SNRdB = current_logSNR;
                    else
                        temp_Pow = multiple_im_struct(i_order).log10power;
                        temp_SNR = multiple_im_struct(i_order).SNRdB;
                        multiple_im_struct(i_order).log10power = cat(1, temp_Pow, current_logpow);
                        multiple_im_struct(i_order).SNRdB = cat(1, temp_SNR, current_logSNR);
                    end
                    counter_im = counter_im + 1;
                end % i_set
            else
                for i_set = start_id_for_coeff:end_id_for_coeff 
                    set_now = coefficience_set(i_set,:);
                    for sign = 1:2
                        if sign == 1
                            im = f1*set_now(1) + f2*set_now(2);
                            ims(counter_im) = im; 
                            % Find closest point to current IM
                            [~, frequency_id_in_array] = findclosest(faxis, im);
                            % Get data
                            current_logpow = squeeze(PowdB(:, frequency_id_in_array, trial_id))'; 
                            current_logpow = (current_logpow / 10); % Change unit of Power
                            current_logSNR= squeeze(SNRdB(:, frequency_id_in_array, trial_id))';
                            
                            if all([trial_id == 1, i_set == start_id_for_coeff, sign == 1])
                                multiple_im_struct(i_order).log10power = current_logpow;
                                multiple_im_struct(i_order).SNRdB = current_logSNR;
                            else
                                temp_Pow = multiple_im_struct(i_order).log10power;
                                temp_SNR = multiple_im_struct(i_order).SNRdB;
                                multiple_im_struct(i_order).log10power = cat(1, temp_Pow, current_logpow);
                                multiple_im_struct(i_order).SNRdB = cat(1, temp_SNR, current_logSNR);
                            end
                        elseif sign == 2
                            im = abs(f1*set_now(1) - f2*set_now(2));
                            ims(counter_im) = im; 
                            % Find closest point to current IM
                            current_logpow = squeeze(PowdB(:, frequency_id_in_array, trial_id))'; 
                            current_logpow = (current_logpow / 10); % Change unit of Power
                            current_SNR= squeeze(SNRdB(:, frequency_id_in_array, trial_id))';
                            temp_Pow = multiple_im_struct(i_order).log10power;
                            temp_SNR = multiple_im_struct(i_order).SNRdB;
                            multiple_im_struct(i_order).log10power = cat(1, temp_Pow, current_logpow);
                            multiple_im_struct(i_order).SNRdB = cat(1, temp_SNR, current_logSNR);
                        end
                        counter_im = counter_im + 1;
                        
                    end % sign
                end % i_set
            end % if i_order 
            
            if trial_id == 1
                multiple_im_struct(i_order).ims = ims; % Remove 50xHz 
            else 
                multiple_im_struct(i_order).ims = cat(1, multiple_im_struct(i_order).ims, ims); % Remove 50xHz
            end 
            
        end % i_order
        
    end % trial_id (The end of IM array part)
    
    
    %%% Remove 50xHz
    for i_order = 1:nOrder
        current_ims = multiple_im_struct(i_order).ims;
        ids_without_50Hz = mod(current_ims, 50) ~= 0;
        multiple_im_struct(i_order).ims = current_ims(ids_without_50Hz);
        multiple_im_struct(i_order).log10power = multiple_im_struct(i_order).log10power(ids_without_50Hz,:); %dim num. of ims x num. of bipolar_ch
        multiple_im_struct(i_order).SNRdB = multiple_im_struct(i_order).SNRdB(ids_without_50Hz,:); %dim num. of ims x num. of bipolar_ch
    end
    
    %%% Plot
    % Data storage
    store_array_mean_for_plot_multiple = zeros(nOrder, size(SNRdB, 1)); % dim = range x num. of bipolar-channel
    store_array_std_for_plot_multiple = zeros(nOrder, size(SNRdB, 1));
    store_array_mean_for_plot_range = zeros(size(frequency_range,1), size(SNRdB, 1)); % dim = range x num. of bipolar-channel
    store_array_std_for_plot_range = zeros(size(frequency_range,1), size(SNRdB, 1));        
    legend_text =strings(1, size(SNRdB, 1));

    for bipolar_ch_id = 1:size(SNRdB,1)
        legend_text(1, bipolar_ch_id) = [num2str(bipolar_ch_id) ,'bp-ch']; 
        for data_type=0:1 % Power or SNR
            if data_type == 0
                logdata_this_ch = squeeze(PowdB(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                fig_file_prefix_multiple = strcat('relation_IM_multiple_PowdB_');
                fig_file_prefix_range = strcat('relation_IM_range_PowdB_');
                plot_txt = 'log10(Power) [-]';
                % Prepare a figure window for all ch plot 
                if bipolar_ch_id == 1
                    f_pow_im_multiple = figure('visible','off');
                    f_pow_im_range = figure('visible','off');
                end
            elseif data_type == 1
                logdata_this_ch = squeeze(SNRdB(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                fig_file_prefix_multiple = strcat('relation_IM_multiple_SNRdB_');
                fig_file_prefix_range = strcat('relation_IM_range_SNRdB_');
                plot_txt = '10log10(SNR) [dB]';
                % Prepare a figure window for all ch plot
                if bipolar_ch_id == 1
                    f_snr_im_multiple = figure('visible','off');
                    f_snr_im_range = figure('visible','off');
                end
            end

            %%% Plot relation frequency and logPow or log SNR
            fig_name_multiple = strcat(file_name_prefix, '_IM_multiple_', string(bipolar_ch_id), 'ch');
            fig_name_range = strcat(file_name_prefix, '_IM_range_', string(bipolar_ch_id), 'ch');

            %%% Multiple
            for ith_order = 1:nOrder
                if data_type == 0
                    store_array_mean_for_plot_multiple(ith_order, bipolar_ch_id) = mean(multiple_im_struct(ith_order).log10power(:,bipolar_ch_id));
                    store_array_std_for_plot_multiple(ith_order, bipolar_ch_id) = std(multiple_im_struct(ith_order).log10power(:,bipolar_ch_id));
                else
                    store_array_mean_for_plot_multiple(ith_order, bipolar_ch_id) = mean(multiple_im_struct(ith_order).SNRdB(:,bipolar_ch_id));
                    store_array_std_for_plot_multiple(ith_order, bipolar_ch_id) = std(multiple_im_struct(ith_order).SNRdB(:,bipolar_ch_id));
                end
            end

            % Figure harmonic for each ch  
            
            f_temp = figure('Name', fig_name_multiple, 'visible','off');
            errorbar(1:nOrder, store_array_mean_for_plot_multiple(:, bipolar_ch_id), store_array_std_for_plot_multiple(:, bipolar_ch_id), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'auto', 'MarkerFaceColor', 'auto', 'LineWidth', 3, 'CapSize', 15);
            %{
            hold on;
            for harmonic_id = 1:nHarmonic
                scatter( harmonic_id * ones(size(store_array_for_mean_harmonic, 1),1), store_array_for_mean_harmonic(:, harmonic_id), 15, 'filled', 'red');
            end
            hold off;
            %}

            set(gca,'XTick', 1:nOrder);
            xlim([1, nOrder]);
            if data_type ==0 
                ylim([-2, 7]);
            elseif data_type == 1
                ylim([-20, 25]);
            end
            xticklabels({'1^{st}','2^{nd}','3^{rd}','4^{th}','5^{th}','6^{th}', '7^{th}', '8^{th}'});   
            xlabel('n^{th}-order IM[Hz]');
            ylabel(plot_txt);
            title(['(f1,f2) = (', finger1, ',', finger2, ')', ', Sweep type : ', sweep_name, ', Target : Intermodulation (1^{st}~8^{th})'], 'FontSize', 20, 'FontWeigh', 'bold');
            format_fig;
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf, strcat(path_save_dir_exp3, fig_file_prefix_multiple, file_name_prefix, '_IM_', num2str(bipolar_ch_id) ,'ch.png'));
            close(f_temp);
            
            % Plot data of all bipolar channels on the same figure 
            % Switch figure window
            if data_type == 0
                figure(f_pow_im_multiple);
            else
                figure(f_snr_im_multiple);
            end

            subplot_id = floor((bipolar_ch_id - 1)/7) + 1; 
            hold on; 
            subplot(1,5,subplot_id);
            errorbar(1:nOrder, store_array_mean_for_plot_multiple(:, bipolar_ch_id), store_array_std_for_plot_multiple(:, bipolar_ch_id), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'auto', 'MarkerFaceColor', 'auto', 'LineWidth', 3, 'CapSize', 15);
            hold off;

            % Description for each subplot 
            if mod(bipolar_ch_id,7) == 0 || bipolar_ch_id == size(SNRdB,1)
                if bipolar_ch_id == size(SNRdB,1)
                    start_ch_id = 29;
                else
                    start_ch_id = bipolar_ch_id - 6;
                end
                title([num2str(start_ch_id), ' bp-ch ~ ', num2str(bipolar_ch_id), ' bp-ch']);
                legend(legend_text(start_ch_id:bipolar_ch_id));
                set(gca,'XTick', 1:nOrder);
                xticklabels({'1^{st}','2^{nd}','3^{rd}','4^{th}','5^{th}','6^{th}', '7^{th}', '8^{th}'});   
                xlim([1, nOrder]);
                if data_type ==0 
                    ylim([-2, 7]);
                elseif data_type == 1
                    ylim([-20, 25]);
                end
                format_fig;
            end

            % Set x-y axis
            if bipolar_ch_id == 1
                xlabel('n^{th}-order IM[Hz]');
                ylabel(plot_txt);
            end

            % Close the window when all ch data have been plotted.
            if bipolar_ch_id == size(SNRdB,1)
                axes_now = findall(gcf, 'type', 'axes');
                linkaxes(axes_now, 'xy');
                sgtitle(['(f1,f2) = (', finger1, ',', finger2, ')', ', Sweep type : ', sweep_name, ', Target : Intermodulation (1^{st}~8^{th})'], 'FontSize', 20, 'FontWeigh', 'bold');
                format_fig;
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf, strcat(path_save_dir_exp3, fig_file_prefix_multiple, file_name_prefix, '_IM_allch.png'));
                close gcf;
            end 

            %%% Range
            %
            % Plot relation between frequency in a range and log data
            for range_id = 1:size(frequency_range,1)
                min_frequency = frequency_range(range_id, 1);
                max_frequency = frequency_range(range_id, 2);
                store_array_for_mean_range = 0;
                for order_id = 1:nOrder
                    frequency_indices = find(multiple_im_struct(order_id).ims > min_frequency & (multiple_im_struct(order_id).ims > max_frequency | multiple_im_struct(order_id).ims == max_frequency));
                    if data_type == 0
                        store_array_for_mean_range= cat(1, store_array_for_mean_range, multiple_im_struct(order_id).log10power(frequency_indices, bipolar_ch_id)); 
                    else 
                        store_array_for_mean_range = cat(1, store_array_for_mean_range, multiple_im_struct(order_id).SNRdB(frequency_indices, bipolar_ch_id)); 
                    end
                end
                % Store mean and std
                store_array_mean_for_plot_range(range_id, bipolar_ch_id) = mean(store_array_for_mean_range);
                store_array_std_for_plot_range(range_id, bipolar_ch_id) = std(store_array_for_mean_range);
            end

            % Plot For each ch
            f_temp = figure('Name', fig_name_range, 'visible','off');
            errorbar(1:range_id, store_array_mean_for_plot_range(:, bipolar_ch_id), store_array_std_for_plot_range(:, bipolar_ch_id), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'auto', 'MarkerFaceColor', 'auto', 'LineWidth', 3, 'CapSize', 15)
            xlim([1,range_id]);
            if data_type ==0 
                ylim([-2, 7]);
            elseif data_type == 1
                ylim([-20, 25]);
            end
            set(gca,'XTick', 1:range_id);
            xticklabels({'0~100','100~200','200~300','300~400'});
            xlabel('Frequency range [Hz]');
            ylabel(plot_txt);
            title(['(f1,f2) = (', finger1, ',', finger2, ')', ', Sweep type : ', sweep_name, ', Target : Intermodulation (1^{st}~8^{th})'], 'FontSize', 20, 'FontWeigh', 'bold');
            format_fig;
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf, strcat(path_save_dir_exp3, fig_file_prefix_range, file_name_prefix, '_IM_', num2str(bipolar_ch_id) ,'ch.png'));
            close(f_temp);
            
            % Plot data of all bipolar channels on the same figure 
            % Switch figure window
            if data_type == 0
                figure(f_pow_im_range);
            else
                figure(f_snr_im_range);
            end

            subplot_id = floor((bipolar_ch_id - 1)/7) + 1; 
            hold on; 
            subplot(1,5,subplot_id);
            errorbar(1:range_id, store_array_mean_for_plot_range(:, bipolar_ch_id), store_array_std_for_plot_range(:, bipolar_ch_id), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'auto', 'MarkerFaceColor', 'auto', 'LineWidth', 3, 'CapSize', 15)
            hold off;

            % Description for each subplot 
            if mod(bipolar_ch_id,7) == 0 || bipolar_ch_id == size(SNRdB,1)
                if bipolar_ch_id == size(SNRdB,1)
                    start_ch_id = 29;
                else
                    start_ch_id = bipolar_ch_id - 6;
                end
                title([num2str(start_ch_id), ' bp-ch ~ ', num2str(bipolar_ch_id), ' bp-ch']);
                legend(legend_text(start_ch_id:bipolar_ch_id));
                xlim([1, range_id]);                    
                if data_type ==0 
                    ylim([-2, 7]);
                elseif data_type == 1
                    ylim([-20, 25]);
                end
                set(gca,'XTick', 1:range_id);
                %{
                l0 = ['0' newline ' ~100'];
                l1 = ['100' newline ' ~200'];
                l2 = ['200' newline ' ~300'];
                l3 = ['300' newline ' ~400'];
                xticklabels({l0, l1, l2, l3});
                %}
                xticklabels({'0~100', '100~200', '200~300', '300~400'});
                xtickangle(gca, 10)
                format_fig;
            end

            % Set x-y axis
            if bipolar_ch_id == 1
                xlabel('Frequency range [Hz]');
                ylabel(plot_txt);
            end

            % Close the window when all ch data have been plotted.
            if bipolar_ch_id == size(SNRdB,1)
                axes_now = findall(gcf, 'type', 'axes');
                linkaxes(axes_now, 'xy');
                sgtitle(['(f1,f2) = (', finger1, ',', finger2, ')', ', Sweep type : ', sweep_name, ', Target : Intermodulation (1^{st}~8^{th})'], 'FontSize', 20, 'FontWeigh', 'bold');
                format_fig;
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf, strcat(path_save_dir_exp3, fig_file_prefix_range, file_name_prefix, '_IM_allch.png'));
                close gcf;
            end 

        end %data_type
    end %bipolar_ch_id
end % nversion

%% Functions
function [coefficience_set, set_start_id] = create_coefficience(nOrder)
    set_start_id = zeros(1, nOrder);

    counter_set = 0;

    for i_order = 1:nOrder
        if i_order == 1
            coefficience_set = [0 1];
            coefficience_set = cat(1, coefficience_set, [1 0]);
            set_start_id(i_order) = 1;
            counter_set = 2;
        else
            for element_1 = 1:i_order-1
                element_2 = i_order - element_1;
                coefficience_set = cat(1, coefficience_set, [element_1 element_2]);
                counter_set = counter_set + 1;
                if element_1 == 1
                    set_start_id(i_order) = counter_set;
                end 
            end % element_1
        end % if i_order
    end % i_order
end