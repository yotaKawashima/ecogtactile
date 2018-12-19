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

%% Looping on the different sweeps, different frequencies, and different fingers
%
%    Take maen of power from two different data which share the same finger, the same
%    frequecy to the finger, and the same sweep type
%    
%    This part shows spectrogram and correlation between spectrogram and
%    stimuli matrix

fingers = {'Thumb', 'Index', 'Middle'};
for sweep_id = 1:length(sweep_types)
    % Get sweep type
    sweep_type = sweep_types{1, sweep_id};
    % Get frequency array for this version
    switch sweep_type 
        case 'A'
            frequency_array = ascending_frequency;
        case 'D'
            frequency_array = descending_frequency;
    end
    
    % Get Intermodulation frequency array
    % Intermodulation array column=1 -> subtraction,  column=2 -> addition
    IM_frequency_array = zeros(size(frequency_array));
    IM_frequency_array(:,1) = abs(frequency_array(:,1) - frequency_array(:,2));
    IM_frequency_array(:,2) = abs(frequency_array(:,1) + frequency_array(:,2));
    
    IM2_frequency_array = 2 * IM_frequency_array;
    
    IM3_frequency_array = 3 * IM_frequency_array;
    
    IM4_frequency_array = 4 * IM_frequency_array;
    
    for frequency_id = 1:2
        % Get target frequency 
        target_frequency_array = nan(size(frequency_array));
        target_frequency_array(:,1) = frequency_array(:, frequency_id);
        
        for finger_id = 1:3
            % Get finger
            finger = fingers{1, finger_id};
            
            % Get file name prefix for the same finger, the same frequency, and the same sweep type
            [file_name1_prefix, file_name2_prefix] = get_two_file_name_prefix(file_name_prefixs, sweep_id, frequency_id, finger);
            % Get file name 
            file_name_prefix = strcat('Exp3_', sweep_type, '_', finger, '_f', num2str(frequency_id));
            % For data1
            file_name1_temp = join([path_exp3data filesep file_name1_prefix, '*.mat'],''); 
            file_name1_info = dir(file_name1_temp{1});  % Here add {1}
            data1_file = file_name1_info.name;
            % For data2 
            file_name2_temp = join([path_exp3data filesep file_name2_prefix, '*.mat'],'');
            file_name2_info = dir(file_name2_temp{1});  % Here add {1}
            data2_file = file_name2_info.name;
            
            % Load data1 and data2
            fprintf('... loading %s\n',data1_file)
            fprintf('... loading %s\n',data2_file)
            data1 = load([path_exp3data filesep data1_file]);
            data2 = load([path_exp3data filesep data2_file]);
            SR = data1.SR;
            
            % Get groupid (frequency condition on each trial)
            fprintf('... finding epochs\n')
            groupid = squeeze(data1.y(cond_channels,1,:));
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
            
            % Preprocess data
            for data_id = 1:2
                
                % Change variable name
                switch data_id
                    case 1
                        y = data1.y;
                    case 2 
                        y = data2.y;
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
                    
                    % This following part would be problematic. %%%%%%%
                    %{
                    f1=100;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f1);
                    f2=150;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f2);
                    f3=200;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f3);
                    f4=250;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f4);
                    f5=300;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f5);
                    f6=350;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f6);
                    f7=400;
                    [dataf,~,~,~]=rmlinesmovingwinc(dataf,movingwin,tau,params,[],[],f7);
                    %}
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    dataf=dataf'; % channel*time
                else
                    dataf=squeeze(refy);
                end
                
                % Change variable name
                switch data_id
                    case 1
                        data1f = dataf;
                    case 2 
                        data2f = dataf;
                end
            end

            % Epoch Data
            fprintf('... epoching data (-1:11s)\n')
            times= 0:1/SR:duration_stim;
            % Initialisation
            sorted_data1=nan(size(data1f,1),length(times),size(epochs,1)); % dim (channels x sample points x epochs)
            sorted_data2=sorted_data1;
            
            for ntr=1:size(epochs,1)
                sorted_data1(:,:,ntr)=squeeze(data1f(:,epochs(ntr,4):epochs(ntr,3))); % epochs(:,4) - start index (sample points), epochs(:, 3) - end index (sample points)
                sorted_data2(:,:,ntr)=squeeze(data2f(:,epochs(ntr,4):epochs(ntr,3))); 
            end

            % Remove erp part for experiment 3
            % erp=sorted_data-repmat(mean(sorted_data(:,times<0,:),2),[1 length(times) 1]); % Computation 

            % Compute log(SNR)
            fprintf('... extracting log(Pow) and log(SNR) by trial\n')
            snrparam=[];
            snrparam.method='fft'; % use FFT
            snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)
            [logSNR, faxis, logPow]=get_logSNR_mean_two_power(sorted_data1, sorted_data2,SR,snrparam); % logSNR/logPow --dim (channels x faxis(frequency) x epochs)
            
            %%% Plot Part

            %%% Get plot parameter
            % Get x-axis 
            group_id_for_plot = 1:21;

            % Get frequency lim for plot 
            [~,idxFmax] = findclosest(faxis,max_frequency);% max frequency = 49
            
            % Get index for 50Hz to remove
            [~, idxF50] = findclosest(faxis, 50);
            [~, idxF100] = findclosest(faxis, 100);
            [~, idxF150] = findclosest(faxis, 150);
            [~, idxF200] = findclosest(faxis, 200);
            [~, idxF250] = findclosest(faxis, 250);
            [~, idxF300] = findclosest(faxis, 300);
            [~, idxF350] = findclosest(faxis, 350);
            [~, idxF400] = findclosest(faxis, 400);
            
            % Get frequency array for correlation
            target_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array);
            target2_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*2);
            target3_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*3);
            target4_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*4);
            target5_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*5);
            target6_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*6);
            target7_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*7);
            target8_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*8);
            target9_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*9);
            target10_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, target_frequency_array*10);
       
        
            fundamental_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array);
            fundamental2_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*2);
            fundamental3_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*3);
            fundamental4_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*4);
            fundamental5_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*5);
            fundamental6_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*6);
            fundamental7_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*7);
            fundamental8_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*8);
            fundamental9_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*9);
            fundamental10_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, frequency_array*10);
            
            IM1_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, IM_frequency_array);
            IM2_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, IM2_frequency_array);
            IM3_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, IM3_frequency_array);
            IM4_frequency_array_for_correlation = get_frequency_array_for_correlation(faxis, IM4_frequency_array);
            
            % Add harmonic frequency
            target_frequency_array_for_correlation = target_frequency_array_for_correlation + target2_frequency_array_for_correlation + target3_frequency_array_for_correlation + target4_frequency_array_for_correlation + target5_frequency_array_for_correlation + target6_frequency_array_for_correlation + target7_frequency_array_for_correlation + target8_frequency_array_for_correlation + target9_frequency_array_for_correlation + target10_frequency_array_for_correlation;
            fundamental_frequency_array_for_correlation = fundamental_frequency_array_for_correlation + fundamental2_frequency_array_for_correlation + fundamental3_frequency_array_for_correlation + fundamental4_frequency_array_for_correlation + fundamental5_frequency_array_for_correlation + fundamental6_frequency_array_for_correlation + fundamental7_frequency_array_for_correlation + fundamental8_frequency_array_for_correlation + fundamental9_frequency_array_for_correlation + fundamental10_frequency_array_for_correlation;
            
            target_frequency_array_for_correlation = target_frequency_array_for_correlation(2:idxFmax, :);          % This avoids assigning 1 to the element corresponding more than 50Hz
            fundamental_frequency_array_for_correlation = fundamental_frequency_array_for_correlation(2:idxFmax, :); 
            IM1_frequency_array_for_correlation = IM1_frequency_array_for_correlation(2:idxFmax, :);
            IM2_frequency_array_for_correlation = IM2_frequency_array_for_correlation(2:idxFmax, :);
            IM3_frequency_array_for_correlation = IM3_frequency_array_for_correlation(2:idxFmax, :);
            IM4_frequency_array_for_correlation = IM4_frequency_array_for_correlation(2:idxFmax, :);
            
            all_IM_frequency_array_for_correlation = IM1_frequency_array_for_correlation + IM2_frequency_array_for_correlation + IM3_frequency_array_for_correlation + IM4_frequency_array_for_correlation;
            fundamental_IM_array_for_correlation = fundamental_frequency_array_for_correlation + all_IM_frequency_array_for_correlation;
            
            % Get frequency arrau for image
            % Background white, fundamental Black, IM 
            fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger(target_frequency_array_for_correlation, fundamental_frequency_array_for_correlation);
            %fundamental_IM_array_for_image = get_frequency_array_for_image_for_one_finger_4IM(target_frequency_array_for_correlation, fundamental_frequency_array_for_correlation, IM1_frequency_array_for_correlation, IM2_frequency_array_for_correlation, IM3_frequency_array_for_correlation, IM4_frequency_array_for_correlation);
            % Initialisation 
            correlation_array_target_fundamental = zeros(size(logSNR, 1), 2); % dim = number of bipolar channel id 
            correlation_array_IM1 = zeros(size(logSNR, 1), 2);
            correlation_array_IM2 = zeros(size(logSNR, 1), 2);
            correlation_array_IM3 = zeros(size(logSNR, 1), 2);
            correlation_array_IM4 = zeros(size(logSNR, 1), 2);
            correlation_array_IM1234 = zeros(size(logSNR, 1), 2);
            correlation_array_fundamental_IM = zeros(size(logSNR, 1), 2);
            
            for bipolar_ch_id = 1:size(logSNR,1)
                for data_type=0:1

                    if data_type == 0
                        logdata_this_ch = squeeze(logPow(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                        fig_file_prefix = strcat('spectrogram_logPower_');
                        %fig_file_prefix = 'spectrogram_logPower_plus1_';
                    elseif data_type == 1
                        logdata_this_ch = squeeze(logSNR(bipolar_ch_id,:,:)); % dim = faxis(frequecny) x epochs
                        fig_file_prefix = strcat('spectrogram_logSNR_');
                        %fig_file_prefix = 'spectrogram_logSNR_from_plus1_power_';
                    end
                    
                    %%% Remove value at 50Hz+-1Hz multiples  (faxis - 1 shift ~0.25Hz)
                    logdata_this_ch(idxF50-2:idxF50+2, :) = nan;
                    logdata_this_ch(idxF100-2:idxF100+2, :) = nan;
                    logdata_this_ch(idxF150-2:idxF150+2, :) = nan;
                    logdata_this_ch(idxF200-2:idxF200+2, :) = nan;
                    logdata_this_ch(idxF250-2:idxF250+2, :) = nan;
                    logdata_this_ch(idxF300-2:idxF300+2, :) = nan;
                    logdata_this_ch(idxF350-2:idxF350+2, :) = nan;
                    logdata_this_ch(idxF400-2:idxF400+2, :) = nan;
                    
                    %%% Plot spectrogram
                    fig_name = strcat(file_name_prefix, '_', string(bipolar_ch_id), 'ch');
                    spectrogram = plot_spectrogram_for_one_finger(fig_name, logdata_this_ch(2:idxFmax, :), faxis(2:idxFmax), group_id_for_plot, fundamental_IM_array_for_image, data_type, finger, sweep_type, frequency_id);
                    saveas(spectrogram, strcat(path_save_dir_exp3, fig_file_prefix, file_name_prefix, '_', num2str(bipolar_ch_id) ,'ch.png'));
                    close gcf;

                    %%% Compute correlation between spectrogram and frequency_array_for_plot
                    correlation_array_target_fundamental(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), target_frequency_array_for_correlation); 
                    correlation_array_IM1(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), IM1_frequency_array_for_correlation);
                    correlation_array_IM2(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), IM2_frequency_array_for_correlation);
                    correlation_array_IM3(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), IM3_frequency_array_for_correlation);
                    correlation_array_IM4(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), IM4_frequency_array_for_correlation);
                    correlation_array_IM1234(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), all_IM_frequency_array_for_correlation);
                    correlation_array_fundamental_IM(bipolar_ch_id, data_type+1) = corr2(logdata_this_ch(2:idxFmax, :), fundamental_IM_array_for_correlation);
                    
                end
            end
            %{
            % Plot correlation for bipolar channel map
            %%% Fundamental %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_target_fundamental(:,1)),  bipolar_reref, 0, 'Fundamental');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_fundamental_grid_logPower_', file_name_prefix, '.png'));
            %grid_logPow = plot_grid_correlation_YV(file_name_prefix, squeeze(correlation_array_target_fundamental(:,1)),  bipolar_reref, 0, 'Fundamental');
            %saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_fundamental_grid_logPower_', file_name_prefix, '_YV.png'));
            close gcf;
            
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_target_fundamental(:,2)),  bipolar_reref, 1, 'Fundamental');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_fundamental_gird_logSNR_', file_name_prefix, '.png'));
            %grid_logSNR = plot_grid_correlation_YV(file_name_prefix, squeeze(correlation_array_target_fundamental(:,2)),  bipolar_reref, 1, 'Fundamental');
            %saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_fundamental_gird_logSNR_', file_name_prefix, '_YV.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% IM1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM1(:,1)),  bipolar_reref, 0, 'IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM1_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM1(:,2)),  bipolar_reref, 1, 'IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM1_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% IM2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM2(:,1)),  bipolar_reref, 0, 'IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM2_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM2(:,2)),  bipolar_reref, 1, 'IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM2_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% IM3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM3(:,1)),  bipolar_reref, 0, 'IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM3_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM3(:,2)),  bipolar_reref, 1, 'IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM3_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% IM4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM4(:,1)),  bipolar_reref, 0, 'IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM4_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM4(:,2)),  bipolar_reref, 1, 'IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM4_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% IM1234 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM1234(:,1)),  bipolar_reref, 0, 'IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_IM1234_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_IM1234(:,2)),  bipolar_reref, 1, 'IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_IM1234_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Fundamental  +IM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Log Power
            grid_logPow = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental_IM(:,1)),  bipolar_reref, 0, 'Fundamental + IM');
            saveas(grid_logPow, strcat(path_save_dir_exp3, 'corr_fundamental_IM1234_grid_logPower_', file_name_prefix, '.png'));
            close gcf;
            % Log SNR
            grid_logSNR = plot_grid_correlation(file_name_prefix, squeeze(correlation_array_fundamental_IM(:,2)),  bipolar_reref, 1, 'Fundamental + IM');
            saveas(grid_logSNR, strcat(path_save_dir_exp3, 'corr_fundamental_IM1234_gird_logSNR_', file_name_prefix, '.png'));
            close gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
        end
    end   
end



