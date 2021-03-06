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
    
    for frequency_id = 1:2
        % Get target frequency 
        target_frequency_array = nan(size(frequency_array));
        target_frequency_array(:,1) = frequency_array(:, frequency_id);
        
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
        fprintf('... extracting log(Pow) and log(SNR) by trial\n')
        snrparam=[];
        snrparam.method='fft'; % use FFT
        snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)
        
        % Legend initialisation
        legend_matrix = frequency_array(:,frequency_id);
        legend_text =strings(1, length(legend_matrix));

        % Create figure
        fig1 = figure('Name','logSNR_RESS_components','NumberTitle','off'); 
        hold on;
        
        for f_id=1:size(target_frequency_array,1)
            target_frequency_now = target_frequency_array(f_id,1);
            paramRESS=[];
            paramRESS.peakfreq1=target_frequency_now;

            paramRESS.peakwidt=0.4; % FWHM at peak frequency
            paramRESS.neighfreq=0.6;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
            paramRESS.neighwidt=0.4; % FWHM of the neighboring frequencies
            paramRESS.fft_res=0.1; % Resolution of FFT (default 0.1 Hz)
            paramRESS.mindist=0.5; % Min distance between peaks (default 0.5Hz)
            paramRESS.lowerfreq=2; % Lower frequency bound (defaul 2Hz)
            sorted_data_now(:,:,1)=sorted_data(:,:,f_id);
            [snrR, snrE, faxis2, maps, components]=get_logSNR_RESS(sorted_data_now,SR,paramRESS);

            snrparam=[];
            snrparam.method='fft'; % use FFT
            snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)
            RESS_comp(1,:,:)=components;
            [logSNR_comp, faxis_comp, logpow_comp]=get_logSNR(RESS_comp,SR,snrparam);

            baseline_logSNR_RESS(f_id,:)=logSNR_comp;
            baseline_Comp_RESS(f_id,:)=components;
            baseline_Maps_RESS(f_id,:)=maps;
            
            %%% Plot 
            plot(faxis_comp, logSNR_comp);
            xlim([0, 50]);

            legend_text(1, f_id) = [num2str(legend_matrix(f_id)) ,'Hz'];
        end % f_id

        title(['(f1,f2) = ', finger1, ',', finger2, ')', ', sweep type : ', sweep_name, ', RESS target : ', 'f', num2str(frequency_id)]);
        legend(legend_text);
        %{
        % Add Stimulus line
        for f_id=1:size(target_frequency_array,1) 
            % Add target line
            xline(legend_matrix(f_id),'--k');
        end % f_id 
        %}
        hold off;
        
        
        
    end % frequency_id
end % nversion


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


