%% Init
clear all;
close all;

% EXPERIMENT 1 on SUBJECT 3
% Stimulation on both hands with attentional focus on left or right hand
% 2 versions (same experiment but counterbalanced frequency mapping)
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
% % Frequency Mapping - Normal:
% % Left Thumb: 8 Hz
% % Left Index: 11 Hz
% % Right Thumb: 17 Hz
% % Right Index: 23 Hz
% %
% % Frequency Mapping - Reverse:
% % Left Thumb: 23 Hz
% % Left Index: 17 Hz
% % Right Thumb: 11 Hz
% % Right Index: 8 Hz
% %
% % Group IDs:
% % 1 ? Stim, Attention Left
% % 2 ? Rest (after Attention Left)
% % 3 ? Stim, Attention Right
% % 4 ? Rest (after Attention Right)

% prepare path
localdef;
path_exp1data=[path_data filesep 'S3' filesep 'Exp1'];
addpath(genpath(path_LSCPtools));
addpath(genpath(path_chronux));

% channels of interest
S1_channels = (2:21);
tot_channels = 100;
cond_channels = 106;
load(path_bipolar_reref);

duration_stim = 11; % in seconds

% preprocessing steps
preproc=[];
preproc.bipolar=1;
preproc.avref=0; % make sure you pick on of the two!
preproc.notch=1;


%% Looping on the different versions of the experiment
for nversion=1:2
    if nversion==1
        file_name='Exp1_FrequencyNormal_20_09_2018_14_26_01_0000.mat';
        LeftF0=[8 11]; % frequencies applied to the left hand
        RightF0=[17 23]; % frequencies applied to the right hand
    elseif nversion==2
        file_name= 'Exp1_FrequencyReverse_20_09_2018_14_51_35_0000.mat';
        LeftF0=[23 17];
        RightF0=[11 8];
    end
    LeftIM=[abs(LeftF0(2)-LeftF0(1)) LeftF0(2)+LeftF0(1)];
    RightIM=[abs(RightF0(2)-RightF0(1)) RightF0(2)+RightF0(1)];
    fprintf('... loading %s\n',file_name)
    load([path_exp1data filesep file_name])
    
    % Identify trials (trial condition and start/end)
    fprintf('... finding epochs\n')
    groupid=squeeze(y(cond_channels,1,:));
    epochs=[];
    for ngroup=[1 3]
        temp_start=(groupid==ngroup);
        temp_end=(groupid==ngroup+1);
        start=find(diff(temp_start)==1)+1;
        ending=find(diff(temp_end)==1)+1;
        duration=ending-start;
        durationsec = (duration/SR);
        fprintf('... ... cond %g ... found %g epochs of length %2.2fs on average\n',ngroup,length(start), mean(durationsec))
        % matrix of epochs: condition; stimulus onset (index); epoch end; epoch start (onset -1s)
        epochs=[epochs ; [ngroup*ones(length(start),1) start start+duration_stim*SR start-SR]];
    end
    
    % Re-reference (bi-polar referencing of the 60 electrodes of the 1st grid)
    if preproc.bipolar
        fprintf('... bi-polar referencing on somatosensory channels\n')
        y=y(S1_channels,:,:);
        refy=nan(size(bipolar_reref,1),1,size(y,3));
        for nch=1:size(bipolar_reref,1)
            refy(nch,1,:)=y(bipolar_reref(nch,1),1,:)-y(bipolar_reref(nch,2),1,:);
        end
    elseif preproc.avref % reference to the average
        fprintf('... average referencing on somatosensory channels\n')
        y=y(S1_channels,:,:);
        refy=y-repmat(mean(y,1),[size(y,1) 1 1]);
    else % no re-referencing
        fprintf('... no referencing and keep all channels\n')
        refy=y(2:tot_channels,:,:);
    end
    
    % filtering
    if preproc.notch
        fprintf('... 50Hz notch-filter using chronux\n')
        data=squeeze(refy)'; % time*channel
        movingwin=[14 7]; % window length and step
        tau=10; %
        params=[];
        params.tapers=[3 4]; %params has the following fields: tapers, Fs, fpass, pad
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
    times=-1:1/SR:duration_stim;
    erp=nan(size(dataf,1),length(times),size(epochs,1));
    for ntr=1:size(epochs,1)
        erp(:,:,ntr)=squeeze(dataf(:,epochs(ntr,4):epochs(ntr,3)));
    end
    erp=erp-repmat(mean(erp(:,times<0,:),2),[1 length(times) 1]);
    
    % Compute log(SNR)
    fprintf('... extracting log(Pow) and log(SNR) by trial\n')
    snrparam=[];
    snrparam.method='fft'; % use FFT
    snrparam.mindist=0.9; % minimal expected distance between peak (in Hz)
    [logSNR, faxis, logPow]=get_logSNR(erp,SR,snrparam);
    
    %%% Plot full spectrum
    figure;
    set(gcf,'Name',file_name)
    subplot(1,2,1); 
    plot(faxis,mean(logPow(bipolar_reref(:,3)==0,:,:),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logPow(bipolar_reref(:,3)==1,:,:),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(Power) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;
    
    subplot(1,2,2); 
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==0,:,:),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==1,:,:),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(SNR) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;

    % mark fundamentals
    subplot(1,2,2); hold on;
    lgd=[];
    for nF=1:length(LeftF0)
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        lgd(1)=scatter(thisF,max(mean(logSNR(:,idxF,:),3),[],1),'Marker','o','MarkerEdgeColor','r','SizeData',36);
    end
    for nF=1:length(RightF0)
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        lgd(2)=scatter(thisF,max(mean(logSNR(:,idxF,:),3),[],1),'Marker','o','MarkerEdgeColor','b','SizeData',36);
    end
        format_fig;

    % mark IM
    subplot(1,2,2); hold on;
    for nF=1:length(LeftIM)
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        lgd(3)=scatter(thisF,max(mean(logSNR(:,idxF,:),3),[],1),'Marker','*','MarkerEdgeColor','r','SizeData',36);
    end
    for nF=1:length(RightIM)
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        lgd(4)=scatter(thisF,max(mean(logSNR(:,idxF,:),3),[],1),'Marker','*','MarkerEdgeColor','b','SizeData',36);
    end
    legend(lgd,{'contraF','ipsiF','contraIM','ipsiIM'})
        format_fig;

    % Plot grid logSNR values
    figure;
    set(gcf,'Name',file_name)
    cmap=colormap('hot');
    cmap=flipud(cmap);
    for nF=1:length(LeftF0)
        subplot(2,4,nF);
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,:)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraF (%g)',LeftF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(LeftIM)
        subplot(2,4,nF+2);
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,:)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraIM (%g)',LeftIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
        for nF=1:length(RightF0)
        subplot(2,4,nF+4);
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,:)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiF (%g)',RightF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightIM)
        subplot(2,4,nF+6);
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,:)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiIM (%g)',RightIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
end