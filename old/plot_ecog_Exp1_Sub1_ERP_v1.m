%% Init
clear all;
close all;

% Organization of the data for  Experiment 1a
% % % CH1: time
% % % CH2-161: ECoG 1-160
% % % CH162: DI (unused)
% % % CH163: StimCode
% % % CH164: GroupId
% % % GroupIds:
% % % 1 ... Left Stim (instruction to attend left hand)
% % % 2 ... Left Rest (rest)
% % % 3 ... Right Stim (instruction to attend right hand)
% % % 4 ... Right Rest (rest)

rootFile    = ['/Users/tand0009/'];
filePath    = [rootFile 'WorkGit/LSCPtools/'];
preDataPath = '/Users/tand0009/Data/ECoG_tactile/S1/data_IM/';
chronuxPath = [rootFile 'Work/local/chronux_2_12/'];

addpath(genpath(filePath));
allversions={'1a','1arev','1b','1brev','2a','2arev','2b'};
load modified_ttestcolortable.mat

%% Looping on the different versions of the experiment
ntrc=0;
    erp_nofilt=[]; %nan(60,5,10,3*SR+1);
    erp_notch=[]; %nan(60,5,10,3*SR+1);
    erp_bp=[]; %nan(60,5,10,3*SR+1);
    erp_bp2=[]; %nan(60,5,10,3*SR+1);
for nversion=1:length(allversions) %was 4:5
    version=allversions{nversion};
    if strcmp(version,'1a') % version 1a: sensors on both hands
        data_path= [preDataPath];
        file_name='1a/Vibrotactile_Regular_16_03_2018_19_05_12_0000.mat';
    elseif strcmp(version,'1arev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath];
        file_name= '1arev/Vibrotactile_Regular_22_03_2018_11_40_37_0000.mat';
    elseif strcmp(version,'1b') % version 1a: sensors on both hands
        data_path= [preDataPath];
        file_name='1b/Vibrotactile_Crossed_16_03_2018_19_14_05_0000.mat';
    elseif strcmp(version,'1brev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath];
        file_name= '1brev/Vibrotactile_Crossed_22_03_2018_11_49_18_0000.mat';
    elseif strcmp(version,'2a') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath];
        file_name= '2a/Vibrotactile_FeetClose_17_03_2018_13_11_09_0000.mat';
    elseif strcmp(version,'2arev') % version 1a: sensors on both hands
        data_path= [preDataPath];
        file_name='2arev/Vibrotactile_FeetClose_22_03_2018_17_39_35_0000.mat';
    elseif strcmp(version,'2b') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath];
        file_name= '2b/Vibrotactile_FeetApart_17_03_2018_13_22_10_0000.mat';
    end
    fprintf('... loading %s\n',file_name) %reverse the order sd
    load([data_path filesep file_name])
    
    %% Identify trials (trial condition and start/end)
    groups=1:4; % 1 to 5 are digits from thumb to pinkie and 6 is rest
    groupid=squeeze(y(164,1,:));
    epochs=[];
    for ngroup=[1 3]
        temp=groupid==ngroup;
        start=find(diff(temp)==1)+1;
        ending=find(diff(temp)==-1);
        duration=ending-start;
        meanSec = ceil(mean(duration/SR));
        fprintf('... ... cond %g ... found %g epochs of length %2.2fs on average\n',ngroup,length(start), meanSec)
        epochs=[epochs ; [ngroup*ones(length(start),1) start start+meanSec*SR start-SR]];
    end
    
    %% Re-reference (elec 70)
    y=y(2:161,:,:);
    elecref=repmat(y(70,:,:),[size(y,1) 1 1]);
    refy=y-elecref;
   
    load([preDataPath 'ch1-60_rerefmat.mat'])
    %endif (useCh_1-60)
    refy2=nan(size(reref_mat,1),1,size(y,3));
    for nch=1:size(reref_mat,1)
        refy2(nch,1,:)=y(reref_mat(nch,1),1,:)-y(reref_mat(nch,2),1,:);
    end
    
    %% Filter out 50Hz noise (chronux)
    addpath(genpath(chronuxPath));
    data=squeeze(refy)'; %
    movingwin=[14 7]; % (4s and 2s step size)
    tau=10; % s
    params=[];
    %         params.tapers=[3 5]; %params has the following fields: tapers, Fs, fpass, pad
    params.tapers=[3 4]; %params has the following fields: tapers, Fs, fpass, pad
    params.Fs=SR;
    f0= 50; %
    [datac,~,~,~]=rmlinesmovingwinc(data,movingwin,tau,params,[],[],f0);
    
    dataf=nan(size(datac));
    for nch=1:size(datac,2)
       dataf(:,nch)=highpass(datac(:,nch),SR,1,5); 
    end
    dataf2=nan(size(datac));
    for nch=1:size(datac,2)
       dataf2(:,nch)=bandpass(datac(:,nch),SR,1,40,3); 
    end
    data_nofilt=data';
    data_notch=datac';
    data_bp=dataf';
    data_bp2=dataf2';
    %% Epoch Data
    times=-1:1/SR:meanSec;

    for ntr=1:size(epochs,1)
        ntrc=ntrc+1;
        erp_nofilt(:,ntrc,:)=squeeze(data_nofilt(:,epochs(ntr,4):epochs(ntr,3)));
        erp_notch(:,ntrc,:)=squeeze(data_notch(:,epochs(ntr,4):epochs(ntr,3)));
        erp_bp(:,ntrc,:)=squeeze(data_bp(:,epochs(ntr,4):epochs(ntr,3)));
        erp_bp2(:,ntrc,:)=squeeze(data_bp2(:,epochs(ntr,4):epochs(ntr,3)));
    end
    

    
end
%% baseline correction
    erp_nofilt=erp_nofilt-repmat(mean(erp_nofilt(:,:,times<0 & times>-0.01),3),[1 1 length(times)]);
    erp_notch=erp_notch-repmat(mean(erp_notch(:,:,times<0 & times>-0.01),3),[1 1 length(times)]);
    erp_bp=erp_bp-repmat(mean(erp_bp(:,:,times<0 & times>-0.01),3),[1 1 length(times)]);
    erp_bp2=erp_bp2-repmat(mean(erp_bp2(:,:,times<0 & times>-0.01),3),[1 1 length(times)]);
    
%%
figure; set(gcf,'Position',[440    13   901   785])
subplot(2,2,1); format_fig
plot(times,squeeze(nanmean(erp_nofilt(1:60,:,:),2))); xlim([-0.01 0.1])
xlabel('Time from onset (s)'); ylabel('ERP (\muV)'); ylim([-100 50])
title({'ERP by Elec (av. 12 trials)','No Filter'});

subplot(2,2,2); format_fig
plot(times,squeeze(nanmean(erp_notch(1:60,:,:),2))); xlim([-0.01 0.1])
xlabel('Time from onset (s)'); ylabel('ERP (\muV)'); ylim([-100 50])
title({'ERP by Elec (av. 12 trials)','50Hz Notch'});

subplot(2,2,3); format_fig
plot(times,squeeze(nanmean(erp_bp(1:60,:,:),2))); xlim([-0.01 0.1])
xlabel('Time from onset (s)'); ylabel('ERP (\muV)'); ylim([-100 50])
title({'ERP by Elec (av. 12 trials)','Notch + 1Hz Highpass'});

subplot(2,2,4); format_fig
plot(times,squeeze(nanmean(erp_bp2(1:60,:,:),2))); xlim([-0.01 0.1])
xlabel('Time from onset (s)'); ylabel('ERP (\muV)'); ylim([-100 50])
title({'ERP by Elec (av. 12 trials)','Notch + [1-40]Hz BP'});

%%
erp_P18=squeeze(nanmean(erp_bp(1:60,1:end,times==findclosest(times,0.015)),3))-(squeeze(nanmean(erp_bp(1:60,1:end,times==findclosest(times,0.0117)),3))+squeeze(nanmean(erp_bp(1:60,1:end,times==findclosest(times,0.0192)),3)))/2;
[~,~,~,stats]=ttest(erp_P18',0);
erp_P18=mean(erp_P18,2);
erp_P18tval=stats.tstat;
erp_P18=(reshape(erp_P18,10,6));
erp_P18tval=(reshape(erp_P18tval,10,6));

erp_N20=squeeze(nanmean(erp_bp(1:60,1:end,times>0.025 & times<0.035),3));
[~,~,~,stats]=ttest(erp_N20',0);
erp_N20=mean(erp_N20,2);
erp_N20tval=stats.tstat;

erp_P30=squeeze(nanmean(erp_bp(1:60,1:end,times>0.045 & times<0.055),3)); %-squeeze(nanmean(erp_bp(1:60,:,times>0.025 & times<0.035),3));
[~,~,~,stats]=ttest(erp_P30',0);
erp_P30=mean(erp_P30,2);
erp_P30tval=stats.tstat;

erp_N20=(reshape(erp_N20,10,6));
erp_P30=(reshape(erp_P30,10,6));
erp_N20tval=(reshape(erp_N20tval,10,6));
erp_P30tval=(reshape(erp_P30tval,10,6));

figure; set(gcf,'Position',[ 100         500        1200         400])
subplot(1,3,1:2); format_fig
plot(times,squeeze(nanmean(erp_bp(1:60,1:end,:),2)),'r'); xlim([-0.01 0.1])
xlabel('Time from onset (s)'); ylabel('ERP (\muV)'); ylim([-100 50])
title({'ERP by Elec (av. 12 trials)','Notch + 1Hz HP'});

subplot(1,3,3); format_fig
imagesc(erp_N20tval); caxis([-1 1]*3); colorbar;
title({'T-value against 0 at ~20ms','Notch + 1Hz HP'});

% subplot(2,2,4); format_fig
% imagesc(erp_P30tval); caxis([-1 1]*3); colorbar;
%% ERP image
figure; set(gca,'Position',[440   300   851   498])
subplot(2,3,1:2); title('Exp 1arev'); format_fig
imagesc(times,1:24,squeeze(erp_bp(6,1:end,:)));
xlim([-0.1 0.5]); xlabel('Time'); ylabel('Trials')
caxis([-300 300]); format_fig

subplot(2,3,3); format_fig
imagesc(times,1:24,squeeze(erp_bp(6,1:end,:)));
xlim([-0.01 0.2]); xlabel('Time'); ylabel('Trials')
caxis([-300 300]); format_fig

subplot(2,3,4:5); title('Exp 1arev'); format_fig
imagesc(times,1:24,squeeze(erp_bp(58,1:end,:)));
xlim([-0.1 0.5]); xlabel('Time'); ylabel('Trials')
caxis([-300 300]); format_fig

subplot(2,3,6); format_fig
imagesc(times,1:24,squeeze(erp_bp(58,1:end,:)));
xlim([-0.01 0.2]); xlabel('Time'); ylabel('Trials')
caxis([-300 300]); format_fig
