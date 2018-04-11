%% Init
clear all;
% close all;

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
% rmpath(genpath('/Users/Thomas/Work/local/toolbox/spm12/'))
addpath(genpath('/Users/Thomas/EEGgit/LSCPtools/'))
allversions={'1a','1arev','1b','2a','2b'};

%% Looping on the different versions of the experiment
for nversion=1:5
    version=allversions{nversion};
    if strcmp(version,'1a') % version 1a: sensors on both hands
        data_path='/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/1a/';
        file_name='Vibrotactile_Regular_16_03_2018_19_05_12_0000.mat';
        
        S1channels=[11 21]; % channels given the strongest frequency-tag responses
        LeftF0=[8 11]; % frequencies applied to the left hand
        RightF0=[17 23]; % frequencies applied to the right hand
        
    elseif strcmp(version,'1arev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path='/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/1arev/';
        file_name='Vibrotactile_Regular_22_03_2018_11_40_37_0000.mat';
        
        S1channels=[41 51];
        LeftF0=[23 17];
        RightF0=[11 8];
        
    elseif strcmp(version,'1b') % version 1a: sensors on both hands, hands crossed
        data_path='/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/1b/';
        file_name='Vibrotactile_Crossed_16_03_2018_19_14_05_0000.mat';
        
        S1channels=[11 21];
        LeftF0=[8 11];
        RightF0=[17 23];
        
    elseif strcmp(version,'2a') % version 1a: sensors on hands and feet kept close
        data_path='/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/2a/';
        file_name='Vibrotactile_FeetClose_17_03_2018_13_11_09_0000.mat';
        
        S1channels=[11 21];
        LeftF0=[8 17];
        RightF0=[11 23];
        
    elseif strcmp(version,'2b') % version 1a: sensors on hands and feet kept apart
        data_path='/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/2b/';
        file_name='Vibrotactile_FeetApart_17_03_2018_13_22_10_0000.mat';
        
        S1channels=[11 21];
        LeftF0=[8 17];
        RightF0=[11 23];
        
    end
    S2channels=[];
    
    load([data_path filesep file_name])
    fprintf('... loading %s\n',file_name)
    
    %% Identify trials (trial condition and start/end)
    groups=1:4; % 1 to 5 are digits from thumb to pinkie and 6 is rest
    groupid=squeeze(y(164,1,:));
    epochs=[];
    for ngroup=[1 3]
        temp=groupid==ngroup;
        start=find(diff(temp)==1)+1;
        ending=find(diff(temp)==-1);
        duration=ending-start;
        fprintf('... ... cond %g ... found %g epochs of length %2.2fs on average\n',ngroup,length(start),mean(duration/SR))
        
        epochs=[epochs ; [ngroup*ones(length(start),1) start start+11*SR start-SR]];
    end
    
    %% Re-reference (bi-polar referencing of the 60 electrodes of the 1st grid)
    load('/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/ch1-60_rerefmat.mat')
    y=y(2:161,:,:);
    %     refy=nan(size(reref_mat,1)+60,1,size(y,3));
    %     for nch=1:size(reref_mat,1)
    %         refy(nch,1,:)=y(reref_mat(nch,1),1,:)-y(reref_mat(nch,2),1,:);
    %     end
    refy=nan(60,1,size(y,3));
    for nch=1:60
        refy(nch,1,:)=y(nch,1,:)-mean(y(1:60,1,:),1);
    end
    %% Filter out 50Hz noise (chronux)
    addpath(genpath('/Users/Thomas/Work/local/toolbox/chronux_2_12/'))
    data=squeeze(refy)'; %
    
    movingwin=[14 7]; % (4s and 2s step size)
    tau=10; % s
    params=[];
    %         params.tapers=[3 5]; %params has the following fields: tapers, Fs, fpass, pad
    params.tapers=[3 4]; %params has the following fields: tapers, Fs, fpass, pad
    params.Fs=SR;
    f0=[50 ]; %
    [datac,datafit,Amps,freqs]=rmlinesmovingwinc(data,movingwin,tau,params,[],[],f0);
    datac=datac';
    
    %% Epoch Data
    times=-1:1/SR:11;
    erp=[]; %nan(60,5,10,3*SR+1);
    erp_raw=[];
    for ntr=1:size(epochs,1)
        %         data(:,ntr,:)=squeeze(y(2:161,1,epochs(ntr,4):epochs(ntr,3)));
        erp(:,ntr,:)=squeeze(datac(:,epochs(ntr,4):epochs(ntr,3)));
        erp_raw(:,ntr,:)=squeeze(refy(:,1,epochs(ntr,4):epochs(ntr,3)));
    end
    
    % baseline correction
    erp=erp-repmat(mean(erp(:,:,times<0),3),[1 1 length(times)]);
    erp2=erp-repmat(mean(erp(:,:,times>-0.05 & times<0),3),[1 1 length(times)]);
    
    
    posX=repmat(1:10,1,6);
    posY=[ones(1,10) 2*ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10) 6*ones(1,10)];
    
    %%
    figure;
    subplot(3,9,1:9)
    plot(times,squeeze(mean(erp2(:,:,:),2))); xlim([-0.05 0.4])
    np=0;
    for t=[-0.025:0.025:0.4]
        np=np+1;
        [~,idx]=findclosest(times,t);
        ampl(np,:,:)=squeeze(erp2(:,:,idx));
        
        
        subplot(3,9,9+np)
        [h, pV, ~, stats]=ttest(squeeze(ampl(np,:,:)),0,'dim',2);
        tp=stats.tstat/7*255; %((tp)/maxtp*511+1);
        tp1=tp; tp1(tp1<0)=0;
        tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
        %         scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
        %         xlim([0 11])
        %         ylim([0 7])
        corrtstat=stats.tstat;
        [corrPV]=fdr(pV,0.05);
        corrtstat(pV>corrPV)=0;
        tp3=reshape(corrtstat,10,6);
        
        tp4=squeeze(mean(ampl(np,:,:),3));
        tp4=reshape(tp4,10,6);
        %         imagesc(tp4); caxis([-120 120]);
        imagesc(tp3); caxis([-6 6]);
        title(sprintf('t=%1.3f',t))
    end
    set(gcf,'Name',version)
end
