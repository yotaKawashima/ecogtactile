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

localdef;
useCh_1_60 = 1; %use channels 1 till 60

if useCh_1_60
    startElct = 1;
    endElct   = 60;
else
    startElct = 85;
    endElct   = 104;
end %useCh_1_60

if(sd)
    rootFile    = '/Users/sdan0007/Documents/MATLAB/';
    filePath    = [rootFile 'Ecog/'];
    preDataPath = [filePath 'Data/data_IM/'];  %prefix of data file
    chronuxPath = [rootFile 'Add-Ons/Toolboxes/chronux_2_12'];
    % Thomas you probably added this path somewhere else
    LSCPPath    = [rootFile 'Add-Ons/Toolboxes/LSCPtools'];
    addpath(genpath(LSCPPath));
else
    rootFile    = ['/Users/Thomas/'];
    filePath    = [rootFile 'EEGgit/LSCPtools/'];
    preDataPath = '/media/tLab_BackUp1/Monash/ECogG_somatosens/data_IM/';
    chronuxPath = [rootFile 'Work/local/toolbox/chronux_2_12/'];
    ttestcolortable;
    load modified_ttestcolortable.mat
end %if sd

addpath(genpath(filePath));
allversions={'1a','1arev','1b','2a','2b'};

%% Looping on the different versions of the experiment
for nversion=1:5
    %     S1channels= (startElct:endElct); %[11 21]; % channels given the strongest frequency-tag responses 1:60
    S1channels = [11 12 13 14 15 16 21 22 23 24 25 26 37 38 39 40 47 48 49 50 59 60];
    
    version=allversions{nversion};
    if strcmp(version,'1a') % version 1a: sensors on both hands
        data_path= [preDataPath '1a'];
        file_name='Vibrotactile_Regular_16_03_2018_19_05_12_0000.mat';
        LeftF0=[8 11]; % frequencies applied to the left hand
        RightF0=[17 23]; % frequencies applied to the right hand
    elseif strcmp(version,'1arev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath '1arev/'];
        file_name= 'Vibrotactile_Regular_22_03_2018_11_40_37_0000.mat';
        LeftF0=[23 17];
        RightF0=[11 8];
    elseif strcmp(version,'1b') % version 1a: sensors on both hands, hands crossed
        data_path=[preDataPath '1b/'];
        file_name='Vibrotactile_Crossed_16_03_2018_19_14_05_0000.mat';
        LeftF0=[8 11];
        RightF0=[17 23];
    elseif strcmp(version,'2a') % version 1a: sensors on hands and feet kept close
        data_path=[preDataPath '2a/'];
        file_name='Vibrotactile_FeetClose_17_03_2018_13_11_09_0000.mat';
        LeftF0=[8 17];
        RightF0=[11 23];
        
    elseif strcmp(version,'2b') % version 1a: sensors on hands and feet kept apart
        data_path=[preDataPath '2b/'];
        file_name='Vibrotactile_FeetApart_17_03_2018_13_22_10_0000.mat';
        LeftF0=[8 17];
        RightF0=[11 23];
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
    
    %% Re-reference (bi-polar referencing of the 60 electrodes of the 1st grid)
    if(useCh_1_60)
        load([preDataPath 'ch1-60_rerefmat.mat'])
    else
        load ([preDataPath 'ch85-104_rerefmat.dat'] ,'-ascii')
        reref_mat = ch85_104_rerefmat;
    end %endif (useCh_1-60)
    y=y(2:161,:,:);
    refy=nan(size(reref_mat,1),1,size(y,3));
    for nch=1:size(reref_mat,1)
        refy(nch,1,:)=y(reref_mat(nch,1),1,:)-y(reref_mat(nch,2),1,:);
    end
    % %
    sizeH = size(find(reref_mat(:,3) == 0),1);  %ref mat size horizontal
    sizeV = size(find(reref_mat(:,3) == 1),1);  %ref mat size verticale
    syms elR elC %number of rows and colums of electrode grid
    [elecR elecC] = solve((elR-1)*elC == sizeH, (elC -1)*elR == sizeV);
    elecR (find (elecR < 0)) = [];
    elecC  (find(elecC < 0))  = [];
    elecR = double(elecR);
    elecC = double(elecC);
    fprintf('number of rows is %d number of columns is %d in Electrode grid\n',elecR,elecC);
    fprintf('... ... cond %g ... found %g epochs of length %2.2fs on average\n',ngroup,length(start), meanSec)
    elecR = elecR -1;
    elecC = elecC -1;
    % %
    % for nch=1:60
    %     refy(size(reref_mat,1)+nch,1,:)=y(nch,1,:);
    % end
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
    f0= 100; %
    [datac,~,~,~]=rmlinesmovingwinc(datac,movingwin,tau,params,[],[],f0);
    datac=datac';
    
    %% Epoch Data
    times=-1:1/SR:meanSec;
    erp=[]; %nan(60,5,10,3*SR+1);
    erp_raw=[];
    for ntr=1:size(epochs,1)
        %         data(:,ntr,:)=squeeze(y(2:161,1,epochs(ntr,4):epochs(ntr,3)));
        erp(:,ntr,:)=squeeze(datac(:,epochs(ntr,4):epochs(ntr,3)));
        erp_raw(:,ntr,:)=squeeze(refy(:,1,epochs(ntr,4):epochs(ntr,3)));
    end
    
    % baseline correction
    erp=erp-repmat(mean(erp(:,:,times<0),3),[1 1 length(times)]);
    erpbp=nan(size(erp));
    pow_gamma=[];
    for ntr=1:size(epochs,1)
        for nch=1:size(erp,1)
            erpbp(nch,ntr,:)=abs(hilbert(bandpass(squeeze(erp(nch,ntr,:)), SR, 60, 90, 4)));
            erpbp(nch,ntr,:)=erpbp(nch,ntr,:)-mean(erpbp(nch,ntr,times<0));
            [faxis_gamma, pow]=get_PowerSpec(squeeze(erpbp(nch,ntr,times>0)),SR,0,0);
            pow_gamma(nch,ntr,:)=pow;
        end
    end
    
%     % % load eeglab
%     addpath(genpath('/Users/Thomas/Work/local/toolbox/eeglab13_2_2b/'));
%     fprintf('... ... tf %3.0f/%3.0f\n',0,size(erp,1))
%     tf_data=[];
%     for nch=1:size(erp,1)
%         fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... ... tf %3.0f/%3.0f\n',nch,size(erp,1))
%         tempdata=squeeze(erp(nch,:,:))';
%         
%         [ersp,~,powbase,times2,freqs,~,~,datatf] =  newtimef(tempdata,size(tempdata,1),...
%             [-1 meanSec]*1000,SR,0,'baseline', [-1000 0] ,'maxfreq', 100, 'plotersp', 'off', 'plotitc', 'off', 'padratio', 1, ...
%             'verbose','off'); %%,'timesout',-1500:40:3000);
%         
%         tf_data(nch,:,:)=(ersp);
%     end
    
    %%
    pow_bych=[];
    pow_bych2=[];
    snr_bych=[];
    powgamma_bych=[];
    gammasnr_bych=[];
    kernel=[-.25 -.25 0 0 1 0 0 -.25 -.25];
    for nch=1:size(erp,1)
        for ntr=1:size(erp,2)
            [faxis, pow_bych(nch,ntr,:)]=get_PowerSpec_new(squeeze(erp(nch,ntr,times>=0 & times<10)),1/SR,sum(times>=0 & times<10)/SR,0,0);
            snr_bych(nch,ntr,:)= conv(squeeze(log(pow_bych(nch,ntr,:))), kernel, 'same');
            
            [faxis, powgamma_bych(nch,ntr,:)]=get_PowerSpec_new(squeeze(erpbp(nch,ntr,times>=0 & times<10)),1/SR,sum(times>=0 & times<10)/SR,0,0);
            gammasnr_bych(nch,ntr,:)= conv(squeeze(log(powgamma_bych(nch,ntr,:))), kernel, 'same');
            
        end
    end
    
    erp2=permute(erp,[1 3 2]);
    logSNR_param=[];
    logSNR_param.method='fft';
    logSNR_param.mindist=1;
    [logSNR, faxis1, logpow]=get_logSNR(erp2(:,times>0 & times<10,:),SR,logSNR_param);
    
    logSNR_param=[];
    logSNR_param.method='taper';
    logSNR_param.mindist=1;
    logSNR_param.numTaper=3;
    [logSNR2, faxis2, logpow2]=get_logSNR(erp2(:,times>0 & times<10,:),SR,logSNR_param);
    
    
    figure;
    subplot(1,3,1); format_fig
    plot(faxis1,squeeze(nanmean(logpow(55,:,:),3)));
    title('log Power')
    xlim([6 13])
    subplot(1,3,2); format_fig
    plot(faxis,squeeze(nanmean(snr_bych(55,:,:),2)));
    title('log SNR - previous kernel')
    xlim([6 13])
       ylim([-2 5])
 subplot(1,3,3); format_fig
    plot(faxis1,squeeze(nanmean(logSNR(55,:,:),3)));
    xlim([6 13])
    ylim([-2 5])
    title('log SNR - new kernel')

    %%
    figure;
    set(gcf,'Name',sprintf('%s - %s',version,'SNR'))
    theseFreqs=[LeftF0 RightF0 abs(LeftF0(2)-LeftF0(1)) abs(RightF0(2)-RightF0(1))];
    colorCode='rrbbrb';
    lineCode={'-','-','-','-','--','--','--','--'};
    fillCode=[ones(1,4) zeros(1,4)];
    
    subplot(1,2,1); format_fig;
    plot(faxis,squeeze(mean(snr_bych(1:sizeH,:,:),2))','Color','k')
    c=0;
    hold on;
    for n=theseFreqs
        c=c+1;
        [~,idxfreqs(c)]=findclosest(faxis,n);
        tp=squeeze(mean(snr_bych(1:sizeH,:,idxfreqs(c)),2));
        if fillCode(c)==1
            scatter(n,max(tp),'MarkerEdgeColor',colorCode(c),'MarkerFaceColor',colorCode(c),'Marker','o','SizeData',64);%,'LineStyle',lineCode{c})
        else
            scatter(n,max(tp),'MarkerEdgeColor',colorCode(c),'Marker','d','SizeData',64);%,'LineStyle',lineCode{c})
        end
    end
        format_fig;
xlim([2 45])
    ylim([-2 5])
    ylabel('SNR')
    xlabel('Frequency (Hz)')
    
    subplot(1,2,2); format_fig
    plot(faxis,squeeze(mean(snr_bych(sizeH+1:sizeH+sizeV,:,:),2))','Color','k')
    c=0;
    hold on;
    for n=theseFreqs
        c=c+1;
        [~,idxfreqs(c)]=findclosest(faxis,n);
        tp=squeeze(mean(snr_bych(sizeH+1:sizeH+sizeV,:,idxfreqs(c)),2));
        if fillCode(c)==1
            scatter(n,max(tp),'MarkerEdgeColor',colorCode(c),'MarkerFaceColor',colorCode(c),'Marker','o','SizeData',64);%,'LineStyle',lineCode{c})
        else
            scatter(n,max(tp),'MarkerEdgeColor',colorCode(c),'Marker','d','SizeData',64);%,'LineStyle',lineCode{c})
        end
    end
    format_fig;
    xlim([2 45])
    ylim([-2 5])
    ylabel('SNR')
    xlabel('Frequency (Hz)')
    
    %%
    % extract fundamentals
    fondFreq=[LeftF0 RightF0];
    pow_tag=[];
    pow_tag2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_tag(:,:,nfreq)=(snr_bych(:,:,findfreq));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_tag2(:,nfreq)=(mean(snr_bych(:,:,findfreq),2));%-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    % extract IM
    fondFreq=[abs(LeftF0(2)-LeftF0(1)) abs(RightF0(2)-RightF0(1))];
    pow_IM=[];
    pow_IM2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_IM(:,:,nfreq)=(snr_bych(:,:,findfreq));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_IM2(:,nfreq)=(mean(snr_bych(:,:,findfreq),2));%-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    
    % plot grids
    %     posX=repmat(1:elecC,1,6);
    posX=repmat(1:elecC,1,elecR+1);
    posY=[];
    for ind =1:elecR +1
        posY = [posY ind * ones(1,elecC)];
    end %for
    %posY=[ones(1,elecC) 2*ones(1,elecC) 3*ones(1,elecC) 4*ones(1,elecC) 5*ones(1,elecC) 6*ones(1,elecC)];
    
    figure;
    set(gcf,'position',[-442        1399        1075         458],'Name',sprintf('%s - %s',version,'Horizontal Bi-Polar'))
    subplot(2,2,1)
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,1:2)),0,'dim',2);
    imagesc(1:elecC,1:elecR+1,reshape(mean(stats.tstat,3),elecC,elecR+1)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('Freq Tag Contra')
    xlim([0.5 elecC + 0.5])
    ylim([0.5 elecR + 1.5])
    
    % end
    
    subplot(2,2,3)
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,1)),0,'dim',2);
    imagesc(1:elecC,1:elecR+1,reshape(mean(stats.tstat,3),elecC,elecR+1)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('IM Tag Contra')
    xlim([0.5 elecC + 0.5])
    ylim([0.5 elecR + 1.5])
    
    subplot(2,2,2)
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,3:4)),0,'dim',2);
    imagesc(1:elecC,1:elecR+1,reshape(mean(stats.tstat,3),elecC,elecR+1)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('Freq Tag Ipsi')
    xlim([0.5 elecC + 0.5])
    ylim([0.5 elecR + 1.5])
    
    subplot(2,2,4)
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,2)),0,'dim',2);
    imagesc(1:elecC,1:elecR+1,reshape(mean(stats.tstat,3),elecC,elecR+1)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('IM Tag Ipsi')
    xlim([0.5 elecC + 0.5])
    ylim([0.5 elecR + 1.5])
    
    figure;
    set(gcf,'position',[-442        1399        1075         458],'Name',sprintf('%s - %s',version,'Vertical Bi-Polar'))
    subplot(2,2,1)
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(sizeH+1:sizeH+sizeV,:,1:2)),0,'dim',2);
    imagesc(1:elecC+1,1:elecR,reshape(mean(stats.tstat,3),elecC+1,elecR)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('Freq Tag Contra')
    xlim([0.5 elecC + 1.5])
    ylim([0.5 elecR + 0.5])
    
    % end
    
    subplot(2,2,2)
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(sizeH+1:sizeH+sizeV,:,1)),0,'dim',2);
    imagesc(1:elecC+1,1:elecR,reshape(mean(stats.tstat,3),elecC+1,elecR)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('IM Tag Contra')
    xlim([0.5 elecC + 1.5])
    ylim([0.5 elecR + 0.5])
    
    subplot(2,2,3)
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(sizeH+1:sizeH+sizeV,:,3:4)),0,'dim',2);
    imagesc(1:elecC+1,1:elecR,reshape(mean(stats.tstat,3),elecC+1,elecR)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('Freq Tag Ipsi')
    xlim([0.5 elecC + 1.5])
    ylim([0.5 elecR + 0.5])
    
    subplot(2,2,4)
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(sizeH+1:sizeH+sizeV,:,2)),0,'dim',2);
    imagesc(1:elecC+1,1:elecR,reshape(mean(stats.tstat,3),elecC+1,elecR)'); set(gca,'ydir','normal');
    caxis([-8 8]); colormap(cmap); colorbar;
    title('IM Tag Ipsi')
    xlim([0.5 elecC + 1.5])
    ylim([0.5 elecR + 0.5])
    
end
