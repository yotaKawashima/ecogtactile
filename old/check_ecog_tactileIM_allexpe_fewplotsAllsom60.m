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

sd = 1; %Thomas change this to 0
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
    rootFile    = [rootFile '/Users/Thomas/'];
    filePath    = [rootFile 'EEGgit/LSCPtools/'];
    preDataPath = [rootFile 'temp_data/ECogG_tapping/data_IM/'];
    chronuxPath = [rootFile 'Work/local/toolbox/chronux_2_12/'];
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
        
% %         S1channels= (1:60); %[11 21]; % channels given the strongest frequency-tag responses 1:60
        LeftF0=[8 11]; % frequencies applied to the left hand
        RightF0=[17 23]; % frequencies applied to the right hand
        
    elseif strcmp(version,'1arev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath '1arev/'];
        file_name= 'Vibrotactile_Regular_22_03_2018_11_40_37_0000.mat';
        
% %         S1channels=[41 51];
        LeftF0=[23 17];
        RightF0=[11 8];
        
    elseif strcmp(version,'1b') % version 1a: sensors on both hands, hands crossed
        data_path=[preDataPath '1b/'];
        file_name='Vibrotactile_Crossed_16_03_2018_19_14_05_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0=[8 11];
        RightF0=[17 23];
        
    elseif strcmp(version,'2a') % version 1a: sensors on hands and feet kept close
        data_path=[preDataPath '2a/'];
        file_name='Vibrotactile_FeetClose_17_03_2018_13_11_09_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0=[8 17];
        RightF0=[11 23];
        
    elseif strcmp(version,'2b') % version 1a: sensors on hands and feet kept apart
        data_path=[preDataPath '2b/'];
        file_name='Vibrotactile_FeetApart_17_03_2018_13_22_10_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0=[8 17];
        RightF0=[11 23];
        
    end
    S2channels=[];
    
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
    %%%%%% Grid search of optimal parameters: num tappers=4 and window=14
% % %     taus=[1:20];
% % %     kappas=[1:8];
% % %     for ntau=1:length(taus)
% % %         for nk=1:length(kappas)
% % %             fprintf('tau %g tappers %g\n',taus(ntau),kappas(nk))
% % %     %         movingwin=[4 2]; % (4s and 2s step size)
% % %             movingwin=[taus(ntau) taus(ntau)/2]; % (4s and 2s step size)
% % %             tau=10; % s
% % %             params=[];
% % %     %         params.tapers=[3 5]; %params has the following fields: tapers, Fs, fpass, pad
% % %             params.tapers=[3 kappas(nk)]; %params has the following fields: tapers, Fs, fpass, pad
% % %             params.Fs=SR;
% % %             f0=50; %
% % %             [datac,datafit,Amps,freqs]=rmlinesmovingwinc(data(:,11),movingwin,tau,params,[],[],f0);
% % %             datac=datac';
% % %     
% % %             [faxisfilt{ntau}(nk,:), powfilt{ntau}(nk,:)]=get_PowerSpec(datac,SR,0,0);
% % %         end
% % %     end
% % %     
% % %     figure;
% % %     for ntau=1:length(taus)
% % %         for nk=1:length(kappas)
% % %            subplot(8,20,(nk-1)*20+ntau)
% % %             plot(faxisfilt{ntau}(nk,:),powfilt{ntau}(nk,:)); xlim([10 60])
% % %         end
% % %     end
    movingwin=[14 7]; % (4s and 2s step size)
    tau=10; % s
    params=[];
    %         params.tapers=[3 5]; %params has the following fields: tapers, Fs, fpass, pad
    params.tapers=[3 4]; %params has the following fields: tapers, Fs, fpass, pad
    params.Fs=SR;
    f0= 50; %
    [datac,datafit,Amps,freqs]=rmlinesmovingwinc(data,movingwin,tau,params,[],[],f0);
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
    
    
    
    %%
    pow_bych=[];
    pow_bych2=[];
    for nch=1:size(erp,1)
        for ntr=1:size(erp,2)
            [faxis, pow_bych(nch,ntr,:)]=get_PowerSpec_new(squeeze(erp(nch,ntr,:)),1/SR,size(erp,3)/SR,0,0);
        end
    end
    
    %%
    figure;
    set(gcf,'position',[-499 801 854 1056],'Name',version)
    subplot(4,2,1)
%     plot(faxis,log(squeeze(mean(mean(pow_bych([find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==0) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==0)],:,:),1),2))),'Color','k','LineWidth',2)
    plot(faxis,log(squeeze(mean(mean(pow_bych([find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==0)],:,:),1),2))),'Color','k','LineWidth',2)
    xlim([1 45])
    xlabel('Freq (Hz)')
    ylabel('logPower')
% %     ylabel(sprintf('SNR (ch %g %g) bipolar I',S1channels(1),S1channels(2)))
    ylabel(sprintf('SNR (ch %g - %g) bipolar I',S1channels(1),S1channels(end)))
    
    subplot(4,2,2)
    %plot(faxis,log(squeeze(mean(mean(pow_bych([find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==1) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==1)],:,:),1),2))),'Color','k','LineWidth',2)
    plot(faxis,log(squeeze(mean(mean(pow_bych([find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==1)],:,:),1),2))),'Color','k','LineWidth',2)
    xlim([1 45])
    xlabel('Freq (Hz)')
    ylabel('logPower')
% %     ylabel(sprintf('SNR (ch %g %g) bipolar II',S1channels(1),S1channels(2)))
    ylabel(sprintf('SNR (ch %g - %g) bipolar I',S1channels(1),S1channels(end)))
    
    %%
    % extract fundamentals
    fondFreq=[LeftF0 RightF0];
    pow_tag=[];
    pow_tag2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_tag(:,:,nfreq)=log(pow_bych(:,:,findfreq))-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_tag2(:,nfreq)=log(mean(pow_bych(:,:,findfreq),2))-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    % extract IM
    fondFreq=[abs(LeftF0(2)-LeftF0(1)) abs(RightF0(2)-RightF0(1))];
    pow_IM=[];
    pow_IM2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_IM(:,:,nfreq)=log(pow_bych(:,:,findfreq))-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_IM2(:,nfreq)=log(mean(pow_bych(:,:,findfreq),2))-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    
    % plot grids
%     posX=repmat(1:elecC,1,6);
    posX=repmat(1:elecC,1,elecR+1);
    posY=[];
    for ind =1:elecR +1
    posY = [posY ind * ones(1,elecC)];
    end %for
    %posY=[ones(1,elecC) 2*ones(1,elecC) 3*ones(1,elecC) 4*ones(1,elecC) 5*ones(1,elecC) 6*ones(1,elecC)];
    
    subplot(4,2,3)
    % titles={sprintf('f0=%g (contra)',LeftF0(1)),sprintf('f0=%g (contra)',LeftF0(2)),sprintf('f0=%g (ipsi)',RightF0(1)),sprintf('f0=%g (ipsi)',RightF0(2))};
    % for nfreq=1:4
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,1:2)),0,'dim',2);
    tp=mean(stats.tstat,3)/3.5*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('Freq Tag Contra')
%     xlim([0 11])
%     ylim([0 7])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])

    % end
    
    subplot(4,2,5)
    % titles={sprintf('IM=%g (contra)',LeftF0(1)+LeftF0(2)),sprintf('IM=%g (ipsi)',RightF0(1)+RightF0(2))};
    % for nfreq=1:2
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,1))',0);
    % tp=((tp)/maxtp*511+1);
    tp=stats.tstat'/2*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('IM Contra')
%     xlim([0 11])
%     ylim([0 7])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    
%     posX=repmat(1:10,1,5);
%     posY=[ones(1,10) 2*ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10)];
    
    subplot(4,2,4)
    % titles={sprintf('f0=%g (contra)',LeftF0(1)),sprintf('f0=%g (contra)',LeftF0(2)),sprintf('f0=%g (ipsi)',RightF0(1)),sprintf('f0=%g (ipsi)',RightF0(2))};
    % for nfreq=1:4
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,3:4)),0,'dim',2);
    tp=mean(stats.tstat,3)/3.5*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('Freq Tag Ipsi')
%     xlim([0 11])
%     ylim([0 6])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    
    subplot(4,2,6)
    % titles={sprintf('IM=%g (contra)',LeftF0(1)+LeftF0(2)),sprintf('IM=%g (ipsi)',RightF0(1)+RightF0(2))};
    % for nfreq=1:2
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,2))',0);
    % tp=((tp)/maxtp*511+1);
    tp=stats.tstat'/2*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('IM Ipsi')
%     xlim([0 11])
%     ylim([0 6])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    %%
    %myChans=[find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==0) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==0)];
    myChans=find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==0);
    subplot(4,2,7)
    simpleBarPlot(1-0.2,squeeze(mean(mean(pow_tag(myChans,:,1),3),1)),'b',0.38,'r',[],2);
    simpleBarPlot(1+0.2,squeeze(mean(mean(pow_tag(myChans,:,2),3),1)),'b',0.38,'r',[],2);
    
    simpleBarPlot(2-0.2,squeeze(mean(mean(pow_tag(myChans,:,3),3),1)),'k',0.38,'r',[],2);
    simpleBarPlot(2+0.2,squeeze(mean(mean(pow_tag(myChans,:,4),3),1)),'k',0.38,'r',[],2);
    
    
    simpleBarPlot(3-0.2,squeeze(mean(pow_IM(myChans,:,1),1)),'b',0.38,'r',[],2);
    simpleBarPlot(3+0.2,squeeze(mean(pow_IM(myChans,:,2),1)),'k',0.38,'r',[],2);
    
    set(gca,'FontSize',18,'FontWeight','bold')
    xlim([0.5 3.5])
    set(gca,'XTick',1:3,'XTickLabel',{'F0-contra','F0-ipsi','IM'})
    ylabel(sprintf('SNR (ch %g - %g)',S1channels(1),S1channels(end)))
    
% %     myChans=[find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==1) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==1)];
    myChans=find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==1);
    subplot(4,2,8)
    simpleBarPlot(1-0.2,squeeze(mean(mean(pow_tag(myChans,:,1),3),1)),'b',0.38,'r',[],2);
    simpleBarPlot(1+0.2,squeeze(mean(mean(pow_tag(myChans,:,2),3),1)),'b',0.38,'r',[],2);
    
    simpleBarPlot(2-0.2,squeeze(mean(mean(pow_tag(myChans,:,3),3),1)),'k',0.38,'r',[],2);
    simpleBarPlot(2+0.2,squeeze(mean(mean(pow_tag(myChans,:,4),3),1)),'k',0.38,'r',[],2);
    
    
    simpleBarPlot(3-0.2,squeeze(mean(pow_IM(myChans,:,1),1)),'b',0.38,'r',[],2);
    simpleBarPlot(3+0.2,squeeze(mean(pow_IM(myChans,:,2),1)),'k',0.38,'r',[],2);
    
    set(gca,'FontSize',18,'FontWeight','bold')
    xlim([0.5 3.5])
    set(gca,'XTick',1:3,'XTickLabel',{'F0-contra','F0-ipsi','IM'})
    ylabel(sprintf('SNR (ch %g - %g)',S1channels(1),S1channels(end)))
    
end
