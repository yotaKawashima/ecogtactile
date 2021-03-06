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
    filePath    = [rootFile 'Ecog_Local/'];
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
end %if sd

    ttestcolortable;
    load modified_ttestcolortable.mat

addpath(genpath(filePath));
allversions={'1a','1arev','1b','2a','2b'};

%% Looping on the different versions of the experiment
for nversion=1:length(allversions) %was 4:5
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
    pow_bych7 = cell(11);
    snr_bych7 = cell(11);
    kernel=[-.25 -.25 0 0 1 0 0 -.25 -.25];
    %%
    logSNR = [];
    faxis = [];%cell(11);
    logpow = [];
    logSNRAll=[];
    faxisAll =[];
    logpowAll = [];
    snrparam=[];
    snrparam.method='fft';
    snrparam.mindist=0.9;

    erp2=permute(erp,[1 3 2]);
    totalTime = (meanSec -1);
    secInd = 1;
    secArr = (totalTime/3:totalTime/3:totalTime);
    for sec=totalTime/3:totalTime/3:totalTime
%         for nch=1:size(erp2,1)
%                 for ntr=1:size(erp2,3)
%     for nch=1:size(erp,1)
%         for ntr=1:size(erp,2)
%             [faxis, pow_bych(nch,ntr,:,sec)]=get_PowerSpec_new(squeeze(erp(nch,ntr,times>=sec -1 & times<sec)),1/SR,sum(times>=sec-1 & times<sec)/SR,0,0);
%             [faxis7{sec}, pow_bych7{sec}(nch,ntr,:,sec)]=get_PowerSpec_new(squeeze(erp(nch,ntr,times>=0 & times<sec)),1/SR,sum(times>=0 & times<sec)/SR,0,0);
%             snr_bych(nch,ntr,:,sec)= conv(squeeze(log(pow_bych(nch,ntr,:,sec))), kernel, 'same');
%             snr_bych7{sec}(nch,ntr,:,sec)= conv(squeeze(log(pow_bych7{sec}(nch,ntr,:,sec))), kernel, 'same');
%             [faxis, powgamma_bych(nch,ntr,:)]=get_PowerSpec_new(squeeze(erpbp(nch,ntr,times>=0 & times<10)),1/SR,sum(times>=0 & times<10)/SR,0,0);
%             gammasnr_bych(nch,ntr,:)= conv(squeeze(log(powgamma_bych(nch,ntr,:))), kernel, 'same');
            
             [logSNR(:,:,:,secInd), faxis(:,secInd), logpow(:,:,:,secInd)]=get_logSNR(erp2(:,times>=(sec-totalTime/3) & times<(sec),:),SR,snrparam); 
%              snr_bych7{secInd} = get_logSNR(erp2(:,times>=(sec-totalTime/3) & times<(sec),:),SR,snrparam); 
             secInd = secInd +1;
%         end
%     end
    end %for sec 1:11
    
    [logSNRAll, faxisAll, logpowAll]=get_logSNR(erp2(:,times>=0 & times<totalTime,:),SR,snrparam); 
%     Debug func
%  figure, plot(faxis(:,1),squeeze(mean(mean(logSNR(:,:,1,1),4),1))) plot
% % [theFreq,findfreq]=findclosest(faxis9,fondFreq(2));
% pow_tag9=(logSNR9(:,111,:));
%  [h, pV, ~, stats]=ttest(squeeze(pow_tag9(:,1,:)),0,'dim',2);
    %%
    
    
   if(0) 
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
end %if 0
    
    %%
    % extract fundamentals
    fondFreq=[LeftF0 RightF0];
    pow_tag=[];
    pow_tag2=[];
    pow_tagA = [];
    pow_tag7 = cell(11);
    freqArr = nan(1,length(fondFreq));
   for sec=1:length(secArr)
    for nfreq=1:length(fondFreq)
        [theFreq,findfreq]=findclosest(faxis(:,sec),fondFreq(nfreq)); %making sure we did find the right freq and not somehing very far way
        freqArr(nfreq) = theFreq;      
%         pow_tag(:,:,nfreq,sec)=(snr_bych(:,:,findfreq,sec));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
%         pow_tag7{sec}(:,:,nfreq,sec)=(snr_bych7{sec}(:,:,findfreq,sec));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
%         pow_tag2(:,nfreq)=(mean(snr_bych(:,:,findfreq),2));%-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
          pow_tag(:,nfreq,:,sec)=(logSNR(:,findfreq,:,sec));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
    end
    end %for sec 
    
    freqArrA = nan(1,length(fondFreq));
    for nfreq=1:length(fondFreq)
    [theFreqA,findfreqA]=findclosest(faxisAll,fondFreq(nfreq)); %making sure we did find the right freq and not somehing very far way
    pow_tagA(:,nfreq,:) = (logSNRAll(:,findfreqA,:));
    freqArrA(nfreq) = theFreqA;   
    end %for nfreq=1:length(fondFreq) 
   
    ipsiContraStr = {'Freq Tag Contra ', 'Freq Tag Contra ', 'Freq Tag Ipsi ', 'Freq Tag Ipsi ' };
    
  
    freqInd = 1;
        for ind = 1:length(freqArr)
        figure;
        subplot(2,2,4);
        k=ind;
        [h, pV, ~, stats]=ttest(squeeze(pow_tagA(:,freqInd,:)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
        title([ipsiContraStr{(k)} num2str(freqArr(freqInd))  'Sec = 0-10']);
        format_fig;
        colormap(cmap);
        
%         subplot(1,3,3);
%         [h, pV, ~, stats]=ttest(squeeze(pow_tag(:,:,k*2 -1:k*2,sec)),0,'dim',2);
%         [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
%         title([ipsiContraStr{(k)} num2str(freqArr(k*2 -1)) ' and ' num2str(freqArr(k*2)) ' mean']);
%         format_fig;
%         colormap(cmap);

            for sec=1:length(secArr)
            set(gcf,'position',[-442        1399         362         458],'Name',sprintf('%s - %s %s ',version,ipsiContraStr{(k)},num2str(freqArrA(freqInd))))
            [h, pV, ~, stats]=ttest(squeeze(pow_tag(:,freqInd,:,sec)),0,'dim',2);
%             [h, pV, ~, stats]=ttest(squeeze(pow_tag7{sec}(:,:,freqInd,sec)),0,'dim',2);
%              [h, pV, ~, stats]=ttest(squeeze(logSNR(:,:,freqInd,sec)),0,'dim',2);
            subplot(2,2,sec);
            [p, FV] = draw_biprref(stats.tstat, reref_mat, [10 6], [-8 8]);
            title([ipsiContraStr{(k)} num2str(freqArr(freqInd)) ' sec ' num2str(secArr(sec))]);
            format_fig;
            colormap(cmap)
        end %for %for sec
        freqInd = freqInd +1;
        end % ind = 1:legnth(freqArr)/2
%     end %for  k = 1:legnth(freqArr)/2
    
    
%   ====================IM Analysis ======================================
    %Make a permutaion of all possibile IM options
% % 
if (0)
    place = 1;
    fondFreqIm = nan (1,300);
    ipsiContrNumRel = zeros(1,300);
    ipsi_ipsi = 1;
    ipsi_lateral = 2;
    lateral_lateral = 3;
    ipsiContaraRelationStr = {'ipsi_ipsi', 'ipsi_lateral', 'lateral_lateral'}; 
    for mult1 = 1:5
        for mult2 = 1:5
            for ind=1:length(fondFreq)
                %define the relaion ipsi lateal  
                for k=ind+1:length(fondFreq)      
                    if( k < 3 & ind < 3)
                        resRel = ipsi_ipsi;
                    elseif (k > 2 & ind < 3 | ind > 2 & k < 3)
                        resRel = ipsi_lateral;
                    else resRel = lateral_lateral;
                    end %if  k < 3 & ind < 3) for relation
                    ipsiContrNumRel(place) = resRel;
                    fondFreqIm(place) = abs(mult1*fondFreq(ind) - fondFreq(k)*mult2);
                    place = place + 1;
                    ipsiContrNumRel(place) = resRel;
                    fondFreqIm(place) = mult1*fondFreq(ind) + fondFreq(k)*mult2;                    
                    place = place + 1;
                end
            end
        end
    end 
    
% % %     [C,IA,IC] = unique(A) also returns index vectors IA and IC such that
% % %     C = A(IA) and A = C(IC) (or A(:) = C(IC), if A is a matrix or array).

    lineNoise = [50,100,150];
    for freq = lineNoise
        fLine = find(freq == fondFreqIm);
        fondFreqIm(fLine) = [];
    end
    
    [fondFreqIm, IA, IC] = unique(fondFreqIm, 'stable'); %Avoid repitions
    ipsiContrNumRel = ipsiContrNumRel(IA);  %delete repition from array
% %     
    % extract IM
    maxDiffValue = 1; %allow less than 1Hz differences between the Im freq and the found freq
%     fondFreq=[abs(LeftF0(2)-LeftF0(1)) abs(RightF0(2)-RightF0(1))]; 
    freqArrIm = nan(1,length(fondFreqIm)); %remember the closest number found
    pow_IM=[];
    pow_IM2=[];
    ignoreNum = 0;
    for nfreq=1:length(fondFreqIm)
        [theImFreq,findfreq]=findclosest(faxis,fondFreqIm(nfreq));
        if(abs(theImFreq - fondFreqIm(nfreq) >= maxDiffValue))
            fprintf('ignoring freq found %d freq search %d\n', theImFreq, fondFreqIm(nfreq));
            ignoreNum = ignoreNum + 1; %if no one pays attention to the printing this will count the num of ingnores
            continue;
        end %if abs
        freqArrIm(nfreq) = theImFreq;
        pow_IM(:,:,nfreq)=(snr_bych(:,:,findfreq));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_IM2(:,nfreq)=(mean(snr_bych(:,:,findfreq),2));%-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end

    ipsiContraStrIm = {'IM Contra ', 'IM Ipsi ' };
    maxNumK = 3;
    [h, pV, ~, stats]=ttest(pow_IM,0,'dim',2);
    tempStat = stats;
    [Y I] = maxk(tempStat.tstat,maxNumK, 3);
    [O P] = max(Y);
    freqArrIm(I(P(1)));
    
    for k = 1:3:length(freqArrIm)
        figure;
        for ind = 1:3
        subplot(1,3,ind); 
        [h, pV, ~, stats]=ttest(squeeze(pow_tag(:,:,1:2)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
        title([ipsiContraStr{(ind)} num2str(freqArr(k*2 -1)) ' and ' num2str(freqArr(k*2)) ' mean']);
        format_fig;
        colormap(cmap);
        for ind = 1:length(freqArr)/2
            set(gcf,'position',[-442        1399         362         458],'Name',sprintf('%s - %s %s %s',version,ipsiContraStr{(k)},num2str(freqArr(k*2-1)), num2str(freqArr(k*2))))
            [h, pV, ~, stats]=ttest(squeeze(pow_tag(:,:,ind)),0,'dim',2);
            subplot(1,3,ind);
            [p, FV] = draw_biprref(stats.tstat, reref_mat, [10 6], [-8 8]);
            title([ipsiContraStr{(k)} num2str(freqArr(ind*k))]);
            format_fig;
            colormap(cmap)
        end %for  ind = 1:legnth(freqArr)/2
    end %for  k = 1:legnth(freqArr)/2
    
    figure;
    set(gcf,'position',[-442   867   362   458],'Name',sprintf('%s - %s',version,'Im Contra'))
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(:,:,1)),0,'dim',2);
    [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
    title('IM Tag Contra')
    format_fig;
    colormap(cmap)
    
% % %     figure;
% % %     set(gcf,'position',[ -79        1399         362         458],'Name',sprintf('%s - %s',version,'Freq Ipsi'))
% % %     [h, pV, ~, stats]=ttest(squeeze(pow_tag(:,:,3:4)),0,'dim',2);
% % %     [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
% % %     title('Freq Tag Ipsi')
% % %     format_fig;
% % %     colormap(cmap)
    
    figure;
    set(gcf,'position',[-79   867   362   458],'Name',sprintf('%s - %s',version,'Im Ipsi'))
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(:,:,2)),0,'dim',2);
    [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
    title('IM Tag Ipsi')
    format_fig;
    colormap(cmap)
    
    
    
    if strcmp(version(1),'2')
        % extract IM
        fondFreq=[abs(8-11) abs(17-23)];
        pow_IM3=[];
        pow_IM4=[];
        for nfreq=1:length(fondFreq)
            [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
            pow_IM3(:,:,nfreq)=(snr_bych(:,:,findfreq));%-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
            pow_IM4(:,nfreq)=(mean(snr_bych(:,:,findfreq),2));%-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
        end
        figure;
        set(gcf,'position',[-442   867   362   458],'Name',sprintf('%s - %s',version,'IM  Hands'))
        [h, pV, ~, stats]=ttest(squeeze(pow_IM3(:,:,1)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
        title('IM Tag Hands')
        format_fig;
        colormap(cmap)
        
        figure;
        set(gcf,'position',[-79   867   362   458],'Name',sprintf('%s - %s',version,'IM Feets'))
        [h, pV, ~, stats]=ttest(squeeze(pow_IM3(:,:,2)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), reref_mat, [10 6], [-8 8]);
        title('IM Tag Feet')
        format_fig;
        colormap(cmap)
        
    end
    end

end %if
end %for