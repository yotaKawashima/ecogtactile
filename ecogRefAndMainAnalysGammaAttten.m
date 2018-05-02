% This script load the epochs data, making the referencing and the main/heavy analysis

% clear all;

% Im at (8 + 11)*3 + 8 = 65
fprintf(' Processing Ref and Main Analysis .......\n');

load ecogPrePros.mat;
load ecogGlobal.mat;

printElc = 0;
%% Re-reference (bi-polar referencing of the 60 electrodes of the 1st grid)
    if(useCh_1_60)
        load([preDataPath 'ch1-60_rerefmat.mat'])
    else
        load ([preDataPath 'ch85-104_rerefmat.dat'] ,'-ascii')
        reref_mat = ch85_104_rerefmat;
    end %endif (useCh_1-60)
     
    
    pow_bych_all=[];
    %pow_bych2=[]; not used
    pow_gamma_all=[];
    %% Looping on the different versions of the experiment
for nversion=1:numOfCondition
    epochs = allEpochs(:,:,nversion);
    y  = rawData{nversion}.y;
    SR = rawData{nversion}.SR;
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
    
    if (printElc)
    fprintf('number of rows is %d number of columns is %d in Electrode grid\n',elecR,elecC);
    fprintf('... found %g epochs of length %2.2fs on average\n',length(start), meanSec)
    end %if printElc
    
    elecR = elecR -1;
    elecC = elecC -1;

    addpath(genpath(chronuxPath));
    data=squeeze(refy)'; %
    
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
    
    att_1 = erp(:,1:6,:);
    att_3 = erp(:,7:12,:);
    att_1 = (squeeze(mean(att_1,2)))'; %Avarage between conditions
    att_3 = (squeeze(mean(att_3,2)))';
    
    att_1_H = erp(1:sizeH,1:6,:);
    att_3_H = erp(1:sizeH,7:12,:);
    att_1_H = (squeeze(mean(att_1_H,2)))'; %Avarage between conditions
    att_3_H = (squeeze(mean(att_3_H,2)))';
    
    att_1_V = erp(1+sizeH:end,1:6,:);
    att_3_V = erp(1+sizeH:end,7:12,:);
    att_1_V = (squeeze(mean(att_1_V,2)))'; %Avarage between conditions
    att_3_V = (squeeze(mean(att_3_V,2)))';
    
    
    p = 0.05; %probability rate
    params.tapers = [ 3 5 ];
    params.fpass = [ 60 90 ];
    params.Fs = 1200; %Fs;
    movingwin = [0.1 0.05];
    params.trialave = 1; % Avarage between electrodes
    params.err = [2 p]; % Jackknife error bars;  [1 p] - Theoretical error bars taken from http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/continuous/mtspecgramc.html
    [S_1,t_1,f_1, Serr_1]=mtspecgramc(att_1,movingwin,params);
    [S_3,t_3,f_3, Serr_3]=mtspecgramc(att_3,movingwin,params);
    figure('Name','att_1 H + V'),  plot_matrix(S_1,t_1,f_1,1,Serr_1);
    figure('Name','att_3 H + V'), plot_matrix(S_3,t_3,f_3,1,Serr_3);
    
    [S_1,t_1,f_1, Serr_1]=mtspecgramc(att_1_H,movingwin,params);
    [S_3,t_3,f_3, Serr_3]=mtspecgramc(att_3_H,movingwin,params);
    figure('Name','att_1 H'),  plot_matrix(S_1,t_1,f_1,1,Serr_1);
    figure('Name','att_3 H'), plot_matrix(S_3,t_3,f_3,1,Serr_3);
    
    [S_1,t_1,f_1, Serr_1]=mtspecgramc(att_1_V,movingwin,params);
    [S_3,t_3,f_3, Serr_3]=mtspecgramc(att_3_V,movingwin,params);
    figure('Name','att_1 V'), plot_matrix(S_1,t_1,f_1,1,Serr_1);
    figure('Name','att_3 V'), plot_matrix(S_3,t_3,f_3,1,Serr_3);
    
    %==================Do not avarage acroos electrodes===========
    posX=repmat(1:elecC,1,elecR+1);
    posY=[];
    for ind =1:elecR +1
    posY = [posY ind * ones(1,elecC)];
    end %for
    
    params.trialave = 0; % DO NOT Avarage between electrodes
    
    [S_1,t_1,f_1, Serr_1]=mtspecgramc(att_1_H,movingwin,params);
           
    figure('Name', 'Horizontal');
    for ind=1:length(f_1)
         subplot(3,1,ind);
         myF = (squeeze(S_1(:,ind, :)))';
         [h, pV, ~, stats]=ttest(myF,0,'dim',2);
         tp=mean(stats.tstat,3)/3.5*255; %((tp)/maxtp*511+1);
         tp1=tp; tp1(tp1<0)=0;
         tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
         scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
         xlim([0 elecC + 2]);
         ylim([0 elecR + 2]);
         title(['freq = ' num2str(f_1(ind))]);
    end %for
    %%
    for nch=1:size(erp,1)
        for ntr=1:size(erp,2)
            [faxis_all(:,nversion), pow_bych_all(nch,ntr,:, nversion)]=get_PowerSpec_new(squeeze(erp(nch,ntr,:)),1/SR,size(erp,3)/SR,0,0);
        end
    end
    
end %for %% Looping on the different versions of the experiment

save ('ecogProDataProcess', 'pow_bych_all', 'pow_gamma_all', 'numOfCondition', 'allversions', 'faxis_all', 'faxis_gamma_all', 'reref_mat', 'elecC', 'elecR', 'sizeH', 'sizeV');

    