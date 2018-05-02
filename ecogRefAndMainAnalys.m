% This script load the epochs data, making the referencing and the main/heavy analysis

% clear all;


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
    erpbp=nan(size(erp));    
    for ntr=1:size(epochs,1)
        for nch=1:size(erp,1)
            erpbp(nch,ntr,:)=abs(hilbert(bandpass(squeeze(erp(nch,ntr,:)), SR, 60, 90, 4)));
            erpbp(nch,ntr,:)=erpbp(nch,ntr,:)-mean(erpbp(nch,ntr,times<0));
            [faxis_gamma_all, pow]=get_PowerSpec(squeeze(erpbp(nch,ntr,times>0)),SR,0,0);
            pow_gamma_all(nch,ntr,:, nversion)=pow;
        end
    end
    
    
    
    %%
    for nch=1:size(erp,1)
        for ntr=1:size(erp,2)
            [faxis_all(:,nversion), pow_bych_all(nch,ntr,:, nversion)]=get_PowerSpec_new(squeeze(erp(nch,ntr,:)),1/SR,size(erp,3)/SR,0,0);
        end
    end
    
end %for %% Looping on the different versions of the experiment

save ('ecogProDataProcess', 'pow_bych_all', 'pow_gamma_all', 'numOfCondition', 'allversions', 'faxis_all', 'faxis_gamma_all', 'reref_mat', 'elecC', 'elecR', 'sizeH', 'sizeV');

    