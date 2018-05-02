%% This Script is doing the Data preprocessing; read the data, divide to allEpochs and save the results.
% clear all;
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

sd = 1; %Thomas change this to 0
useCh_1_60 = 1; %use channels 1 till 60
saveData   = 1; % save the data results 

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
numOfCondition = length(allversions);

rawData = cell(numOfCondition, 1); % pre allocate memory
allEpochs=[];
%% Looping on the different versions of the experiment
for nversion=1:numOfCondition
%     S1channels= (startElct:endElct); %[11 21]; % channels given the strongest frequency-tag responses 1:60
    S1channels = [11 12 13 14 15 16 21 22 23 24 25 26 37 38 39 40 47 48 49 50 59 60];
    version=allversions{nversion};
    if strcmp(version,'1a') % version 1a: sensors on both hands
        data_path= [preDataPath '1a'];
        file_name='Vibrotactile_Regular_16_03_2018_19_05_12_0000.mat';
        
% %         S1channels= (1:60); %[11 21]; % channels given the strongest frequency-tag responses 1:60
        LeftF0_all (:, nversion)=[8 11]; % frequencies applied to the left hand
        RightF0_all(:, nversion)=[17 23]; % frequencies applied to the right hand
        
    elseif strcmp(version,'1arev') % version 1a: sensors on both hands, frequencies flipped compared to 1a
        data_path= [preDataPath '1arev/'];
        file_name= 'Vibrotactile_Regular_22_03_2018_11_40_37_0000.mat';
        
% %         S1channels=[41 51];
        LeftF0_all(:, nversion)=[23 17];
        RightF0_all(:, nversion)=[11 8];
        
    elseif strcmp(version,'1b') % version 1a: sensors on both hands, hands crossed
        data_path=[preDataPath '1b/'];
        file_name='Vibrotactile_Crossed_16_03_2018_19_14_05_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0_all(:, nversion)=[8 11];
        RightF0_all(:, nversion)=[17 23];
        
    elseif strcmp(version,'2a') % version 1a: sensors on hands and feet kept close
        data_path=[preDataPath '2a/'];
        file_name='Vibrotactile_FeetClose_17_03_2018_13_11_09_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0_all(:, nversion)=[8 17];
        RightF0_all(:, nversion)=[11 23];
        
    elseif strcmp(version,'2b') % version 1a: sensors on hands and feet kept apart
        data_path=[preDataPath '2b/'];
        file_name='Vibrotactile_FeetApart_17_03_2018_13_22_10_0000.mat';
        
% %         S1channels=[11 21];
        LeftF0_all(:, nversion)=[8 17];
        RightF0_all(:, nversion)=[11 23];
        
    end
    S2channels=[];
    
    fprintf('... loading %s\n',file_name) %reverse the order sd
    rawData {nversion} =  load([data_path filesep file_name]);
    y = rawData{nversion}.y;
    SR = rawData{nversion}.SR;
    %% Identify trials (trial condition and start/end)
    groups=1:4; % 1 to 5 are digits from thumb to pinkie and 6 is rest
    groupid=squeeze(y(164,1,:));    
    epocTemp = [];
    for ngroup=[1 3]
        temp=groupid==ngroup;
        start=find(diff(temp)==1)+1;
        ending=find(diff(temp)==-1);
        duration=ending-start;
        meanSec = ceil(mean(duration/SR));
        fprintf('... ... cond %g ... found %g allEpochs of length %2.2fs on average\n',ngroup,length(start), meanSec)
        
        epocTemp =[epocTemp ; [ngroup*ones(length(start),1) start start+meanSec*SR start-SR]]; 
    end
    
    allEpochs(:,:, nversion) =  epocTemp;
    
end %for  nversion=1:5

if(saveData)
    save('ecogGLobal.mat','numOfCondition', 'useCh_1_60', 'preDataPath', 'meanSec', 'chronuxPath', 'S1channels', 'LeftF0_all', 'RightF0_all' );
    save ('ecogPrePros.mat', 'allEpochs', 'allversions', 'rawData');
end %(if saveData)
    