
clear;
clc;

%%
load Vigordata_AK.mat % PeakTime is saved at index = 12

% Please load these cell arrays: dataHCTask  dataPD1Task   dataPD2Task
load ("C:\Elham\EEG-PatientIdentification\SOOJIN_DATA\PT_TimeLocked_1000Back.mat")
%% Parameters
Fs = 500;
ChanNum = 27;
Phi = 1;
Rho = 0;
nLag = 150; 
nFFT = 512;

%% Variables for the final table
SID     = []; % Subject ID
Health  = []; % HC, PD1, PD2
Stim    = []; % Sham, GVS7, GVS8
Vigor     = [];
Feature = []; % feature values
Channel = []; % Channel number
OddEven = []; %Odd or Even samples
Trial   = []; %Trial number
stimString = ["Sham","GVS7","GVS8"];

%% Features for dataHCTask
healthStr = "HC";
bData = B_HC;
% dat = dataHCTask2([1,7,8],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 
dat = dataHCTask2([1,2,3],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 
for subIdx = 1:size(dat,2)   % loop on subjects
    for stimIdx = 1:size(dat,1)  % loop on stimuli
        signal = dat{stimIdx,subIdx}; % signal should be a 27x1000x10
        if(isempty(signal))
            continue;
        end
        signalOdd = signal(:,1:2:end,:); % signalOdd should be a 27x500x10
        signalEven = signal(:,2:2:end,:); % signalEven should be a 27x500x10
        
        for chanIdx = 1:size(signal,1) % loop on channels
            for trialIdx = 1:size(signal,3) % loop on channels
                VigorValue = 1/bData{stimIdx,subIdx}(trialIdx,12);
                Sig = signalOdd(chanIdx,:,trialIdx);
                oddEvenStr = "Odd";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
                Sig = signalEven(chanIdx,:,trialIdx);
                oddEvenStr = "Even";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
            end
        end
    end
end
DataHC = table(SID,Health,Stim,Vigor, Channel,Trial,OddEven,Feature);
writetable(DataHC,"FeatureDataHC.csv")

%% Features for dataPD1Task
healthStr = "PD1";
bData = B_PD1;
%dat = dataPD1Task2([1,7,8],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 
dat = dataPD1Task2([1,2,3],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 

for subIdx = 1:size(dat,2)   % loop on subjects
    for stimIdx = 1:size(dat,1)  % loop on stimuli
        signal = dat{stimIdx,subIdx}; % signal should be a 27x1000x10
        if(isempty(signal))
            continue;
        end
        signalOdd = signal(:,1:2:end,:); % signalOdd should be a 27x500x10
        signalEven = signal(:,2:2:end,:); % signalEven should be a 27x500x10
        
        for chanIdx = 1:size(signal,1) % loop on channels
            for trialIdx = 1:size(signal,3) % loop on channels
                VigorValue = 1/bData{stimIdx,subIdx}(trialIdx,12);
                Sig = signalOdd(chanIdx,:,trialIdx);
                oddEvenStr = "Odd";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
                Sig = signalEven(chanIdx,:,trialIdx);
                oddEvenStr = "Even";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
            end
        end
    end
end
DataPD1 = table(SID,Health,Stim,Vigor, Channel,Trial,OddEven,Feature);
writetable(DataPD1,"FeatureDataPD1.csv")

%% Features for dataPD1Task
healthStr = "PD2";
bData = B_PD2;
%dat = dataPD2Task2([1,7,8],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 
dat = dataPD2Task2([1,2,3],:); % Assuming that dataHCTask is a Cell array (Stim X Subjects) 

for subIdx = 1:size(dat,2)   % loop on subjects
    for stimIdx = 1:size(dat,1)  % loop on stimuli
        signal = dat{stimIdx,subIdx}; % signal should be a 27x1000x10
        if(isempty(signal))
            continue;
        end
        signalOdd = signal(:,1:2:end,:); % signalOdd should be a 27x500x10
        signalEven = signal(:,2:2:end,:); % signalEven should be a 27x500x10
        
        for chanIdx = 1:size(signal,1) % loop on channels
            for trialIdx = 1:size(signal,3) % loop on channels
                VigorValue = 1/bData{stimIdx,subIdx}(trialIdx,12);
                Sig = signalOdd(chanIdx,:,trialIdx);
                oddEvenStr = "Odd";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
                Sig = signalEven(chanIdx,:,trialIdx);
                oddEvenStr = "Even";
                FeatEx;  %-------> Extract Features
                CatAllVariables;%-----------> Concatenate variables
            end
        end
    end
end
DataPD2 = table(SID,Health,Stim,Vigor,Channel,Trial,OddEven,Feature);
writetable(DataPD2,"FeatureDataPD2.csv")