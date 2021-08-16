clear;
clc
channelNum = 27;

SubBand_str= ["Delta","Theta","Alpha","Sigma","Beta","Gamma"];
for subbandIdx = 3:length(SubBand_str)
%%
Health = [];
Stim  = [];
Medication = [];
Run = [];
Channel = [];
MAE = [];
CorVal = [];
Normalization = [];
Subband = [];

load FeatLabels.mat
desiredFeats = find(contains([CH_Lab{:}],SubBand_str(subbandIdx)));
CH_Lab = CH_Lab(desiredFeats);
featLabels = [];
for chIdx = 1:channelNum
    for featIdx = 1:length(CH_Lab)
        featLabels = cat(1,featLabels,"Channel:"+chIdx+"  "+CH_Lab{featIdx});
    end
end

desiredFeats = desiredFeats+7;
%% PD MinMax
load TestOutIDs_PD.mat


Health_str = "PD";
Normalization_str = "MinMax";

Stim_strS = ["Sham","GVS7","GVS8"];
MedS = [0,1];

for stimIdx = 1:length(Stim_strS)
    for medIdx = 1:length(MedS)
        load NormalizedFeatures.mat
        
        Stim_str = Stim_strS(stimIdx);
        Med = MedS(medIdx);
        disp(Health_str+" -- "+Stim_str+" -- Med: "+Med+" -- "+Normalization_str)
        
        normalizedFeats = normalizedFeats(normalizedFeats.Stim==Stim_str &...
                                          normalizedFeats.Health==Health_str &...
                                          normalizedFeats.Medication==Med,:);
        [CorrLasso,ACCs_Lasso,~,FeatIdx2] = LassoRegression(normalizedFeats,testOutIDs,desiredFeats);
        FeatIdx2 = squeeze(sum(reshape(sign(abs(FeatIdx2)),size(testOutIDs,1),[],channelNum),2));
        Run = cat(1,Run,(1:size(testOutIDs,1))');
        Health = cat(1,Health,repmat(Health_str,size(testOutIDs,1),1));
        Stim  = cat(1,Stim,repmat(Stim_str,size(testOutIDs,1),1));
        Medication = cat(1,Medication,repmat(Med,size(testOutIDs,1),1));
        Normalization = cat(1,Normalization,repmat(Normalization_str,size(testOutIDs,1),1));
        MAE = cat(1,MAE,ACCs_Lasso(:,2));
        CorVal = cat(1,CorVal,CorrLasso(:,2));
        Channel = cat(1,Channel,FeatIdx2);
        Subband = cat(1,Subband,repmat(SubBand_str(subbandIdx),size(testOutIDs,1),1));
    end
end


%% PD ZScore
load TestOutIDs_PD.mat


Health_str = "PD";
Normalization_str = "ZScore";

Stim_strS = ["Sham","GVS7","GVS8"];
MedS = [0,1];

for stimIdx = 1:length(Stim_strS)
    for medIdx = 1:length(MedS)
        load ZScoredFeatures.mat
        
        Stim_str = Stim_strS(stimIdx);
        Med = MedS(medIdx);
        disp(Health_str+" -- "+Stim_str+" -- Med: "+Med+" -- "+Normalization_str)

        normalizedFeats = normalizedFeats(normalizedFeats.Stim==Stim_str &...
                                          normalizedFeats.Health==Health_str &...
                                          normalizedFeats.Medication==Med,:);
        [CorrLasso,ACCs_Lasso,~,FeatIdx2] = LassoRegression(normalizedFeats,testOutIDs,desiredFeats);
        FeatIdx2 = squeeze(sum(reshape(sign(abs(FeatIdx2)),size(testOutIDs,1),[],channelNum),2));
        Run = cat(1,Run,(1:size(testOutIDs,1))');
        Health = cat(1,Health,repmat(Health_str,size(testOutIDs,1),1));
        Stim  = cat(1,Stim,repmat(Stim_str,size(testOutIDs,1),1));
        Medication = cat(1,Medication,repmat(Med,size(testOutIDs,1),1));
        Normalization = cat(1,Normalization,repmat(Normalization_str,size(testOutIDs,1),1));
        MAE = cat(1,MAE,ACCs_Lasso(:,2));
        CorVal = cat(1,CorVal,CorrLasso(:,2));
        Channel = cat(1,Channel,FeatIdx2);
        Subband = cat(1,Subband,repmat(SubBand_str(subbandIdx),size(testOutIDs,1),1));
    end
end

%% HC MinMax
load TestOutIDs_HC.mat


Health_str = "HC";
Normalization_str = "MinMax";

Stim_strS = ["Sham","GVS7","GVS8"];
MedS = 0;

for stimIdx = 1:length(Stim_strS)
    for medIdx = 1:length(MedS)
        load NormalizedFeatures.mat
        
        Stim_str = Stim_strS(stimIdx);
        Med = MedS(medIdx);
        disp(Health_str+" -- "+Stim_str+" -- Med: "+Med+" -- "+Normalization_str)

        normalizedFeats = normalizedFeats(normalizedFeats.Stim==Stim_str &...
                                          normalizedFeats.Health==Health_str &...
                                          normalizedFeats.Medication==Med,:);
        [CorrLasso,ACCs_Lasso,~,FeatIdx2] = LassoRegression(normalizedFeats,testOutIDs,desiredFeats);
        FeatIdx2 = squeeze(sum(reshape(sign(abs(FeatIdx2)),size(testOutIDs,1),[],channelNum),2));
        Run = cat(1,Run,(1:size(testOutIDs,1))');
        Health = cat(1,Health,repmat(Health_str,size(testOutIDs,1),1));
        Stim  = cat(1,Stim,repmat(Stim_str,size(testOutIDs,1),1));
        Medication = cat(1,Medication,repmat(Med,size(testOutIDs,1),1));
        Normalization = cat(1,Normalization,repmat(Normalization_str,size(testOutIDs,1),1));
        MAE = cat(1,MAE,ACCs_Lasso(:,2));
        CorVal = cat(1,CorVal,CorrLasso(:,2));
        Channel = cat(1,Channel,FeatIdx2);
        Subband = cat(1,Subband,repmat(SubBand_str(subbandIdx),size(testOutIDs,1),1));
    end
end


%% HC ZScore
load TestOutIDs_HC.mat


Health_str = "HC";
Normalization_str = "ZScore";

Stim_strS = ["Sham","GVS7","GVS8"];
MedS = 0;

for stimIdx = 1:length(Stim_strS)
    for medIdx = 1:length(MedS)
        load ZScoredFeatures.mat
        
        Stim_str = Stim_strS(stimIdx);
        Med = MedS(medIdx);
        disp(Health_str+" -- "+Stim_str+" -- Med: "+Med+" -- "+Normalization_str)

        normalizedFeats = normalizedFeats(normalizedFeats.Stim==Stim_str &...
                                          normalizedFeats.Health==Health_str &...
                                          normalizedFeats.Medication==Med,:);
        [CorrLasso,ACCs_Lasso,~,FeatIdx2] = LassoRegression(normalizedFeats,testOutIDs,desiredFeats);
        FeatIdx2 = squeeze(sum(reshape(sign(abs(FeatIdx2)),size(testOutIDs,1),[],channelNum),2));
        Run = cat(1,Run,(1:size(testOutIDs,1))');
        Health = cat(1,Health,repmat(Health_str,size(testOutIDs,1),1));
        Stim  = cat(1,Stim,repmat(Stim_str,size(testOutIDs,1),1));
        Medication = cat(1,Medication,repmat(Med,size(testOutIDs,1),1));
        Normalization = cat(1,Normalization,repmat(Normalization_str,size(testOutIDs,1),1));
        MAE = cat(1,MAE,ACCs_Lasso(:,2));
        CorVal = cat(1,CorVal,CorrLasso(:,2));
        Channel = cat(1,Channel,FeatIdx2);
        Subband = cat(1,Subband,repmat(SubBand_str(subbandIdx),size(testOutIDs,1),1));
    end
end


T = table(Health,Medication,Stim,Normalization,Run,Subband,MAE,CorVal);
T1 = array2table(Channel);
T = cat(2,T,T1);
writetable(T,"LassoRegResults"+SubBand_str(subbandIdx)+".csv")
end
