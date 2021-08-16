clear
clc
load ForMatlabNormalization.mat

variableIndexes = 7:76;
for sid = unique(dat.SID)'
    for med = unique(dat.Medication)'
        for stim = unique(dat.Stim)'
            for ch = unique(dat.Channel)'
                temp = table2array(dat(dat.SID==sid & dat.Medication == med & dat.Stim==stim & dat.Channel==ch,variableIndexes));
                temp = (temp-repmat(min(temp),size(temp,1),1))./(repmat(max(temp),size(temp,1),1)-repmat(min(temp),size(temp,1),1));
                dat(dat.SID==sid & dat.Medication == med & dat.Stim==stim & dat.Channel==ch,variableIndexes) = array2table(temp);
            end
        end
    end
end

writetable(dat,"MMnormalizedFeats.csv");
%%
clear
clc
load ForMatlabNormalization.mat

variableIndexes = 7:76;
for sid = unique(dat.SID)'
    for med = unique(dat.Medication)'
        for stim = unique(dat.Stim)'
            for ch = unique(dat.Channel)'
                temp = table2array(dat(dat.SID==sid & dat.Medication == med & dat.Stim==stim & dat.Channel==ch,variableIndexes));
                temp = (temp-repmat(mean(temp),size(temp,1),1))./(repmat(std(temp),size(temp,1),1));
                dat(dat.SID==sid & dat.Medication == med & dat.Stim==stim & dat.Channel==ch,variableIndexes) = array2table(temp);
            end
        end
    end
end

writetable(dat,"zScoreNormalizedFeats.csv");