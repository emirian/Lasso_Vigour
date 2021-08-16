########################### Initialization ##########################
library(pacman)
p_load(reshape2,
       ez,
       lme4,
       lmerTest,
       ggplot2,
       grid,
       tidyr,
       plyr,
       dplyr,
       effects,
       gridExtra,
       DescTools,
       Cairo, #alternate image writing package with superior performance.
       corrplot,
       knitr,
       PerformanceAnalytics,
       afex,
       ggpubr,
       readxl,
       officer,
       psych,
       rstatix,
       emmeans,
       export,
       standardize,
       performance,stringr)


###################################################### Load data -----
d1 = read.table("LassoRegResultsAlpha.csv", header=TRUE, sep=",", strip.white = TRUE)
d2 = read.table("LassoRegResultsBeta.csv", header=TRUE, sep=",", strip.white = TRUE)
d3 = read.table("LassoRegResultsDelta.csv", header=TRUE, sep=",", strip.white = TRUE)
d4 = read.table("LassoRegResultsGamma.csv", header=TRUE, sep=",", strip.white = TRUE)
d5 = read.table("LassoRegResultsSigma.csv", header=TRUE, sep=",", strip.white = TRUE)
d6 = read.table("LassoRegResultsTheta.csv", header=TRUE, sep=",", strip.white = TRUE)
d7 = read.table("LassoRegResults.csv", header=TRUE, sep=",", strip.white = TRUE)
d7$Subband = unique("Full")
DatAll = rbind(d1,d2,d3,d4,d5,d6)

d7 = d7[,c("Health", "Medication","Stim","Normalization", "Run","Subband", "MAE", "CorVal")]
DatAll = DatAll[,c("Health", "Medication","Stim","Normalization", "Run","Subband", "MAE", "CorVal")]
DatAll = rbind(d7,DatAll)

DatAll$CorVal[DatAll$CorVal=="NaN"] = NA
DatAll = DatAll[complete.cases(DatAll),]
DatAll$CorVal = as.numeric(DatAll$CorVal)

DatAll$Subband = factor(DatAll$Subband,levels = c("Full","Delta","Theta","Alpha","Sigma","Beta","Gamma"),
                        labels = c("BB","Delta","Theta","AlphaL","AlphaH","Beta","Gamma"))
DatAll$Stim = factor(DatAll$Stim, levels = c("Sham","GVS7","GVS8"))

NameMap  = read.table("FeatNameMap.csv", header=TRUE, sep=",", strip.white = TRUE)
dFeat = read.table("zScoreNormalizedFeats.csv", header=TRUE, sep=",", strip.white = TRUE)
dFeat = melt(dFeat, id.vars = c("SID","Health", "Medication","Stim","Channel", "Trial","Vigor"),
                variable.name = "Feature")

dFeat = as.data.frame(summarise(group_by(dFeat,SID,Health,Medication,Stim,Trial,Feature),value = mean(value,na.rm=T)))
dFeat_NoTrial = as.data.frame(summarise(group_by(dFeat,SID,Health,Medication,Stim,Feature),value = mean(value,na.rm=T)))

dFeat_NoTrial$FeatureOrg = dFeat_NoTrial$Feature 
dFeat_NoTrial$Feature = factor(dFeat_NoTrial$Feature, levels = NameMap$Feats,labels = NameMap$Band)
dFeat_NoTrial$Stim = factor(dFeat_NoTrial$Stim, levels = c("Sham","GVS7","GVS8"))
dFeat_NoTrial = dFeat_NoTrial[dFeat_NoTrial$Feature %in% c("Delta","Theta","AlphaL","AlphaH","Beta","Gamma"),]


###################################################### Accuracies -----
graphdat = DatAll[DatAll$Medication==0 & DatAll$Normalization =="ZScore",c("Health","Stim","Run","Subband", "CorVal")]

graphdat = as.data.frame(summarise(group_by(graphdat,Stim,Health,Run,Subband),N=n(),CorVal = mean(CorVal,na.rm=T)))


Fig1 = ggplot(graphdat, aes(x=Health, y=CorVal, fill=Subband)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Subband))+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size = .5)+
  facet_grid(~Stim)+
  ggtitle("Lasso Performance")

plot(Fig1)

graph2ppt(file="Fig1.pptx",width = 9, height = 5)


ggplot(graphdat, aes(x=Subband, y=CorVal, fill=Stim)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
  #                                             jitter.height = 0,
  #                                             dodge.width = .9),shape = 21,fill="grey",aes(colour = Subband))+
  facet_grid(~Health)+
  ggtitle("Lasso Performance")


as.data.frame(summarise(group_by(graphdat,Stim,Health,Subband),N=n(),CorVal = mean(CorVal,na.rm=T)))

Sdat = reshape2::dcast(graphdat,Subband+Stim+Health ~ Run, value.var="CorVal")
Sdat = melt(Sdat, id.vars = c("Health","Stim","Subband"),
                variable.name = "Run")

Sdat$value[is.na(Sdat$value)]=unique(0)
results=as.data.frame(ezANOVA(data=Sdat, dv="value", wid=.("Run"), within=.("Subband"),
                              between = c("Stim","Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results

Report = as.data.frame(summarise(group_by(Sda,Stim,Health,Subband),M = round(mean(value,na.rm=T),2),SD = round(sd(value,na.rm=T),2)))
write.csv(Report,"Fig1CorrelationFullReportTable.csv",row.names = F)

#----------------------- main effect Stim
SdatStim = as.data.frame(summarise(group_by(Sdat,Stim,Health,Run),value = mean(value,na.rm=T)))
SdatStim = as.data.frame(summarise(group_by(SdatStim,Stim,Run),value = mean(value,na.rm=T)))

as.data.frame(summarise(group_by(SdatStim,Stim),M = round(mean(value,na.rm=T),2),SD = round(sd(value,na.rm=T),2) ))
pairwise_t_test(SdatStim,value~Stim,paired = T)

ggplot(SdatStim, aes(x=Stim, y=value, fill=Stim)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")
#----------------------- main effect Band
SdatBand = as.data.frame(summarise(group_by(Sdat,Subband,Health,Run),value = mean(value,na.rm=T)))
SdatBand = as.data.frame(summarise(group_by(SdatBand,Subband,Run),value = mean(value,na.rm=T)))

as.data.frame(summarise(group_by(SdatBand,Subband),M = round(mean(value,na.rm=T),2),SD = round(sd(value,na.rm=T),2) ))
Compares = as.data.frame(pairwise_t_test(SdatBand,value~Subband,paired = T))
write.csv(Compares, "tempdat.csv")

ggplot(SdatBand, aes(x=Subband, y=value, fill=Subband)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")

#----------------------- interaction Health:subband

SdatHealthSubband = as.data.frame(summarise(group_by(Sdat,Health,Subband,Run),value = mean(value,na.rm=T)))


ggplot(SdatHealthSubband, aes(x=Subband, y=value, fill=Health)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")

SdatHealthSubband$Test = paste(SdatHealthSubband$Health,SdatHealthSubband$Subband,sep = "_")

Compares = pairwise_t_test(SdatHealthSubband,value~Test,paired = T)
write.csv(Compares, "tempdat.csv")

#----------------------- interaction Stim:Health:subband


Sdat$Test = paste(Sdat$Health,Sdat$Subband,Sdat$Stim,sep = "_")

Compares = pairwise_t_test(Sdat,value~Test,paired = T)
write.csv(Compares, "tempdat.csv")

ggplot(Sdat, aes(x=Subband, y=value, fill=Health)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Health))+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size = .5)+
  facet_grid(~Stim)+
  stat_compare_means(method = "t.test",paired = FALSE )
##----------------------------- Test stat ---------
sDat = graphdat[complete.cases(graphdat),]
sDat$Test = paste(sDat$Health,sDat$Stim,sDat$Subband,sep = "_")
pwc <- sDat  %>%
  pairwise_t_test(CorVal ~ Test, pool.sd = T, paired = F,
                  p.adjust.method = "bonferroni")
pwc$p = round(pwc$p,3)
pwc$p.adj = round(pwc$p.adj,3)
Result = as.data.frame(pwc)

write.csv(Result,"Ttest_BandResults.csv",row.names = F)




###################################################### Feature Values -----
graphdat = dFeat_NoTrial[dFeat_NoTrial$Medication==0 ,c("SID","Health","Stim","Feature","FeatureOrg", "value")]
graphdat$Feature = as.character(graphdat$Feature)

graphdat$FeatureOrg = paste(graphdat$FeatureOrg,graphdat$Feature,sep = "_")
Sdat = reshape2::dcast(graphdat,SID+Stim+Health ~ FeatureOrg, value.var="value")
Sdat[!complete.cases(Sdat),]
Sdat = melt(Sdat, id.vars = c("Health","SID","Stim"),
            variable.name = "Feature")

results=as.data.frame(ezANOVA(data=Sdat, dv="value", wid=.("SID"), within=.("Stim","Feature"),
                              between = c("Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results

graphdat = as.data.frame(summarise(group_by(graphdat,SID,Health,Stim,Feature),value = mean(value,na.rm=T)))

results=as.data.frame(ezANOVA(data=graphdat, dv="value", wid=.("SID"), within=.("Stim","Feature"),
                              between = c("Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results


ggplot(graphdat, aes(x=Feature, y=value, fill=Stim)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
  #                                             jitter.height = 0,
  #                                             dodge.width = .9),shape = 21,fill="grey",aes(colour = Subband))+
  facet_grid(~Health)+
  ggtitle("Features")

results=as.data.frame(ezANOVA(data=graphdat, dv="value", wid=.("SID"), within=.("Stim","Feature"),
                              between = c("Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results

##----------------------------- Test stat ---------
sDat = graphdat[complete.cases(graphdat),]
sDat$Test = paste(sDat$Health,sDat$Stim,sDat$Feature,sep = "_")
pwc <- sDat  %>%
  pairwise_t_test(value ~ Test, pool.sd = T, paired = F,
                  p.adjust.method = "bonferroni")
pwc$p = round(pwc$p,3)
pwc$p.adj = round(pwc$p.adj,3)
Result = as.data.frame(pwc)

write.csv(Result,"Ttest_FeaturesResults.csv",row.names = F)

######################################################
graphdat = d[d$Health=="PD",c("Health", "Medication","Stim",
                              "Normalization", "Run", "MAE", "CorVal")]

graphdat = melt(graphdat, id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
                variable.name = "Performance")

ggplot(graphdat, aes(x=Normalization, y=value, fill=Stim)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Stim))+
  ggtitle("Med On=1 and off=0")+
  facet_grid(Performance~Medication)


graphdat = d[(d$Health=="PD" & d$Medication==1)|d$Health=="HC",c("Health", "Medication","Stim",
                                                                 "Normalization", "Run", "MAE", "CorVal")]

graphdat = melt(graphdat, id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
                variable.name = "Performance")

ggplot(graphdat, aes(x=Normalization, y=value, fill=Stim)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Stim))+
  ggtitle("HC vs PD Med on")+
  facet_grid(Performance~Health)


graphdat = d[,c("Health", "Medication","Stim","Normalization", "Run", "MAE", "CorVal")]
graphdat$Health = paste(graphdat$Health,graphdat$Medication,sep = "")

graphdat = melt(graphdat, id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
                variable.name = "Performance")

ggplot(graphdat, aes(x=Stim, y=value, fill=Health)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Health))+
  ggtitle("HC vs PD Med on")+
  facet_grid(Performance~Normalization)






###################################################### Features -----

Sdat = melt(d, id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
            variable.name = "Feature")
Sdat = Sdat[!(Sdat$Feature %in% c("MAE","CorVal")),]

Sdat$value[abs(Sdat$value)<0.001]=NA
Sdat = Sdat[complete.cases(Sdat),]

Sdat$Band = case_when(grepl("Delta",Sdat$Feature)~"Delta",
                      grepl("Theta",Sdat$Feature)~"Theta",
                      grepl("Alpha",Sdat$Feature)~"Alpha",
                      grepl("Sigma",Sdat$Feature)~"Sigma",
                      grepl("Beta",Sdat$Feature)~"Beta",
                      grepl("Gamma",Sdat$Feature)~"Gamma")

Sdat$Channel = str_extract(Sdat$Feature,"Channel.[0-9]+")



FeatDat = as.data.frame(summarise(group_by(Sdat,Health,Medication,Stim,Normalization,Feature,Band,Channel),N=n(),value = mean(value)))
write.csv(FeatDat,"DetailedFeaturesBetavals.csv")

HeadplotDat = as.data.frame(summarise(group_by(Sdat,Health,Medication,Stim,Normalization,Feature,Band,Channel),N=n(),value = sum(value)))
HeadplotDat$value = HeadplotDat$value*HeadplotDat$N
HeadplotDat = as.data.frame(summarise(group_by(HeadplotDat,Health,Medication,Stim,Normalization,Band,Channel),N=sum(N),value = sum(value)))

write.csv(HeadplotDat,"BetaValsHeadPlot.csv")


FeatDat = FeatDat[FeatDat$N>20,]
FeatDat$value = sign(FeatDat$value)

ggplot(FeatDat, aes(x=Band, y=value, fill=Health)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                              jitter.height = 0,
                                              dodge.width = .9),shape = 21,fill="grey",aes(colour = Health))+
  facet_grid(Stim~Normalization+Medication)

BandDat = as.data.frame(summarise(group_by(FeatDat,Health,Medication,Stim,Normalization,Band),value=mean(value)))
write.csv(BandDat,"Band_DatMorethan37.Csv")

SpatialDat = as.data.frame(summarise(group_by(FeatDat,Health,Medication,Stim,Normalization,Band,Channel),value=n()))
SpatialDat$Channel = gsub("Channel.","",SpatialDat$Channel)

SpatialDat = reshape2::dcast(SpatialDat,Normalization+Stim+Health+Medication+Band ~ Channel, value.var="value")
write.csv(SpatialDat,"Spatial_DatMorethan37.Csv")
# 
# ggplot(SpatialDat, aes(x=Band, y=value, fill=Channel)) + 
#     geom_bar(stat="summary",fun="mean",position="dodge")+
#     # stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
#     geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
#                                                 jitter.height = 0,
#                                                 dodge.width = .9),shape = 21,fill="grey",aes(colour = Channel))+
#     ggtitle("HC vs PD Med on")+
#     facet_grid(Stim+Health~Normalization+Medication)

################################################



testDat = graphdat[graphdat$Subspace %in% c("AllperCh","LSTM","LSTM_Reg","AllChPCA")
                   & graphdat$RES==3 & graphdat$Method == "RF",]        

testDat$Subspace = factor(testDat$Subspace, levels = c("AllperCh","AllChPCA","LSTM","LSTM_Reg"),
                          labels = c("AllPerCh_Feats","PCA_Feats","LSTM_Feats","LSTM_LMM_Feats"))
ggplot(testDat, aes(x=Subspace, y=ACC, fill=Subspace)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
  labs(x="Feature Subspace",y="3-level RT Classification Accuracy", size=16)

graph2ppt(file="Fig4.pptx",width = 9, height = 5)

t.test(testDat$ACC[testDat$Subspace=="AllPerCh_Feats"],testDat$ACC[testDat$Subspace=="PCA_Feats"],paired = T)
t.test(testDat$ACC[testDat$Subspace=="AllPerCh_Feats"],testDat$ACC[testDat$Subspace=="LSTM_Feats"],paired = T)
t.test(testDat$ACC[testDat$Subspace=="AllPerCh_Feats"],testDat$ACC[testDat$Subspace=="LSTM_LMM_Feats"],paired = T)

t.test(testDat$ACC[testDat$Subspace=="PCA_Feats"],testDat$ACC[testDat$Subspace=="LSTM_Feats"],paired = T)
t.test(testDat$ACC[testDat$Subspace=="PCA_Feats"],testDat$ACC[testDat$Subspace=="LSTM_LMM_Feats"],paired = T)


testDat = graphdat[graphdat$Method=="RF",]
t.test(testDat$ACC[testDat$Subspace=="LSTM"],testDat$ACC[testDat$Subspace=="AllChLasso"],paired = T)

