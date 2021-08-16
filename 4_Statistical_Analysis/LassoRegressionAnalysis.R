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
       standardize,
       performance,stringr)


###################################################### Load data -----
d=read.table("LassoRegResults.csv", header=TRUE, sep=",", strip.white = TRUE)

d$CorVal[d$CorVal=="NaN"] = NA
d = d[complete.cases(d),]
d$CorVal = as.numeric(d$CorVal)



###################################################### Accuracies -----
graphdat = d[d$Medication==0,c("Health", "Medication","Stim",
             "Normalization", "Run", "MAE", "CorVal")]

graphdat = melt(graphdat, id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
                variable.name = "Performance")

ggplot(graphdat, aes(x=Normalization, y=value, fill=Stim)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Stim))+
    ggtitle("Med Off")+
    facet_grid(Performance~Health)



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







###################################################### Features Average Betas-----

Sdat = melt(d[d$Normalization=="ZScore" & d$Medication==0,], id.vars = c("Health", "Medication","Stim","Normalization", "Run"),
            variable.name = "Feature")

Sdat = Sdat[!(Sdat$Feature %in% c("MAE","CorVal")),]

Sdat$value[Sdat$value==0]=NA
Sdat = Sdat[complete.cases(Sdat),]

Sdat$Band = case_when(grepl("Delta",Sdat$Feature)~"Delta",
                      grepl("Theta",Sdat$Feature)~"Theta",
                      grepl("Alpha",Sdat$Feature)~"Alpha",
                      grepl("Sigma",Sdat$Feature)~"Sigma",
                      grepl("Beta",Sdat$Feature)~"Beta",
                      grepl("Gamma",Sdat$Feature)~"Gamma")

Sdat$FType = case_when(grepl("RSP",Sdat$Feature)~"RSP",
                       grepl("HP",Sdat$Feature)~"HP",
                       grepl("BTS_Amplitude",Sdat$Feature)~"BTS_Amplitude",
                       grepl("BTS_Angle",Sdat$Feature)~"BTS_Phase")

Sdat$Band = factor(Sdat$Band,levels = c("Delta","Theta","Alpha","Sigma","Beta","Gamma"),
                        labels = c("Delta","Theta","AlphaL","AlphaH","Beta","Gamma"))
Sdat$Stim = factor(Sdat$Stim, levels = c("Sham","GVS7","GVS8"))

Sdat$Channel = str_extract(Sdat$Feature,"Channel.[0-9]+")

FeatDat = as.data.frame(summarise(group_by(Sdat,Health,Stim,FType,Run,Band,Channel),N=n(),value = mean(value,na.rm=T)))
FeatDat = as.data.frame(summarise(group_by(FeatDat,Health,Stim,FType,Run,Band),N=n(),value = mean(value,na.rm=T)))
FeatDat$value = FeatDat$value*FeatDat$N

ggplot(FeatDat[FeatDat$FType == "RSP",], aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(~Stim)+
    ggtitle("Only RSP features")

ggplot(FeatDat[FeatDat$FType == "HP",], aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(~Stim)+
    ggtitle("Only HP features")

ggplot(FeatDat[FeatDat$FType == "BTS_Amplitude",], aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(~Stim)+
    ggtitle("Only BTS_Amplitude features")

ggplot(FeatDat[FeatDat$FType == "BTS_Phase",], aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(~Stim)+
    ggtitle("Only BTS_Phase features")


ggplot(FeatDat, aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(Stim~FType)+
    ggtitle("Only RSP features")









HeadplotDat = as.data.frame(summarise(group_by(Sdat,Health,Stim,Run,Band,Channel),N=n(),value = mean(value)))
HeadplotDat = as.data.frame(summarise(group_by(HeadplotDat,Health,Stim,Band,Channel),N=n(),value = mean(value)))
write.csv(HeadplotDat,"BetaValsHeadPlot2.csv")


FeatDat = as.data.frame(summarise(group_by(Sdat,Health,Stim,Run,Band,Channel),N=n(),value = mean(value)))
FeatDat = as.data.frame(summarise(group_by(FeatDat,Health,Stim,Run,Band),N=n(),value = mean(value)))

Fig2 = ggplot(FeatDat, aes(x=Health, y=value, fill=Band)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Band))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=.5)+
    facet_grid(~Stim)

plot(Fig2)

graph2ppt(file="Fig2.pptx",width = 9, height = 5)

AVG = as.data.frame(summarise(group_by(FeatDat,Stim,Health,Band),M = round(mean(value),2), SD = round(sd(value),2)))
write.csv(AVG,"Fig2Average_Betavals.csv",row.names = F)

ggplot(FeatDat, aes(x=Band, y=value, fill=Stim)) + 
    geom_bar(stat="summary",fun="mean",position="dodge")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1,
                                                jitter.height = 0,
                                                dodge.width = .9),shape = 21,fill="grey",aes(colour = Stim))+
    stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=1)+
    facet_grid(~Health)

Sdat = reshape2::dcast(FeatDat,Band+Stim+Health ~ Run, value.var="value")
Sdat = melt(Sdat, id.vars = c("Health","Stim","Band"),
            variable.name = "Run")
Sdat$value[is.na(Sdat$value)]=0

results=as.data.frame(ezANOVA(data=Sdat, dv="value", wid=.("Run"), within=.("Band"),
                              between = c("Stim","Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results

as.data.frame(summarise(group_by(FeatDat,Health,Stim,Band),N=n(),M = mean(value),SD = sd(value)))
##----------------------------- Test stat ---------
sDat = FeatDat
sDat$Test = paste(sDat$Health,sDat$Stim,sDat$Band,sep = "_")
Result <- sDat %>%
    group_by(Test) %>%                       
    summarise(res = list(tidy(t.test(value, mu=0)))) %>%
    unnest()
Result$p.value = round(Result$p.value,3)
Result$Sig = unique("")
# Result$p.value = Result$p.value*36
Result$Sig[Result$p.value<0.1] =unique(".") 
Result$Sig[Result$p.value<0.05] =unique("***") 
Result = as.data.frame(Result)
write.csv(Result,"TtestBands_BB.csv",row.names = F)


SpatialDat = as.data.frame(summarise(group_by(FeatDat,Health,Medication,Stim,Normalization,Band,Channel),value=n()))
SpatialDat$Channel = gsub("Channel.","",SpatialDat$Channel)

SpatialDat = reshape2::dcast(SpatialDat,Normalization+Stim+Health+Medication+Band ~ Channel, value.var="value")
write.csv(SpatialDat,"Spatial_DatMorethan37.Csv")

