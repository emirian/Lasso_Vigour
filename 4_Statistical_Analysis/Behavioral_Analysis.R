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
       psych,
       rstatix,
       emmeans,
       standardize,
       performance,
       stringr,
       ggeffects,
       sjPlot,
       sjmisc)


############################################# Normal Models  Behavioral Results section-----------
DatAll = read.csv("BehaviorData.csv",header = T)

ID_Ignore = c("HC8","HC9","PD3","PD7")


VigorDat = DatAll[DatAll$Behavior=="PeakTime",]
VigorDat = VigorDat[!VigorDat$SID%in%ID_Ignore,]


VigorDat$Value = 1/VigorDat$Value

VigorDat = VigorDat[complete.cases(VigorDat),]
VigorDat$PDDiag = case_when(VigorDat$Health=="HC"~0,
                            VigorDat$Health=="PD"~1)
VigorDat$Severity[VigorDat$Health=="HC"]=unique(1)


VigorDat$Health = as.factor(VigorDat$Health)
VigorDat$SID = as.factor(VigorDat$SID)
VigorDat$Stim = as.factor(VigorDat$Stim)

Graphdat = as.data.frame(summarise(group_by(VigorDat,SID,Age,Sex,Severity,Health,Stim,Behavior),Vigor = mean(Value,na.rm=T)))

p1 = ggplot(Graphdat,aes(x=Stim , y=Vigor, fill = Health)) +
        # geom_bar(stat="summary",fun="mean",position="dodge")+
        geom_jitter(position = position_jitterdodge(jitter.width = NULL,
                                                    jitter.height = 0,
                                                    dodge.width = .75),shape = 21,fill="grey",aes(colour = Health))+
        # geom_boxplot(outlier.colour="black",outlier.shape=8,
        #              size=.5,fill = NA,aes(colour = condition))+
        stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=1,
                     aes(colour = Health))+
        theme_bw(base_family = "serif")+
        stat_compare_means(method = "t.test",paired = F)

plot(p1)


Sdat = reshape2::dcast(Graphdat,SID+Health ~ Stim, value.var="Vigor")
Sdat = melt(Sdat, id.vars = c("SID","Health"),
            variable.name = "Stim")

results=as.data.frame(ezANOVA(data=Sdat, dv="value", wid=.("SID"), within=.("Stim"),
                              between = c("Health"), type=3,detailed=T)$ANOVA)
results$pareta=results$SSn/(results$SSn+results$SSd)
is.num=sapply(results, is.numeric)
results[is.num] =lapply(results[is.num], round, 3)
results

Sdat = as.data.frame(summarise(group_by(Sdat,SID,Stim),value = mean(value,na.rm=T)))
as.data.frame(summarise(group_by(Sdat,Stim),M = round(mean(value,na.rm=T),5), sd  =round(sd(value,na.rm=T),5)))

ggplot(Sdat,aes(x=Stim , y=value, fill = Stim)) +
        geom_bar(stat="summary",fun="mean",position="dodge")+
        geom_jitter(position = position_jitterdodge(jitter.width = NULL,
                                                    jitter.height = 0,
                                                    dodge.width = .75),shape = 21,fill="grey",aes(colour = Stim))+
        stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",size=1)

t.test(Sdat$value[Sdat$Stim=="Sham"],Sdat$value[Sdat$Stim=="GVS7"],paired = T)
t.test(Sdat$value[Sdat$Stim=="Sham"],Sdat$value[Sdat$Stim=="GVS8"],paired = T)
t.test(Sdat$value[Sdat$Stim=="GVS7"],Sdat$value[Sdat$Stim=="GVS8"],paired = T)




