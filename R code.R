##################################################################################################################
# Code used to conduct analyses for "Redlines and power lines: historical redlining,                             #
# the siting of U.S. fossil fuel power plants, and present-day power plant emissions"                            #
# Code by: Sherlock Li and Lara J. Cushing                                                                       #
# Last Update: 5/22/2022                                                                                         #
##################################################################################################################
#Import Library
library("readxl")
library(aod)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(tidyverse)
library(tidyr)
library("sandwich")
library("lmtest")
library(pscl)
library(MASS)
library(boot)
library(yhat)
library(AER)
library(bestNormalize)
library(lme4)
library(forestplot)
library(reshape2)

#Change the path to where the data is
#setwd("Change here")

#Change the dataset to generate estimates for 5km or 10km (indicated by change here)
#Create stratified data
EIA1939_time <-read.csv("PP_1939_5km_Wind.csv")  #change here
EIA1939_time$time<-1939
EIA1969_time <-read.csv("PP_1969_5km_Wind.csv")  #change here
EIA1969_time$time<-1969
EIA1999_time <-read.csv("PP_1999_5km_Wind.csv")  #change here
EIA1999_time$time<-1999
EIA2019_time <-read.csv("PP_2019_5km_Wind.csv")  #change here
EIA2019_time$time<-2019

spatial_time<-rbind(EIA1939_time, EIA1969_time, EIA1999_time, EIA2019_time)

spatial_time$presence_plant<-ifelse(spatial_time$Join_Count>0,1,0)
spatial_time$presence_coplant<-ifelse(spatial_time$new_type_coal_oil>0,1,0)
spatial_time$presence_peaker<-ifelse(spatial_time$Peaker_plant>0,1,0)
spatial_time$presence_retired<-ifelse(spatial_time$Year_Of_Operation>0,1,0)

colnames(spatial_time)[2]<-"numberofplant"
colnames(spatial_time)[1]<-"OBJECTID"

#Create non-stratified emission data
EIA1939_notime <-read.csv("PP_1939_5km_Wind_Emission.csv") #change here
EIA1939_notime$time<-1939
EIA1969_notime <-read.csv("PP_1969_5km_Wind_Emission.csv") #change here
EIA1969_notime$time<-1969
EIA1999_notime <-read.csv("PP_1999_5km_Wind_Emission.csv") #change here
EIA1999_notime$time<-1999
EIA2019_notime <-read.csv("PP_2019_5km_Wind_Emission.csv") #change here
EIA2019_notime$time<-2019

spatial_notime<-rbind(EIA1939_notime, EIA1969_notime, EIA1999_notime, EIA2019_notime)
colnames(spatial_notime)[1]<-"OBJECTID"

spatial_notime_sumUNNOX<-spatial_notime %>%
  group_by(OBJECTID) %>%
  summarise(UNNOX = if(all(is.na(UNNOX))) NA_real_ else sum(UNNOX, na.rm = TRUE))
spatial_notime_sumUNSO2<-spatial_notime %>%
  group_by(OBJECTID) %>%
  summarise(UNSO2 = if(all(is.na(UNSO2))) NA_real_ else sum(UNSO2, na.rm = TRUE))
spatial_notime_sumUNPM25<-spatial_notime %>%
  group_by(OBJECTID) %>%
  summarise(UNPM25 = if(all(is.na(UNPM25))) NA_real_ else sum(UNPM25, na.rm = TRUE))
spatial_notime_sumJoin_Count<-spatial_notime %>%
  group_by(OBJECTID) %>%
  summarise(numberofplant = if(all(is.na(Join_Count))) NA_real_ else sum(Join_Count, na.rm = TRUE))

spatial_notime_1<-EIA1939_notime[,c(1,3:6)]
colnames(spatial_notime_1)[1]<-"OBJECTID"

spatial_notime_1<-merge(spatial_notime_1, spatial_notime_sumUNNOX,by="OBJECTID")
spatial_notime_1<-merge(spatial_notime_1, spatial_notime_sumUNSO2,by="OBJECTID")
spatial_notime_1<-merge(spatial_notime_1, spatial_notime_sumUNPM25,by="OBJECTID")
spatial_notime_1<-merge(spatial_notime_1, spatial_notime_sumJoin_Count,by="OBJECTID")

spatial_notime_1$presence_plant<-ifelse(spatial_notime_1$numberofplant>0,1,0)
spatial_notime_1<-spatial_notime_1[spatial_notime_1$presence_plant==1,]

######### Figure 2 and 4 Distribution of power plant characteristics and its emissions by HOLC grades ###########
#note: figure 4 was created in excel
spatial_time %>%
  dplyr::select(time, holc_grade, presence_plant, numberofplant, presence_coplant, presence_peaker, nameplate, presence_retired, Year_Of_Operation) %>%
  mutate(time = paste("time", time)) %>%
  tbl_strata(
    strata = time,
    ~.x %>%
      tbl_summary(
        type=list(c(nameplate, numberofplant, Year_Of_Operation)~"continuous", c(presence_plant, presence_coplant, presence_peaker, presence_retired)~"categorical"),
        by = holc_grade,
        statistic = all_continuous() ~ c("{mean} {sd}"),
        digits = all_continuous() ~ 4,
        missing="no"
      ) %>%
      add_p(pvalue_fun = ~style_pvalue(.x, digits = 3))%>%
      bold_labels()
  )

spatial_notime_1 %>%
  dplyr::select(holc_grade, UNNOX, UNSO2, UNPM25) %>%
      tbl_summary(
        type=list(c(UNNOX, UNSO2, UNPM25)~"continuous"),
        by = holc_grade,
        statistic = all_continuous() ~ c("{mean} ({sd})"),
        missing="no"
      ) %>%
      add_p(pvalue_fun = ~style_pvalue(.x, digits = 3))%>%
      bold_labels()

#Create Box Plot
bp.vals <- function(x, probs=c(0.1, 0.25, 0.75, .9)) {
  r <- quantile(x, probs=probs , na.rm=TRUE)
  r = c(r[1:2], exp(mean(log(x+1))), r[3:4])
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#Emission bar graphs for figure 2
ggplot(spatial_notime_1, aes(x=holc_grade, y=UNNOX, fill=holc_grade), outline=FALSE) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("green", "blue","yellow","red"))+
  coord_cartesian(ylim=c(0, 800))+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+
  labs(x = "HOLC Grade", y = "2019 NOx (tons)")+
  stat_summary(fun.data=bp.vals, geom="boxplot")

ggplot(spatial_notime_1, aes(x=holc_grade, y=UNSO2, fill=holc_grade), outline=FALSE) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("green", "blue","yellow","red"))+
  coord_cartesian(ylim=c(0, 25))+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+
  labs(x = "HOLC Grade", y = "2019 SO2 (tons)")+
  stat_summary(fun.data=bp.vals, geom="boxplot")

ggplot(spatial_notime_1, aes(x=holc_grade, y=UNPM25, fill=holc_grade), outline=FALSE) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("green", "blue","yellow","red"))+
  coord_cartesian(ylim=c(0, 60))+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+
  labs(x = "HOLC Grade", y = "2018 PM2.5 (tons)")+
  stat_summary(fun.data=bp.vals, geom="boxplot")

############################################################################################################
######### Table S1-S7 Model Estimates ###########

#################################getting the estimate for binary outcome###################################
#import the region, city level variable and HOLC polygon population variable
peter <- read.csv("holc_peter_region.csv")
spatial_time <- merge(spatial_time, peter, by.x='OBJECTID', by.y='peter_OBJECTID', all.x=TRUE)
spatial_time$regions<-as.factor(spatial_time$regions)

#Get the baseline presence of plants
spatial_time_1939 <- spatial_time[ which(spatial_time$time==1939), c("OBJECTID", "presence_plant")]
names(spatial_time_1939)[2]<- "baseline_plant"
spatial_time <- merge(spatial_time, spatial_time_1939, by="OBJECTID", all.x=TRUE)
spatial_time$baseline_plant <- as.factor(spatial_time$baseline_plant)
spatial_time$citypop_cat <- ifelse(spatial_time$citypop_cat=="Large",1,ifelse(spatial_time$citypop_cat=="Mediu",2,ifelse(spatial_time$citypop_cat=="Small",3,NA)))
spatial_time$citypop_cat<-as.factor(spatial_time$citypop_cat)

#subset into different analysis dataset
spatial_time_1969ab <- spatial_time[ which(spatial_time$time==1969 & (spatial_time$holc_grade== 'A' | spatial_time$holc_grade== 'B')),]
spatial_time_1969bc <- spatial_time[ which(spatial_time$time==1969 & (spatial_time$holc_grade== 'B' | spatial_time$holc_grade== 'C')),]
spatial_time_1969cd <- spatial_time[ which(spatial_time$time==1969 & (spatial_time$holc_grade== 'C' | spatial_time$holc_grade== 'D')),]

spatial_time_1999ab <- spatial_time[ which(spatial_time$time==1999 & (spatial_time$holc_grade== 'A' | spatial_time$holc_grade== 'B')),]
spatial_time_1999bc <- spatial_time[ which(spatial_time$time==1999 & (spatial_time$holc_grade== 'B' | spatial_time$holc_grade== 'C')),]
spatial_time_1999cd <- spatial_time[ which(spatial_time$time==1999 & (spatial_time$holc_grade== 'C' | spatial_time$holc_grade== 'D')),]

spatial_time_2019ab <- spatial_time[ which(spatial_time$time==2019 & (spatial_time$holc_grade== 'A' | spatial_time$holc_grade== 'B')),]
spatial_time_2019bc <- spatial_time[ which(spatial_time$time==2019 & (spatial_time$holc_grade== 'B' | spatial_time$holc_grade== 'C')),]
spatial_time_2019cd <- spatial_time[ which(spatial_time$time==2019 & (spatial_time$holc_grade== 'C' | spatial_time$holc_grade== 'D')),]

#model 1
#get the outcome variable
dvList <- names(spatial_time)[c(23:25)]#make sure this is the number of columns
rowlist1 <-c("intercept_pp",
             "grade_pp",
             "region2_pp",
             "region3_pp",
             "region4_pp",
             "baseline_plant_pp",
             "intercept_co",
             "grade_co",
             "region2_co",
             "region3_co",
             "region4_co",
             "baseline_plant_co",
             "intercept_peaker",
             "grade_peaker",
             "region2_peaker",
             "region3_peaker",
             "region4_peaker",
             "baseline_plant_peaker")
#set the variables to keep after regression
keep <- c("PR", "CL","Pr(>|z|)")

# run n regressions
my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab), vcov = sandwich)[, 1:4]})
results_1969ab_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1969ab_poisson) = rowlist1
results_1969ab_poisson$PR=round(exp(results_1969ab_poisson$Estimate), digits=2)
results_1969ab_poisson$LCL_exp=round(exp(results_1969ab_poisson$Estimate-1.96*results_1969ab_poisson$"Std. Erro"), digits=2)
results_1969ab_poisson$UCL_exp=round(exp(results_1969ab_poisson$Estimate+1.96*results_1969ab_poisson$"Std. Erro"), digits=2)
results_1969ab_poisson$CL <- paste(results_1969ab_poisson$LCL_exp, results_1969ab_poisson$UCL_exp, sep=",")
results_1969ab_poisson<-results_1969ab_poisson[c(-1,-7,-13),]
results_1969ab_poisson<-results_1969ab_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1969bc), vcov = sandwich)[, 1:4]})

results_1969bc_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1969bc_poisson) = rowlist1
results_1969bc_poisson$PR=round(exp(results_1969bc_poisson$Estimate), digits=2)
results_1969bc_poisson$LCL_exp=round(exp(results_1969bc_poisson$Estimate-1.96*results_1969bc_poisson$"Std. Erro"), digits=2)
results_1969bc_poisson$UCL_exp=round(exp(results_1969bc_poisson$Estimate+1.96*results_1969bc_poisson$"Std. Erro"), digits=2)
results_1969bc_poisson$CL <- paste(results_1969bc_poisson$LCL_exp, results_1969bc_poisson$UCL_exp, sep=",")
results_1969bc_poisson<-results_1969bc_poisson[c(-1,-7,-13),]
results_1969bc_poisson<-results_1969bc_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant  , list(i = as.name(x))), family="poisson", data=spatial_time_1969cd), vcov = sandwich)[, 1:4]})
results_1969cd_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1969cd_poisson) = rowlist1
results_1969cd_poisson$PR=round(exp(results_1969cd_poisson$Estimate), digits=2)
results_1969cd_poisson$LCL_exp=round(exp(results_1969cd_poisson$Estimate-1.96*results_1969cd_poisson$"Std. Erro"), digits=2)
results_1969cd_poisson$UCL_exp=round(exp(results_1969cd_poisson$Estimate+1.96*results_1969cd_poisson$"Std. Erro"), digits=2)
results_1969cd_poisson$CL <- paste(results_1969cd_poisson$LCL_exp, results_1969cd_poisson$UCL_exp, sep=",")
results_1969cd_poisson<-results_1969cd_poisson[c(-1,-7,-13),]
results_1969cd_poisson<-results_1969cd_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1999ab), vcov = sandwich)[, 1:4]})
results_1999ab_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1999ab_poisson) = rowlist1
results_1999ab_poisson$PR=round(exp(results_1999ab_poisson$Estimate), digits=2)
results_1999ab_poisson$LCL_exp=round(exp(results_1999ab_poisson$Estimate-1.96*results_1999ab_poisson$"Std. Erro"), digits=2)
results_1999ab_poisson$UCL_exp=round(exp(results_1999ab_poisson$Estimate+1.96*results_1999ab_poisson$"Std. Erro"), digits=2)
results_1999ab_poisson$CL <- paste(results_1999ab_poisson$LCL_exp, results_1999ab_poisson$UCL_exp, sep=",")
results_1999ab_poisson<-results_1999ab_poisson[c(-1,-7,-13),]
results_1999ab_poisson<-results_1999ab_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant  , list(i = as.name(x))), family="poisson", data=spatial_time_1999bc), vcov = sandwich)[, 1:4]})
results_1999bc_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1999bc_poisson) = rowlist1
results_1999bc_poisson$PR=round(exp(results_1999bc_poisson$Estimate), digits=2)
results_1999bc_poisson$LCL_exp=round(exp(results_1999bc_poisson$Estimate-1.96*results_1999bc_poisson$"Std. Erro"), digits=2)
results_1999bc_poisson$UCL_exp=round(exp(results_1999bc_poisson$Estimate+1.96*results_1999bc_poisson$"Std. Erro"), digits=2)
results_1999bc_poisson$CL <- paste(results_1999bc_poisson$LCL_exp, results_1999bc_poisson$UCL_exp, sep=",")
results_1999bc_poisson<-results_1999bc_poisson[c(-1,-7,-13),]
results_1999bc_poisson<-results_1999bc_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1999cd), vcov = sandwich)[, 1:4]})
results_1999cd_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_1999cd_poisson) = rowlist1
results_1999cd_poisson$PR=round(exp(results_1999cd_poisson$Estimate), digits=2)
results_1999cd_poisson$LCL_exp=round(exp(results_1999cd_poisson$Estimate-1.96*results_1999cd_poisson$"Std. Erro"), digits=2)
results_1999cd_poisson$UCL_exp=round(exp(results_1999cd_poisson$Estimate+1.96*results_1999cd_poisson$"Std. Erro"), digits=2)
results_1999cd_poisson$CL <- paste(results_1999cd_poisson$LCL_exp, results_1999cd_poisson$UCL_exp, sep=",")
results_1999cd_poisson<-results_1999cd_poisson[c(-1,-7,-13),]
results_1999cd_poisson<-results_1999cd_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019ab), vcov = sandwich)[, 1:4]})
results_2019ab_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_2019ab_poisson) = rowlist1
results_2019ab_poisson$PR=round(exp(results_2019ab_poisson$Estimate), digits=2)
results_2019ab_poisson$LCL_exp=round(exp(results_2019ab_poisson$Estimate-1.96*results_2019ab_poisson$"Std. Erro"), digits=2)
results_2019ab_poisson$UCL_exp=round(exp(results_2019ab_poisson$Estimate+1.96*results_2019ab_poisson$"Std. Erro"), digits=2)
results_2019ab_poisson$CL <- paste(results_2019ab_poisson$LCL_exp, results_2019ab_poisson$UCL_exp, sep=",")
results_2019ab_poisson<-results_2019ab_poisson[c(-1,-7,-13),]
results_2019ab_poisson<-results_2019ab_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019bc), vcov = sandwich)[, 1:4]})
results_2019bc_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_2019bc_poisson) = rowlist1
results_2019bc_poisson$PR=round(exp(results_2019bc_poisson$Estimate), digits=2)
results_2019bc_poisson$LCL_exp=round(exp(results_2019bc_poisson$Estimate-1.96*results_2019bc_poisson$"Std. Erro"), digits=2)
results_2019bc_poisson$UCL_exp=round(exp(results_2019bc_poisson$Estimate+1.96*results_2019bc_poisson$"Std. Erro"), digits=2)
results_2019bc_poisson$CL <- paste(results_2019bc_poisson$LCL_exp, results_2019bc_poisson$UCL_exp, sep=",")
results_2019bc_poisson<-results_2019bc_poisson[c(-1,-7,-13),]
results_2019bc_poisson<-results_2019bc_poisson[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019cd), vcov = sandwich)[, 1:4]})
results_2019cd_poisson <-do.call(rbind.data.frame, my_lms)
rownames(results_2019cd_poisson) = rowlist1
results_2019cd_poisson$PR=round(exp(results_2019cd_poisson$Estimate), digits=2)
results_2019cd_poisson$LCL_exp=round(exp(results_2019cd_poisson$Estimate-1.96*results_2019cd_poisson$"Std. Erro"), digits=2)
results_2019cd_poisson$UCL_exp=round(exp(results_2019cd_poisson$Estimate+1.96*results_2019cd_poisson$"Std. Erro"), digits=2)
results_2019cd_poisson$CL <- paste(results_2019cd_poisson$LCL_exp, results_2019cd_poisson$UCL_exp, sep=",")
results_2019cd_poisson<-results_2019cd_poisson[c(-1,-7,-13),]
results_2019cd_poisson<-results_2019cd_poisson[keep]


result_poisson <-cbind(results_1969ab_poisson, results_1999ab_poisson, results_2019ab_poisson,
                       results_1969bc_poisson,  results_1999bc_poisson,  results_2019bc_poisson,
                       results_1969cd_poisson, results_1999cd_poisson,  results_2019cd_poisson)
result_poisson[,2]<-paste0("[", result_poisson[,2], "]")
result_poisson[,5]<-paste0("[", result_poisson[,5], "]")
result_poisson[,8]<-paste0("[", result_poisson[,8], "]")
result_poisson[,11]<-paste0("[", result_poisson[,11], "]")
result_poisson[,14]<-paste0("[", result_poisson[,14], "]")
result_poisson[,17]<-paste0("[", result_poisson[,17], "]")
result_poisson[,20]<-paste0("[", result_poisson[,20], "]")
result_poisson[,23]<-paste0("[", result_poisson[,23], "]")
result_poisson[,26]<-paste0("[", result_poisson[,26], "]")

colnames(result_poisson) <- c("1969ab ", "1969ab CL","1969ab p-value", "1999ab ", "1999ab CL","1999ab p-value", "2019ab ", "2019ab CL","2019ab p-value",
                              "1969bc ", "1969bc CL","1969bc p-value", "1999bc ", "1999bc CL","1999bc p-value", "2019bc ", "2019bc CL","2019bc p-value",
                              "1969cd ", "1969cd CL","1969cd p-value", "1999cd ", "1999cd CL","1999cd p-value", "2019cd ", "2019cd CL","2019cd p-value")

write.csv(result_poisson,"result_all_poisson_model1.csv", row.names = TRUE)

#check for observations used for analysis
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd)))})
# 
# 
# #check for dispersion
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd))})
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd))})
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd))})
# 
# 
# #check for model fitness
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1969ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1969bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1969cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1999ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1999bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_1999cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant , list(i = as.name(x))), family="poisson", data=spatial_time_2019cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# 
# #get AIC
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd))})
# 
# AIC1969AB <-do.call(rbind.data.frame, AIC1)
# AIC1969BC <-do.call(rbind.data.frame, AIC2)
# AIC1969CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd))})
# 
# AIC1999AB <-do.call(rbind.data.frame, AIC1)
# AIC1999BC <-do.call(rbind.data.frame, AIC2)
# AIC1999CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd))})
# 
# AIC2019AB <-do.call(rbind.data.frame, AIC1)
# AIC2019BC <-do.call(rbind.data.frame, AIC2)
# AIC2019CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC_all <- cbind(AIC1969AB, AIC1999AB, AIC2019AB, AIC1969BC, AIC1999BC,AIC2019BC, AIC1969CD,AIC1999CD ,AIC2019CD)
# 
# write.csv (AIC_all, "EIA2019/AIC.csv")

#model 2
dvList <- names(spatial_time)[c(23:25)]#make sure this is the number of columns
rowlist2 <-c("intercept_pp",	
             "grade_pp",
             "region2_pp",
             "region3_pp",
             "region4_pp",
             'baseline_plant_pp',
             'median_pp',
             'small_pp',
             'nowhite_pp',
             "intercept_co",	
             "grade_co",
             "region2_co",
             "region3_co",
             'region4_co',
             'baseline_plant_co',
             'median_co',
             'small_co',	
             'nowhite_co',
             "intercept_peaker",
             "grade_peaker",
             "region2_peaker",
             "region3_peaker",
             'region4_peaker',
             'baseline_plant_peaker',
             'median_peaker',
             'small_peaker',
             'nowhite_peaker'
)
keep <- c("PR", "CL","Pr(>|z|)")

# run n regressions
my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab), vcov = sandwich)[, 1:4]})
results_1969ab_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969ab_poisson2) = rowlist2
results_1969ab_poisson2$PR=round(exp(results_1969ab_poisson2$Estimate), digits=2)
results_1969ab_poisson2$LCL_exp=round(exp(results_1969ab_poisson2$Estimate-1.96*results_1969ab_poisson2$"Std. Erro"), digits=2)
results_1969ab_poisson2$UCL_exp=round(exp(results_1969ab_poisson2$Estimate+1.96*results_1969ab_poisson2$"Std. Erro"), digits=2)
results_1969ab_poisson2$CL <- paste(results_1969ab_poisson2$LCL_exp, results_1969ab_poisson2$UCL_exp, sep=",")
results_1969ab_poisson2<-results_1969ab_poisson2[c(-1,-10,-19),]
results_1969ab_poisson2<-results_1969ab_poisson2[keep]


my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc), vcov = sandwich)[, 1:4]})
results_1969bc_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969bc_poisson2) = rowlist2
results_1969bc_poisson2$PR=round(exp(results_1969bc_poisson2$Estimate), digits=2)
results_1969bc_poisson2$LCL_exp=round(exp(results_1969bc_poisson2$Estimate-1.96*results_1969bc_poisson2$"Std. Erro"), digits=2)
results_1969bc_poisson2$UCL_exp=round(exp(results_1969bc_poisson2$Estimate+1.96*results_1969bc_poisson2$"Std. Erro"), digits=2)
results_1969bc_poisson2$CL <- paste(results_1969bc_poisson2$LCL_exp, results_1969bc_poisson2$UCL_exp, sep=",")
results_1969bc_poisson2<-results_1969bc_poisson2[c(-1,-10,-19),]
results_1969bc_poisson2<-results_1969bc_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd), vcov = sandwich)[, 1:4]})
results_1969cd_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969cd_poisson2) = rowlist2
results_1969cd_poisson2$PR=round(exp(results_1969cd_poisson2$Estimate), digits=2)
results_1969cd_poisson2$LCL_exp=round(exp(results_1969cd_poisson2$Estimate-1.96*results_1969cd_poisson2$"Std. Erro"), digits=2)
results_1969cd_poisson2$UCL_exp=round(exp(results_1969cd_poisson2$Estimate+1.96*results_1969cd_poisson2$"Std. Erro"), digits=2)
results_1969cd_poisson2$CL <- paste(results_1969cd_poisson2$LCL_exp, results_1969cd_poisson2$UCL_exp, sep=",")
results_1969cd_poisson2<-results_1969cd_poisson2[c(-1,-10,-19),]
results_1969cd_poisson2<-results_1969cd_poisson2[keep]


# run n regressions - spatial_time_1999ab
my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab), vcov = sandwich)[, 1:4]})
results_1999ab_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1999ab_poisson2) = rowlist2
results_1999ab_poisson2$PR=round(exp(results_1999ab_poisson2$Estimate), digits=2)
results_1999ab_poisson2$LCL_exp=round(exp(results_1999ab_poisson2$Estimate-1.96*results_1999ab_poisson2$"Std. Erro"), digits=2)
results_1999ab_poisson2$UCL_exp=round(exp(results_1999ab_poisson2$Estimate+1.96*results_1999ab_poisson2$"Std. Erro"), digits=2)
results_1999ab_poisson2$CL <- paste(results_1999ab_poisson2$LCL_exp, results_1999ab_poisson2$UCL_exp, sep=",")
results_1999ab_poisson2<-results_1999ab_poisson2[c(-1,-10,-19),]
results_1999ab_poisson2<-results_1999ab_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc), vcov = sandwich)[, 1:4]})
results_1999bc_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1999bc_poisson2) = rowlist2
results_1999bc_poisson2$PR=round(exp(results_1999bc_poisson2$Estimate), digits=2)
results_1999bc_poisson2$LCL_exp=round(exp(results_1999bc_poisson2$Estimate-1.96*results_1999bc_poisson2$"Std. Erro"), digits=2)
results_1999bc_poisson2$UCL_exp=round(exp(results_1999bc_poisson2$Estimate+1.96*results_1999bc_poisson2$"Std. Erro"), digits=2)
results_1999bc_poisson2$CL <- paste(results_1999bc_poisson2$LCL_exp, results_1999bc_poisson2$UCL_exp, sep=",")
results_1999bc_poisson2<-results_1999bc_poisson2[c(-1,-10,-19),]
results_1999bc_poisson2<-results_1999bc_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd), vcov = sandwich)[, 1:4]})
results_1999cd_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1999cd_poisson2) = rowlist2
results_1999cd_poisson2$PR=round(exp(results_1999cd_poisson2$Estimate), digits=2)
results_1999cd_poisson2$LCL_exp=round(exp(results_1999cd_poisson2$Estimate-1.96*results_1999cd_poisson2$"Std. Erro"), digits=2)
results_1999cd_poisson2$UCL_exp=round(exp(results_1999cd_poisson2$Estimate+1.96*results_1999cd_poisson2$"Std. Erro"), digits=2)
results_1999cd_poisson2$CL <- paste(results_1999cd_poisson2$LCL_exp, results_1999cd_poisson2$UCL_exp, sep=",")
results_1999cd_poisson2<-results_1999cd_poisson2[c(-1,-10,-19),]
results_1999cd_poisson2<-results_1999cd_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab), vcov = sandwich)[, 1:4]})
results_2019ab_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_2019ab_poisson2) = rowlist2
results_2019ab_poisson2$PR=round(exp(results_2019ab_poisson2$Estimate), digits=2)
results_2019ab_poisson2$LCL_exp=round(exp(results_2019ab_poisson2$Estimate-1.96*results_2019ab_poisson2$"Std. Erro"), digits=2)
results_2019ab_poisson2$UCL_exp=round(exp(results_2019ab_poisson2$Estimate+1.96*results_2019ab_poisson2$"Std. Erro"), digits=2)
results_2019ab_poisson2$CL <- paste(results_2019ab_poisson2$LCL_exp, results_2019ab_poisson2$UCL_exp, sep=",")
results_2019ab_poisson2<-results_2019ab_poisson2[c(-1,-10,-19),]
results_2019ab_poisson2<-results_2019ab_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc), vcov = sandwich)[, 1:4]})
results_2019bc_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_2019bc_poisson2) = rowlist2
results_2019bc_poisson2$PR=round(exp(results_2019bc_poisson2$Estimate), digits=2)
results_2019bc_poisson2$LCL_exp=round(exp(results_2019bc_poisson2$Estimate-1.96*results_2019bc_poisson2$"Std. Erro"), digits=2)
results_2019bc_poisson2$UCL_exp=round(exp(results_2019bc_poisson2$Estimate+1.96*results_2019bc_poisson2$"Std. Erro"), digits=2)
results_2019bc_poisson2$CL <- paste(results_2019bc_poisson2$LCL_exp, results_2019bc_poisson2$UCL_exp, sep=",")
results_2019bc_poisson2<-results_2019bc_poisson2[c(-1,-10,-19),]
results_2019bc_poisson2<-results_2019bc_poisson2[keep]

my_lms <- lapply(dvList, function(x) {
  coeftest(glm(substitute(i ~ holc_grade  +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd), vcov = sandwich)[, 1:4]})
results_2019cd_poisson2 <-do.call(rbind.data.frame, my_lms)
rownames(results_2019cd_poisson2) = rowlist2
results_2019cd_poisson2$PR=round(exp(results_2019cd_poisson2$Estimate), digits=2)
results_2019cd_poisson2$LCL_exp=round(exp(results_2019cd_poisson2$Estimate-1.96*results_2019cd_poisson2$"Std. Erro"), digits=2)
results_2019cd_poisson2$UCL_exp=round(exp(results_2019cd_poisson2$Estimate+1.96*results_2019cd_poisson2$"Std. Erro"), digits=2)
results_2019cd_poisson2$CL <- paste(results_2019cd_poisson2$LCL_exp, results_2019cd_poisson2$UCL_exp, sep=",")
results_2019cd_poisson2<-results_2019cd_poisson2[c(-1,-10,-19),]
results_2019cd_poisson2<-results_2019cd_poisson2[keep]


result_poisson2 <-cbind(results_1969ab_poisson2, results_1999ab_poisson2, results_2019ab_poisson2,
                        results_1969bc_poisson2,  results_1999bc_poisson2,  results_2019bc_poisson2,
                        results_1969cd_poisson2, results_1999cd_poisson2,  results_2019cd_poisson2)

colnames(result_poisson2) <- c("1969ab ", "1969ab CL","1969ab p-value", "1999ab ", "1999ab CL","1999ab p-value", "2019ab ", "2019ab CL","2019ab p-value",
                               "1969bc ", "1969bc CL","1969bc p-value", "1999bc ", "1999bc CL","1999bc p-value", "2019bc ", "2019bc CL","2019bc p-value",
                               "1969cd ", "1969cd CL","1969cd p-value", "1999cd ", "1999cd CL","1999cd p-value", "2019cd ", "2019cd CL","2019cd p-value")
result_poisson2[,2]<-paste0("[", result_poisson2[,2], "]")
result_poisson2[,5]<-paste0("[", result_poisson2[,5], "]")
result_poisson2[,8]<-paste0("[", result_poisson2[,8], "]")
result_poisson2[,11]<-paste0("[", result_poisson2[,11], "]")
result_poisson2[,14]<-paste0("[", result_poisson2[,14], "]")
result_poisson2[,17]<-paste0("[", result_poisson2[,17], "]")
result_poisson2[,20]<-paste0("[", result_poisson2[,20], "]")
result_poisson2[,23]<-paste0("[", result_poisson2[,23], "]")
result_poisson2[,26]<-paste0("[", result_poisson2[,26], "]")
write.csv(result_poisson2,"result_all_poisson_model2.csv", row.names = TRUE)

# #Check number of observations used and diagnostics of models
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc)))})
# obs <- lapply(dvList, function(x) {
#   nobs((glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd)))})
# 
# 
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd))})
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd))})
# disperse1 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab))})
# disperse2 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc))})
# disperse3 <- lapply(dvList, function(x) {
#   dispersiontest(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd))})
# 
# 
# 
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1969ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1969bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1969cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1999ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1999bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_1999cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD1 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_2019ab), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD2 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_2019bc), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# GOOD3 <- lapply(dvList, function(x) {
#   with(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite , list(i = as.name(x))), family="poisson", data=spatial_time_2019cd), cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))})
# 
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1969cd))})
# 
# AIC1969AB <-do.call(rbind.data.frame, AIC1)
# AIC1969BC <-do.call(rbind.data.frame, AIC2)
# AIC1969CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_1999cd))})
# 
# AIC1999AB <-do.call(rbind.data.frame, AIC1)
# AIC1999BC <-do.call(rbind.data.frame, AIC2)
# AIC1999CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC1 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019ab))})
# AIC2 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019bc))})
# AIC3 <- lapply(dvList, function(x) {
#   AIC(glm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), family="poisson", data=spatial_time_2019cd))})
# 
# AIC2019AB <-do.call(rbind.data.frame, AIC1)
# AIC2019BC <-do.call(rbind.data.frame, AIC2)
# AIC2019CD <-do.call(rbind.data.frame, AIC3)
# 
# AIC_all2 <- cbind(AIC1969AB, AIC1999AB, AIC2019AB, AIC1969BC, AIC1999BC,AIC2019BC, AIC1969CD,AIC1999CD ,AIC2019CD)
# write.csv (AIC_all2, "EIA2019/AIC.csv")


############################################################################################################
##########number of plants   negative binomial###############################
#model1
keep<-c("IRR", "CL","Pr(>|z|)")
nb_1969ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969ab))$coefficients[, 1:4])
nb_1969bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969bc))$coefficients[, 1:4])
nb_1969cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969cd))$coefficients[, 1:4])

nb_1969ab$IRR=round(exp(nb_1969ab$Estimate), digits=2)
nb_1969ab$LCL_exp=round(exp(nb_1969ab$Estimate-1.96*nb_1969ab$"Std. Erro"), digits=2)
nb_1969ab$UCL_exp=round(exp(nb_1969ab$Estimate+1.96*nb_1969ab$"Std. Erro"), digits=2)
nb_1969ab$CL <- paste(nb_1969ab$LCL_exp, nb_1969ab$UCL_exp, sep=",")
nb_1969ab<-nb_1969ab[c(-1),]
nb_1969ab<-nb_1969ab[keep]

nb_1969bc$IRR=round(exp(nb_1969bc$Estimate), digits=2)
nb_1969bc$LCL_exp=round(exp(nb_1969bc$Estimate-1.96*nb_1969bc$"Std. Erro"), digits=2)
nb_1969bc$UCL_exp=round(exp(nb_1969bc$Estimate+1.96*nb_1969bc$"Std. Erro"), digits=2)
nb_1969bc$CL <- paste(nb_1969bc$LCL_exp, nb_1969bc$UCL_exp, sep=",")
nb_1969bc<-nb_1969bc[c(-1),]
nb_1969bc<-nb_1969bc[keep]

nb_1969cd$IRR=round(exp(nb_1969cd$Estimate), digits=2)
nb_1969cd$LCL_exp=round(exp(nb_1969cd$Estimate-1.96*nb_1969cd$"Std. Erro"), digits=2)
nb_1969cd$UCL_exp=round(exp(nb_1969cd$Estimate+1.96*nb_1969cd$"Std. Erro"), digits=2)
nb_1969cd$CL <- paste(nb_1969cd$LCL_exp, nb_1969cd$UCL_exp, sep=",")
nb_1969cd<-nb_1969cd[c(-1),]
nb_1969cd<-nb_1969cd[keep]

nb_1999ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999ab))$coefficients[, 1:4])
nb_1999bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999bc))$coefficients[, 1:4])
nb_1999cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999cd))$coefficients[, 1:4])

nb_1999ab$IRR=round(exp(nb_1999ab$Estimate), digits=2)
nb_1999ab$LCL_exp=round(exp(nb_1999ab$Estimate-1.96*nb_1999ab$"Std. Erro"), digits=2)
nb_1999ab$UCL_exp=round(exp(nb_1999ab$Estimate+1.96*nb_1999ab$"Std. Erro"), digits=2)
nb_1999ab$CL <- paste(nb_1999ab$LCL_exp, nb_1999ab$UCL_exp, sep=",")
nb_1999ab<-nb_1999ab[c(-1),]
nb_1999ab<-nb_1999ab[keep]

nb_1999bc$IRR=round(exp(nb_1999bc$Estimate), digits=2)
nb_1999bc$LCL_exp=round(exp(nb_1999bc$Estimate-1.96*nb_1999bc$"Std. Erro"), digits=2)
nb_1999bc$UCL_exp=round(exp(nb_1999bc$Estimate+1.96*nb_1999bc$"Std. Erro"), digits=2)
nb_1999bc$CL <- paste(nb_1999bc$LCL_exp, nb_1999bc$UCL_exp, sep=",")
nb_1999bc<-nb_1999bc[c(-1),]
nb_1999bc<-nb_1999bc[keep]

nb_1999cd$IRR=round(exp(nb_1999cd$Estimate), digits=2)
nb_1999cd$LCL_exp=round(exp(nb_1999cd$Estimate-1.96*nb_1999cd$"Std. Erro"), digits=2)
nb_1999cd$UCL_exp=round(exp(nb_1999cd$Estimate+1.96*nb_1999cd$"Std. Erro"), digits=2)
nb_1999cd$CL <- paste(nb_1999cd$LCL_exp, nb_1999cd$UCL_exp, sep=",")
nb_1999cd<-nb_1999cd[c(-1),]
nb_1999cd<-nb_1999cd[keep]

nb_2019ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019ab))$coefficients[, 1:4])
nb_2019bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019bc))$coefficients[, 1:4])
nb_2019cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019cd))$coefficients[, 1:4])

nb_2019ab$IRR=round(exp(nb_2019ab$Estimate), digits=2)
nb_2019ab$LCL_exp=round(exp(nb_2019ab$Estimate-1.96*nb_2019ab$"Std. Erro"), digits=2)
nb_2019ab$UCL_exp=round(exp(nb_2019ab$Estimate+1.96*nb_2019ab$"Std. Erro"), digits=2)
nb_2019ab$CL <- paste(nb_2019ab$LCL_exp, nb_2019ab$UCL_exp, sep=",")
nb_2019ab<-nb_2019ab[c(-1),]
nb_2019ab<-nb_2019ab[keep]

nb_2019bc$IRR=round(exp(nb_2019bc$Estimate), digits=2)
nb_2019bc$LCL_exp=round(exp(nb_2019bc$Estimate-1.96*nb_2019bc$"Std. Erro"), digits=2)
nb_2019bc$UCL_exp=round(exp(nb_2019bc$Estimate+1.96*nb_2019bc$"Std. Erro"), digits=2)
nb_2019bc$CL <- paste(nb_2019bc$LCL_exp, nb_2019bc$UCL_exp, sep=",")
nb_2019bc<-nb_2019bc[c(-1),]
nb_2019bc<-nb_2019bc[keep]

nb_2019cd$IRR=round(exp(nb_2019cd$Estimate), digits=2)
nb_2019cd$LCL_exp=round(exp(nb_2019cd$Estimate-1.96*nb_2019cd$"Std. Erro"), digits=2)
nb_2019cd$UCL_exp=round(exp(nb_2019cd$Estimate+1.96*nb_2019cd$"Std. Erro"), digits=2)
nb_2019cd$CL <- paste(nb_2019cd$LCL_exp, nb_2019cd$UCL_exp, sep=",")
nb_2019cd<-nb_2019cd[c(-1),]
nb_2019cd<-nb_2019cd[keep]

result_nb <-cbind(nb_1969ab, nb_1999ab, nb_2019ab,
                  nb_1969bc,  nb_1999bc,  nb_2019bc,
                  nb_1969cd, nb_1999cd,  nb_2019cd)
result_nb[,2]<-paste0("[", result_nb[,2], "]")
result_nb[,5]<-paste0("[", result_nb[,5], "]")
result_nb[,8]<-paste0("[", result_nb[,8], "]")
result_nb[,11]<-paste0("[", result_nb[,11], "]")
result_nb[,14]<-paste0("[", result_nb[,14], "]")
result_nb[,17]<-paste0("[", result_nb[,17], "]")
result_nb[,20]<-paste0("[", result_nb[,20], "]")
result_nb[,23]<-paste0("[", result_nb[,23], "]")
result_nb[,26]<-paste0("[", result_nb[,26], "]")
colnames(result_nb) <- c("1969ab ", "1969ab CL","1969ab p-value", "1999ab ", "1999ab CL","1999ab p-value", "2019ab ", "2019ab CL","2019ab p-value",
                         "1969bc ", "1969bc CL","1969bc p-value", "1999bc ", "1999bc CL","1999bc p-value", "2019bc ", "2019bc CL","2019bc p-value",
                         "1969cd ", "1969cd CL","1969cd p-value", "1999cd ", "1999cd CL","1999cd p-value", "2019cd ", "2019cd CL","2019cd p-value")

write.csv(result_nb,"result_nb_model1.csv", row.names = TRUE)

# #Check number of observations used and diagnostics of models
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969cd)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999cd)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019cd)))
# 
# #Dispersion in Poisson Model
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1969ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1969bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1969cd))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1999ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1999bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_1999cd))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_2019ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_2019bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant, family="poisson", data=spatial_time_2019cd))
# 
# AIC1969AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969ab))
# AIC1969BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969bc))
# AIC1969CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1969cd))
# 
# AIC1999AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999ab))
# AIC1999BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999bc))
# AIC1999CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_1999cd))
# 
# AIC2019AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019ab))
# AIC2019BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019bc))
# AIC2019CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant, data=spatial_time_2019cd))
# 
# AIC_all <- cbind(AIC1969AB, AIC1999AB, AIC2019AB, AIC1969BC, AIC1999BC,AIC2019BC, AIC1969CD,AIC1999CD ,AIC2019CD)
# write.csv (AIC_all, "EIA2019/AIC.csv")
###############################################################################################################################
#model2
nb2_1969ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite+citypop_cat+pct_nowhite, data=spatial_time_1969ab))$coefficients[, 1:4])
nb2_1969bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969bc))$coefficients[, 1:4])
nb2_1969cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969cd))$coefficients[, 1:4])

nb2_1969ab$IRR=round(exp(nb2_1969ab$Estimate), digits=2)
nb2_1969ab$LCL_exp=round(exp(nb2_1969ab$Estimate-1.96*nb2_1969ab$"Std. Erro"), digits=2)
nb2_1969ab$UCL_exp=round(exp(nb2_1969ab$Estimate+1.96*nb2_1969ab$"Std. Erro"), digits=2)
nb2_1969ab$CL <- paste(nb2_1969ab$LCL_exp, nb2_1969ab$UCL_exp, sep=",")
nb2_1969ab<-nb2_1969ab[c(-1),]
nb2_1969ab<-nb2_1969ab[keep]

nb2_1969bc$IRR=round(exp(nb2_1969bc$Estimate), digits=2)
nb2_1969bc$LCL_exp=round(exp(nb2_1969bc$Estimate-1.96*nb2_1969bc$"Std. Erro"), digits=2)
nb2_1969bc$UCL_exp=round(exp(nb2_1969bc$Estimate+1.96*nb2_1969bc$"Std. Erro"), digits=2)
nb2_1969bc$CL <- paste(nb2_1969bc$LCL_exp, nb2_1969bc$UCL_exp, sep=",")
nb2_1969bc<-nb2_1969bc[c(-1),]
nb2_1969bc<-nb2_1969bc[keep]

nb2_1969cd$IRR=round(exp(nb2_1969cd$Estimate), digits=2)
nb2_1969cd$LCL_exp=round(exp(nb2_1969cd$Estimate-1.96*nb2_1969cd$"Std. Erro"), digits=2)
nb2_1969cd$UCL_exp=round(exp(nb2_1969cd$Estimate+1.96*nb2_1969cd$"Std. Erro"), digits=2)
nb2_1969cd$CL <- paste(nb2_1969cd$LCL_exp, nb2_1969cd$UCL_exp, sep=",")
nb2_1969cd<-nb2_1969cd[c(-1),]
nb2_1969cd<-nb2_1969cd[keep]

nb2_1999ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999ab))$coefficients[, 1:4])
nb2_1999bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999bc))$coefficients[, 1:4])
nb2_1999cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999cd))$coefficients[, 1:4])

nb2_1999ab$IRR=round(exp(nb2_1999ab$Estimate), digits=2)
nb2_1999ab$LCL_exp=round(exp(nb2_1999ab$Estimate-1.96*nb2_1999ab$"Std. Erro"), digits=2)
nb2_1999ab$UCL_exp=round(exp(nb2_1999ab$Estimate+1.96*nb2_1999ab$"Std. Erro"), digits=2)
nb2_1999ab$CL <- paste(nb2_1999ab$LCL_exp, nb2_1999ab$UCL_exp, sep=",")
nb2_1999ab<-nb2_1999ab[c(-1),]
nb2_1999ab<-nb2_1999ab[keep]

nb2_1999bc$IRR=round(exp(nb2_1999bc$Estimate), digits=2)
nb2_1999bc$LCL_exp=round(exp(nb2_1999bc$Estimate-1.96*nb2_1999bc$"Std. Erro"), digits=2)
nb2_1999bc$UCL_exp=round(exp(nb2_1999bc$Estimate+1.96*nb2_1999bc$"Std. Erro"), digits=2)
nb2_1999bc$CL <- paste(nb2_1999bc$LCL_exp, nb2_1999bc$UCL_exp, sep=",")
nb2_1999bc<-nb2_1999bc[c(-1),]
nb2_1999bc<-nb2_1999bc[keep]

nb2_1999cd$IRR=round(exp(nb2_1999cd$Estimate), digits=2)
nb2_1999cd$LCL_exp=round(exp(nb2_1999cd$Estimate-1.96*nb2_1999cd$"Std. Erro"), digits=2)
nb2_1999cd$UCL_exp=round(exp(nb2_1999cd$Estimate+1.96*nb2_1999cd$"Std. Erro"), digits=2)
nb2_1999cd$CL <- paste(nb2_1999cd$LCL_exp, nb2_1999cd$UCL_exp, sep=",")
nb2_1999cd<-nb2_1999cd[c(-1),]
nb2_1999cd<-nb2_1999cd[keep]

nb2_2019ab <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019ab))$coefficients[, 1:4])
nb2_2019bc <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019bc))$coefficients[, 1:4])
nb2_2019cd <-as.data.frame(summary(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019cd))$coefficients[, 1:4])

nb2_2019ab$IRR=round(exp(nb2_2019ab$Estimate), digits=2)
nb2_2019ab$LCL_exp=round(exp(nb2_2019ab$Estimate-1.96*nb2_2019ab$"Std. Erro"), digits=2)
nb2_2019ab$UCL_exp=round(exp(nb2_2019ab$Estimate+1.96*nb2_2019ab$"Std. Erro"), digits=2)
nb2_2019ab$CL <- paste(nb2_2019ab$LCL_exp, nb2_2019ab$UCL_exp, sep=",")
nb2_2019ab<-nb2_2019ab[c(-1),]
nb2_2019ab<-nb2_2019ab[keep]

nb2_2019bc$IRR=round(exp(nb2_2019bc$Estimate), digits=2)
nb2_2019bc$LCL_exp=round(exp(nb2_2019bc$Estimate-1.96*nb2_2019bc$"Std. Erro"), digits=2)
nb2_2019bc$UCL_exp=round(exp(nb2_2019bc$Estimate+1.96*nb2_2019bc$"Std. Erro"), digits=2)
nb2_2019bc$CL <- paste(nb2_2019bc$LCL_exp, nb2_2019bc$UCL_exp, sep=",")
nb2_2019bc<-nb2_2019bc[c(-1),]
nb2_2019bc<-nb2_2019bc[keep]

nb2_2019cd$IRR=round(exp(nb2_2019cd$Estimate), digits=2)
nb2_2019cd$LCL_exp=round(exp(nb2_2019cd$Estimate-1.96*nb2_2019cd$"Std. Erro"), digits=2)
nb2_2019cd$UCL_exp=round(exp(nb2_2019cd$Estimate+1.96*nb2_2019cd$"Std. Erro"), digits=2)
nb2_2019cd$CL <- paste(nb2_2019cd$LCL_exp, nb2_2019cd$UCL_exp, sep=",")
nb2_2019cd<-nb2_2019cd[c(-1),]
nb2_2019cd<-nb2_2019cd[keep]

result_nb2 <-cbind(nb2_1969ab, nb2_1999ab, nb2_2019ab,
                   nb2_1969bc,  nb2_1999bc,  nb2_2019bc,
                   nb2_1969cd, nb2_1999cd,  nb2_2019cd)

result_nb2[,2]<-paste0("[", result_nb2[,2], "]")
result_nb2[,5]<-paste0("[", result_nb2[,5], "]")
result_nb2[,8]<-paste0("[", result_nb2[,8], "]")
result_nb2[,11]<-paste0("[", result_nb2[,11], "]")
result_nb2[,14]<-paste0("[", result_nb2[,14], "]")
result_nb2[,17]<-paste0("[", result_nb2[,17], "]")
result_nb2[,20]<-paste0("[", result_nb2[,20], "]")
result_nb2[,23]<-paste0("[", result_nb2[,23], "]")
result_nb2[,26]<-paste0("[", result_nb2[,26], "]")
colnames(result_nb2) <- c("1969ab ", "1969ab CL","1969ab p-value", "1999ab ", "1999ab CL","1999ab p-value", "2019ab ", "2019ab CL","2019ab p-value",
                          "1969bc ", "1969bc CL","1969bc p-value", "1999bc ", "1999bc CL","1999bc p-value", "2019bc ", "2019bc CL","2019bc p-value",
                          "1969cd ", "1969cd CL","1969cd p-value", "1999cd ", "1999cd CL","1999cd p-value", "2019cd ", "2019cd CL","2019cd p-value")

write.csv(result_nb2,"result_nb2_model2.csv", row.names = TRUE)

# #Check number of observations used and diagnostics of models
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969cd)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999cd)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019ab)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019bc)))
# nobs((glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019cd)))
# 
# #Dispersion in Poisson Model
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1969ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1969bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1969cd))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1999ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1999bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_1999cd))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_2019ab))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_2019bc))
# dispersiontest(glm(numberofplant ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, family="poisson", data=spatial_time_2019cd))
# 
# AIC1969AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969ab))
# AIC1969BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969bc))
# AIC1969CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1969cd))
# 
# AIC1999AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999ab))
# AIC1999BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999bc))
# AIC1999CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_1999cd))
# 
# AIC2019AB<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019ab))
# AIC2019BC<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019bc))
# AIC2019CD<-AIC(glm.nb(numberofplant ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, data=spatial_time_2019cd))
# 
# AIC_all2 <- cbind(AIC1969AB, AIC1999AB, AIC2019AB, AIC1969BC, AIC1999BC,AIC2019BC, AIC1969CD,AIC1999CD ,AIC2019CD)
# write.csv (AIC_all2, "EIA2019/AIC.csv")
############################################################################################################
############################################################################################################
##########Linear Model for Emission Data###############################
spatial_notime_1 <- merge(spatial_notime_1, peter, by.x='OBJECTID', by.y='peter_OBJECTID', all.x=TRUE)
spatial_notime_1$regions<-as.factor(spatial_notime_1$regions)

#process other variables
spatial_time_1939 <- spatial_time[ which(spatial_time$time==1939), c("OBJECTID", "presence_plant")]
names(spatial_time_1939)[2]<- "baseline_plant"
spatial_notime_1 <- merge(spatial_notime_1, spatial_time_1939, by="OBJECTID", all.x=TRUE)
spatial_notime_1$baseline_plant <- factor(spatial_notime_1$baseline_plant)

spatial_notime_1$UNNOX   <- log(spatial_notime_1$UNNOX+1)
spatial_notime_1$UNSO2   <- log(spatial_notime_1$UNSO2+1)
spatial_notime_1$UNPM25  <- log(spatial_notime_1$UNPM25+1)

spatial_notime_1$citypop_cat <- ifelse(spatial_notime_1$citypop_cat=="Large",1,ifelse(spatial_notime_1$citypop_cat=="Mediu",2,ifelse(spatial_notime_1$citypop_cat=="Small",3,NA)))
spatial_notime_1$citypop_cat<-as.factor(spatial_notime_1$citypop_cat)


#subset into different analysis dataset
spatial_notime_ab <- spatial_notime_1[ which( (spatial_notime_1$holc_grade== 'A' | spatial_notime_1$holc_grade== 'B')),]
spatial_notime_bc <- spatial_notime_1[ which( (spatial_notime_1$holc_grade== 'B' | spatial_notime_1$holc_grade== 'C')),]
spatial_notime_cd <- spatial_notime_1[ which((spatial_notime_1$holc_grade== 'C' | spatial_notime_1$holc_grade== 'D')),]

dvList3 <- names(spatial_notime)[c(10,12,16)]#make sure this is the number of columns
rowlist1 <-c("intercept_nox",
             "grade_nox",
             "region2_nox",
             "region3_nox",
             "region4_nox",
             "baseline_plant_nox",
             "intercept_so2",
             "grade_so2",
             "region2_so2",
             "region3_so2",
             "region4_so2",
             "baseline_plant_so2",
             "intercept_pm",
             "grade_pm",
             "region2_pm",
             "region3_pm",
             "region4_pm",
             "baseline_plant_pm")
keep <- c("GMR", "CL","Pr(>|t|)")

# model 1
# run n regressions
my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_ab))$coefficient[, 1:4]})
results_1969ab_linear <-do.call(rbind.data.frame, my_lms)
rownames(results_1969ab_linear) = rowlist1
results_1969ab_linear$GMR=round(exp(results_1969ab_linear$Estimate), digits=2)
results_1969ab_linear$LCL_exp=round(exp(results_1969ab_linear$Estimate-1.96*results_1969ab_linear$"Std. Erro"), digits=2)
results_1969ab_linear$UCL_exp=round(exp(results_1969ab_linear$Estimate+1.96*results_1969ab_linear$"Std. Erro"), digits=2)
results_1969ab_linear$CL <- paste(results_1969ab_linear$LCL_exp, results_1969ab_linear$UCL_exp, sep=",")
results_1969ab_linear<-results_1969ab_linear[c(-1,-7,-13),]
results_1969ab_linear<-results_1969ab_linear[keep]

my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_bc))$coefficient[, 1:4]})
results_1969bc_linear <-do.call(rbind.data.frame, my_lms)
rownames(results_1969bc_linear) = rowlist1
results_1969bc_linear$GMR=round(exp(results_1969bc_linear$Estimate), digits=2)
results_1969bc_linear$LCL_exp=round(exp(results_1969bc_linear$Estimate-1.96*results_1969bc_linear$"Std. Erro"), digits=2)
results_1969bc_linear$UCL_exp=round(exp(results_1969bc_linear$Estimate+1.96*results_1969bc_linear$"Std. Erro"), digits=2)
results_1969bc_linear$CL <- paste(results_1969bc_linear$LCL_exp, results_1969bc_linear$UCL_exp, sep=",")
results_1969bc_linear<-results_1969bc_linear[c(-1,-7,-13),]
results_1969bc_linear<-results_1969bc_linear[keep]

my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_cd))$coefficient[, 1:4]})
results_1969cd_linear <-do.call(rbind.data.frame, my_lms)
rownames(results_1969cd_linear) = rowlist1
results_1969cd_linear$GMR=round(exp(results_1969cd_linear$Estimate), digits=2)
results_1969cd_linear$LCL_exp=round(exp(results_1969cd_linear$Estimate-1.96*results_1969cd_linear$"Std. Erro"), digits=2)
results_1969cd_linear$UCL_exp=round(exp(results_1969cd_linear$Estimate+1.96*results_1969cd_linear$"Std. Erro"), digits=2)
results_1969cd_linear$CL <- paste(results_1969cd_linear$LCL_exp, results_1969cd_linear$UCL_exp, sep=",")
results_1969cd_linear<-results_1969cd_linear[c(-1,-7,-13),]
results_1969cd_linear<-results_1969cd_linear[keep]


result_linear <-cbind(results_1969ab_linear,
                      results_1969bc_linear,
                      results_1969cd_linear)
result_linear[,2]<-paste0("[", result_linear[,2], "]")
result_linear[,5]<-paste0("[", result_linear[,5], "]")
result_linear[,8]<-paste0("[", result_linear[,8], "]")

colnames(result_linear) <- c("ab GMR", "ab CL", "ab p-value", 
                             "bc GMR", "bc CL", "bc p-value",
                             "cd GMR", "cd CL", "cd p-value")

write.csv(result_linear,"result_all_linear_emission_model1.csv", row.names = TRUE)

#check for observations used for analysis
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_ab)))})
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_bc)))})
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_cd)))})
#R^2
# r1 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_ab))$r.squared})
# r2 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_bc))$r.squared})
# r3 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant, list(i = as.name(x))), data=spatial_notime_cd))$r.squared})
# r1 <-do.call(rbind.data.frame, r1)
# r2 <-do.call(rbind.data.frame, r2)
# r3 <-do.call(rbind.data.frame, r3)
# r_model1<-cbind(r1, r2, r3)

# model 2
rowlist2 <-c("intercept_nox",	
             "grade_nox",
             "region2_nox",
             "region3_nox",
             "region4_nox",
             'baseline_plant_nox',
             'medium_nox',
             'small_nox',
             'nowhite_nox',
             "intercept_so2",	
             "grade_so2",
             "region2_so2",
             "region3_so2",
             'region4_so2',
             'baseline_plant_so2',
             'medium_so2',
             'small_so2',	
             'nowhite_so2',
             "intercept_pm",	
             "grade_pm",
             "region2_pm",
             "region3_pm",
             'region4_pm',
             'baseline_plant_pm',
             'medium_pm',
             'small_pm',	
             'nowhite_pm'
)

# run n regressions
my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_ab))$coefficient[, 1:4]})
results_1969ab_linear2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969ab_linear2) = rowlist2
results_1969ab_linear2$GMR=round(exp(results_1969ab_linear2$Estimate), digits=2)
results_1969ab_linear2$LCL_exp=round(exp(results_1969ab_linear2$Estimate-1.96*results_1969ab_linear2$"Std. Erro"), digits=2)
results_1969ab_linear2$UCL_exp=round(exp(results_1969ab_linear2$Estimate+1.96*results_1969ab_linear2$"Std. Erro"), digits=2)
results_1969ab_linear2$CL <- paste(results_1969ab_linear2$LCL_exp, results_1969ab_linear2$UCL_exp, sep=",")
results_1969ab_linear2<-results_1969ab_linear2[c(-1,-10,-19),]
results_1969ab_linear2<-results_1969ab_linear2[keep]

my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_bc))$coefficient[, 1:4]})
results_1969bc_linear2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969bc_linear2) = rowlist2
results_1969bc_linear2$GMR=round(exp(results_1969bc_linear2$Estimate), digits=2)
results_1969bc_linear2$LCL_exp=round(exp(results_1969bc_linear2$Estimate-1.96*results_1969bc_linear2$"Std. Erro"), digits=2)
results_1969bc_linear2$UCL_exp=round(exp(results_1969bc_linear2$Estimate+1.96*results_1969bc_linear2$"Std. Erro"), digits=2)
results_1969bc_linear2$CL <- paste(results_1969bc_linear2$LCL_exp, results_1969bc_linear2$UCL_exp, sep=",")
results_1969bc_linear2<-results_1969bc_linear2[c(-1,-10,-19),]
results_1969bc_linear2<-results_1969bc_linear2[keep]

my_lms <- lapply(dvList3, function(x) {
  summary(lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_cd))$coefficient[, 1:4]})
results_1969cd_linear2 <-do.call(rbind.data.frame, my_lms)
rownames(results_1969cd_linear2) = rowlist2
results_1969cd_linear2$GMR=round(exp(results_1969cd_linear2$Estimate), digits=2)
results_1969cd_linear2$LCL_exp=round(exp(results_1969cd_linear2$Estimate-1.96*results_1969cd_linear2$"Std. Erro"), digits=2)
results_1969cd_linear2$UCL_exp=round(exp(results_1969cd_linear2$Estimate+1.96*results_1969cd_linear2$"Std. Erro"), digits=2)
results_1969cd_linear2$CL <- paste(results_1969cd_linear2$LCL_exp, results_1969cd_linear2$UCL_exp, sep=",")
results_1969cd_linear2<-results_1969cd_linear2[c(-1,-10,-19),]
results_1969cd_linear2<-results_1969cd_linear2[keep]


result_linear2 <-cbind(results_1969ab_linear2, 
                       results_1969bc_linear2,
                       results_1969cd_linear2)

result_linear2[,2]<-paste0("[", result_linear2[,2], "]")
result_linear2[,5]<-paste0("[", result_linear2[,5], "]")
result_linear2[,8]<-paste0("[", result_linear2[,8], "]")

colnames(result_linear2) <- c("ab GMR", "ab CL", "ab p-value", 
                              "bc GMR", "bc CL", "bc p-value",
                              "cd GMR", "cd CL", "cd p-value")

write.csv(result_linear2,"result_all_linear_emission_model2.csv", row.names = TRUE)

#check for observations used for analysis
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_ab)))})
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_bc)))})
# obs <- lapply(dvList3, function(x) {
#   nobs((lm(substitute(i ~ holc_grade +regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_cd)))})
#R^2
# r1 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_ab))$r.squared})
# r2 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_bc))$r.squared})
# r3 <- lapply(dvList3, function(x) {summary(lm(substitute(i ~ holc_grade + regions+baseline_plant+citypop_cat+pct_nowhite, list(i = as.name(x))), data=spatial_notime_cd))$r.squared})
# r1 <-do.call(rbind.data.frame, r1)
# r2 <-do.call(rbind.data.frame, r2)
# r3 <-do.call(rbind.data.frame, r3)
# r_model2<-cbind(r1, r2, r3)


############################################################################################################
############################################################################################################

######### Figure 3 and 5 Adjusted Effect Estimates for 5km #########

#Create forest plot for 5km
############################################################################################################
#################################getting the estimate for binary outcome###################################
#import data
result_spatial_time<-rbind(result_poisson, result_nb)
result_spatial_time<-result_spatial_time[,c(-3,-6,-9,-12,-15,-18,-21,-24,-27)]
result_spatial_time<-result_spatial_time[c(1,6,11,16),]
result_spatial_time$Outcome<-c("Presence of power plant within 5km (PR)", "Presence of coal or oil plant within 5km (PR)", "Presence of peaker plant within 5km (PR)", "Number of power plants within 5km (IRR)")
rownames(result_spatial_time)<-c(1:4)

result_spatial_time_1 <- melt(result_spatial_time, id.vars = c("Outcome"),
                              measure.vars = c(1,3,5,7,9,11,13,15,17))

colnames(result_spatial_time_1)[3]<-c("estimate")

result_spatial_time_2 <- melt(result_spatial_time, id.vars = c("Outcome"),
                              measure.vars = c(2,4,6,8,10,12,14,16,18))
result_spatial_time_2[c('lower_CL', 'upper_CL')] <- str_split_fixed(result_spatial_time_2$value, ',', 2)
result_spatial_time_2<-result_spatial_time_2[,c(4,5)]
result_spatial_time_3<-cbind(result_spatial_time_1, result_spatial_time_2)

result_spatial_time_3$year<-substr(result_spatial_time_3$variable, start = 1, stop = 4)
result_spatial_time_3$comparison<-substr(result_spatial_time_3$variable, start = 5, stop = 6)
result_spatial_time_3$outcome_new <- paste(result_spatial_time_3$Outcome, result_spatial_time_3$year)


result_spatial_time_4 <- result_spatial_time_3[
  order( result_spatial_time_3[,1], result_spatial_time_3[,2] ),
]

result_spatial_time_4$index<-c('3', 
                               '2', 
                               '1', 
                               '3', 
                               '2', 
                               '1', 
                               '3', 
                               '2', 
                               '1', 
                               '9', 
                               '8', 
                               '7', 
                               '9', 
                               '8', 
                               '7', 
                               '9', 
                               '8', 
                               '7', 
                               '6', 
                               '5', 
                               '4', 
                               '6', 
                               '5', 
                               '4', 
                               '6', 
                               '5', 
                               '4', 
                               '12', 
                               '11', 
                               '10', 
                               '12', 
                               '11', 
                               '10', 
                               '12', 
                               '11', 
                               '10' 
)


result_spatial_time_4$lower_CL<-as.numeric(result_spatial_time_4$lower_CL)
result_spatial_time_4$upper_CL<-as.numeric(result_spatial_time_4$upper_CL)
result_spatial_time_4$index<-as.numeric(result_spatial_time_4$index)
result_spatial_time_4$variable<-as.character(result_spatial_time_4$variable)

result_spatial_time_4 <- result_spatial_time_4[
  order( result_spatial_time_4[,6] ),
]

library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
p <-ggplot(data=result_spatial_time_4, aes(y=index, x=estimate, xmin=lower_CL, xmax=upper_CL)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(name = "", breaks=1:nrow(result_spatial_time_4))+
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()+
  coord_trans(x ="log10", xlim =c(0.5,3))
p + facet_grid(cols= vars(comparison), scales = "free")

#########################################################################################
#########################################################################################
result_spatial_no_time<-result_linear
result_spatial_no_time<-result_spatial_no_time[,c(-3,-6,-9)]
result_spatial_no_time<-result_spatial_no_time[c(1,6,11),]
result_spatial_no_time$Outcome<-c("2019 annual NOx (tons)", "2019 annual SO2 (tons)", "2018 annual PM2.5 (tons)")
rownames(result_spatial_no_time)<-c(1:3)

result_spatial_no_time_1 <- melt(result_spatial_no_time, id.vars = c("Outcome"),
                                 measure.vars = c(1,3,5))

colnames(result_spatial_no_time_1)[3]<-c("estimate")

result_spatial_no_time_2 <- melt(result_spatial_no_time, id.vars = c("Outcome"),
                                 measure.vars = c(2,4,6))
result_spatial_no_time_2[c('lower_CL', 'upper_CL')] <- str_split_fixed(result_spatial_no_time_2$value, ',', 2)

result_spatial_no_time_2<-result_spatial_no_time_2[,c(4,5)]
result_spatial_no_time_3<-cbind(result_spatial_no_time_1, result_spatial_no_time_2)

result_spatial_no_time_3$index<-c('3', 
                                  '2', 
                                  '1', 
                                  '3', 
                                  '2', 
                                  '1', 
                                  '3', 
                                  '2',
                                  '1')

result_spatial_no_time_3$lower_CL<-as.numeric(result_spatial_no_time_3$lower_CL)
result_spatial_no_time_3$upper_CL<-as.numeric(result_spatial_no_time_3$upper_CL)
result_spatial_no_time_3$index<-as.numeric(result_spatial_no_time_3$index)
result_spatial_no_time_3$variable<-as.character(result_spatial_no_time_3$variable)

p <-ggplot(data=result_spatial_no_time_3, aes(y=index, x=estimate, xmin=lower_CL, xmax=upper_CL)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(name = "", breaks=1:nrow(result_spatial_no_time_3))+
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()+
  coord_trans(x ="log10", xlim =c(0.5,2.5))
p + facet_grid(cols= vars(variable), scales = "free")



############################################################################################################
############################################################################################################
######### Figure S2 and S3 Adjusted Effect Estimates for 10km #########

#Create forest plot for 10km
############################################################################################################
#################################getting the estimate for binary outcome###################################
#import data
result_spatial_time<-rbind(result_poisson, result_nb)
result_spatial_time<-result_spatial_time[,c(-3,-6,-9,-12,-15,-18,-21,-24,-27)]
result_spatial_time<-result_spatial_time[c(1,6,11,16),]
result_spatial_time$Outcome<-c("Presence of power plant within 5km (PR)", "Presence of coal or oil plant within 5km (PR)", "Presence of peaker plant within 5km (PR)", "Number of power plants within 5km (IRR)")
rownames(result_spatial_time)<-c(1:4)

result_spatial_time_1 <- melt(result_spatial_time, id.vars = c("Outcome"),
                              measure.vars = c(1,3,5,7,9,11,13,15,17))

colnames(result_spatial_time_1)[3]<-c("estimate")

result_spatial_time_2 <- melt(result_spatial_time, id.vars = c("Outcome"),
                              measure.vars = c(2,4,6,8,10,12,14,16,18))
result_spatial_time_2[c('lower_CL', 'upper_CL')] <- str_split_fixed(result_spatial_time_2$value, ',', 2)
result_spatial_time_2<-result_spatial_time_2[,c(4,5)]
result_spatial_time_3<-cbind(result_spatial_time_1, result_spatial_time_2)

result_spatial_time_3$year<-substr(result_spatial_time_3$variable, start = 1, stop = 4)
result_spatial_time_3$comparison<-substr(result_spatial_time_3$variable, start = 5, stop = 6)
result_spatial_time_3$outcome_new <- paste(result_spatial_time_3$Outcome, result_spatial_time_3$year)


result_spatial_time_4 <- result_spatial_time_3[
  order( result_spatial_time_3[,1], result_spatial_time_3[,2] ),
]

result_spatial_time_4$index<-c('3', 
                               '2', 
                               '1', 
                               '3', 
                               '2', 
                               '1', 
                               '3', 
                               '2', 
                               '1', 
                               '9', 
                               '8', 
                               '7', 
                               '9', 
                               '8', 
                               '7', 
                               '9', 
                               '8', 
                               '7', 
                               '6', 
                               '5', 
                               '4', 
                               '6', 
                               '5', 
                               '4', 
                               '6', 
                               '5', 
                               '4', 
                               '12', 
                               '11', 
                               '10', 
                               '12', 
                               '11', 
                               '10', 
                               '12', 
                               '11', 
                               '10' 
)


result_spatial_time_4$lower_CL<-as.numeric(result_spatial_time_4$lower_CL)
result_spatial_time_4$upper_CL<-as.numeric(result_spatial_time_4$upper_CL)
result_spatial_time_4$index<-as.numeric(result_spatial_time_4$index)
result_spatial_time_4$variable<-as.character(result_spatial_time_4$variable)

result_spatial_time_4 <- result_spatial_time_4[
  order( result_spatial_time_4[,6] ),
]

library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
p <-ggplot(data=result_spatial_time_4, aes(y=index, x=estimate, xmin=lower_CL, xmax=upper_CL)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(name = "", breaks=1:nrow(result_spatial_time_4))+
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()+
  coord_trans(x ="log10", xlim =c(0.5,3))
p + facet_grid(cols= vars(comparison), scales = "free")

#########################################################################################
#########################################################################################
result_spatial_no_time<-result_linear
result_spatial_no_time<-result_spatial_no_time[,c(-3,-6,-9)]
result_spatial_no_time<-result_spatial_no_time[c(1,6,11),]
result_spatial_no_time$Outcome<-c("2019 annual NOx (tons)", "2019 annual SO2 (tons)", "2018 annual PM2.5 (tons)")
rownames(result_spatial_no_time)<-c(1:3)

result_spatial_no_time_1 <- melt(result_spatial_no_time, id.vars = c("Outcome"),
                                 measure.vars = c(1,3,5))

colnames(result_spatial_no_time_1)[3]<-c("estimate")

result_spatial_no_time_2 <- melt(result_spatial_no_time, id.vars = c("Outcome"),
                                 measure.vars = c(2,4,6))
result_spatial_no_time_2[c('lower_CL', 'upper_CL')] <- str_split_fixed(result_spatial_no_time_2$value, ',', 2)

result_spatial_no_time_2<-result_spatial_no_time_2[,c(4,5)]
result_spatial_no_time_3<-cbind(result_spatial_no_time_1, result_spatial_no_time_2)

result_spatial_no_time_3$index<-c('3', 
                                  '2', 
                                  '1', 
                                  '3', 
                                  '2', 
                                  '1', 
                                  '3', 
                                  '2',
                                  '1')

result_spatial_no_time_3$lower_CL<-as.numeric(result_spatial_no_time_3$lower_CL)
result_spatial_no_time_3$upper_CL<-as.numeric(result_spatial_no_time_3$upper_CL)
result_spatial_no_time_3$index<-as.numeric(result_spatial_no_time_3$index)
result_spatial_no_time_3$variable<-as.character(result_spatial_no_time_3$variable)

p <-ggplot(data=result_spatial_no_time_3, aes(y=index, x=estimate, xmin=lower_CL, xmax=upper_CL)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(name = "", breaks=1:nrow(result_spatial_no_time_3))+
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()+
  coord_trans(x ="log10", xlim =c(0.5,2.5))
p + facet_grid(cols= vars(variable), scales = "free")


