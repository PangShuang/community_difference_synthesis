################################################################################
##  results_figures_all_years: Compiles Bayesian output and makes figures for the primary analysis of richness and compositonal differences between treatment and control plots.
##
##  Author: Kimberly Komatsu
##  Date created: January 17, 2018
##  See https://github.com/klapierre/Converge_Diverge/blob/master/core%20data%20paper_bayesian%20results_figures_sig%20test_expinteractions_20yr.R for full history.
################################################################################

library(grid)
library(tidyverse)

#kim laptop
setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm")

#kim desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}

#function to get standard deviations of columns in a dataframe
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

##################################################################################
##################################################################################
#import experiment information --------------------------------------------------------
expRaw <- read.csv('ExperimentInformation_March2019.csv')


expInfo <- expRaw%>%
  #remove any pre-treatment data for the few experiments that have it -- pre-treatment data for experiments is awesome and we should all strive to collect it!
  filter(treatment_year!=0)%>%
  #make columns for irrigation and drought from precip column
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  #calcualte minumum years for each project
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought), experiment_length=max(treatment_year))

#import treatment data
trtInfo1 <- read.csv('ExperimentInformation_March2019.csv')

#import diversity metrics that went into Bayesian analysis
rawData <- read.csv('ForAnalysis_allAnalysisAllDatasets_04082019.csv')

#calculate means and standard deviations across all data for richness and compositonal differences to backtransform
rawData2<- rawData%>%
  left_join(trtInfo1)%>%
  filter(anpp!='NA', treatment_year!=0)%>%
  summarise(mean_mean=mean(composition_diff), std_mean=sd(composition_diff), mean_rich=mean(S_lnRR), std_rich=sd(S_lnRR)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  left_join(trtInfo1)%>%
  filter(anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length))

#for table of experiment summarizing various factors
expInfoSummary <- rawData%>%
  left_join(trtInfo1)%>%
  filter(anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_mean=mean(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=mean(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_mean=mean(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_mean=mean(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_mean=mean(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_mean=mean(MAT), MAT_min=min(MAT), MAT_max=max(MAT))%>%
  gather(variable, estimate)

#treatment info
trtInfo2 <- trtInfo1%>%
  select(site_code, project_name, community_type, treatment, plot_mani, trt_type)%>%
  unique()
  
trtInfo <- rawData%>%
  select(site_code, project_name, community_type, treatment, trt_type, experiment_length, rrich, anpp, MAT, MAP)%>%
  unique()%>%
  left_join(expInfo)


################################################################################
trtShape <- read.csv('stdtimebytrt_shape_classification_04072019.csv')

#do the parabolic shapes cluster in certain trt types?
para <- trtShape%>%filter(N01_shape>6)%>%group_by(trt_type, variable)%>%summarise(n=length(variable))%>%ungroup()


###main figure (Figure 1)
# compositional response panels------------ 
#------------------------
meanPlot0 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(j) 64.4%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.545480733206665 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.504045957021865 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.50367841925215 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6259143004259 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.9271318519415 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6577366494802 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.54717614780275 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.7691894927915 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.508050154259935 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.760856177924 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.34300590478995 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.532765724114475 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.98770558265 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.9791406781 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.2376585794 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.702553974635 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.5359916485545 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.55221453564825 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.394852928070305 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.8871011677 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.12093266155 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.1529636751 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.0999206158 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.84519298515 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.74702547968655 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.501572583840015 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.52092184493615 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.49392402759095 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.670362081235515 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.62918121790425 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.4366389011433 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(-1.11472587315 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-12.4)/8.13428956127495) + 0*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(-0.6947335607 + 0*((x-12.4)/8.13428956127495) + 0*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(-0.45419539490155 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.5413393330107 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.595589913364665 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.510148863962825 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.6046016951045 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.8765448348 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-1.46236688885 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.48132610355 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.4282105185 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.3485251765 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.41247707715 + 0*((x-2)/1.82574185835055) + 0*((x-2)/1.82574185835055)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.3557881426 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.12521574425 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.98397474633 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6122132334114 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.812490673855 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6350909278199 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.820090585135 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.873271353405 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.58052938308255 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.66263617918535 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.7112730361855 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.933594453465 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8196480894635 + 0*((x-0.5)/0.707106781186548) + 0*((x-0.5)/0.707106781186548)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.54688832083035 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.945971467235 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8145671214 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.4244541049696 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.5800785026429 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.68802146450225 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.61826160934515 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.905688863525 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.136039694 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.851569504715 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.2595024512 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.0010112708 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.778125626249 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.980487691255 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.00594504958 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8356944125405 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(-0.566422788495 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(-0.612674536965 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.7990143139 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.478779171577017 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.5371836977597 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.899392652775 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.909070571015 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-0.767175357895 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.50071538698345 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.494416543905975 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.562203525027955 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.54718912006985 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.48604732128265 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6666443147044 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.64911421345625 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.845947630602 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.94346259944 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6453278423 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.74712925535 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.6867478263 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.35700259785915 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.7179803302 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-1.09128317765 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.59140875908 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.8554793132 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.6115479184417 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.6606535244515 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.53899831885705 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.61607761700525 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8184050907725 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.10314257985 + 0*((x-3.85714285714286)/3.02371578407382) + 0*((x-3.85714285714286)/3.02371578407382)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5245327834018 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.563549018135 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.4929689143508 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.47989270748365 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.51641315044845 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.87410271805 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.80478865422 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.46546882971055 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.49510630436085 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.657745614635 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.665669995855 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-1.0289388977 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.612667394291 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.06171917175 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.8388129089 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.8855030433 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.780139908805 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.6890367846175 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.01663735065 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/3.1091263510296) + 0*((x-2.5)/3.1091263510296)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-1.19165092655 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.652933286087 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.576399093452 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.655513949034 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.50829436093355 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-1.38904248705 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.09022906045 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.526756346224735 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.38752314319431 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.526728356572 + 0*((x-4)/2.73861278752583) + 0*((x-4)/2.73861278752583)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.6871925186435 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.5671588196935 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.52526170972715 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.571114582735 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.725071650938 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.64037959996015 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.752402493865 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.837331703075 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.44087567540385 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5122297855131 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.89569963532395 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.61531488178865 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6217351624599 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5844351171141 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2))
  



#------------------------
meanPlot1 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(k) 21.2%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.6803039682 + 0.430476176647*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.74415501475 + 0.56593681165*((x-5.63636363636364)/3.93122696553448) + 0*((x-5.63636363636364)/3.93122696553448)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.5077651173095 + 0.70707648185*((x-4.55555555555556)/3.39525813124389) + 0*((x-4.55555555555556)/3.39525813124389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.939197666785 + 0.58603055848825*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.52800332965 + 0.5511530404556*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.7231580312 + 0.610295458185*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.822079113545 + 0.43057015238*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.97960949465 + 0.32053652899245*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.85736697415 + 0.38657382636377*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.3939441859328*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.680654619855 + 0.334713303344242*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.677419098281385 + 0.408112037863755*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.853619161826 + 0.40618255251205*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8764962578 + 0.29837308447585*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.0068238421 + 0.3032016775261*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.48003943519555*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.25432929144525 + 0.292701064411*((x-14.5)/8.8034084308295) + 0*((x-14.5)/8.8034084308295)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(0.44765436463285 + 0.49452518667*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.67625446435 + 0.47569696294*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.3580524789886 + 0.22192015433835*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0.44287111605097*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.28796826685875*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.60073782851*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.5068574942194*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.613789725541*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.36753638431227*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.32592826664705*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.4301502386059*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.259430428669975*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.30987355187474*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.25480458089676*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0.37683853432689*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.348481005974*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.328623359513835*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(2.6116390945 + 0.22219075814475*((x-10.5454545454545)/6.57359271723775) + 0*((x-10.5454545454545)/6.57359271723775)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,22)) +
  stat_function(fun=function(x){(0 + 0.20472525251375*((x-11.5)/7.07106781186548) + 0*((x-11.5)/7.07106781186548)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,23)) +
  stat_function(fun=function(x){(2.4674400965 + 0.2355444930078*((x-10.5454545454545)/6.57359271723775) + 0*((x-10.5454545454545)/6.57359271723775)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,22)) +
  stat_function(fun=function(x){(0.8542267644 + 0.51872717605*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0.277726192532513*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.376528574015*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(-0.50335395754 + 0.1894311045816*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(0 + 0.4287243420425*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.380082103255495 + 0.500107560609*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.4005678004 + 0.6876418731*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.473297558201 + 0.39006194079*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.26045037652748*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.8717840402 + 0.446592661285*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.68878186085 + 0.3069201674818*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.722688624114265 + 0.582563623957*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.67019377830385 + 0.638830470767*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.555495613511405 + 0.397701666187*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.37094548691397*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.580788401945*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.8433867876 + 0.3472551039599*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.7803558616*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.918215782*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.83434552115*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.87364900275*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.9031644986 + 0.29817841354959*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.476933982582 + 0.9081486985*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.91729567145*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.39374122870315 + 0.80254220255*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.67069945685*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.83393073465 + 0.7622925016*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.81807536985*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.3567618750544 + 0.6713265659*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.40682767780303*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.445918082501455*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.41027941420365*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.777079564869 + 0.3425142535101*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.4806705429667*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.15340105225 + 0.5335746106655*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.626901996367 + 0.43233617108778*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.866875849055 + 0.52176735504*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.268837415479145*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.5210401485158 + 0.476075763562*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.5821444142905 + 0.57561196332*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.1556762707 + 0.58952983175*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.4276471874298 + 0.3028414700992*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.53852158321*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.325820528784135 + 0.4105799283045*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.4960894732815 + 0.4385016266227*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.56717341321465 + 0.36265423492823*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.6975386681 + 0.2173522062522*((x-8)/5.04975246918104) + 0*((x-8)/5.04975246918104)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.6283455843156 + 0.9831555087*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.45046067925 + 0.4271376600328*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.548984280387 + 0.396880747016*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.46902919659885 + 0.4576081689205*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.37631119970735*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.4029262177136*((x-1)/1) + 0*((x-1)/1)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.41158721058955*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.53122533002615*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0.549224870555*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4))

  
  
#------------------------
meanPlot2 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(o) 1.1%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.3586980483978 + 0.6703595036*((x-5)/3.3166247903554) + 0.186086150542465*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.592662817399 + 0.2644996450147*((x-5)/3.3166247903554) + 0.25718613591405*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.371873978805955*((x-5)/3.3166247903554) + 0.266917185904*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.29906943651425*((x-5)/3.3166247903554) + 0.2886446814222*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.73156942115 + 0.59735170005*((x-5)/3.3166247903554) + 0.160181017669155*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10))

  
  
#------------------------
meanPlot3 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(p) 3.2%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(1.9018044975 + 0.6200037609*((x-12.4)/8.13428956127495) + -0.1480241259277*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(2.4725529915 + 0.6659252855*((x-12.4)/8.13428956127495) + -0.253791414825*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(1.8511288835 + 0.319969817045*((x-14.5)/8.8034084308295) + -0.1579025123209*((x-14.5)/8.8034084308295)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(2.296333018 + 1.00170299275*((x-4.5)/3.02765035409749) + -0.4094681206415*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.00292234445 + 0.76663043295*((x-4.5)/3.02765035409749) + -0.24063731252835*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.9216300635 + 0.916410754*((x-4.5)/3.02765035409749) + -0.34339574293675*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0.3071684168334*((x-6)/3.89444048184931) + -0.214452960182675*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(1.1530448811 + 0.52754066004*((x-4.5)/3.02765035409749) + -0.27180220948376*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.2874751585 + 0.45517560019535*((x-10.4)/8.38450952650183) + -0.3894945451535*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.8558449948 + 0.375836336980095*((x-10.4)/8.38450952650183) + -0.323885182469815*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.6347570236 + 0.433408250725*((x-5.5)/3.60555127546399) + -0.18948855723609*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.9872661903 + 0.517159120775*((x-5)/3.3166247903554) + -0.16900730386952*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.7760028282 + 0.9193847094*((x-4.33333333333333)/3.278719262151) + -0.23044859642712*((x-4.33333333333333)/3.278719262151)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(1.1198145972 + 0.535434036975*((x-2)/1.58113883008419) + -0.312310815499*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4))

  
#------------------------
meanPlot4 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(l) 0%', size=15, hjust='left')
  #below are the individual treatment lines
  

#------------------------ 
meanPlot5 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(q) 0%', size=15, hjust='left')
  #below are the individual treatment lines
  

#------------------------
meanPlot6 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('Year of Experiment') +
  ylab('') +
  annotate('text', x=0, y=1, label='(r) 0%', size=15, hjust='left')
  #below are the individual treatment lines
  

#------------------------
meanPlot7 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(m) 9.6%', size=15, hjust='left') +
  #below are the individual treatment lines  
  stat_function(fun=function(x){(-0.395011920057185 + 0*((x-2)/1.58113883008419) + -0.22781205201671*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0.23487945604354*((x-6.30769230769231)/4.44193305113542) + -0.16609392205037*((x-6.30769230769231)/4.44193305113542)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + -0.2115421428786*((x-6.30769230769231)/4.44193305113542)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0.4318474776365 + 0.31148480736855*((x-6.30769230769231)/4.44193305113542) + -0.327798032475*((x-6.30769230769231)/4.44193305113542)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-0.313710325810765 + 0*((x-6.30769230769231)/4.44193305113542) + -0.228210351897*((x-6.30769230769231)/4.44193305113542)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-0.56006330782319 + 0*((x-6.30769230769231)/4.44193305113542) + -0.16721394665819*((x-6.30769230769231)/4.44193305113542)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(1.554994847 + 0*((x-14.5)/8.8034084308295) + -0.16334472050735*((x-14.5)/8.8034084308295)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(0.216378395276745 + 0*((x-11.625)/7.30581959810123) + -0.1989212860055*((x-11.625)/7.30581959810123)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,26)) +
  stat_function(fun=function(x){(1.04828739115 + -0.16772800586285*((x-11.625)/7.30581959810123) + -0.268499872966*((x-11.625)/7.30581959810123)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,26)) +
  stat_function(fun=function(x){(1.856767805 + 0*((x-11.625)/7.30581959810123) + -0.4170446933*((x-11.625)/7.30581959810123)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,26)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.16704347309252*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.19570036236166*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.164808191494804*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0.200749649105064*((x-6)/3.89444048184931) + -0.182624542295875*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.254158438392*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(-0.538210080575 + 0*((x-6)/3.89444048184931) + -0.170116008252525*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(1.0540430299 + 0*((x-10.4)/8.38450952650183) + -0.405653450785*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.45463343011396 + 0*((x-10.4)/8.38450952650183) + -0.368831228927315*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.744766359425 + 0*((x-10.4)/8.38450952650183) + -0.3113249493185*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.8606966598 + 0*((x-10.4)/8.38450952650183) + -0.29602125631555*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.72821377093 + 0*((x-10.4)/8.38450952650183) + -0.3809400054321*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.5054699825002 + 0*((x-10.4)/8.38450952650183) + -0.34246588549645*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.6953080298 + 0*((x-10.4)/8.38450952650183) + -0.24988960017935*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(1.3670729974 + 0*((x-10.4)/8.38450952650183) + -0.3292465913525*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.60771110163385 + 0*((x-10.4)/8.38450952650183) + -0.23805132390306*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + -0.335872903728396*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + -0.27788609160225*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-7.9375)/5.20856666144023) + -0.18632041834299*((x-7.9375)/5.20856666144023)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(-0.53851840823 + 0*((x-8)/5.04975246918104) + -0.364407863438*((x-8)/5.04975246918104)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.5612626151754 + 0*((x-2)/1.58113883008419) + -0.2342103005654*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + -0.25337391111503*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.45195946500037 + 0*((x-2)/1.58113883008419) + -0.203767693344245*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(2.4899950905 + 1.0031510469*((x-12.4)/8.13428956127495) + -0.4587720423*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(3.027837951 + 1.1991097208*((x-12.4)/8.13428956127495) + -0.5488815412*((x-12.4)/8.13428956127495)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + -0.26915183488465*((x-1.5)/1.29099444873581)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(2.9038057825 + 0.83233612835*((x-4.5)/3.02765035409749) + -0.49141921535*((x-4.5)/3.02765035409749)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.159697735040825*((x-6)/3.89444048184931)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(1.13205143975 + 0*((x-10.4)/8.38450952650183) + -0.3515405161109*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.80510660301 + 0*((x-10.4)/8.38450952650183) + -0.3242980867898*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.49055266415735 + 0*((x-10.4)/8.38450952650183) + -0.2201309364034*((x-10.4)/8.38450952650183)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0.62353208324 + 0.28457301718525*((x-5.5)/3.60555127546399) + -0.302824828935595*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11)) +
    stat_function(fun=function(x){(0.38975824907325 + 0.32296697133*((x-5.5)/3.60555127546399) + -0.3387971598615*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,11))

  
  
  
#------------------------
meanPlot8 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1.1))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(n) 0.5%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.64363426886 + 0*((x-5)/3.3166247903554) + 0.2695130746677*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.210021988195509*((x-5)/3.3166247903554) + 0.23105613000231*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=2, xlim=c(0,10))
  
  
  
  


#richness response panels------------ 
#------------------------
richnessPlot0 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(a) 77.7%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.650592009021 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(1.13083171091 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5940877836788 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0.487666091431 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.6264150098267 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.78283818560395 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(-0.3154749897207 + 0*((x-12.4)/8.13428956127495) + 0*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(0 + 0*((x-14.5)/8.8034084308295) + 0*((x-14.5)/8.8034084308295)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(0 + 0*((x-11.625)/7.30581959810123) + 0*((x-11.625)/7.30581959810123)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,26)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.60376731417895 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.5892379315675 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.82574185835055) + 0*((x-2)/1.82574185835055)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.98577706229 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.5979885167908 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.61438531256556 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.147419271538 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.86955055998205 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.7473322990956 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.980329199903 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.94938506524862 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.67068346526686 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.771467231321 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-0.5)/0.707106781186548) + 0*((x-0.5)/0.707106781186548)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.9467636269765 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.72286218476515 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.8015945126825 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.4635595511531 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.577034705875 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.34976046843965 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.33775528346405 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.35243635715063 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.55425884024 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.661491584654355 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.607088633674705 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.81187943663315 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8572512809 + 0*((x-10.5454545454545)/6.57359271723775) + 0*((x-10.5454545454545)/6.57359271723775)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,22)) +
  stat_function(fun=function(x){(0 + 0*((x-10.5454545454545)/6.57359271723775) + 0*((x-10.5454545454545)/6.57359271723775)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,22)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.8870235110185 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.89026285699165 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.5681757804485 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.7404078798 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.5846494120425 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.6866898944915 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.67943265546995 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.66057632253987 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.81530198322 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.6924297226635 + 0*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.64905947435968 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.670169755362536 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.8131765112405 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.09731714428865 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.176236140415 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.888716422522 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.69644855413035 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.485535105425995 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.49676541645428 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.6977978951752 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.548586606793415 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.73951201977 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.61699411939735 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.59740243286725 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/3.1091263510296) + 0*((x-2.5)/3.1091263510296)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.63081826336785 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.8012934940341 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.56008776692305 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.20179414335 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.10636148805 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.80828250729 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.5991666940875 + 0*((x-8)/5.04975246918104) + 0*((x-8)/5.04975246918104)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.629166080485 + 0*((x-7.9375)/5.20856666144023) + 0*((x-7.9375)/5.20856666144023)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.617926324233 + 0*((x-8)/5.04975246918104) + 0*((x-8)/5.04975246918104)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(1.0975945177 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.87643474341 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-4)/2.73861278752583) + 0*((x-4)/2.73861278752583)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.61475320033424 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.594343003960425 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.63794110580525 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.593413538308 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.78435135393325 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.4825988670629 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.8718872198175 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.75617965438965 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.815003046618 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(1.56324844085 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.63541840624099 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.59554868748475 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.5372180921315 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.9146431695025 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.2445415543 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2))

  
  

#------------------------
richnessPlot1 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(b) 0.2%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.67000616894983 + 0.370184221889863*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2))

  
  
#------------------------
richnessPlot2 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(f) 0%', size=15, hjust='left')
  #below are the individual treatment lines

    
  
#------------------------
richnessPlot3 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(g) 0%', size=15, hjust='left')
  #below are the individual treatment lines


  
#-------------------
richnessPlot4 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(c) 9.1%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + -0.41770228673315*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.3458910584223*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.36191836206965*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.4728528144592*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.38681950137475*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.5296404601618*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.345244279420645*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.36208579645053*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.3628618292547*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.8014051261375 + -0.43332442983147*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.42804269951375*((x-1)/1) + 0*((x-1)/1)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.405148175399975*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + -0.39812718185805*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + -0.40829995855865*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.610978884333 + -0.4915479178925*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.89789083342665 + -0.52850750349*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.4088449975983*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.50086932714235*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.51551411386325 + -0.66640669035*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-1.0532788899 + -0.9182304598*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + -0.526015470192995*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.4901730209704*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.54790015785 + -0.7552618819*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.54408938945 + -0.73410215085*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.817664147712 + -0.412015052128*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.79785025200305 + -0.42471979736001*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + -0.299194472988724*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + -0.286261179528314*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + -0.36744386644474*((x-4.55555555555556)/3.39525813124389) + 0*((x-4.55555555555556)/3.39525813124389)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.701411512495 + -0.3335981879547*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.7187275681175 + -0.80764256825*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + -0.31025200882605*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + -0.320782996153*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.41548348893615 + -0.478514308212*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + -0.43514721184*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + -0.225665205217257*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(-0.4234832688944 + -0.525843455245*((x-5.63636363636364)/3.93122696553448) + 0*((x-5.63636363636364)/3.93122696553448)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(-0.41957730933265 + -0.240937156125566*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(0 + -0.452761853155*((x-11.5)/7.07106781186548) + 0*((x-11.5)/7.07106781186548)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,23)) +
  stat_function(fun=function(x){(-0.542277303835 + -0.3609012922445*((x-11.625)/7.30581959810123) + 0*((x-11.625)/7.30581959810123)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,26)) 
  
  
#------------------------
richnessPlot5 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(h) 0.5%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + -0.30220571316688*((x-3.5)/2.44948974278318) + -0.2990267342657*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + -0.48135776152*((x-3.5)/2.44948974278318) + -0.24315536556735*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7))
  
  
  
  
#------------------------
richnessPlot6 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('Year of Experiment') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(i) 2.7%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-3.102091284 + -1.01257748045*((x-4.5)/3.02765035409749) + 0.8036600226*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-4.9400637535 + -1.0788990159*((x-4.5)/3.02765035409749) + 0.9792800206*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.8261056362 + -0.6873250233*((x-4.5)/3.02765035409749) + 0.4722011099245*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-3.288962019 + -1.2791187157*((x-4.5)/3.02765035409749) + 0.6702237996*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.52262975845 + -0.40972010141752*((x-4.5)/3.02765035409749) + 0.3429761638577*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-2.643309637 + -0.6435519235*((x-4.5)/3.02765035409749) + 0.4491099024845*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.2356541086 + -0.6040997*((x-5)/3.3166247903554) + 0.32559142799155*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-1.712944763 + -0.56279400535*((x-11.625)/7.30581959810123) + 0.240710489948*((x-11.625)/7.30581959810123)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,26)) +
  stat_function(fun=function(x){(-2.449392705 + -0.303558463065*((x-14.5)/8.8034084308295) + 0.2220554894946*((x-14.5)/8.8034084308295)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(-2.174090084 + -0.9503249484*((x-12.4)/8.13428956127495) + 0.47342528445*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(-2.95632 + -1.0834205519*((x-12.4)/8.13428956127495) + 0.50340868775*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30))
  

  
  
#------------------------
richnessPlot7 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(d) 5.0%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.620495545878 + 0*((x-3.5)/2.44948974278318) + -0.333058498201*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(1.1291262488 + 0*((x-3.75)/2.81577190634672) + -0.3049487342346*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.1341533837 + 0*((x-3.75)/2.81577190634672) + -0.46311407707*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.2008446141 + 0*((x-3.75)/2.81577190634672) + -0.50278835439*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.9033541678 + 0*((x-3.75)/2.81577190634672) + -0.438582481895*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.55771477553655 + 0*((x-3.85714285714286)/3.02371578407382) + -0.399389785947345*((x-3.85714285714286)/3.02371578407382)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.221092227 + 0*((x-3.75)/2.81577190634672) + -0.361044008534*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.46462606745 + 0*((x-3.75)/2.81577190634672) + -0.387152006506*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.67603380335 + 0*((x-3.75)/2.81577190634672) + -0.501107442505*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.16610973725 + 0*((x-3.75)/2.81577190634672) + -0.28717878147051*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.38036159873984 + -0.34096970849605*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.600865077982 + 0*((x-5)/3.3166247903554) + -0.18796292398005*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.3649931757596 + -0.274782907407785*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.40175809367805 + 0*((x-6)/3.89444048184931) + -0.16669812802297*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.5551974686874 + 0*((x-6)/3.89444048184931) + -0.20253955152085*((x-6)/3.89444048184931)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.5263556761875 + -0.21189837623905*((x-12.4)/8.13428956127495) + 0*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(1.17543346265 + 0.324511862397185*((x-2)/1.58113883008419) + -0.249972157410805*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.84609538322 + 0*((x-3.5)/2.44948974278318) + -0.358515586745*((x-3.5)/2.44948974278318)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0.95820603055 + 0.30566118621644*((x-3.75)/2.81577190634672) + -0.26699188723695*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.1668083223 + 0.2894939463382*((x-3.75)/2.81577190634672) + -0.32152655082853*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.25181614665 + 0.28629832972905*((x-3.75)/2.81577190634672) + -0.34072782501858*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.30280411645 + 0*((x-3.75)/2.81577190634672) + -0.43453353346765*((x-3.75)/2.81577190634672)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,8))
  
  

  
  
  
#------------------------
richnessPlot8 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.3))  +
  scale_x_continuous(limits=c(0,30), breaks=seq(4,30,5), labels=seq(5,30,5)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Response') +
  annotate('text', x=0, y=2.0, label='(e) 5.5%', size=15, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0.301447826221775*((x-2)/1.58113883008419)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.7487398982095 + 0*((x-10.4)/8.38450952650183) + 0.364892731240955*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-0.687413417266 + 0*((x-10.4)/8.38450952650183) + 0.3739524085455*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0.295446859626955*((x-10.4)/8.38450952650183) + 0.3942686475703*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.43411013887*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.31015400581755*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.3764619159309*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-2.3950727465 + 0.49201426165*((x-12.4)/8.13428956127495) + 0*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0.466488219709*((x-1.5)/1.29099444873581)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.5770625952574 + 0*((x-4.5)/3.02765035409749) + 0.4561795620635*((x-4.5)/3.02765035409749)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0.18432388414583*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0.206056753143395*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.8330653986 + -0.29290575240145*((x-5)/3.3166247903554) + 0.310636515261*((x-5)/3.3166247903554)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0.192183396638396*((x-6.30769230769231)/4.44193305113542)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-0.7034541959145 + 0.42048706490945*((x-10.4)/8.38450952650183) + 0.525221783935*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-0.61216208494903 + 0*((x-10.4)/8.38450952650183) + 0.466763756885*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.442683225505*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-0.52469706403585 + 0*((x-10.4)/8.38450952650183) + 0.3585465121225*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.292080211981935*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-0.85725215779 + 0*((x-10.4)/8.38450952650183) + 0.451452531407*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.352080979313215*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(0 + 0*((x-10.4)/8.38450952650183) + 0.3633458976162*((x-10.4)/8.38450952650183)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,21)) +
  stat_function(fun=function(x){(-2.020313481 + 0*((x-14.5)/8.8034084308295) + 0.2240675035869*((x-14.5)/8.8034084308295)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,29)) +
  stat_function(fun=function(x){(-1.9003253585 + 0.239539404698525*((x-12.4)/8.13428956127495) + 0.147712848168115*((x-12.4)/8.13428956127495)^2)*(0.340217)+(-0.1183477)}, size=2, xlim=c(0,30)) 
  
  
  

  
#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(9,2)))
print(richnessPlot0, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(richnessPlot1, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(richnessPlot4, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(richnessPlot7, vp=viewport(layout.pos.row=4, layout.pos.col=1))
print(richnessPlot8, vp=viewport(layout.pos.row=5, layout.pos.col=1))
print(richnessPlot2, vp=viewport(layout.pos.row=6, layout.pos.col=1))
print(richnessPlot3, vp=viewport(layout.pos.row=7, layout.pos.col=1))
print(richnessPlot5, vp=viewport(layout.pos.row=8, layout.pos.col=1))
print(richnessPlot6, vp=viewport(layout.pos.row=9, layout.pos.col=1))
print(meanPlot0, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanPlot1, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(meanPlot4, vp=viewport(layout.pos.row=3, layout.pos.col=2))
print(meanPlot7, vp=viewport(layout.pos.row=4, layout.pos.col=2))
print(meanPlot8, vp=viewport(layout.pos.row=5, layout.pos.col=2))
print(meanPlot2, vp=viewport(layout.pos.row=6, layout.pos.col=2))
print(meanPlot3, vp=viewport(layout.pos.row=7, layout.pos.col=2))
print(meanPlot5, vp=viewport(layout.pos.row=8, layout.pos.col=2))
print(meanPlot6, vp=viewport(layout.pos.row=9, layout.pos.col=2))
#export at 2400 x 4800



###by magnitude of resource manipulated---------------------------------
#N addition
nData <- read.csv('ForAnalysis_allAnalysisNmag.csv')

nDataMean <- nData%>%
  summarise(mean_mean_change=mean(composition_diff), sd_mean_change=sd(composition_diff), mean_S_PC=mean(S_PC), sd_S_PC=sd(S_PC))

#mean change
Nmean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_N_MeanChange.csv', comment.char='#')
NmeanMean <- as.data.frame(colMeans(Nmean))%>%
  add_rownames('parameter')
names(NmeanMean)[names(NmeanMean) == 'colMeans(Nmean)'] <- 'mean'
NmeanSD <- as.data.frame(colSd(Nmean))%>%
  add_rownames('parameter')
names(NmeanSD)[names(NmeanSD) == 'colSd(Nmean)'] <- 'sd'
NmeanOverall <- NmeanMean%>%
  left_join(NmeanSD)

#richness difference
Nrichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_N_Richness.csv', comment.char='#')
NrichnessMean <- as.data.frame(colMeans(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessMean)[names(NrichnessMean) == 'colMeans(Nrichness)'] <- 'mean'
NrichnessSD <- as.data.frame(colSd(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessSD)[names(NrichnessSD) == 'colSd(Nrichness)'] <- 'sd'
NrichnessOverall <- NrichnessMean%>%
  left_join(NrichnessSD)

###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=NmeanOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.n.', 'std_func.MAP.', 'std_func.n..std_func.MAP.'),
                   labels=c('intercept', 'N', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

richnessOverallPlot <- ggplot(data=NrichnessOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.n.', 'std_func.MAP.', 'std_func.n..std_func.MAP.'),
                   labels=c('intercept', 'N', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 1600x1000


#get mean and sd to transform
nDataSummary <- nData%>%
  summarise(mean_change_mean=mean(composition_diff), mean_change_sd=sd(composition_diff), richness_mean=mean(S_PC), richness_sd=sd(S_PC), n_mean=mean(n), n_sd=sd(n), MAP_mean=mean(MAP), MAP_sd=sd(MAP))

nDataTransform <- nData%>%
  #transform mean change
  mutate(mean_change_transform=((composition_diff-mean(composition_diff))/sd(composition_diff)))%>%
  #transform proportional richness change
  mutate(S_lnRR_transform=((S_lnRR-mean(S_lnRR))/sd(S_lnRR)))%>%
  #transform N treatment magnitude
  mutate(n_transform=((n-mean(n))/sd(n)))%>%
  #transform MAP
  mutate(MAP_transform=((MAP-mean(MAP))/sd(MAP)))

meanNPlotFinal <- ggplot(data=subset(nData), aes(x=n, y=composition_diff, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='Composition Response') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.0, label='(d)', size=12, hjust='left') +
  theme(legend.position='none')

richnessNPlotFinal <- ggplot(data=nData, aes(x=n, y=S_PC, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(-0.8,1)) +
  scale_y_continuous(name='Richness Response') +
  stat_function(fun=function(x){(-0.1436441 + -0.4789404*((1000-709.6087)/313.6992) + 0.3777062*(x-15.65652)/18.28532 + 0.28058497*((1000-709.6087)/313.6992)*(x-9.992142)/9.108662)*0.2170204-0.1416691}, size=5, color='#4793CF')  +
  stat_function(fun=function(x){(-0.1436441 + -0.4789404*((600-709.6087)/313.6992) + 0.3777062*(x-15.65652)/18.28532 + 0.28058497*((600-709.6087)/313.6992)*(x-9.992142)/9.108662)*0.2170204-0.1416691}, size=5, color='#2D5E88') +
  stat_function(fun=function(x){(-0.1436441 + -0.4789404*((200-709.6087)/313.6992) + 0.3777062*(x-15.65652)/18.28532 + 0.28058497*((200-709.6087)/313.6992)*(x-9.992142)/9.108662)*0.2170204-0.1416691}, size=5, color='#153049') +
  xlab('') +
  annotate('text', x=1.0, y=1.0, label='(a)', size=12, hjust='left') +
  theme(legend.position=c(0.75,0.70), legend.justification=c(0,0), legend.title=element_text(size=24)) 


#drought change
droData <- read.csv('ForAnalysis_allAnalysisH2Omag_drought.csv')

#mean change
dromean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_Drought_MeanChange.csv', comment.char='#')
dromeanMean <- as.data.frame(colMeans(dromean))%>%
  add_rownames('parameter')
names(dromeanMean)[names(dromeanMean) == 'colMeans(dromean)'] <- 'mean'
dromeanSD <- as.data.frame(colSd(dromean))%>%
  add_rownames('parameter')
names(dromeanSD)[names(dromeanSD) == 'colSd(dromean)'] <- 'sd'
dromeanOverall <- dromeanMean%>%
  left_join(dromeanSD)

#richness difference
drorichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_Drought_Richness.csv', comment.char='#')
drorichnessMean <- as.data.frame(colMeans(drorichness))%>%
  add_rownames('parameter')
names(drorichnessMean)[names(drorichnessMean) == 'colMeans(drorichness)'] <- 'mean'
drorichnessSD <- as.data.frame(colSd(drorichness))%>%
  add_rownames('parameter')
names(drorichnessSD)[names(drorichnessSD) == 'colSd(drorichness)'] <- 'sd'
drorichnessOverall <- drorichnessMean%>%
  left_join(drorichnessSD)

###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=dromeanOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.precip.', 'std_func.MAP.', 'std_func.precip..std_func.MAP.'),
                   labels=c('intercept', 'Drought', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

richnessOverallPlot <- ggplot(data=drorichnessOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.precip.', 'std_func.MAP.', 'std_func.precip..std_func.MAP.'),
                   labels=c('intercept', 'Drought', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 1600x1000

meanDroPlotFinal <- ggplot(data=droData, aes(x=precip, y=composition_diff)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-100, y=1, label='(e)', size=12, hjust='left')

richnessDroPlotFinal <- ggplot(data=droData, aes(x=precip, y=S_PC)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(-0.8,1)) +
  scale_y_continuous(name='') +
  xlab('') +
  annotate('text', x=-100, y=1, label='(b)', size=12, hjust='left') +
  theme(legend.position='none')


#irrigation change
irrData <- read.csv('ForAnalysis_allAnalysisH2Omag_irr.csv')

#mean change
irrmean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_Irr_MeanChange.csv', comment.char='#')
irrmeanMean <- as.data.frame(colMeans(irrmean))%>%
  add_rownames('parameter')
names(irrmeanMean)[names(irrmeanMean) == 'colMeans(irrmean)'] <- 'mean'
irrmeanSD <- as.data.frame(colSd(irrmean))%>%
  add_rownames('parameter')
names(irrmeanSD)[names(irrmeanSD) == 'colSd(irrmean)'] <- 'sd'
irrmeanOverall <- irrmeanMean%>%
  left_join(irrmeanSD)

#richness difference
irrrichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_Irr_Richness.csv', comment.char='#')
irrrichnessMean <- as.data.frame(colMeans(irrrichness))%>%
  add_rownames('parameter')
names(irrrichnessMean)[names(irrrichnessMean) == 'colMeans(irrrichness)'] <- 'mean'
irrrichnessSD <- as.data.frame(colSd(irrrichness))%>%
  add_rownames('parameter')
names(irrrichnessSD)[names(irrrichnessSD) == 'colSd(irrrichness)'] <- 'sd'
irrrichnessOverall <- irrrichnessMean%>%
  left_join(irrrichnessSD)

###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=irrmeanOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.precip.', 'std_func.MAP.', 'std_func.precip..std_func.MAP.'),
                   labels=c('intercept', 'Irr', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

richnessOverallPlot <- ggplot(data=irrrichnessOverall, aes(x=parameter, y=mean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.4)) +
  scale_x_discrete(limits=c('Intercept', 'std_func.precip.', 'std_func.MAP.', 'std_func.precip..std_func.MAP.'),
                   labels=c('intercept', 'Irr', 'MAP', 'N*MAP')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0))

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 1600x1000

meanIrrPlotFinal <- ggplot(data=irrData, aes(x=precip, y=composition_diff)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=0, y=1, label='(f)', size=12, hjust='left')

richnessIrrPlotFinal <- ggplot(data=irrData, aes(x=precip, y=S_PC, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(-0.8,1)) +
  scale_y_continuous(name='') +
  xlab('') +
  annotate('text', x=0, y=1, label='(c)', size=12, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,3)))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(richnessDroPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanDroPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessIrrPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(meanIrrPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
#export at 2700 x 1600











###for reviewer 1--------------------
#comparing the fitted lines with non-significant parameters (i.e., intercept, linear slope, quadratic slope) included vs lines where non-significant parameters are set to 0
plot18 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(1,19), breaks=seq(5,20,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.7, label='(a) ANG Watering - Spring Water Addition', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='watering'&treatment=='S'), aes(x=treatment_year-1, y=composition_diff), size=5) +
#below are the individual treatment lines
  stat_function(fun=function(x){(0.6803039682 + 0.430476176647*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.2770874)}, size=3, xlim=c(1,12), color='black') +
  #fix line below
  stat_function(fun=function(x){(0.6803039682 + 0.430476176647*((x-5.5)/3.60555127546399) + -0.09874*((x-5.5)/3.60555127546399)^2)*(0.1860342)+(0.2770874)}, size=3, xlim=c(1,12), color='grey')

plot79 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.5, label='(b) CAU RMAPC - Water Addition', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='RMAPC'&treatment=='H2O'), aes(x=treatment_year-2, y=composition_diff), size=5) +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.670362081235515 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3270874)}, size=3, xlim=c(1,7), color='black') +
  #fix line below
  stat_function(fun=function(x){(-0.670362081235515 + 0.146227*((x-2.33333333333333)/3.21455025366432) + -0.01875*((x-2.33333333333333)/3.21455025366432)^2)*(0.1860342)+(0.3270874)}, size=3, xlim=c(1,7), color='grey')

plot117 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.1, label='(c) CUL Culardoch - N Addition and Clipping', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='Culardoch'&treatment=='N10clip'), aes(x=treatment_year-1, y=composition_diff), size=5) +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-1.4282105185 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,5), color='black') +
  #fix line below
  stat_function(fun=function(x){(-1.4282105185 + -0.00574*((x-2)/1.58113883008419) + 0.053445*((x-2)/1.58113883008419)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,5), color='grey')

plot147 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.25, label='(d) IMGERS Yu - N Addition', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='Yu'&treatment=='N1'), aes(x=treatment_year, y=composition_diff), size=5) +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.8145671214 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,8), color='black') +
  #fix line below
  stat_function(fun=function(x){(-0.8145671214 + 0.102746*((x-3.5)/2.44948974278318) + 0.07383*((x-3.5)/2.44948974278318)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,8), color='grey')

plot192 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.92, label='(e) KBS T7 - N Addition and Tilling', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='T7'&treatment=='T1F1'), aes(x=treatment_year, y=composition_diff), size=5) +
  #below are the individual treatment lines
  stat_function(fun=function(x){(2.6116390945 + 0.22219075814475*((x-10.5454545454545)/6.57359271723775) + 0*((x-10.5454545454545)/6.57359271723775)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,24), color='black') +
  #fix line below
  stat_function(fun=function(x){(2.6116390945 + 0.22219075814475*((x-10.5454545454545)/6.57359271723775) + -0.05851*((x-10.5454545454545)/6.57359271723775)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,24), color='grey')

plot254 <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=1, y=0.6, label='(f) KUFS E6 - N and P Addition', size=10, hjust='left') +
  geom_point(data=subset(rawData, project_name=='E6'&treatment=='N4P8S0'), aes(x=treatment_year, y=composition_diff), size=5) +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.26045037652748*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,12), color='black') +
  #fix line below
  stat_function(fun=function(x){(0.18753 + 0.26045037652748*((x-5)/3.3166247903554) + -0.04082*((x-5)/3.3166247903554)^2)*(0.1860342)+(0.3070874)}, size=3, xlim=c(1,12), color='grey')

pushViewport(viewport(layout=grid.layout(2,3)))
print(plot18, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(plot79, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(plot117, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(plot147, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(plot192, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(plot254, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 3600x2400













