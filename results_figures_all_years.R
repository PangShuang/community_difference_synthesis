################################################################################
##  results_figures_all_years: Compiles Bayesian output and makes figures for the primary analysis of richness and compositonal differences between treatment and control plots.
##
##  Author: Kimberly La Pierre
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


###main figure (Figure 1)
# compositional response panels 
#------------------------
meanPlot0 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(b) 63.0%', size=10, hjust='left') +
  #below are the individual treatment lines
  



#------------------------
meanPlot1 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(b) 63.0%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
#------------------------
meanPlot2 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(f) 0.9%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
#------------------------
meanPlot3 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(h) 3.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  
#------------------------
meanPlot4 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(j) 0%', size=10, hjust='left')
  #below are the individual treatment lines
  

#------------------------ 
meanPlot5 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(l) 0%', size=10, hjust='left')
  #below are the individual treatment lines
  

#------------------------
meanPlot6 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(n) 0%', size=10, hjust='left')
  #below are the individual treatment lines
  

#------------------------
meanPlot7 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(p) 9.4%', size=10, hjust='left') +
  #below are the individual treatment lines  
  
#------------------------
meanPlot8 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(r) 0.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  
  


#richness response panels 
#------------------------
richnessPlot0 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(a) 76.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  

#------------------------
richnessPlot1 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(c) 0.9%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
#------------------------
richnessPlot2 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(e) 0.2%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
#------------------------
richnessPlot3 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(g)', size=10, hjust='left')
  #below are the individual treatment lines
  
#-------------------
richnessPlot4 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(i) 0%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
#------------------------
richnessPlot5 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(k) 0.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  
#------------------------
richnessPlot6 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(m) 2.5%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  
#------------------------
richnessPlot7 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(o) 4.1%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  
  
#------------------------
richnessPlot8 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(q) 5.5%', size=10, hjust='left') +
  #below are the individual treatment lines
  
  

  
#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(9,2)))
print(richnessPlot0, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(richnessPlot1, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(richnessPlot2, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(richnessPlot3, vp=viewport(layout.pos.row=4, layout.pos.col=1))
print(richnessPlot4, vp=viewport(layout.pos.row=5, layout.pos.col=1))
print(richnessPlot5, vp=viewport(layout.pos.row=6, layout.pos.col=1))
print(richnessPlot6, vp=viewport(layout.pos.row=7, layout.pos.col=1))
print(richnessPlot7, vp=viewport(layout.pos.row=8, layout.pos.col=1))
print(richnessPlot8, vp=viewport(layout.pos.row=9, layout.pos.col=1))
print(meanPlot0, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanPlot1, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(meanPlot2, vp=viewport(layout.pos.row=3, layout.pos.col=2))
print(meanPlot3, vp=viewport(layout.pos.row=4, layout.pos.col=2))
print(meanPlot4, vp=viewport(layout.pos.row=5, layout.pos.col=2))
print(meanPlot5, vp=viewport(layout.pos.row=6, layout.pos.col=2))
print(meanPlot6, vp=viewport(layout.pos.row=7, layout.pos.col=2))
print(meanPlot7, vp=viewport(layout.pos.row=8, layout.pos.col=2))
print(meanPlot8, vp=viewport(layout.pos.row=9, layout.pos.col=2))
#export at 1200 x 3600



###by magnitude of resource manipulated---------------------------------
#N addition
nData <- read.csv('ForAnalysis_allAnalysisNmag.csv')

nDataMean <- nData%>%
  summarise(mean_mean_change=mean(composition_diff), sd_mean_change=sd(composition_diff), mean_S_PC=mean(S_PC), sd_S_PC=sd(S_PC))

#mean change
Nmean <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\magnitude_042019\\posteriors_N_MeanChange.csv', comment.char='#')
NmeanMean <- as.data.frame(colMeans(Nmean))%>%
  add_rownames('parameter')
names(NmeanMean)[names(NmeanMean) == 'colMeans(Nmean)'] <- 'mean'
NmeanSD <- as.data.frame(colSd(Nmean))%>%
  add_rownames('parameter')
names(NmeanSD)[names(NmeanSD) == 'colSd(Nmean)'] <- 'sd'
NmeanOverall <- NmeanMean%>%
  left_join(NmeanSD)


# #get mean and sd to transform
# nDataSummary <- nData%>%
#   summarise(mean_change_mean=mean(composition_diff), mean_change_sd=sd(composition_diff), richness_mean=mean(S_PC), richness_sd=sd(S_PC), n_mean=mean(n), n_sd=sd(n), MAP_mean=mean(MAP), MAP_sd=sd(MAP))

nDataTransform <- nData%>%
  #transform mean change
  mutate(mean_change_transform=((composition_diff-mean(composition_diff))/sd(composition_diff)))%>%
  #transform proportional richness change
  mutate(S_PC_transform=((S_PC-mean(S_PC))/sd(S_PC)))%>%
  #transform N treatment magnitude
  mutate(n_transform=((n-mean(n))/sd(n)))%>%
  #transform MAP
  mutate(MAP_transform=((MAP-mean(MAP))/sd(MAP)))

meanNPlotFinal <- ggplot(data=subset(nData), aes(x=n, y=composition_diff, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='Composition Response') +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((1000-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((1000-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#4793CF')  +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((600-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((600-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#2D5E88') +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((200-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((200-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#153049') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.0, label='(d)', size=12, hjust='left') +
  theme(legend.position=c(0.8,0.05), legend.justification=c(0,0), legend.title=element_text(size=24))


#richness difference
Nrichness <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\nate_results\\manipulation\\posteriors_N_Richness.csv', comment.char='#')
NrichnessMean <- as.data.frame(colMeans(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessMean)[names(NrichnessMean) == 'colMeans(Nrichness)'] <- 'mean'
NrichnessSD <- as.data.frame(colSd(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessSD)[names(NrichnessSD) == 'colSd(Nrichness)'] <- 'sd'
NrichnessOverall <- NrichnessMean%>%
  left_join(NrichnessSD)

richnessNPlotFinal <- ggplot(data=nData, aes(x=n, y=S_PC, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(-0.8,1)) +
  scale_y_continuous(name='Richness Response') +
  stat_function(fun=function(x){(-0.005589416 + -0.562241618*(x-9.992142)/9.108662)*0.2548196-0.1338463}, size=5) +
  xlab('') +
  annotate('text', x=1.0, y=1.0, label='(a)', size=12, hjust='left') +
  theme(legend.position='none')


#drought change
droData <- read.csv('ForAnalysis_allAnalysisH2Omag_drought.csv')

meanDroPlotFinal <- ggplot(data=droData, aes(x=precip, y=composition_diff)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=1, label='(e)', size=12, hjust='left')

richnessDroPlotFinal <- ggplot(data=droData, aes(x=precip, y=S_PC, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(-0.8,1)) +
  scale_y_continuous(name='') +
  xlab('') +
  annotate('text', x=-80, y=1, label='(b)', size=12, hjust='left') +
  theme(legend.position='none')


#irrigation change
irrData <- read.csv('ForAnalysis_allAnalysisH2Omag_irr.csv')

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













