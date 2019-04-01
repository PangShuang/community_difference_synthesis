################################################################################
##  results_figures_10yr.R: Compiles Bayesian output and makes figures for the primary analysis of richness and compositonal differences between treatment and control plots for datasets cut off at 10 years.
##
##  Author: Kimberly La Pierre
##  Date created: December 19, 2018
################################################################################

library(grid)
library(tidyverse)

#kim
setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

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
expRaw <- read.csv('ExperimentInformation_Nov2017.csv')


expInfo <- expRaw%>%
  #remove any pre-treatment data for the few experiments that have it -- pre-treatment data for experiments is awesome and we should all strive to collect it!
  filter(treatment_year!=0)%>%
  #make columns for irrigation and drought from precip column
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  #calcualte minumum years for each project
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

#import treatment data
trtInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  select(-X)

#import diversity metrics that went into Bayesian analysis
rawData <- read.csv('ForAnalysis_allAnalysis10yr_pairwise_12142018.csv')

#calculate means and standard deviations across all data for richness and compositonal differences to backtransform
rawData2<- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_rich=mean(S_PC), std_rich=sd(S_PC)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length))

#for table of experiment summarizing various factors
expInfoSummary <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
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
trtInfo <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, trt_type)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich),
            anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(resource_mani=(nutrients+carbon+irrigation+drought), id=1:length(treatment))

#list of all treatments
studyInfo <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  select(site_code, community_type, project_name, treatment)%>%
  unique()


################################################################################
################################################################################
###Bayesian output processing

#only run to generate initial chains files
#raw chains data --------------------------------------------------------
memory.limit(size=50000)
chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10yr models_12172018\\noninf_lnRR_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10yr models_12172018\\noninf_lnRR_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10yr models_12172018\\noninf_lnRR_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10yr models_12172018\\noninf_lnRR_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)


#density plot of chains --------------------------------------------------------
plot(density(chainsCommunity$D.1.1.1))
plot(density(chainsCommunity$D.1.1.2))
plot(density(chainsCommunity$D.1.1.3))


#get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
#mean change are the 1's, richness are the 2's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__,
         # #trt_type intercepts: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.1, U.2.1.1, U.1.2.1, U.2.2.1, U.1.3.1, U.2.3.1, U.1.4.1, U.2.4.1, U.1.5.1, U.2.5.1,
         # U.1.6.1, U.2.6.1, U.1.7.1, U.2.7.1, U.1.8.1, U.2.8.1, U.1.9.1, U.2.9.1, U.1.10.1, U.2.10.1,
         # U.1.11.1, U.2.11.1, U.1.12.1, U.2.12.1, U.1.13.1, U.2.13.1, U.1.14.1, U.2.14.1, U.1.15.1, U.2.15.1,
         # U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1, U.1.19.1, U.2.19.1, U.1.20.1, U.2.20.1,
         # U.1.21.1, U.2.21.1, U.1.22.1, U.2.22.1, U.1.23.1, U.2.23.1, U.1.24.1, U.2.24.1, U.1.25.1, U.2.25.1,
         # U.1.26.1, U.2.26.1, U.1.27.1, U.2.27.1, U.1.28.1, U.2.28.1, U.1.29.1, U.2.29.1, U.1.30.1, U.2.30.1,
         # U.1.31.1, U.2.31.1, U.1.32.1, U.2.32.1, U.1.33.1, U.2.33.1, U.1.34.1, U.2.34.1, U.1.35.1, U.2.35.1,
         # U.1.36.1, U.2.36.1, U.1.37.1, U.2.37.1, U.1.38.1, U.2.38.1, U.1.39.1, U.2.39.1, U.1.40.1, U.2.40.1,
         # U.1.41.1, U.2.41.1, U.1.42.1, U.2.42.1, U.1.43.1, U.2.43.1, U.1.44.1, U.2.44.1, U.1.45.1, U.2.45.1,
         # U.1.46.1, U.2.46.1, U.1.47.1, U.2.47.1, U.1.48.1, U.2.48.1, U.1.49.1, U.2.49.1, U.1.50.1, U.2.50.1,
         # U.1.51.1, U.2.51.1, U.1.52.1, U.2.52.1, U.1.53.1, U.2.53.1, U.1.54.1, U.2.54.1,
         # #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
         # U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
         # U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
         # U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2, U.1.19.2, U.2.19.2, U.1.20.2, U.2.20.2,
         # U.1.21.2, U.2.21.2, U.1.22.2, U.2.22.2, U.1.23.2, U.2.23.2, U.1.24.2, U.2.24.2, U.1.25.2, U.2.25.2,
         # U.1.26.2, U.2.26.2, U.1.27.2, U.2.27.2, U.1.28.2, U.2.28.2, U.1.29.2, U.2.29.2, U.1.30.2, U.2.30.2,
         # U.1.31.2, U.2.31.2, U.1.32.2, U.2.32.2, U.1.33.2, U.2.33.2, U.1.34.2, U.2.34.2, U.1.35.2, U.2.35.2,
         # U.1.36.2, U.2.36.2, U.1.37.2, U.2.37.2, U.1.38.2, U.2.38.2, U.1.39.2, U.2.39.2, U.1.40.2, U.2.40.2,
         # U.1.41.2, U.2.41.2, U.1.42.2, U.2.42.2, U.1.43.2, U.2.43.2, U.1.44.2, U.2.44.2, U.1.45.2, U.2.45.2,
         # U.1.46.2, U.2.46.2, U.1.47.2, U.2.47.2, U.1.48.2, U.2.48.2, U.1.49.2, U.2.49.2, U.1.50.2, U.2.50.2,
         # U.1.51.2, U.2.51.2, U.1.52.2, U.2.52.2, U.1.53.2, U.2.53.2, U.1.54.2, U.2.54.2,
         # #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
         # U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
         # U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
         # U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3, U.1.19.3, U.2.19.3, U.1.20.3, U.2.20.3,
         # U.1.21.3, U.2.21.3, U.1.22.3, U.2.22.3, U.1.23.3, U.2.23.3, U.1.24.3, U.2.24.3, U.1.25.3, U.2.25.3,
         # U.1.26.3, U.2.26.3, U.1.27.3, U.2.27.3, U.1.28.3, U.2.28.3, U.1.29.3, U.2.29.3, U.1.30.3, U.2.30.3,
         # U.1.31.3, U.2.31.3, U.1.32.3, U.2.32.3, U.1.33.3, U.2.33.3, U.1.34.3, U.2.34.3, U.1.35.3, U.2.35.3,
         # U.1.36.3, U.2.36.3, U.1.37.3, U.2.37.3, U.1.38.3, U.2.38.3, U.1.39.3, U.2.39.3, U.1.40.3, U.2.40.3,
         # U.1.41.3, U.2.41.3, U.1.42.3, U.2.42.3, U.1.43.3, U.2.43.3, U.1.44.3, U.2.44.3, U.1.45.3, U.2.45.3,
         # U.1.46.3, U.2.46.3, U.1.47.3, U.2.47.3, U.1.48.3, U.2.48.3, U.1.49.3, U.2.49.3, U.1.50.3, U.2.50.3,
         # U.1.51.3, U.2.51.3, U.1.52.3, U.2.52.3, U.1.53.3, U.2.53.3, U.1.54.3, U.2.54.3,
         #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
         D.1.2.1, D.2.2.1,
         D.1.2.2, D.2.2.2,
         D.1.2.3, D.2.2.3,
         #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
         D.1.3.1, D.2.3.1,
         D.1.3.2, D.2.3.2,
         D.1.3.3, D.2.3.3,
         #overall intercept, linear, and quad slopes (center digit): 1=overall
         D.1.1.1, D.2.1.1,
         D.1.1.2, D.2.1.2,
         D.1.1.3, D.2.1.3)%>%
  gather(key=parameter, value=value, D.1.2.1:D.2.1.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))

# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteractions_10yr_noninf_12192018.csv')

chainsCommunity2 <- read.csv('ayesian_output_summary_expinteractions_10yr_noninf_12192018.csv')

#gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
#numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
#variable (second place): 1=mean change, 2=richness change
#parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
#set any that are not significant (CI overlaps 0) as 0

#get mean parameter values across all runs for each experiment, treatment, etc
chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,7482:9935]))%>% #may need to delete original four chains dataframes to get this to work
  add_rownames('parameter')
names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 7482:9935])'] <- 'mean'
#get sd of parameter values across all runs for each experiment, treatment, etc
chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,7482:9935]))
names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 7482:9935])'] <- 'sd'

chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
  #split names into parts
  separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
  select(-B)%>%
  #rename parts to be more clear
  mutate(variable=ifelse(variable==1, 'mean', 'richness'),
         parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
         id=as.integer(id))%>%
  #if 95% confidence interval overlaps 0, then set mean to 0
  mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
  #spread by variable
  select(variable, id, parameter, mean)%>%
  spread(key=parameter, value=mean)

# write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_10yr_noninf_12192018.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_10yr_noninf_12192018.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

#generate equations for main figure of richness and compositional responses through time
chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=20, 19, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1527948)+(0.3026725), (intercept+linear*7+quadratic*7^2)*(0.2331644)+(-0.04520116)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1527948)+(0.3026725),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2331644)+(-0.04520116)))%>%
  mutate(color=ifelse(rrich<31, '#1104DC44', ifelse(rrich<51&rrich>30, '#4403AE55', ifelse(rrich<71&rrich>50, '#77038166', ifelse(rrich>70, '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1527948)+(0.3026725)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2331644)+(-0.04520116)}, size=2, xlim=c(0,'),
         curve5='), colour=',
         curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_20yr_11302018_noninf.csv', row.names=F)

#summary lines
chainsEquationsSummary <- chainsEquations%>%
  group_by(variable, trt_overall)%>%
  summarise(intercept_mean=mean(intercept), intercept_sd=sd(intercept), linear_mean=mean(linear), linear_sd=sd(linear), quadratic_mean=mean(quadratic), quadratic_sd=sd(quadratic))%>%
  ungroup()%>%
  mutate(intercept_high=intercept_mean+1.96*intercept_sd, intercept_low=intercept_mean-1.96*intercept_sd, linear_high=linear_mean+1.96*linear_sd, linear_low=linear_mean-1.96*linear_sd, quadratic_high=quadratic_mean+1.96*quadratic_sd, quadratic_low=quadratic_mean-1.96*quadratic_sd)%>%
  mutate(intercept_high2=ifelse(intercept_high>0, 1, 0), intercept_low2=ifelse(intercept_low<0, 1, 0), linear_high2=ifelse(linear_high>0, 1, 0), linear_low2=ifelse(linear_low<0, 1, 0), quadratic_high2=ifelse(quadratic_high>0, 1, 0), quadratic_low2=ifelse(quadratic_low<0, 1, 0))%>%
  mutate(intercept_cross=intercept_high2+intercept_low2, linear_cross=linear_high2+linear_low2, quadratic_cross=quadratic_high2+quadratic_low2)%>%
  mutate(intercept_final=ifelse(intercept_cross==2, 0, intercept), linear_final=ifelse(linear_cross==2, 0, linear), quadratic_final=ifelse(quadratic_cross==2, 0, quadratic))

###main figure (Figure 1)
# compositional response panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,9), breaks=seq(4,9,5), labels=seq(5,10,5)) +
  ylim(-10,10) +
xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines

  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.4346700 + 0.2557790*x + -0.0166317*x^2)*(0.1527948)+(0.3026725)}, size=5, xlim=c(0,19), colour='black')

# print(meanPlot) #export at 1200x1000



#richness response panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Richness Difference') +
  annotate('text', x=0, y=2.0, label='(a)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.38809837413322*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.549865645952*x + 0.04155431570449*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6797793152601 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36437437647475*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.7665116302726 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.5829194632377 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6937042106611 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2929943482313*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.12228230946005 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.851838790537 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7854440563095 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.917212876386 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.7073021996115 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26995861021*x + -0.02290064001761*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.60273060558085 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.8898756921111 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.260304269470896*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0.69796085051145 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.74180954418 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.9481443055035 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.9571713994465 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.7450114646917 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.52612479178705 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.7580112992 + 0.365141400253755*x + -0.04441614850448*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.02686448527 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.65749743265*x + 0.0509804921357*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.9389478576385 + -0.4645262989799*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.632252538675*x + 0.0508405379941*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.184652007662 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.50861047764355*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.648483685759795 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.54649226143295 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.6802061809763 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.90800883790375 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.312898902338887*x + 0.040109741399055*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.84616617712 + 0.24704903004215*x + -0.0273150110606*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.61324013882109 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26027937424825*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8022844560655 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.716593936017 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.511290830274695 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52009887474995*x + -0.06468771994736*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.412342424164215*x + -0.05612744537087*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3822651674508*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.1835035735 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.57547128714925 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.880584483162 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8548644983225 + 0.522302252995*x + -0.05383626433115*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.8778739585024 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.63830493912725 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.8075485714415 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.7601291066727 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5655918889325 + 0.3715093547343*x + -0.043411200137875*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7107029681655 + 0.377613232431*x + -0.039382787140074*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.61801542558885 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.70187147401435 + 0*x + 0.054447526770485*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.164603634 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.85574838786906 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.702657371845 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.5898678342303 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.9722089927835 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.568977569272945 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3844686079597*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.278247670895*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.54606502296996 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17705417525154*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.59449760674877 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.687631404636295 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.3555563211 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.240245467102 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.24426893265 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.29923271915 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.96796071262 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.86070898424275 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.71750157264875 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0167479959699725*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.61866131308145 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0161108543478081*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.45640945514*x + 0.0329172347135*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.208257193364532*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8023619746556 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.73130550131755 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.372916984451*x + 0.025774441626798*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3357364713795*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4407078654795*x + -0.0431912439939615*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4011735277753*x + -0.047569938912425*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.468070830991*x + -0.04705722184215*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52102630291*x + -0.05313916178065*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.10899622371 + 0.54292919182*x + -0.068551284375*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.7249241395103 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.87442188269535 + -0.50515541592393*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -1.0474519769627*x + 0.69393793125*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.29645892275 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.0796769840882 + 0*x + 0.66951665795*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.8198706866345 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.5493467930334 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.63950911699275 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.868012831363 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.6152130319396 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.6030455180305 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.89304101384935 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.65385323725 + 0.7825907042402*x + -0.0783151407005*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.70174102150455*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.55315499693454 + -0.33188621933*x + 0.017441332685*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5769346181*x + 0.0316451064*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.674939381848 + -0.55670096705*x + 0.032300423165*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.01757441605 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.93670374874285 + -0.48567069585*x + 0.0174565922675*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.387327975*x + 0.01347276128305*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.40440615465*x + 0.0195696882575*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.584934419864005 + -0.40815973515*x + 0.0196372802875*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.140469259864091*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.86629442345*x + 0.06610571839*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.7067941104675 + -0.8986110425*x + 0.07248182809*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.7394203402*x + 0.06188394378*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.82655624085*x + 0.05416132770905*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4208851340895*x + 0.028474236609865*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.59994385518635 + -0.583366453775*x + 0.04226565804765*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.365765788542555*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3469637705682*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.30256507572641*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0216078253124015*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.251148943893405*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3134956139725*x + 0.0209631028047902*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.50218825483*x + 0.032120278347*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.558229314820735 + 0*x + -0.035320708732155*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.9095460015*x + -0.0860449631766*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.2462924078 + 0.86132571885*x + -0.12069858703*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.041074775495 + 0.3968408563853*x + -0.05708563290022*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7979759551624 + 0.73701668575*x + -0.10780230049*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31244117882025*x + -0.0463892129097*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.264307503733153*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6135223897838 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.024615088735 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.1007688899435*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.26837650 + 0*x + 0*x^2)*0.2407859 + -0.06393594}, size=5, xlim=c(0,19), colour='black')

# print(richnessPlot) #export at 1200x1000

#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,1)))
print(richnessPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1200 x 2400



#summary stats from bayesian output --------------------------------------------------------
#gather summary stats needed and relabel them
chainsCommunitySummary <- chainsCommunity%>%
  select(
         # #trt_type intercepts: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.1, U.2.1.1, U.1.2.1, U.2.2.1, U.1.3.1, U.2.3.1, U.1.4.1, U.2.4.1, U.1.5.1, U.2.5.1,
         # U.1.6.1, U.2.6.1, U.1.7.1, U.2.7.1, U.1.8.1, U.2.8.1, U.1.9.1, U.2.9.1, U.1.10.1, U.2.10.1,
         # U.1.11.1, U.2.11.1, U.1.12.1, U.2.12.1, U.1.13.1, U.2.13.1, U.1.14.1, U.2.14.1, U.1.15.1, U.2.15.1,
         # U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1, U.1.19.1, U.2.19.1, U.1.20.1, U.2.20.1,
         # U.1.21.1, U.2.21.1, U.1.22.1, U.2.22.1, U.1.23.1, U.2.23.1, U.1.24.1, U.2.24.1, U.1.25.1, U.2.25.1,
         # U.1.26.1, U.2.26.1, U.1.27.1, U.2.27.1, U.1.28.1, U.2.28.1, U.1.29.1, U.2.29.1, U.1.30.1, U.2.30.1,
         # U.1.31.1, U.2.31.1, U.1.32.1, U.2.32.1, U.1.33.1, U.2.33.1, U.1.34.1, U.2.34.1, U.1.35.1, U.2.35.1,
         # U.1.36.1, U.2.36.1, U.1.37.1, U.2.37.1, U.1.38.1, U.2.38.1, U.1.39.1, U.2.39.1, U.1.40.1, U.2.40.1,
         # U.1.41.1, U.2.41.1, U.1.42.1, U.2.42.1, U.1.43.1, U.2.43.1, U.1.44.1, U.2.44.1, U.1.45.1, U.2.45.1,
         # U.1.46.1, U.2.46.1, U.1.47.1, U.2.47.1, U.1.48.1, U.2.48.1, U.1.49.1, U.2.49.1, U.1.50.1, U.2.50.1,
         # U.1.51.1, U.2.51.1, U.1.52.1, U.2.52.1, U.1.53.1, U.2.53.1, U.1.54.1, U.2.54.1,
         # #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
         # U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
         # U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
         # U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2, U.1.19.2, U.2.19.2, U.1.20.2, U.2.20.2,
         # U.1.21.2, U.2.21.2, U.1.22.2, U.2.22.2, U.1.23.2, U.2.23.2, U.1.24.2, U.2.24.2, U.1.25.2, U.2.25.2,
         # U.1.26.2, U.2.26.2, U.1.27.2, U.2.27.2, U.1.28.2, U.2.28.2, U.1.29.2, U.2.29.2, U.1.30.2, U.2.30.2,
         # U.1.31.2, U.2.31.2, U.1.32.2, U.2.32.2, U.1.33.2, U.2.33.2, U.1.34.2, U.2.34.2, U.1.35.2, U.2.35.2,
         # U.1.36.2, U.2.36.2, U.1.37.2, U.2.37.2, U.1.38.2, U.2.38.2, U.1.39.2, U.2.39.2, U.1.40.2, U.2.40.2,
         # U.1.41.2, U.2.41.2, U.1.42.2, U.2.42.2, U.1.43.2, U.2.43.2, U.1.44.2, U.2.44.2, U.1.45.2, U.2.45.2,
         # U.1.46.2, U.2.46.2, U.1.47.2, U.2.47.2, U.1.48.2, U.2.48.2, U.1.49.2, U.2.49.2, U.1.50.2, U.2.50.2,
         # U.1.51.2, U.2.51.2, U.1.52.2, U.2.52.2, U.1.53.2, U.2.53.2, U.1.54.2, U.2.54.2,
         # #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
         # U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
         # U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
         # U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
         # U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3, U.1.19.3, U.2.19.3, U.1.20.3, U.2.20.3,
         # U.1.21.3, U.2.21.3, U.1.22.3, U.2.22.3, U.1.23.3, U.2.23.3, U.1.24.3, U.2.24.3, U.1.25.3, U.2.25.3,
         # U.1.26.3, U.2.26.3, U.1.27.3, U.2.27.3, U.1.28.3, U.2.28.3, U.1.29.3, U.2.29.3, U.1.30.3, U.2.30.3,
         # U.1.31.3, U.2.31.3, U.1.32.3, U.2.32.3, U.1.33.3, U.2.33.3, U.1.34.3, U.2.34.3, U.1.35.3, U.2.35.3,
         # U.1.36.3, U.2.36.3, U.1.37.3, U.2.37.3, U.1.38.3, U.2.38.3, U.1.39.3, U.2.39.3, U.1.40.3, U.2.40.3,
         # U.1.41.3, U.2.41.3, U.1.42.3, U.2.42.3, U.1.43.3, U.2.43.3, U.1.44.3, U.2.44.3, U.1.45.3, U.2.45.3,
         # U.1.46.3, U.2.46.3, U.1.47.3, U.2.47.3, U.1.48.3, U.2.48.3, U.1.49.3, U.2.49.3, U.1.50.3, U.2.50.3,
         # U.1.51.3, U.2.51.3, U.1.52.3, U.2.52.3, U.1.53.3, U.2.53.3, U.1.54.3, U.2.54.3,
         #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
         D.1.2.1, D.2.2.1,
         D.1.2.2, D.2.2.2,
         D.1.2.3, D.2.2.3,
         #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
         D.1.3.1, D.2.3.1,
         D.1.3.2, D.2.3.2,
         D.1.3.3, D.2.3.3,
         #overall intercept, linear, and quad slopes (center digit): 1=overall
         D.1.1.1, D.2.1.1,
         D.1.1.2, D.2.1.2,
         D.1.1.3, D.2.1.3)

chainsCommunitySummary <- chainsCommunitySummary%>%
  gather(key=parameter, value=value, D.1.2.1:D.2.1.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  ungroup()%>%
  mutate(CI=sd*2)%>%
  separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
  #rename parts to be more clear
  mutate(variable=ifelse(variable==1, 'mean', 'richness'),
         parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
         predictor2=ifelse(predictor==2, 'ANPP', ifelse(predictor==3, 'rrich', 'overall')))%>%
  select(variable, parameter, predictor2, median, sd, CI)

# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_20yr_11302018_noninf.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_20yr_11302018_noninf.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  mutate(type=paste(predictor2, parameter, sep='_'))



###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear', 'rrich_intercept', 'ANPP_intercept', 'overall_intercept'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  geom_vline(aes(xintercept=6.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Community Difference') +
  annotate('text', x=9.2, y=-0.8, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='richness' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear', 'rrich_intercept', 'ANPP_intercept', 'overall_intercept'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  geom_vline(aes(xintercept=6.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Richness Difference') +
  annotate('text', x=9.2, y=-0.8, label='(a)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 1600x1000




###by magnitude of resource manipulated---------------------------------
#N addition
nData <- read.csv('ForAnalysis_allAnalysisNmag.csv')

nDataMean <- nData%>%
  summarise(mean_mean_change=mean(mean_change), sd_mean_change=sd(mean_change), mean_S_PC=mean(S_PC), sd_S_PC=sd(S_PC))

#mean change
Nmean <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\posteriors_N_MeanChange.csv', comment.char='#')
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
#   summarise(mean_change_mean=mean(mean_change), mean_change_sd=sd(mean_change), richness_mean=mean(S_PC), richness_sd=sd(S_PC), n_mean=mean(n), n_sd=sd(n), MAP_mean=mean(MAP), MAP_sd=sd(MAP))

nDataTransform <- nData%>%
  #transform mean change
  mutate(mean_change_transform=((mean_change-mean(mean_change))/sd(mean_change)))%>%
  #transform proportional richness change
  mutate(S_PC_transform=((S_PC-mean(S_PC))/sd(S_PC)))%>%
  #transform N treatment magnitude
  mutate(n_transform=((n-mean(n))/sd(n)))%>%
  #transform MAP
  mutate(MAP_transform=((MAP-mean(MAP))/sd(MAP)))

meanNPlotFinal <- ggplot(data=subset(nData), aes(x=n, y=mean_change, color=MAP)) +
  geom_point(size=5) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name='Compositional Response') +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((1000-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((1000-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#4793CF')  +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((600-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((600-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#2D5E88') +
  stat_function(fun=function(x){(0.02512656 + 0.40341207*((200-661.9362)/298.3696) + 0.54133077*(x-9.992142)/9.108662 + 0.28058497*((200-661.9362)/298.3696)*(x-9.992142)/9.108662)*0.1658319+0.3699378}, size=5, color='#153049') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.0, label='(d)', size=12, hjust='left') +
  theme(legend.position=c(0.8,0.05), legend.justification=c(0,0), legend.title=element_text(size=24))


#richness difference
Nrichness <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\posteriors_N_Richness.csv', comment.char='#')
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

meanDroPlotFinal <- ggplot(data=droData, aes(x=precip, y=mean_change)) +
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

meanIrrPlotFinal <- ggplot(data=irrData, aes(x=precip, y=mean_change)) +
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