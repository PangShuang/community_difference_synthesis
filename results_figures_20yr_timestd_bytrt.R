################################################################################
##  results_figures_20yr.R: Compiles Bayesian output and makes figures for the primary analysis of richness and compositonal differences between treatment and control plots for datasets cut off at 20 years.
##
##  Author: Kimberly La Pierre
##  Date created: January 17, 2018
##  See https://github.com/klapierre/Converge_Diverge/blob/master/core%20data%20paper_bayesian%20results_figures_sig%20test_expinteractions_20yr.R for full history.
################################################################################

library(grid)
library(tidyverse)

#kim laptop
setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

#kim desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

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
rawData <- read.csv('ForAnalysis_allAnalysis20yr_pairwise_04032019.csv')

#calculate means and standard deviations across all data for richness and compositonal differences to backtransform
rawData2<- rawData%>%
  left_join(trtInfo1)%>%
  filter(anpp!='NA', treatment_year!=0)%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_rich=mean(S_lnRR), std_rich=sd(S_lnRR)) #to backtransform

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
################################################################################
# ###Bayesian output processing
# 
# only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\noninf_lnRR\\noninf_timestdbytrt_20yr_lnRR_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\noninf_lnRR\\noninf_timestdbytrt_20yr_lnRR_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\noninf_lnRR\\noninf_timestdbytrt_20yr_lnRR_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\noninf_lnRR\\noninf_timestdbytrt_20yr_lnRR_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$D.1.1.1))
# plot(density(chainsCommunity$D.1.1.2))
# plot(density(chainsCommunity$D.1.1.3))
# 
# 
# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# #mean change are the 1's, richness are the 2's
# chainsCommunity2 <- chainsCommunity%>%
#   select(lp__,
#          #trt_type intercepts: center digit refers to trts
#          E.1.1.1, E.2.1.1, E.1.2.1, E.2.2.1, E.1.3.1, E.2.3.1, E.1.4.1, E.2.4.1, E.1.5.1, E.2.5.1,
#          E.1.6.1, E.2.6.1, E.1.7.1, E.2.7.1, E.1.8.1, E.2.8.1, E.1.9.1, E.2.9.1, E.1.10.1, E.2.10.1,
#          E.1.11.1, E.2.11.1, E.1.12.1, E.2.12.1, E.1.13.1, E.2.13.1, E.1.14.1, E.2.14.1, E.1.15.1, E.2.15.1,
#          E.1.16.1, E.2.16.1,
#          #trt_type linear slopes: center digit refers to trts
#          E.1.1.2, E.2.1.2, E.1.2.2, E.2.2.2, E.1.3.2, E.2.3.2, E.1.4.2, E.2.4.2, E.1.5.2, E.2.5.2,
#          E.1.6.2, E.2.6.2, E.1.7.2, E.2.7.2, E.1.8.2, E.2.8.2, E.1.9.2, E.2.9.2, E.1.10.2, E.2.10.2,
#          E.1.11.2, E.2.11.2, E.1.12.2, E.2.12.2, E.1.13.2, E.2.13.2, E.1.14.2, E.2.14.2, E.1.15.2, E.2.15.2,
#          E.1.16.2, E.2.16.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          E.1.1.3, E.2.1.3, E.1.2.3, E.2.2.3, E.1.3.3, E.2.3.3, E.1.4.3, E.2.4.3, E.1.5.3, E.2.5.3,
#          E.1.6.3, E.2.6.3, E.1.7.3, E.2.7.3, E.1.8.3, E.2.8.3, E.1.9.3, E.2.9.3, E.1.10.3, E.2.10.3,
#          E.1.11.3, E.2.11.3, E.1.12.3, E.2.12.3, E.1.13.3, E.2.13.3, E.1.14.3, E.2.14.3, E.1.15.3, E.2.15.3,
#          E.1.16.3, E.2.16.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
#          D.1.2.1, D.2.2.1,
#          D.1.2.2, D.2.2.2,
#          D.1.2.3, D.2.2.3,
#          #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
#          D.1.3.1, D.2.3.1,
#          D.1.3.2, D.2.3.2,
#          D.1.3.3, D.2.3.3,
#          #overall intercept, linear, and quad slopes (center digit): 1=overall
#          D.1.1.1, D.2.1.1,
#          D.1.1.2, D.2.1.2,
#          D.1.1.3, D.2.1.3)%>%
#   gather(key=parameter, value=value, E.1.1.1:D.2.1.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# # write.csv(chainsCommunity2, 'bayesian_output_summary_expinteraction_20yr_stdtimebytrt_04072019.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteraction_20yr_stdtimebytrt_04072019.csv')

#gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
#numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
#variable (second place): 1=mean change, 2=richness change
#parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
#set any that are not significant (CI overlaps 0) as 0

# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,8696:11323]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 8696:11323])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,8696:11323]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 8696:11323])'] <- 'sd'
# 
# chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
#   #spread by variable
#   select(variable, id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# # write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_04072019_noninf.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_04072019_noninf.csv')


#merge together with experiment list
trtID <- read.csv('bayesian_trt_index.csv')%>%
  select(site_code, project_name, community_type, treatment, treat_INT)%>%
  unique()%>%
  rename(id=treat_INT)
timeStd <- read.csv('bayesian_trt_index.csv')%>%
  group_by(comm_INT, treat_INT)%>%
  summarise(time_mean=mean(time), time_std=sd(time))%>%
  ungroup()%>%
  rename(id=treat_INT)
chainsExperiment <- chainsFinal%>%
  left_join(trtID)%>%
  left_join(trtInfo)%>%
  left_join(timeStd)

#generate equations for main figure of richness and compositional responses through time
chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=20, 19, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1725573)+(0.3215148), (intercept+linear*7+quadratic*7^2)*(0.3294397)+(-0.1075646)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1725573)+(0.3215148),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.3294397)+(-0.1075646)))%>%
  mutate(color=ifelse(rrich<31, '#1104DC44', ifelse(rrich<51&rrich>30, '#4403AE55', ifelse(rrich<71&rrich>50, '#77038166', ifelse(rrich>70, '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*((x-',
         curve4=')/',
         curve5=') + ',
         curve6='*((x-',
         curve7=')/',
         curve8=ifelse(variable=='mean', ')^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,',
                       ')^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,'),
         curve9=')) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, time_mean, curve4, time_std, curve5, quadratic, curve6, time_mean, curve7, time_std, curve8, alt_length, curve9, sep=''))
  # mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))

#export, group by shape type, and paste lines below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_20yr_stdtimebytrt_04072019.csv', row.names=F)


trtShape <- read.csv('treatment_response_shape_classification_stdtimebytrt_04072019.csv')

mean(trtShape$experiment_length) #8.075342
median(trtShape$experiment_length) #6

#plot of shape type by length
ggplot(barGraphStats(data=trtShape, variable="experiment_length", byFactorNames=c("shape_category","variable")), aes(x=shape_category, y=mean)) +
  geom_bar(stat='identity') +
  facet_grid(~variable)


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
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6401008315651 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.4709532317155 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.54089440738607 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.3691754406512 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(1.14471942605 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.14935338765 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.59218103903 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.6670398764655 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.74517583654405 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.51024117363095 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.481013538228 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.86695832445 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.6454678548 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-1.2006242818 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.25939339475 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.23877545645 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.9472963638 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.4950797928252 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.544712450949425 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.50822456845342 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.642342220266635 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.6636239351495 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.474060350137 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(-0.28186474175064 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(-1.1457828895 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.47439841149425 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.65130379063355 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.686425937831 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.6046795152315 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.86405550363 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-1.49985274445 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.4307299435 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.43275046585 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.3489631471 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.47789658645 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.3510168526 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.11524441245 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-1.19736416625 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.68542849555255 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.62992865215155 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.742271309868 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.52384602882725 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8030773985505 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8500707111182 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5469972727325 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5357380495301 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.55377901050675 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.699244639513 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.8051895177 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.389713770587245 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.322938759481455 + 0*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.48634290153365 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.736343040663 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.87421853279 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.7556074724866 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-1.1089358432 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.0094152744 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.2294250567 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.076016852375 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.061843125 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.13205751525 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.934205493975 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.09618201375 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.472677304087 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.85499895505 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.755929450413 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.541923989273 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.843333315945 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.947167469930615 + 0*((x-7.75)/6.84957419601151) + 0*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.592063782538 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.849215014145 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5423125806732 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.61184075564005 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5602856434643 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5094238544233 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.77288450305685 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6601455386855 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.90049744484 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.05087181815 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6555821337 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.80516677045 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.69189432275 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-1.21200712045 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.8462310307 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.9865182035 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.640313556396 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.7892365990405 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.609355438662 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.64460857400375 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.96520561545 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.02877071425 + 0*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.559047531857975 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.45405023881115 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.5408435498817 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.425595296066475 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.4466422569062 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.60018395094 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.7172894933855 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.4139737015519 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.4005145893168 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.742643770155 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.78401590256 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.75117609831 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.73932803658 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.01165978295 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.00132111095 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.76316509182 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.6961553882025 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.0499546583 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/3.1091263510296) + 0*((x-2.5)/3.1091263510296)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-1.36029149515 + 0*((x-3.14285714285714)/2.41029537806548) + 0*((x-3.14285714285714)/2.41029537806548)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.54167496260135 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.476553870787 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.6061341863038 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-1.5071642933 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.31483237225 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.90624733605 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(-0.56625444731615 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.51065151312525 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.524516709192 + 0*((x-4)/2.73861278752583) + 0*((x-4)/2.73861278752583)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.56431890010335 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.567682489938 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.54126033346695 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.8663086466865 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.63329920765335 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.817628941785 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.91451470215 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.647551194977 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.4969145202982 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.76652581053 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.64339563799145 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.5434659686049 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.52698049275652 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) 



#------------------------
meanPlot1 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(b) 63.0%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.648460224905 + 0.466609335325*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.73537874155 + 0.6085903779*((x-5.63636363636364)/3.93122696553448) + 0*((x-5.63636363636364)/3.93122696553448)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.4485647350565 + 0.7643202309*((x-4.55555555555556)/3.39525813124389) + 0*((x-4.55555555555556)/3.39525813124389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(1.029155089385 + 0.65550679989*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.67500792955 + 0.6419503636*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.896664874 + 0.69785951025*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.87926623346 + 0.432226902973*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(1.12990685415 + 0.39773178986667*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.9394478342 + 0.4060182535675*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.42386840635945*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.706636425547 + 0.44401527942585*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.895724100085 + 0.4471879522911*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.9496811273 + 0.34937672950285*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-1.16380539685 + 0.33856363415695*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.5218400209935*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.3933319282128 + 0.15916001006591*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.4499046703288 + 0.50714120245*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.67742118305 + 0.49577817065*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.7072990372148 + 0.33809600252165*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.4711223846705*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.2884047032719*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.2734400970362*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.67511187235*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0.6009743335345*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.705254464275*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.37198017885925*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.42560587076766*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.3974763085367*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.228103702536315*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0.277422997930855*((x-3.5)/2.44948974278318) + 0*((x-3.5)/2.44948974278318)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.928743556095 + 0.326213787892015*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.24155125729625*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.60125587186245 + 0.31380862585691*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.337942027481535*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.34142469611215*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(2.568941679 + 0.17680566297945*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0.2803716516834*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.38969752672*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(-0.643095562655 + 0.19674427619475*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(-0.70653887585 + 0.227510433596885*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.414177368934735 + 0.4748088536245*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.50234271881995 + 0.556604542155*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.9240168 + 0.476112352*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(1.40496754115 + 0.7400813754*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.4767489929415 + 0.401911352345*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.248371058767849*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.72502102075 + 0.29153577291595*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(1.06825644295 + 0.5294026741*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.6718394385134 + 0.73333577104*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.46036571372795*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.60485726794075 + 0.430249763322755*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0.6763052094*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.8360653928*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.90895637595*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.87416079995*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.684079106025 + 0.3935584868145*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(-0.704294368315 + 0.30197021927885*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.569639643427605 + 0.93262101915*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.93822871595*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.98244815685 + 0.75090990165*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.8603454414*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.37517404670615 + 0.8031353379*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.91801384975*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.374572134941494 + 0.80746168145*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.531351889765*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0.40592845578535*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.38862270846275*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.36231963709115*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.514569075692*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.9585038089 + 0.34166578984935*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.541755086495*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.6834788965625 + 0.36074731774405*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.38426427269705*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.181327652 + 0.591989395045*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.676334048663 + 0.555972115076*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.6917406479357 + 0.49208590756*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.6618047464185 + 0.4956726193905*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.9448282535 + 0.50610364586*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.2833880132 + 0.610198020365*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.56822849906165 + 0.26482193352851*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.35119867118685 + 0.5434063679*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.40600681013116 + 0.45936887173*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-1.1161897188 + 0.26755970548834*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.630122764947 + 0.49902246252*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.6114769839685 + 0.294499907653*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.1040156104 + 0.4818117421615*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.55900464926085 + 1.0450543834*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(-0.5539495581988 + 0.477531073671*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.44009576534335 + 0.542286036005*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.630211600678 + 0.31800687202845*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.5493514534052 + 0.41511963461515*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.6213942110515 + 0.43498566311285*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0.5056801409*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.52369465349795 + 0.380725642828055*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(-0.46037215443201 + 0.40532913347485*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0.590896647655*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0.5672213702685*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0.341444883264445*((x-1)/1) + 0*((x-1)/1)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,2))
  
#------------------------
meanPlot2 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(f) 0.9%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.593661917435 + 0.2553582834994*((x-5)/3.3166247903554) + 0.24367603430132*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.37348114333365 + 0.66492015455*((x-5)/3.3166247903554) + 0.2066100937079*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.7540852598 + 0.63824753185*((x-5)/3.3166247903554) + 0.19784419664395*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.3561763521215*((x-5)/3.3166247903554) + 0.29172652755585*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) 
  
#------------------------
meanPlot3 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(h) 3.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.6911753343455 + 0*((x-2.33333333333333)/2.51661147842358) + -0.283119187932195*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.9343722095 + 0.55711504585*((x-9.5)/5.91607978309962) + -0.2581652549465*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(2.5071341835 + 0.75241550155*((x-9.5)/5.91607978309962) + -0.2712518293587*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(3.1167620145 + 1.28924450285*((x-9.5)/5.91607978309962) + -0.69208632615*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(2.509327407 + 1.00629284395*((x-4.5)/3.02765035409749) + -0.459278441845*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(3.1852865085 + 0.87421238435*((x-4.5)/3.02765035409749) + -0.59276925525*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.09499768735 + 0.7762416523*((x-4.5)/3.02765035409749) + -0.246188968841005*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(2.0641595805 + 0.95614083735*((x-4.5)/3.02765035409749) + -0.365101186069*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0.2324320355198*((x-6)/3.89444048184931) + -0.16722764969448*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.95349945275 + 0.61603318035*((x-4.5)/3.02765035409749) + -0.215524884471519*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.1790456512 + 0.56074669355*((x-4.5)/3.02765035409749) + -0.252644285540525*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.348275541123525 + 0.3036686604655*((x-5.5)/3.60555127546399) + -0.2463428108864*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.6824642105 + 0.42998601908*((x-5.5)/3.60555127546399) + -0.2536188667354*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.85776759785 + 0.6192598810655*((x-3.33333333333333)/3.05505046330389) + -0.333610313837625*((x-3.33333333333333)/3.05505046330389)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.74474589395 + 0.9619131735*((x-4.33333333333333)/3.278719262151) + -0.22598530394343*((x-4.33333333333333)/3.278719262151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(1.1057840104 + 0.625586158595*((x-2)/1.58113883008419) + -0.292368737867015*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4))
  
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
  stat_function(fun=function(x){(1.42399368825 + 0*((x-2.33333333333333)/2.51661147842358) + -0.3994702040472*((x-2.33333333333333)/2.51661147842358)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.41716202934925 + 0*((x-2)/1.58113883008419) + -0.2502824073133*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0.241361169910025*((x-6.30769230769231)/4.44193305113542) + -0.2104470107111*((x-6.30769230769231)/4.44193305113542)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + -0.236483860804*((x-6.30769230769231)/4.44193305113542)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0.4196127584335 + 0.3052211171406*((x-6.30769230769231)/4.44193305113542) + -0.340175218935*((x-6.30769230769231)/4.44193305113542)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-0.3811998793665 + 0*((x-6.30769230769231)/4.44193305113542) + -0.2171542896657*((x-6.30769230769231)/4.44193305113542)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-0.51122882649 + 0*((x-9.5)/5.91607978309962) + -0.2229047547752*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(2.733298538 + 0.99260452775*((x-9.5)/5.91607978309962) + -0.7513991392*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.886468257 + 0*((x-9.5)/5.91607978309962) + -0.321778843865*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.931675426 + 0.254464071887*((x-9.5)/5.91607978309962) + -0.325155154075*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.45785583995 + 0*((x-9.5)/5.91607978309962) + -0.262320649931*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.494252662 + 0*((x-9.5)/5.91607978309962) + -0.556604522*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(2.1916718575 + 0.2563068784935*((x-9.5)/5.91607978309962) + -0.5203419805*((x-9.5)/5.91607978309962)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + -0.217540206318055*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0.22250410273955*((x-6)/3.89444048184931) + -0.19410127786708*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.2025647310746*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + -0.192759097655555*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.28907880475315 + 0.201322486840305*((x-6)/3.89444048184931) + -0.2147877167529*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.4493614828235 + 0*((x-6)/3.89444048184931) + -0.30508495602*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.265763637126124 + 0.20067863224878*((x-6)/3.89444048184931) + -0.18502391885129*((x-6)/3.89444048184931)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(2.81264697 + 0.169313197729875*((x-9)/5.62731433871138) + -0.1997391874648*((x-9)/5.62731433871138)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.57525233229345 + 0.3223192109833*((x-2.66666666666667)/2.16024689946929) + -0.26149622558815*((x-2.66666666666667)/2.16024689946929)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(1.19890958825 + 0*((x-7.75)/6.84957419601151) + -0.57518905525*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.7217976979 + 0*((x-7.75)/6.84957419601151) + -0.65755744445*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.004336517745 + 0*((x-7.75)/6.84957419601151) + -0.545922164785*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.51538115765 + 0.3544353269439*((x-7.75)/6.84957419601151) + -0.62048209265*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.51253681515 + 0*((x-7.75)/6.84957419601151) + -0.6523646746*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.53165079635 + 0*((x-7.75)/6.84957419601151) + -0.62714663665*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.0801634156 + 0.29274236285193*((x-7.75)/6.84957419601151) + -0.473480664849*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.94487593055 + 0*((x-7.75)/6.84957419601151) + -0.59398211015*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.47368064175 + 0*((x-7.75)/6.84957419601151) + -0.79395326735*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.46159971185 + 0*((x-7.75)/6.84957419601151) + -0.633683998795*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(1.5552715191 + 0*((x-7.75)/6.84957419601151) + -0.7070844257*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.88426814121 + 0*((x-7.75)/6.84957419601151) + -0.491065428035*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.723621832526 + 0*((x-7.75)/6.84957419601151) + -0.430137598808*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.43081551046951 + 0*((x-7.75)/6.84957419601151) + -0.4270724637559*((x-7.75)/6.84957419601151)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.70626079822 + 0.25907806773245*((x-5.5)/3.60555127546399) + -0.347565306437*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.49287332307 + 0*((x-5.5)/3.60555127546399) + -0.21393782526305*((x-5.5)/3.60555127546399)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + -0.378221651254*((x-4.5)/3.02765035409749)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + -0.2925299288263*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + -0.279793431078316*((x-2)/1.58113883008419)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,4))
  
#------------------------
meanPlot8 <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=1, label='(r) 0.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.612594226855 + 0*((x-5)/3.3166247903554) + 0.256146976411*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.201856177582915*((x-5)/3.3166247903554) + 0.2531590546857*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0.3413042733755*((x-5)/3.3166247903554) + 0.367617449785*((x-5)/3.3166247903554)^2)*(0.1725573)+(0.3215148)}, size=2, xlim=c(0,10)) 
  
  
  


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
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.638629026395 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.6593714811254 + 0*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(1.2062405682 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.6)/2.40831891575846) + 0*((x-2.6)/2.40831891575846)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0.5270674399205 + 0*((x-6.30769230769231)/4.44193305113542) + 0*((x-6.30769230769231)/4.44193305113542)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.5983117137335 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.8122021034055 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0 + 0*((x-6.5)/4.18330013267038) + 0*((x-6.5)/4.18330013267038)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,13)) +
  stat_function(fun=function(x){(0.632679211855 + 0*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0.349959489226 + 0*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.65700748175825 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.54081606492855 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(1.0109321058 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.655072142106 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0.731972712426 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.59626868737095 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.8582671076338 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.7592761634865 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(1.1552530148265 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.986088985815 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.68556931230535 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(1.022799489115 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.97965927109185 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.803511432861 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.79344460459175 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.72159696350755 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.4814433628784 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.60521234486 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.3537450782117 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.37327263227335 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.377628686065944 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.5437516236842 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0.504423500693375 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.80648721753 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.7913683010295 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.648945935484336 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.66666666666667)/2.16024689946929) + 0*((x-2.66666666666667)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.9489028629335 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.8672464011055 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.9196388284695 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(0.404605727390785 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0.560041003595 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.76325318845 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.6254828221695 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.7026579292775 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-3.83333333333333)/3.31159578853861) + 0*((x-3.83333333333333)/3.31159578853861)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.5984250184945 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.78454079931975 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.6450010890139 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.820802279640465 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.710819521829125 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.62095066979663 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.636464377586266 + 0*((x-3.75)/2.81577190634672) + 0*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.650119418499025 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.655265056879635 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.728854690464225 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.33333333333333)/1.52752523165195) + 0*((x-1.33333333333333)/1.52752523165195)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.915302275897 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.200074881 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.316926175 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.98469447145 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.677649315057 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.5294942625926 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.473296830502247 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-4.25)/3.37003603202441) + 0*((x-4.25)/3.37003603202441)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.454664785767535 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.72328309763 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.5580855949988 + 0*((x-4.33333333333333)/3.278719262151) + 0*((x-4.33333333333333)/3.278719262151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.767692033 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.56490095048035 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.54574996832919 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0.5999509596228 + 0*((x-2.5)/3.1091263510296) + 0*((x-2.5)/3.1091263510296)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + 0*((x-3.14285714285714)/2.41029537806548) + 0*((x-3.14285714285714)/2.41029537806548)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0.63022365484735 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.6740871452329 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.6398150331664 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.80549033045465 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.68714984592705 + 0*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(1.201980756405 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0.88545121095 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0.5533615423922 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(0.5091896997848 + 0*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(1.14079278196 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-4)/2.73861278752583) + 0*((x-4)/2.73861278752583)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.36643191323985) + 0*((x-3)/2.36643191323985)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.639763269731185 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.70585616707015 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.595055200228685 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.513235067671825 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.50705390927465 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.8595159857035 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.809539969636 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0.7283989734076 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.81378356889533 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.58842212902895 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.59828139604185 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.48790044105748 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.73032363300804 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-1.53679094205 + 0*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2))

#------------------------
richnessPlot1 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(c) 0.9%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.44181191409335*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0.39616670639991 + 0.290830427392585*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,16)) +
  stat_function(fun=function(x){(1.022318756037 + 0.4117991038616*((x-2)/1.58113883008419) + 0*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) +
  stat_function(fun=function(x){(1.60835336305 + 0.42401222999135*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2))
  
#------------------------
richnessPlot2 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(e) 0.2%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0.30913779666705*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19))
  
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
  stat_function(fun=function(x){(0 + -0.49914040582*((x-5.5)/3.60555127546399) + 0*((x-5.5)/3.60555127546399)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,11)) +
  stat_function(fun=function(x){(-0.5016173969134 + -0.52904972436*((x-5.63636363636364)/3.93122696553448) + 0*((x-5.63636363636364)/3.93122696553448)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + -0.3985731688194*((x-4.55555555555556)/3.39525813124389) + 0*((x-4.55555555555556)/3.39525813124389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.6799896493018 + -0.5685201383569*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.6836207002 + -1.03031594695*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-1.68661937945 + -0.9893542251*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.9647520531495 + -0.453976046013965*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.981598805713 + -0.46609825721715*((x-2.33333333333333)/2.51661147842358) + 0*((x-2.33333333333333)/2.51661147842358)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + -0.37692491154944*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.6938334216581 + -0.39008620205657*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.641797242278*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.38198149206934*((x-2.33333333333333)/3.21455025366432) + 0*((x-2.33333333333333)/3.21455025366432)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(-0.31850231964505 + -0.400667321885*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-1.5229561034 + -0.7919432698*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + -0.6219897490857*((x-3.33333333333333)/3.51188458428425) + 0*((x-3.33333333333333)/3.51188458428425)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-1.012632767455 + -0.51005056398355*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(-0.64538445091145 + -0.459474467253*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0.3765325433606 + -0.2601960535371*((x-6)/3.89444048184931) + 0*((x-6)/3.89444048184931)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,12)) +
  stat_function(fun=function(x){(0 + -0.400686828045*((x-9.5)/5.91607978309962) + 0*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + -0.28230266293315*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.47875854389625 + -0.26701596375974*((x-9)/5.62731433871138) + 0*((x-9)/5.62731433871138)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,18)) +
  stat_function(fun=function(x){(-0.4421715212544 + -0.55558134545*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.6850095014485 + -0.9524255893*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + -0.288246920735535*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0.3723023802952 + -0.32460902555635*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + -0.3636965609855*((x-5)/3.3166247903554) + 0*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.9947703694865 + -0.6480181671225*((x-3.33333333333333)/3.05505046330389) + 0*((x-3.33333333333333)/3.05505046330389)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.42973725003101*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.4113923563682*((x-1)/1) + 0*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) +
  stat_function(fun=function(x){(0 + -0.36831846878358*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.547215056013918*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.47795707230105*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.45121444385405*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.36555878037525*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.39801070576235*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.54461740339925 + -0.4411218379335*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.5050888534076 + -0.3828643416076*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(-0.7027819350575 + -0.586156930725*((x-2.5)/1.87082869338697) + 0*((x-2.5)/1.87082869338697)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,5)) +
  stat_function(fun=function(x){(0 + -0.28588225567995*((x-3)/2.16024689946929) + 0*((x-3)/2.16024689946929)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,6)) +
  stat_function(fun=function(x){(0 + -0.42267954417805*((x-1.5)/1.29099444873581) + 0*((x-1.5)/1.29099444873581)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,3)) +
  stat_function(fun=function(x){(0 + -0.36812705989132*((x-4.5)/3.02765035409749) + 0*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) 
  
#------------------------
richnessPlot5 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(k) 0.7%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + -0.530878891905*((x-3.5)/2.44948974278318) + -0.287721361557*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.45950199277256 + -0.8059016976*((x-3.5)/2.44948974278318) + -0.2778404319041*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(-0.95648348015 + -1.16335031075*((x-3.5)/2.44948974278318) + -0.27147101485848*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) 
  
#------------------------
richnessPlot6 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(m) 2.5%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-2.0939772975 + -0.95962104105*((x-9.5)/5.91607978309962) + 0.4619298023*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-2.8331772775 + -1.07016848605*((x-9.5)/5.91607978309962) + 0.53571449275*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-3.2987774655 + -1.09201081485*((x-4.5)/3.02765035409749) + 0.8577669581*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-5.5241242285 + -1.019568862*((x-4.5)/3.02765035409749) + 1.35809478305*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-2.0324486655 + -0.68145109765*((x-4.5)/3.02765035409749) + 0.5793310788*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-3.488116406 + -1.42622611105*((x-4.5)/3.02765035409749) + 0.6967644382*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.6693239423 + -0.4212131046075*((x-4.5)/3.02765035409749) + 0.4172076921975*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-2.8742205255 + -0.72351667055*((x-4.5)/3.02765035409749) + 0.5722398146*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-0.61393220605132 + 0*((x-7.75)/6.84957419601151) + 0.303487926440495*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.52292354969942 + 0*((x-7.75)/6.84957419601151) + 0.2881972127081*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-1.2649985157 + -0.7065344681*((x-5)/3.3166247903554) + 0.27770669447673*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) 
  
#------------------------
richnessPlot7 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(o) 4.1%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.829185616815 + 0*((x-3.5)/2.44948974278318) + -0.308305109611965*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0.641172248081 + 0*((x-3.5)/2.44948974278318) + -0.331081780903*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(0 + -0.283735142101505*((x-3.5)/2.44948974278318) + -0.30319786719905*((x-3.5)/2.44948974278318)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,7)) +
  stat_function(fun=function(x){(1.04675907795 + 0.32190340621385*((x-3.75)/2.81577190634672) + -0.31777817122565*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.18326023845 + 0*((x-3.75)/2.81577190634672) + -0.35045500417095*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.0725616444 + 0*((x-3.75)/2.81577190634672) + -0.37782777309975*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.13349721865 + 0*((x-3.75)/2.81577190634672) + -0.407504888705*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.49089506684225 + 0*((x-3.75)/2.81577190634672) + -0.302244075123*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.29897727345 + 0*((x-3.75)/2.81577190634672) + -0.429341620217*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(0.805117022935 + 0*((x-3.75)/2.81577190634672) + -0.3207876707255*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.74250032365 + 0*((x-3.75)/2.81577190634672) + -0.52838150332*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.21906051815 + 0*((x-3.75)/2.81577190634672) + -0.3331415715565*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.343092555 + 0*((x-3.75)/2.81577190634672) + -0.440030659518*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.18415247325 + 0*((x-3.75)/2.81577190634672) + -0.352000053709*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.58006527065 + 0*((x-3.75)/2.81577190634672) + -0.50197641578*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.24694652735 + 0*((x-3.75)/2.81577190634672) + -0.33600256542555*((x-3.75)/2.81577190634672)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,8)) +
  stat_function(fun=function(x){(1.16463387175 + 0*((x-4.5)/3.02765035409749) + -0.278579694842705*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(1.033945351775 + 0*((x-2)/1.58113883008419) + -0.2858358888419*((x-2)/1.58113883008419)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,4)) 
  
#------------------------
richnessPlot8 <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('') +
  annotate('text', x=0, y=2.0, label='(q) 5.5%', size=10, hjust='left') +
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.3958865882397 + 0*((x-6.30769230769231)/4.44193305113542) + 0.198739653553624*((x-6.30769230769231)/4.44193305113542)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,15)) +
  stat_function(fun=function(x){(-1.02862283025 + 0*((x-9.5)/5.91607978309962) + 0.7251356974*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-3.3471712875 + 0.3267829701715*((x-9.5)/5.91607978309962) + 1.349283108*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-4.153622019 + 0.6764299824*((x-9.5)/5.91607978309962) + 1.5614351265*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-9.5)/5.91607978309962) + 0.230395758030415*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-2.729288661 + -0.250968531383845*((x-9.5)/5.91607978309962) + 0.78380787565*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-3.06269179 + -0.260621205556815*((x-9.5)/5.91607978309962) + 0.84092354635*((x-9.5)/5.91607978309962)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.5816567103435 + 0*((x-4.5)/3.02765035409749) + 0.4543186206718*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(0 + 0*((x-4.5)/3.02765035409749) + 0.27785170266791*((x-4.5)/3.02765035409749)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,9)) +
  stat_function(fun=function(x){(-1.13480577335 + 0*((x-9)/5.62731433871138) + 0.31210782598045*((x-9)/5.62731433871138)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-1.1487965346 + -0.33450547466097*((x-7.75)/6.84957419601151) + 0.484316575829*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.7814852596725 + 0*((x-7.75)/6.84957419601151) + 0.49445645693*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0.39474695175635*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-1.3854939163 + 0*((x-7.75)/6.84957419601151) + 0.62869052435*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-1.148691291 + 0*((x-7.75)/6.84957419601151) + 0.59001968075*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.5800057757827 + 0*((x-7.75)/6.84957419601151) + 0.3753986771225*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0.327870395070475*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0.36502134736995*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.59346385122955 + 0*((x-7.75)/6.84957419601151) + 0.3558651766508*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(0 + 0*((x-7.75)/6.84957419601151) + 0.288271603737875*((x-7.75)/6.84957419601151)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,19)) +
  stat_function(fun=function(x){(-0.84401680855 + -0.28473236147263*((x-5)/3.3166247903554) + 0.30733364305005*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(0 + 0*((x-5)/3.3166247903554) + 0.2292657911235*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-0.9475037901 + -0.274701452995*((x-5)/3.3166247903554) + 0.349366397897*((x-5)/3.3166247903554)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,10)) +
  stat_function(fun=function(x){(-1.16758772225 + 0*((x-1)/1) + 0.333685543724983*((x-1)/1)^2)*(0.3294397)+(-0.1075646)}, size=2, xlim=c(0,2)) 

  
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



# #summary stats from bayesian output --------------------------------------------------------
# #gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(
#         #trt_type intercepts: center digit refers to trts
#         E.1.1.1, E.2.1.1, E.1.2.1, E.2.2.1, E.1.3.1, E.2.3.1, E.1.4.1, E.2.4.1, E.1.5.1, E.2.5.1,
#         E.1.6.1, E.2.6.1, E.1.7.1, E.2.7.1, E.1.8.1, E.2.8.1, E.1.9.1, E.2.9.1, E.1.10.1, E.2.10.1,
#         E.1.11.1, E.2.11.1, E.1.12.1, E.2.12.1, E.1.13.1, E.2.13.1, E.1.14.1, E.2.14.1, E.1.15.1, E.2.15.1,
#         E.1.16.1, E.2.16.1,
#         #trt_type linear slopes: center digit refers to trts
#         E.1.1.2, E.2.1.2, E.1.2.2, E.2.2.2, E.1.3.2, E.2.3.2, E.1.4.2, E.2.4.2, E.1.5.2, E.2.5.2,
#         E.1.6.2, E.2.6.2, E.1.7.2, E.2.7.2, E.1.8.2, E.2.8.2, E.1.9.2, E.2.9.2, E.1.10.2, E.2.10.2,
#         E.1.11.2, E.2.11.2, E.1.12.2, E.2.12.2, E.1.13.2, E.2.13.2, E.1.14.2, E.2.14.2, E.1.15.2, E.2.15.2,
#         E.1.16.2, E.2.16.2,
#         #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#         E.1.1.3, E.2.1.3, E.1.2.3, E.2.2.3, E.1.3.3, E.2.3.3, E.1.4.3, E.2.4.3, E.1.5.3, E.2.5.3,
#         E.1.6.3, E.2.6.3, E.1.7.3, E.2.7.3, E.1.8.3, E.2.8.3, E.1.9.3, E.2.9.3, E.1.10.3, E.2.10.3,
#         E.1.11.3, E.2.11.3, E.1.12.3, E.2.12.3, E.1.13.3, E.2.13.3, E.1.14.3, E.2.14.3, E.1.15.3, E.2.15.3,
#         E.1.16.3, E.2.16.3,
#         #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
#         D.1.2.1, D.2.2.1,
#         D.1.2.2, D.2.2.2,
#         D.1.2.3, D.2.2.3,
#         #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
#         D.1.3.1, D.2.3.1,
#         D.1.3.2, D.2.3.2,
#         D.1.3.3, D.2.3.3,
#         #overall intercept, linear, and quad slopes (center digit): 1=overall
#         D.1.1.1, D.2.1.1,
#         D.1.1.2, D.2.1.2,
#         D.1.1.3, D.2.1.3)
# 
# chainsCommunitySummary <- chainsCommunitySummary%>%
#   gather(key=parameter, value=value, D.1.2.1:D.2.1.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   ungroup()%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor2=ifelse(predictor==2, 'ANPP', ifelse(predictor==3, 'rrich', 'overall')))%>%
#   select(variable, parameter, predictor2, median, sd, CI)

# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_04072019.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_04072019.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  mutate(type=paste(predictor2, parameter, sep='_'))



###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  # scale_y_continuous(limits=c(-0.15, 0.25), breaks=seq(-0.1, 0.2, 0.1)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Compositional Difference') +
  annotate('text', x=6.3, y=-0.15, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='richness' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  # scale_y_continuous(limits=c(-0.15, 0.25), breaks=seq(-0.1, 0.2, 0.1)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Richness Difference') +
  annotate('text', x=6.3, y=-0.15, label='(a)', size=10, hjust='left')

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