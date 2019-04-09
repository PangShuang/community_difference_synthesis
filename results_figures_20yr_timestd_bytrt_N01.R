################################################################################
##  results_figures_20yr.R: Compiles Bayesian output and makes figures for the primary analysis of richness and compositonal differences between treatment and control plots for datasets cut off at 20 years.
##
##  Author: Kimberly La Pierre
##  Date created: January 17, 2018
##  See https://github.com/klapierre/Converge_Diverge/blob/master/core%20data%20paper_bayesian%20results_figures_sig%20test_expinteractions_20yr.R for full history.
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
expRaw <- read.csv('ExperimentInformation_March2019.csv')


expInfo <- expRaw%>%
  #remove any pre-treatment data for the few experiments that have it -- pre-treatment data for experiments is awesome and we should all strive to collect it!
  filter(treatment_year!=0)%>%
  #make columns for irrigation and drought from precip column
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  #calcualte minumum years for each project
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

#import treatment data
trtInfo <- read.csv('ExperimentInformation_March2019.csv')

#import diversity metrics that went into Bayesian analysis
rawData <- read.csv('ForAnalysis_allAnalysis20yr_pairwise_04032019.csv')

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
  mutate(resource_mani=(nutrients+carbon+irrigation+drought))

#list of all treatments
studyInfo <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  select(site_code, community_type, project_name, treatment)%>%
  unique()


################################################################################
################################################################################
###Bayesian output processing

# only run to generate initial chains files
#raw chains data --------------------------------------------------------
memory.limit(size=50000)
chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\N01_lnRR\\N01_timestdbytrt_20yr_lnRR_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\N01_lnRR\\N01_timestdbytrt_20yr_lnRR_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\N01_lnRR\\N01_timestdbytrt_20yr_lnRR_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_04062019\\N01_lnRR\\N01_timestdbytrt_20yr_lnRR_3.csv', comment.char='#')
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
         #trt_type intercepts: center digit refers to trts
         E.1.1.1, E.2.1.1, E.1.2.1, E.2.2.1, E.1.3.1, E.2.3.1, E.1.4.1, E.2.4.1, E.1.5.1, E.2.5.1,
         E.1.6.1, E.2.6.1, E.1.7.1, E.2.7.1, E.1.8.1, E.2.8.1, E.1.9.1, E.2.9.1, E.1.10.1, E.2.10.1,
         E.1.11.1, E.2.11.1, E.1.12.1, E.2.12.1, E.1.13.1, E.2.13.1, E.1.14.1, E.2.14.1, E.1.15.1, E.2.15.1,
         E.1.16.1, E.2.16.1,
         #trt_type linear slopes: center digit refers to trts
         E.1.1.2, E.2.1.2, E.1.2.2, E.2.2.2, E.1.3.2, E.2.3.2, E.1.4.2, E.2.4.2, E.1.5.2, E.2.5.2,
         E.1.6.2, E.2.6.2, E.1.7.2, E.2.7.2, E.1.8.2, E.2.8.2, E.1.9.2, E.2.9.2, E.1.10.2, E.2.10.2,
         E.1.11.2, E.2.11.2, E.1.12.2, E.2.12.2, E.1.13.2, E.2.13.2, E.1.14.2, E.2.14.2, E.1.15.2, E.2.15.2,
         E.1.16.2, E.2.16.2,
         #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
         E.1.1.3, E.2.1.3, E.1.2.3, E.2.2.3, E.1.3.3, E.2.3.3, E.1.4.3, E.2.4.3, E.1.5.3, E.2.5.3,
         E.1.6.3, E.2.6.3, E.1.7.3, E.2.7.3, E.1.8.3, E.2.8.3, E.1.9.3, E.2.9.3, E.1.10.3, E.2.10.3,
         E.1.11.3, E.2.11.3, E.1.12.3, E.2.12.3, E.1.13.3, E.2.13.3, E.1.14.3, E.2.14.3, E.1.15.3, E.2.15.3,
         E.1.16.3, E.2.16.3,
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
  gather(key=parameter, value=value, E.1.1.1:D.2.1.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))

# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteraction_20yr_stdtimebytrt_N01_04072019.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteraction_20yr_stdtimebytrt_N01_04072019.csv')

#gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
#numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
#variable (second place): 1=mean change, 2=richness change
#parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
#set any that are not significant (CI overlaps 0) as 0

# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,8470:11025]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 8470:11025])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,8470:11025]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 8470:11025])'] <- 'sd'
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
# # write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_N01_04072019_noninf.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_N01_04072019_noninf.csv')


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
  left_join(trtID, by='id')%>%
  left_join(trtInfo)%>%
  left_join(timeStd)

#generate equations for main figure of richness and compositional responses through time
chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=20, 19, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1218158)+(0.2366542), (intercept+linear*7+quadratic*7^2)*(0.1821148)+(-0.02788138)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1218158)+(0.2366542),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1821148)+(-0.02788138)))%>%
  mutate(color=ifelse(rrich<31, '#1104DC44', ifelse(rrich<51&rrich>30, '#4403AE55', ifelse(rrich<71&rrich>50, '#77038166', ifelse(rrich>70, '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*(x*',
         curve4='+',
         curve5=') + ',
         curve6='*(x*',
         curve7='+',
         curve8=ifelse(variable=='mean', ')^2)*(0.1218158)+(0.2366542)}, size=2, xlim=c(0,',
                       ')^2)*(0.1821148)+(-0.02788138)}, size=2, xlim=c(0,'),
         curve9=')) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, time_std, curve4, time_mean, curve5, quadratic, curve6, time_std, curve7, time_mean, curve8, alt_length, curve9, sep=''))

#export, group by shape type, and paste lines below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_20yr_stdtimebytrt__N01_04072019.csv', row.names=F)

#summary lines by treatment type
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
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines

  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(0 + 0.2301080*x + -0.1112255*x^2)*(0.165778)+(0.3160571)}, size=5, xlim=c(0,19), colour='black')

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

  #overall line
  stat_function(fun=function(x){(0.24736400 + 0*x + 0*x^2)*0.2407859 + -0.06393594}, size=5, xlim=c(0,19), colour='black')

# print(richnessPlot) #export at 1200x1000

#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,1)))
print(richnessPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1200 x 2400



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

# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_N01_04072019.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_N01_04072019.csv')

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