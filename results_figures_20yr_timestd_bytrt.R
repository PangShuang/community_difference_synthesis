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
rawData <- read.csv('ForAnalysis_allAnalysis20yr_pairwise_03132019.csv')

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

# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_03162019\\noninformative_timestdbytrt\\noninformative_timestdbytrt_20yr_lnRR_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_03162019\\noninformative_timestdbytrt\\noninformative_timestdbytrt_20yr_lnRR_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_03162019\\noninformative_timestdbytrt\\noninformative_timestdbytrt_20yr_lnRR_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\final models_03162019\\noninformative_timestdbytrt\\noninformative_timestdbytrt_20yr_lnRR_3.csv', comment.char='#')
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
#          E.1.16.1, E.2.16.1, E.1.17.1, E.2.17.1, E.1.18.1, E.2.18.1,
#          #trt_type linear slopes: center digit refers to trts
#          E.1.1.2, E.2.1.2, E.1.2.2, E.2.2.2, E.1.3.2, E.2.3.2, E.1.4.2, E.2.4.2, E.1.5.2, E.2.5.2,
#          E.1.6.2, E.2.6.2, E.1.7.2, E.2.7.2, E.1.8.2, E.2.8.2, E.1.9.2, E.2.9.2, E.1.10.2, E.2.10.2,
#          E.1.11.2, E.2.11.2, E.1.12.2, E.2.12.2, E.1.13.2, E.2.13.2, E.1.14.2, E.2.14.2, E.1.15.2, E.2.15.2,
#          E.1.16.2, E.2.16.2, E.1.17.2, E.2.17.2, E.1.18.2, E.2.18.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          E.1.1.3, E.2.1.3, E.1.2.3, E.2.2.3, E.1.3.3, E.2.3.3, E.1.4.3, E.2.4.3, E.1.5.3, E.2.5.3,
#          E.1.6.3, E.2.6.3, E.1.7.3, E.2.7.3, E.1.8.3, E.2.8.3, E.1.9.3, E.2.9.3, E.1.10.3, E.2.10.3,
#          E.1.11.3, E.2.11.3, E.1.12.3, E.2.12.3, E.1.13.3, E.2.13.3, E.1.14.3, E.2.14.3, E.1.15.3, E.2.15.3,
#          E.1.16.3, E.2.16.3, E.1.17.3, E.2.17.3, E.1.18.3, E.2.18.3,
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
# # write.csv(chainsCommunity2, 'bayesian_output_summary_expinteraction_20yr_stdtimebytrt_03192019.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteraction_20yr_stdtimebytrt_03192019.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,8220:10691]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 8220:10691])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,8220:10691]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 8220:10691])'] <- 'sd'
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
# # write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_03192019_noninf.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_20yr_stdtimebytrt_03192019_noninf.csv')

# chainsFinal2 <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #spread by variable
#   select(variable, id, parameter, mean)%>%
#   spread(key=parameter, value=mean)

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
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.165778)+(0.3160571), (intercept+linear*7+quadratic*7^2)*(0.2407859)+(-0.06393594)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.165778)+(0.3160571),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2407859)+(-0.06393594)))%>%
  mutate(color=ifelse(rrich<31, '#1104DC44', ifelse(rrich<51&rrich>30, '#4403AE55', ifelse(rrich<71&rrich>50, '#77038166', ifelse(rrich>70, '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*(x*',
         curve4='+',
         curve5=') + ',
         curve6='*(x*',
         curve7='+',
         curve8=ifelse(variable=='mean', ')^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,',
                       ')^2)*(0.2407859)+(-0.06393594)}, size=2, xlim=c(0,'),
         curve9=')) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, time_std, curve4, time_mean, curve5, quadratic, curve6, time_std, curve7, time_mean, curve8, alt_length, curve9, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))%>%
  left_join(read.csv('treatment_response_shape_classification_stdtimebytrt_03192019.csv'))
#export, group by shape type, and paste lines below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_20yr_stdtimebytrt_03192019.csv', row.names=F)

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



#summary stats from bayesian output --------------------------------------------------------
# #gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(
#         #trt_type intercepts: center digit refers to trts
#         E.1.1.1, E.2.1.1, E.1.2.1, E.2.2.1, E.1.3.1, E.2.3.1, E.1.4.1, E.2.4.1, E.1.5.1, E.2.5.1,
#         E.1.6.1, E.2.6.1, E.1.7.1, E.2.7.1, E.1.8.1, E.2.8.1, E.1.9.1, E.2.9.1, E.1.10.1, E.2.10.1,
#         E.1.11.1, E.2.11.1, E.1.12.1, E.2.12.1, E.1.13.1, E.2.13.1, E.1.14.1, E.2.14.1, E.1.15.1, E.2.15.1,
#         E.1.16.1, E.2.16.1, E.1.17.1, E.2.17.1, E.1.18.1, E.2.18.1,
#         #trt_type linear slopes: center digit refers to trts
#         E.1.1.2, E.2.1.2, E.1.2.2, E.2.2.2, E.1.3.2, E.2.3.2, E.1.4.2, E.2.4.2, E.1.5.2, E.2.5.2,
#         E.1.6.2, E.2.6.2, E.1.7.2, E.2.7.2, E.1.8.2, E.2.8.2, E.1.9.2, E.2.9.2, E.1.10.2, E.2.10.2,
#         E.1.11.2, E.2.11.2, E.1.12.2, E.2.12.2, E.1.13.2, E.2.13.2, E.1.14.2, E.2.14.2, E.1.15.2, E.2.15.2,
#         E.1.16.2, E.2.16.2, E.1.17.2, E.2.17.2, E.1.18.2, E.2.18.2,
#         #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#         E.1.1.3, E.2.1.3, E.1.2.3, E.2.2.3, E.1.3.3, E.2.3.3, E.1.4.3, E.2.4.3, E.1.5.3, E.2.5.3,
#         E.1.6.3, E.2.6.3, E.1.7.3, E.2.7.3, E.1.8.3, E.2.8.3, E.1.9.3, E.2.9.3, E.1.10.3, E.2.10.3,
#         E.1.11.3, E.2.11.3, E.1.12.3, E.2.12.3, E.1.13.3, E.2.13.3, E.1.14.3, E.2.14.3, E.1.15.3, E.2.15.3,
#         E.1.16.3, E.2.16.3, E.1.17.3, E.2.17.3, E.1.18.3, E.2.18.3,
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

# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_03192019.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_20yr_stdtimebytrt_03192019.csv')

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









####testing lines  ---  mean change
BGP <- ggplot(data=subset(rawData, project_name=='BGP'&treatment=='b_u_n'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.449065262*x + -0.065372924*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

BGP2 <- ggplot(data=subset(rawData, project_name=='BGP'&treatment=='b_u_p'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-1.04788779 + 0.480665506*x + -0.06235107*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

nfert <- ggplot(data=subset(rawData, project_name=='Nfert'&treatment=='F'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(0.822832418 + 0.699235591*x + -0.04209497*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,18), colour='black')

wenndex <- ggplot(data=subset(rawData, project_name=='WENNDEx'&treatment=='PN'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.693766164*x + -0.038141653*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,7), colour='black')

fireplots <- ggplot(data=subset(rawData, project_name=='fireplots'&treatment=='wnug'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.520600303*x + -0.0378649*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,9), colour='black')

tmece <- ggplot(data=subset(rawData, project_name=='TMECE'&community_type=='SP'&treatment=='E'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.711655574*x + -0.036138435*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

pushViewport(viewport(layout=grid.layout(2,3)))
print(BGP, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(BGP2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(nfert, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(wenndex, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(fireplots, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(tmece, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))



####testing lines  --  richness change
BGP <- ggplot(data=subset(rawData, project_name=='BGP'&treatment=='b_u_n'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.449065262*x + -0.065372924*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

BGP2 <- ggplot(data=subset(rawData, project_name=='BGP'&treatment=='b_u_p'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-1.04788779 + 0.480665506*x + -0.06235107*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

nfert <- ggplot(data=subset(rawData, project_name=='Nfert'&treatment=='F'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(0.822832418 + 0.699235591*x + -0.04209497*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,18), colour='black')

wenndex <- ggplot(data=subset(rawData, project_name=='WENNDEx'&treatment=='PN'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.693766164*x + -0.038141653*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,7), colour='black')

fireplots <- ggplot(data=subset(rawData, project_name=='fireplots'&treatment=='wnug'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.520600303*x + -0.0378649*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,9), colour='black')

tmece <- ggplot(data=subset(rawData, project_name=='TMECE'&community_type=='SP'&treatment=='E'), aes(x=treatment_year, y=mean_change)) + 
  # coord_cartesian(ylim=c(0,1))  +
  # scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  geom_point(color='red', size=5) +
  stat_function(fun=function(x){(-0 + 0.711655574*x + -0.036138435*x^2)*(0.1718192)+(0.2935109)}, size=1, xlim=c(1,20), colour='black')

pushViewport(viewport(layout=grid.layout(2,3)))
print(BGP, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(BGP2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(nfert, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(wenndex, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(fireplots, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(tmece, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
