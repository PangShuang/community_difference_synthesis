################################################################################
##  time_comparisons.R: Calculating changes and differences between treatment and control plots between first and last years of treatments for only the treatments that exhibited significant curvature.
##
##  Author: Kimberly La Pierre
##  Date created: March 18, 2019
################################################################################

library(tidyverse)
library(codyn)


#kim
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')


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


#import relative species abundance data
alldata <- read.csv("SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))

maxMin <- alldata%>%
  select(site_code, project_name, community_type, treatment, calendar_year)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(max_year=max(calendar_year), min_year=min(calendar_year))%>%
  ungroup()%>%
  select(-calendar_year)%>%
  unique()

#import treatment information and subset to get number of factors manipulated (i.e., plot_mani) for each plot
expinfo <- read.csv("ExperimentInformation_Nov2017.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

#merge treatment information with species relative abundancecs
alldata2 <- merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  #merge first and last years of data for each experiment
  left_join(maxMin)%>%
  #filter to only keep first and last years
  mutate(keep=ifelse(max_year==calendar_year, 'last', ifelse(min_year==calendar_year, 'first', 'drop')))%>%
  filter(keep!='drop')%>%
  #kbs needs special treatment, see below
  filter(site_code!='KBS')




###calculating bray-curtis dissimilarities within and among treatments to get distances between treatment centroids and dispersion among replicate plots within a treatment
#make a new dataframe with just the label
exp_year=alldata2%>%
  select(site_project_comm)%>%
  unique()

#makes an empty dataframe
for.analysis.difference=data.frame(row.names=1) 
for.analysis.change=data.frame(row.names=1) 

###first: composition_diff is the distance between trt and controls
####second: dispersion is the average dispersion of plots within a treatment to treatment centriod
####third: richness and exp_H also calculated
for(i in 1:length(exp_year$site_project_comm)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- alldata2[alldata2$site_project_comm==as.character(exp_year$site_project_comm[i]),]%>%
    select(site_project_comm, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
    mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))
  
  #need this to keep track of treatment
  treatments <- subset%>%
    select(plot_id, treatment)%>%
    unique()
  
  #calculating composition difference and abs(dispersion difference)
  difference <- multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
    rename(treatment=treatment22)%>%
    select(-treatment2, -trt_greater_disp)%>%
    mutate(site_project_comm=exp_year$site_project_comm[i])
  
  #calculating composition difference and abs(dispersion difference)
  change <- multivariate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=min(subset$calendar_year))%>%
    rename(treatment=treatment2)%>%
    mutate(site_project_comm=exp_year$site_project_comm[i])
  
  #pasting dispersions into the dataframe made for this analysis
  for.analysis.difference=rbind(difference, for.analysis.difference)  
  for.analysis.change=rbind(change, for.analysis.change)  
}

#kbs T7 special treatment because tilling began one year after fertilization
kbsUntilled <- merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)%>%
  filter(site_code=='KBS', calendar_year==1989|calendar_year==2012, treatment=='T0F0'|treatment=='T0F1')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  #filter to only keep first and last years
  mutate(keep=ifelse(calendar_year==1989, 'first', 'last'))%>%
  group_by(treatment)%>%
  mutate(max_year=max(calendar_year), min_year=min(calendar_year))%>%
  ungroup()%>%
  select(site_project_comm, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
  mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))

differenceKBSuntilled <- multivariate_difference(kbsUntilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
  rename(treatment=treatment22)%>%
  select(-treatment2, -trt_greater_disp)%>%
  mutate(site_project_comm='KBS::T7::0')
  

changeKBSuntilled <- multivariate_change(kbsUntilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=1989)%>%
  rename(treatment=treatment2)%>%
  mutate(site_project_comm='KBS::T7::0')


kbsTilled <- merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)%>%
  filter(site_code=='KBS', calendar_year==1990|calendar_year==2012, treatment!='T0F1')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  #filter to only keep first and last years
  mutate(keep=ifelse(calendar_year==1990, 'first', 'last'))%>%
  group_by(treatment)%>%
  mutate(max_year=max(calendar_year), min_year=min(calendar_year))%>%
  ungroup()%>%
  select(site_project_comm, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
  mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))

#calculating composition difference and abs(dispersion difference)
differenceKBStilled <- multivariate_difference(kbsTilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
  rename(treatment=treatment22)%>%
  select(-treatment2, -trt_greater_disp)%>%
  mutate(site_project_comm='KBS::T7::0')

#calculating composition difference and abs(dispersion difference)
changeKBStilled <- multivariate_change(kbsTilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=1990)%>%
  rename(treatment=treatment2)%>%
  mutate(site_project_comm='KBS::T7::0')

for.analysis.difference=rbind(for.analysis.difference, differenceKBStilled, differenceKBSuntilled)  
for.analysis.change=rbind(for.analysis.change, changeKBStilled, changeKBSuntilled) 


###final numbers for difference and change
curveShape <- read.csv('treatment_response_shape_classification_03192019.csv')%>%
  filter(variable=='mean')%>%
  select(site_code, project_name, community_type, treatment, shape_category)

differenceAll <- for.analysis.difference%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(shape_category==1|shape_category==2)%>%
  group_by(site_project_comm)%>%
  mutate(comparison=ifelse(calendar_year==min(calendar_year), 'trt-ctl_first', 'trt-ctl_last'))%>%
  ungroup()%>%
  rename(value=composition_diff)%>%
  select(comparison, site_project_comm, treatment, calendar_year, value)

changeTrt <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  filter(shape_category==1|shape_category==2)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

expToUse <- changeTrt%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(site_project_comm)%>%
  unique()

changeCtl <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  full_join(curveShape)%>%
  mutate(shape_category=ifelse(treatment=='TRUECONTROL', 999, shape_category))%>%
  filter(shape_category==999)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(site_project_comm %in% expToUse$site_project_comm)

changeAll <- rbind(changeTrt, changeCtl)%>%
  mutate(comparison=ifelse(treatment=='TRUECONTROL', 'ctl_first-last', 'trt_first-last'))%>%
  rename(value=composition_change)%>%
  select(site_project_comm, treatment, comparison, calendar_year, value)

comparisons <- rbind(differenceAll, changeAll)


###figures
ggplot(data=barGraphStats(data=comparisons, variable="value", byFactorNames=c("comparison")), aes(x=comparison, y=mean)) +
  geom_bar(stat='identity', color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  ylab('Response Magnitude') +
  scale_x_discrete(limits=c('ctl_first-last', 'trt_first-last', 'trt-ctl_first', 'trt-ctl_last'),
                   labels=c('ctl change', 'trt change', 'initial diff', 'final diff')) +
  theme(axis.title.x=element_blank())

ggplot(data=comparisons, aes(x=comparison, y=value)) +
  geom_boxplot() +
  ylab('Response Magnitude') +
  scale_x_discrete(limits=c('ctl_first-last', 'trt_first-last', 'trt-ctl_first', 'trt-ctl_last'),
                   labels=c('ctl change', 'trt change', 'initial diff', 'final diff')) +
  theme(axis.title.x=element_blank())
