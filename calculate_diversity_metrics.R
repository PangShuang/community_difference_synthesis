################################################################################
##  calculate_diversity_metrics.R: Calculating diversity metrics (richness, evenness, Bray-Curtis dissimilarity). Formatting data for Bayesian analysis in python.
##
##  Author: Kimberly La Pierre, Meghan Avolio, Emily Grman, Forest Isbell
##  Date created: November 19, 2015
##  See https://github.com/klapierre/Converge_Diverge/blob/master/Merging%20to%20a%20single%20file_11172015.R for full history.
################################################################################

library(vegan)
library(tidyverse)
library(gridExtra)
library(grid)
library(codyn)


#kim laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")



#import relative species abundance data
alldata<-read.csv("SpeciesRelativeAbundance_March2019.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))

#import treatment information and subset to get number of factors manipulated (i.e., plot_mani) for each plot
expinfo<-read.csv("ExperimentInformation_March2019.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

#merge treatment information with species relative abundancecs
alldata2<-merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))



###calculating bray-curtis dissimilarities within and among treatments to get distances between treatment centroids and dispersion among replicate plots within a treatment
#make a new dataframe with just the label
exp_year=alldata2%>%
  select(site_project_comm)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1) 

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
  multivariate <- multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
    rename(treatment=treatment22)%>%
    select(-treatment2, -trt_greater_disp)

  #calculating univariate community metrics for each plot
  univariate <- community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')%>%
    left_join(treatments)%>%
    group_by(calendar_year, treatment)%>%
    summarise(S=mean(richness))%>%
    ungroup()
  
  #calculate e^H
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  H <- diversity(species[,7:ncol(species)])%>%
    cbind(species[,1:6])%>%
    mutate(expH=exp(.))%>%
    group_by(calendar_year, treatment)%>%
    summarise(exp_H=mean(expH))%>%
    ungroup()

  #merge multivariate and univariate metrics
  all <- univariate%>%
    full_join(H)%>%
    full_join(multivariate)%>%
    mutate(site_project_comm=exp_year$site_project_comm[i])

  #pasting dispersions into the dataframe made for this analysis
  for.analysis=rbind(all, for.analysis)  
}

rm(list=setdiff(ls(), "for.analysis"))


###formatting data for Bayesian analysis in python

#import treatment information
expInfo <- read.csv('ExperimentInformation_March2019.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))

#diversity data
div <- for.analysis%>%
  left_join(expInfo)%>%
  filter(treatment_year!=0)

#import site and project level data (MAP, MAT, ANPP of controls, and rarefied gamma diversity)
SiteExp<-read.csv("SiteExperimentDetails_March2019.csv")%>%
  select(-X)

###calculate difference in dispersion, H, S, and evenness between treatment and control plots for each year

#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, exp_H, S, calendar_year, treatment_year)%>%
  rename(ctl_expH=exp_H, ctl_S=S)

#filtering to get only non-e002, e001 treatments, later will remove a subset of the treatments from these two projects because they have so many levels and make up a disproportionate amount of the entire dataset
divTrt1 <- div%>%
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

#removing a subset of CDR e001 and e002 treatments to prevent the majority of data being from CDR; keeping lowest, highest, and 10 gm-2 N additions (levels most comparable to other studies)
divCDRe001<-div%>%
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)
divCDRe002<-div%>%
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)

#combine the two CDR experiments
divTrt<-rbind(divTrt1, divCDRe002, divCDRe001)

##16% of our data is from CDR and 10% is from KNZ
#merge controls and treatments
divCompare <- divControls%>%
  left_join(divTrt)%>%
  #calculate the proportional difference and log response ratio of S and eH
  mutate(expH_PC=(exp_H-ctl_expH)/ctl_expH, 
         S_PC=(S-ctl_S)/ctl_S,  
         expH_lnRR=log(exp_H/ctl_expH), 
         S_lnRR=log(S/ctl_S),
         S_lnRR_abs=abs(S_lnRR))%>%
  select(exp_year, treatment_year, treatment, plot_mani, composition_diff, expH_PC, expH_lnRR, S_PC, S_lnRR, S_lnRR_abs, site_code, project_name, community_type, calendar_year)

#some preliminary histograms
theme_set(theme_bw(16))
m<-qplot(composition_diff, data=divCompare, geom="histogram")+
  ggtitle("Among Treatment Change")+
  geom_vline(xintercept = 0, size=2)

s1<-qplot(S_PC, data=divCompare, geom="histogram")+
  ggtitle("Richness Percent Change")+
  geom_vline(xintercept = 0, size=2)

e2<-qplot(expH_PC, data=divCompare, geom="histogram")+
  ggtitle("Effective Diversity Percent Change")+
  geom_vline(xintercept = 0, size=2)

s2<-qplot(S_lnRR, data=divCompare, geom="histogram")+
  ggtitle("Richness ln Response Ratio")+
  geom_vline(xintercept = 0, size=2)

e3<-qplot(expH_lnRR, data=divCompare, geom="histogram")+
  ggtitle("Effective Diversity ln Response Ratio")+
  geom_vline(xintercept = 0, size=2)

grid.arrange(m, s1, e2, s2, e3, ncol=3)


###merging with treatment information
SiteExp<-read.csv("SiteExperimentDetails_March2019.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))%>%
  left_join(read.csv('ExperimentInformation_March2019.csv'))

# full dataset
# write.csv(ForAnalysis, "ForBayesianAnalysis_March2019.csv")


###generating treatment categories (resource, non-resource, and interactions)
allAnalysis <- ForAnalysis%>%
  select(-plot_mani)%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #modify categorical treatment column
  mutate(trt_type2=ifelse(trt_type=='CO2', 'CO2', ifelse(trt_type=='drought', 'drought', ifelse(trt_type=='irr', 'irr', ifelse(trt_type=='N', 'N', ifelse(trt_type=='P', 'P', ifelse(trt_type=='precip_vari', 'precip_vari', ifelse(trt_type %in% c('light','lime'), 'other_resource', ifelse(trt_type=='herb_removal', 'herb_removal', ifelse(trt_type=='plant_mani', 'plant_mani', ifelse(trt_type=='mow_clip', 'mow_clip', ifelse(trt_type=='temp', 'temp', ifelse(trt_type %in% c('burn','fungicide','disturbance','stone','till'), 'other_nonresource', ifelse(trt_type %in% c('irr*CO2','N*CO2','N*drought','N*irr','N*P'), 'RxR', ifelse(trt_type %in% c('burn*graze','burn*mow_clip','herb_removal*mow_clip','plant_mani*herb_removal','plant_mani*other','temp*mow_clip'), 'NxN', ifelse(trt_type %in% c('CO2*temp','drought*mow_clip','drought*temp','irr*mow_clip','irr*plant_mani','irr*temp','N*burn','N*mow_clip','N*other','N*plant_mani','N*stone','N*temp','N*till','P*burn','P*mow_clip','precip_vari*temp'), 'RxN', ifelse(trt_type %in% c('mult_nutrient','N*irr*CO2'), 'RxRxR', 'threeway')))))))))))))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type2, composition_diff, expH_PC, expH_lnRR, S_PC, S_lnRR, S_lnRR_abs, experiment_length, rrich, anpp, MAT, MAP)%>%
    rename(trt_type=trt_type2)


#subset out datasets with less than 3 temporal data points;  drops GVN FACE
numPoints <- allAnalysis%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
allAnalysisAllDatasets <- allAnalysis%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)%>%
  #filter datasets with NA for comp diff (imaginary axis problem)
  filter(!is.na(composition_diff))
# write.csv(allAnalysisAllDatasets, 'ForAnalysis_allAnalysisAllDatasets_04082019.csv')


# #subset out treatment years 20 or less (i.e., cut off datasets at 20 years)
# allAnalysis20yr <- allAnalysisAllDatasets%>%
#   filter(treatment_year<21)%>%
#   select(-num_datapoints)%>%
#   rename(mean_change=composition_diff)
# numPoints <- allAnalysis20yr%>%
#   select(site_code, project_name, community_type, treatment, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   summarise(num_datapoints=length(treatment_year))%>%
#   ungroup()
# # write.csv(allAnalysis20yr, 'ForAnalysis_allAnalysis20yr_pairwise_04082019.csv')
# 
# 
# #subset out treatment years 10 or less (i.e., cut off datasets at 10 years)
# allAnalysis10yr <- allAnalysisAllDatasets%>%
#   filter(treatment_year<11)%>%
#   select(-num_datapoints)%>%
#   filter(!is.na(composition_diff))%>%
#   rename(mean_change=composition_diff)
# numPoints <- allAnalysis10yr%>%
#   select(site_code, project_name, community_type, treatment, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   summarise(num_datapoints=length(treatment_year))
# # write.csv(allAnalysis10yr, 'ForAnalysis_allAnalysis10yr_pairwise_04082019.csv')

#subset out final year of all data
allAnalysisFinalYear <- allAnalysis%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  filter(treatment_year==max(treatment_year))%>%
  ungroup()
# write.csv(allAnalysisFinalYear, 'ForAnalysis_allAnalysisFinalYear_04082019.csv')


###treatment magnitude data
nMag <- allAnalysisFinalYear%>%
  filter(trt_type=='N')%>%
  left_join(expInfo)%>%
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, composition_diff, expH_PC, expH_lnRR, S_PC, S_lnRR, S_lnRR_abs, experiment_length, rrich, anpp, MAT, MAP, n)
# write.csv(nMag, 'ForAnalysis_allAnalysisNmag.csv', row.names=F)

irrMag <- allAnalysisFinalYear%>%
  filter(trt_type=='irr')%>%
  left_join(expInfo)%>%
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, composition_diff, expH_PC, expH_lnRR, S_PC, S_lnRR, S_lnRR_abs, experiment_length, rrich, anpp, MAT, MAP, precip)
# write.csv(irrMag, 'ForAnalysis_allAnalysisH2Omag_irr.csv', row.names=F)

droughtMag <- allAnalysisFinalYear%>%
  filter(trt_type=='drought')%>%
  left_join(expInfo)%>%
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, composition_diff, expH_PC, expH_lnRR, S_PC, S_lnRR, S_lnRR_abs, experiment_length, rrich, anpp, MAT, MAP, precip)
# write.csv(droughtMag, 'ForAnalysis_allAnalysisH2Omag_drought.csv', row.names=F)
