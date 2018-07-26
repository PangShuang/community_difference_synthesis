################################################################################
##  site_information.R: Calculates gamma diversity and compiles site-level ANPP (control plots) for each experiment. Compiles MAP and MAT for each site.
##
##  Author: Kimberly La Pierre, Meghan Avolio
##  Date created: November 20, 2015
##  See https://github.com/klapierre/Converge_Diverge/blob/master/site_information_11202015.R for full history.
################################################################################

library(tidyverse)
library(vegan)

#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")


#import the list of all experiment's site information from datafile
expInfo <- read.csv("SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)

#generate a list of each site, project (i.e., experiment), and community type
expList<-expInfo%>%
  select(site_code, project_name, community_type)%>%
  unique()


###generate list of project-level ANPP under ambient conditions (i.e., no treatments)
#import ANPP data
ANPP<-read.csv("ANPP_Oct2017.csv")

#import data to get number of treatment manipulations (i.e., plot_mani) for each treatment
plotMani<-read.csv("ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment, plot_mani)%>%
  unique()

#generate a list of each site, project, and community type
expList2<-plotMani%>%
  select(site_code, project_name, community_type)%>%
  unique()

#for projects that collected ANPP, galculate mean ANPP across only control plots (i.e., where plot_mani==0)
controlANPP<-merge(ANPP, plotMani, by=c("site_code","project_name","community_type","treatment"))%>%
  filter(plot_mani==0)%>%
  na.omit%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize(anpp=mean(anpp))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(anpp=mean(anpp))%>%
  ungroup()

#for projects that did not collect ANPP, import ANPP values provided by experiment PIs
ANPPnocont<-read.csv("ANPP_noControls.csv")

#bind ANPP from projects with and without ANPP data
allANPP<-rbind(ANPPnocont, controlANPP)

#merge project ANPP values with project list
expANPP<-merge(allANPP, expList2, by=c("site_code","project_name","community_type"), all=T)



###calculate chao richness and rarefied richness for each site

#import species abundance data
species <- read.csv("SpeciesRawAbundance_Oct2017.csv")%>%
  select(site_code, project_name, community_type, plot_id, calendar_year, genus_species, abundance)%>%
  mutate(exp=paste(site_code, project_name, community_type, sep='::'))%>%
  tbl_df()

#determine sampling intensity for each project
sampleIntensity<-species%>%
  group_by(exp, plot_id, calendar_year)%>%
  summarize(sample_intensity=length(abundance))%>%
  ungroup()%>%
  group_by(exp)%>%
  summarize(sample_intensity=length(sample_intensity))%>% #how many plots were sampled over the course of the experiment
  ungroup()

#generate list of projects to calculate rarefied richness for
exp<-sampleIntensity%>%
  select(exp)


#create empty dataframe for loop
estimatedRichness=data.frame(row.names=1) 

for(i in 1:length(exp$exp)) {

  #creates a dataset for each unique experiment
  subset <- species%>%
    filter(exp==exp$exp[i])%>%
    select(exp, plot_id, calendar_year, genus_species, abundance)
  
  #transpose data into wide form
  speciesData <- subset%>%
    spread(genus_species, abundance, fill=0)
  
  #calculate species accumulation curves
  pool <- poolaccum(speciesData[,4:ncol(speciesData)], permutations=100)
  chao <- as.data.frame(as.matrix(pool$chao))#this gives us estimated richness from 1-X samples
  chao$aveChao<-rowMeans(chao)
  chao$n<-row.names(chao)
  chao$exp<-exp$exp[i]
  chao2<-chao%>%
    select(exp,n, aveChao)
  
  #rbind back
  estimatedRichness<-rbind(chao2, estimatedRichness)
}

#generate a list of estimated richness for each project
expRichness<-estimatedRichness%>%
  filter(n==34)%>%#the lowest sampling intensity
  separate(exp, c("site_code", "project_name", "community_type"), sep="::")%>%
  mutate(rrich=aveChao)%>%
  select(-n, -aveChao)

#generate a list of project lengths
expLength<-expInfo%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))

#merge estimated richness with project length and ANPP
expDetails<-expRichness%>%
  left_join(expLength)%>%
  left_join(expANPP)


#import climate data for each site by coordinates
siteClimate<-read.csv("siteList_climate_Feb2016.csv")%>%
  mutate(MAP=ifelse(site_code=="Finse", 1030, MAP))%>%
  select(site_code, MAP, MAT)
#for Finse_WarmNut there is a big differnce between this and what they published, and their coordinates were VERY vauge. Replacing with their value: 1030 mm.


#merge site and project dataframes
siteExpDetails<-merge(siteClimate, expDetails, by="site_code")

# #kim
# write.csv(siteExpDetails, 'C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\SiteExperimentDetails_Dec2016.csv')

rm(list=setdiff(ls(), "siteExpDetails"))
