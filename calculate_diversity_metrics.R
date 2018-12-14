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


#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")



#import relative species abundance data
alldata<-read.csv("SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))

#import treatment information and subset to get number of factors manipulated (i.e., plot_mani) for each plot
expinfo<-read.csv("ExperimentInformation_Nov2017.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

#merge treatment information with species relative abundancecs
alldata2<-merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)



###calculating bray-curtis dissimilarities within and among treatments to get distances between treatment centroids and dispersion among replicate plots within a treatment
#make a new dataframe with just the label
exp_year=alldata2%>%
  select(exp_year)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1) 

###first, get bray curtis dissimilarity values for each year within each experiment between all combinations of plots
###second, get distance of each plot within a trt to the trt centroid 
###third: mean_change is the distance between trt and control centriods
####fourth: dispersion is the average dispersion of plots within a treatment to treatment centriod
for(i in 1:length(exp_year$exp_year)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=alldata2[alldata2$exp_year==as.character(exp_year$exp_year[i]),]%>%
    select(exp_year, treatment, plot_mani, genus_species, relcov, plot_id)
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
   spread(genus_species, relcov, fill=0)

  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid", bias.adjust=T) #bias.adjust takes sqrt
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  #merging back with labels to get back plot_mani
  centroid=merge(cent_C_T, labels, by="treatment")
  
  #collecting and labeling distances to centroid from betadisper
  trt_disp=data.frame(data.frame(exp_year=exp_year$exp_year[i], 
                                 plot_id=species$plot_id,
                                 treatment=species$treatment,
                                 dist=disp$distances))%>%
    tbl_df()%>%
    group_by(exp_year, treatment)%>%
    summarize(dispersion=mean(dist))
  
  #merge together change in mean and dispersion data
  distances<-merge(centroid, trt_disp, by=c("exp_year","treatment"))
  
  #getting diversity indixes
    H<-diversity(species[,5:ncol(species)])
    S<-specnumber(species[,5:ncol(species)])
    InvD<-diversity(species[,5:ncol(species)],"inv")
    SimpEven<-InvD/S
    out1<-cbind(H, S)
    output<-cbind(out1, SimpEven)
    divmeasure<-cbind(species, output)%>%
      select(exp_year, treatment, H, S, SimpEven)%>%
      tbl_df()%>%
      group_by(exp_year, treatment)%>%
      summarise(H=mean(H), S=mean(S), SimpEven=mean(SimpEven))
    
    ##merging all measures of diversity
    alldiv<-merge(distances, divmeasure, by=c("exp_year","treatment"))

    #pasting dispersions into the dataframe made for this analysis
    for.analysis=rbind(alldiv, for.analysis)  
}

# write.csv(for.analysis, 'DiversityMetrics_Nov2018_2.csv')

rm(list=setdiff(ls(), "for.analysis"))



###formatting data for Bayesian analysis in python

#import treatment information
expInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))

#diversity data
div <- for.analysis%>%
  separate(exp_year, c('site_code', 'project_name','community_type', 'calendar_year'), sep='::')%>%
  mutate(calendar_year=as.integer(calendar_year))%>%
  left_join(expInfo)%>%
  filter(treatment_year!=0)%>%
  #calculate e^H metric
  mutate(expH=exp(H))

#import site and project level data (MAP, MAT, ANPP of controls, and rarefied gamma diversity)
SiteExp<-read.csv("SiteExperimentDetails_09062018.csv")%>%
  select(-X)

###calculate difference in dispersion, H, S, and evenness between treatment and control plots for each year

#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, expH, S, SimpEven, calendar_year, treatment_year)
names(divControls)[names(divControls)=='dispersion'] <- 'ctl_dispersion'
names(divControls)[names(divControls)=='expH'] <- 'ctl_expH'
names(divControls)[names(divControls)=='S'] <- 'ctl_S'
names(divControls)[names(divControls)=='SimpEven'] <- 'ctl_SimpEven'

#filtering to get only non-e002, e001 treatments, later will remove a subset of the treatments from these two projects because they have so many levels and make up a disproportionate amount of the entire dataset
divTrt1 <- div%>%
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

#calculate average dispersion among just control plots (this is an estimate of average community dissimilarity, to use as a baseline for dissimilarity change between treatment and control plots)
summary(divControls$ctl_dispersion) #mean=0.2811, median=0.2778

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
  mutate(expH_PC=(expH-ctl_expH)/ctl_expH, 
         S_PC=(S-ctl_S)/ctl_S,  
         expH_lnRR=log(expH/ctl_expH), 
         S_lnRR=log(S/ctl_S))%>%
  select(exp_year, treatment_year, treatment, plot_mani, mean_change, expH_PC, expH_lnRR, S_PC, S_lnRR, site_code, project_name, community_type, calendar_year)

#some preliminary histograms
theme_set(theme_bw(16))
m<-qplot(mean_change, data=divCompare, geom="histogram")+
  ggtitle("Among Treatment Change")+
  xlab(" Distance between Centriods")+
  geom_vline(xintercept = 0, size=2)

s1<-qplot(S_PC, data=divCompare, geom="histogram")+
  ggtitle("Richness Percent Change")+
  xlab("Percent Change in Richness")+
  geom_vline(xintercept = 0, size=2)

e2<-qplot(expH_PC, data=divCompare, geom="histogram")+
  ggtitle("Effective Diversity Percent Change")+
  xlab("Trt Evenness - Cont Evenness")+
  geom_vline(xintercept = 0, size=2)

s2<-qplot(S_lnRR, data=divCompare, geom="histogram")+
  ggtitle("Richness ln Response Ratio")+
  xlab("Percent Change in Richness")+
  geom_vline(xintercept = 0, size=2)

e3<-qplot(expH_lnRR, data=divCompare, geom="histogram")+
  ggtitle("Effective Diversity ln Response Ratio")+
  xlab("Trt Evenness - Cont Evenness")+
  geom_vline(xintercept = 0, size=2)

grid.arrange(m, s1, e2, s2, e3, ncol=3)


###merging with treatment information
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

# full dataset
# write.csv(ForAnalysis, "ForBayesianAnalysis_Nov2018_2.csv")


###generating treatment categories (resource, non-resource, and interactions)
###4 steps: (1) single resource, (2) single non-resource, (3) 2-way interactions, (4) 3+ way interactions
#step 1: single resource
singleResource <- ForAnalysis%>%
  select(-plot_mani)%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==1, plot_mani==1)%>%
  #set CEH Megarich nutrient values to 0 (added to all megaliths, not a treatment)
  mutate(n2=ifelse(site_code=='CEH', 0, n), p2=ifelse(site_code=='CEH', 0, p), k2=ifelse(site_code=='CEH', 0, k))%>%
  #drop lime added, as only one trt does this
  filter(other_trt!='lime added')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(n2>0, 'N', ifelse(p2>0, 'P', ifelse(k2>0, 'K', ifelse(precip<0, 'drought', ifelse(precip>0, 'irr', ifelse(CO2>0, 'CO2', 'precip_vari')))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, expH_PC, expH_lnRR, S_PC, S_lnRR, experiment_length, rrich, anpp, MAT, MAP)

#step 2: single non-resource
singleNonresource <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==0, plot_mani==1)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==1, 'burn', ifelse(mow_clip==1, 'mow_clip', ifelse(herb_removal==1, 'herb_rem', ifelse(temp>0, 'temp', ifelse(plant_trt==1, 'plant_mani', 'other'))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, expH_PC, expH_lnRR, S_PC, S_lnRR, experiment_length, rrich, anpp, MAT, MAP)

#step 3: 2-way interactions
twoWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani==2)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani==1&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&herb_removal==1, 'R*herb_rem', ifelse(resource_mani==1&temp>0, 'R*temp', ifelse(resource_mani==1&plant_trt==1, 'R*plant_mani', ifelse(resource_mani==1&other_trt!=0, 'R*other', ifelse(n>0&p>0, 'R*R', ifelse(n>0&CO2>0, 'R*R', ifelse(n>0&precip!=0, 'R*R', ifelse(p>0&k>0, 'R*R', ifelse(CO2>0&precip!=0, 'R*R', 'N*N'))))))))))))%>%
  #drop R*herb_removal (single rep)
  filter(trt_type!='R*herb_rem')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, expH_PC, expH_lnRR, S_PC, S_lnRR, experiment_length, rrich, anpp, MAT, MAP)

#step 4: 3+ way interactions
threeWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani>2, plot_mani<6)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==0&mow_clip==0&herb_removal==0&temp==0&plant_trt==0, 'all_resource', ifelse(n==0&p==0&k==0&CO2==0&precip==0, 'all_nonresource', 'both')))%>%
  #drop single all-nonresource treatment (NIN herbdiv 5NF)
  filter(trt_type!='all_nonresource')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, expH_PC, expH_lnRR, S_PC, S_lnRR, experiment_length, rrich, anpp, MAT, MAP)

#combine for analysis - one big model, 19 trt types
allAnalysis <- rbind(singleResource, singleNonresource, twoWay, threeWay)

#subset out datasets with less than 3 temporal data points;  drops GVN FACE
numPoints <- allAnalysis%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
allAnalysisAllDatasets <- allAnalysis%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)
# write.csv(allAnalysisAllDatasets, 'ForAnalysis_allAnalysisAllDatasets.csv')


#subset out treatment years 20 or less (i.e., cut off datasets at 20 years)
allAnalysis20yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<21)%>%
  select(-num_datapoints)
numPoints <- allAnalysis20yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))


#subset out treatment years 10 or less (i.e., cut off datasets at 20 years)
allAnalysis10yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<11)%>%
  select(-num_datapoints)
numPoints <- allAnalysis10yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))

# write.csv(allAnalysis10yr, 'ForAnalysis_allAnalysis10yr_pairwise_12142018.csv')

#subset out 20th or final year of all data
allAnalysisFinalYear <- allAnalysis20yr%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  filter(treatment_year==max(treatment_year))
# write.csv(allAnalysisFinalYear, 'ForAnalysis_allAnalysisFinalYear_09062018.csv')