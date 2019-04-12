################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among significant responses by treatment or treatment category for richness and compositional responses.
##
##  Author: Kimberly Komatsu
##  Date created: March 22, 2018
################################################################################
library(RVAideMemoire)
library(tidyverse)


#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#significant lines
correTotal <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
  group_by(variable, trt_type)%>%
  summarise(total_possible=length(trt_type))%>%
  ungroup()

correProp <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
  mutate(sig_line=ifelse(N01_shape==0, 0, 1))%>%
  group_by(variable, trt_type)%>%
  summarise(count_sig=sum(sig_line))%>%
  ungroup()%>%
  left_join(correTotal)%>%
  mutate(count_nonsig=total_possible-count_sig)%>%
  filter(trt_type %in% c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp'))%>%
  mutate(proportion=count_sig/total_possible)

###treatment type -- single factor treatments only
#richness responses
correPropRich <- correProp%>%filter(variable=='richness')
prop.test(x=as.matrix(correPropRich[,c('count_sig', 'total_possible')]), alternative='two.sided')
#compositional responses
correPropComp <- correProp%>%filter(variable=='mean')
prop.test(x=as.matrix(correPropComp[,c('count_sig', 'total_possible')]), alternative='two.sided')


###treatment type -- single factor treatments only
#richness responses
correPropRich <- correProp%>%filter(variable=='richness')%>%select(trt_type,count_sig,count_nonsig)
trtType <- as.character(correPropRich$trt_type)
correPropRichMatrix <- as.matrix(correPropRich[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp')), colnames(c('count_sig','count_nonsig'))), byrow=T)
fisher.test(correPropRichMatrix, workspace = 2e8)
fisher.multcomp(correPropRichMatrix)
#compositional responses
correPropComp <- correProp%>%filter(variable=='mean')%>%select(trt_type,count_sig,count_nonsig)
trtType <- as.character(correPropComp$trt_type)
correPropCompMatrix <- as.matrix(correPropComp[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp')), colnames(c('count_sig','count_nonsig'))), byrow=T)
fisher.test(correPropCompMatrix, workspace = 2e8)
fisher.multcomp(correPropCompMatrix)



###treatment category
#significant lines
correTotalCategory <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
  mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('temp','herb_removal','mow_clip','plant_mani','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
  group_by(variable, trt_category)%>%
  summarise(total_possible=length(trt_category))%>%
  ungroup()

correPropCategory <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
  mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('temp','herb_removal','mow_clip','plant_mani','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
  mutate(sig_line=ifelse(N01_shape==0, 0, 1))%>%
  group_by(variable, trt_category)%>%
  summarise(count_sig=sum(sig_line))%>%
  ungroup()%>%
  left_join(correTotalCategory)%>%
  mutate(count_nonsig=total_possible-count_sig)%>%
  mutate(proportion=count_sig/total_possible)

#richness responses
correPropCategoryRich <- correPropCategory%>%filter(variable=='richness')
prop.test(x=as.matrix(correPropCategoryRich[,c('count_sig', 'total_possible')]), alternative='two.sided')
#compositional responses
correPropCategoryComp <- correPropCategory%>%filter(variable=='mean')
prop.test(x=as.matrix(correPropCategoryComp[,c('count_sig', 'total_possible')]), alternative='two.sided')


###treatment type -- single factor treatments only
#richness responses
correPropRich <- correPropCategory%>%filter(variable=='richness')%>%select(trt_category,count_sig,count_nonsig)
trtType <- as.character(correPropRich$trt_category)
correPropRichMatrix <- as.matrix(correPropRich[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('NxN','RxN','RxR','RxRxR','single_nonresource','single_resource','threeway')), colnames(c('count_sig','count_nonsig'))), byrow=T)
fisher.test(correPropRichMatrix, workspace = 2e8)
fisher.multcomp(correPropRichMatrix)
#compositional responses
correPropComp <- correPropCategory%>%filter(variable=='mean')%>%select(trt_category,count_sig,count_nonsig)
trtType <- as.character(correPropComp$trt_category)
correPropCompMatrix <- as.matrix(correPropComp[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('NxN','RxN','RxR','RxRxR','single_nonresource','single_resource','threeway')), colnames(c('count_sig','count_nonsig'))), byrow=T)
fisher.test(correPropCompMatrix, workspace = 2e8)
fisher.multcomp(correPropCompMatrix)
