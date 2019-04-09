################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among significant responses by treatment or treatment category for richness and compositional responses.
##
##  Author: Kimberly La Pierre
##  Date created: March 22, 2018
################################################################################


library(tidyverse)


#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#significant lines
correTotal <- read.csv('treatment_response_shape_classification_stdtimebytrt_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, shape_category)%>%
  group_by(variable, trt_type)%>%
  summarise(total_possible=length(trt_type))%>%
  ungroup()

correProp <- read.csv('treatment_response_shape_classification_stdtimebytrt_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, shape_category)%>%
  mutate(sig_line=ifelse(shape_category==0, 0, 1))%>%
  group_by(variable, trt_type)%>%
  summarise(count_sig=sum(sig_line))%>%
  ungroup()%>%
  left_join(correTotal)%>%
  mutate(count_nonsig=total_possible-count_sig)

###treatment type
#richness responses
correPropRich <- correProp%>%filter(variable=='richness')
prop.test(x=as.matrix(correPropRich[,c('count_sig', 'count_nonsig')]), alternative='two.sided')
#compositional responses
correPropComp <- correProp%>%filter(variable=='mean')
prop.test(x=as.matrix(correPropComp[,c('count_sig', 'count_nonsig')]), alternative='two.sided')


###treatment category
#significant lines
correTotalCategory <- read.csv('treatment_response_shape_classification_stdtimebytrt_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, shape_category)%>%
  mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('burn','herb_removal','mow_clip','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
  group_by(variable, trt_category)%>%
  summarise(total_possible=length(trt_category))%>%
  ungroup()

correPropCategory <- read.csv('treatment_response_shape_classification_stdtimebytrt_04072019.csv')%>%
  select(variable, trt_type, site_code, project_name, community_type, treatment, shape_category)%>%
  mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('burn','herb_removal','mow_clip','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
  mutate(sig_line=ifelse(shape_category==0, 0, 1))%>%
  group_by(variable, trt_category)%>%
  summarise(count_sig=sum(sig_line))%>%
  ungroup()%>%
  left_join(correTotalCategory)%>%
  mutate(count_nonsig=total_possible-count_sig)

#richness responses
correPropCategoryRich <- correPropCategory%>%filter(variable=='richness')
prop.test(x=as.matrix(correPropCategoryRich[,c('count_sig', 'count_nonsig')]), alternative='two.sided')
#compositional responses
correPropCategoryComp <- correPropCategory%>%filter(variable=='mean')
prop.test(x=as.matrix(correPropCategoryComp[,c('count_sig', 'count_nonsig')]), alternative='two.sided')
