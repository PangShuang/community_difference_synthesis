################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among significant responses by treatment or treatment category for richness and compositional responses.
##
##  Author: Kimberly La Pierre
##  Date created: March 22, 2018
################################################################################


library(tidyverse)


#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#without intercept-only significance
correProp <- read.csv('proportions significant comp difference_by trt type.csv')%>%
  filter(total_possible>4)%>%
  mutate(no_mean_change=total_possible-mean_change_nointercept, no_richness_change=total_possible-richness_change_nointercept)

###treatment type
#richness responses
prop.test(x=as.matrix(correProp[,c('richness_change_nointercept', 'no_richness_change')]), alternative='two.sided')
#compositional responses
prop.test(x=as.matrix(correProp[,c('mean_change_nointercept', 'no_mean_change')]), alternative='two.sided')


###treatment category
correPropCategory <- correProp%>%
  group_by(resource_category)%>%
  summarise(sum_richness_change_nointercept=sum(richness_change_nointercept), sum_no_richness_change=sum(no_richness_change), sum_mean_change_nointercept=sum(mean_change_nointercept), sum_no_mean_change=sum(no_mean_change))%>%
  ungroup()

#richness responses
prop.test(x=as.matrix(correPropCategory[,c('sum_richness_change_nointercept', 'sum_no_richness_change')]), alternative='two.sided')
#compositional responses
prop.test(x=as.matrix(correPropCategory[,c('sum_mean_change_nointercept', 'sum_no_mean_change')]), alternative='two.sided')
