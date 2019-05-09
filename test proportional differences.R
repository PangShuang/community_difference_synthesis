################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among significant responses by treatment or treatment category for richness and compositional responses.
##
##  Author: Kimberly Komatsu
##  Date created: March 22, 2018
################################################################################

library(grid)
library(RVAideMemoire)
library(tidyverse)


#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')



theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))




####all data---------------
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



###figure 2
rich <- ggplot(data=subset(correPropRich, trt_type!='other_resource'&trt_type!='other_nonresource'), aes(x=trt_type, y=proportion)) +
  geom_bar(stat='identity')  +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'precip_vari', 'N', 'P', 'temp', 'mow_clip', 'herb_removal', 'plant_mani'),
                   labels=c('CO2', 'drought', 'irrigation', 'precip. vari.', 'nitrogen', 'phosphorus', 'temperature', 'mow', 'herbivore rem.', 'plant manip.')) +
  xlab('Treatment Type') + ylab('Proportion') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  annotate('text', x=0.5, y=0.4, label='(a) Richness Response', size=12, hjust='left')

comp <- ggplot(data=subset(correPropComp, trt_type!='other_resource'&trt_type!='other_nonresource'), aes(x=trt_type, y=proportion)) +
  geom_bar(stat='identity')  +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'precip_vari', 'N', 'P', 'temp', 'mow_clip', 'herb_removal', 'plant_mani'),
                   labels=c('CO2', 'drought', 'irrigation', 'precip. vari.', 'nitrogen', 'phosphorus', 'temperature', 'mow', 'herbivore rem.', 'plant manip.')) +
  xlab('Treatment Type') + ylab('') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  annotate('text', x=0.5, y=0.4, label='(b) Composition Response', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(rich, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(comp, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400x1000





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



###figure 3
rich <- ggplot(data=correPropCategoryRich, aes(x=trt_category, y=proportion)) +
  geom_bar(stat='identity')  +
  scale_x_discrete(limits=c('single_resource', 'single_nonresource', 'RxR', 'NxN', 'RxN', 'RxRxR', 'threeway'),
                   labels=c('R', 'N', 'R*R', 'N*N', 'R*N', 'R*R*R', '3+')) +
  xlab('Treatment Category') + ylab('Proportion') +
  annotate('text', x=0.5, y=0.7, label='(a) Richness Response', size=12, hjust='left') +
  annotate('text', x=1, y=0.20, label='a', size=12, hjust='center') +
  annotate('text', x=2, y=0.14, label='a', size=12, hjust='center') +
  annotate('text', x=3, y=0.29, label='ab', size=12, hjust='center') +
  annotate('text', x=4, y=0.18, label='ab', size=12, hjust='center') +
  annotate('text', x=5, y=0.20, label='ab', size=12, hjust='center') +
  annotate('text', x=6, y=0.59, label='c', size=12, hjust='center') +
  annotate('text', x=7, y=0.38, label='b', size=12, hjust='center')

comp <- ggplot(data=correPropCategoryComp, aes(x=trt_category, y=proportion)) +
  geom_bar(stat='identity')  +
  scale_x_discrete(limits=c('single_resource', 'single_nonresource', 'RxR', 'NxN', 'RxN', 'RxRxR', 'threeway'),
                   labels=c('R', 'N', 'R*R', 'N*N', 'R*N', 'R*R*R', '3+')) +
  xlab('Treatment Category') + ylab('') +
  annotate('text', x=0.5, y=0.7, label='(b) Composition Response', size=12, hjust='left') +
  annotate('text', x=1, y=0.32, label='a', size=12, hjust='center') +
  annotate('text', x=2, y=0.22, label='a', size=12, hjust='center') +
  annotate('text', x=3, y=0.55, label='bc', size=12, hjust='center') +
  annotate('text', x=4, y=0.26, label='abc', size=12, hjust='center') +
  annotate('text', x=5, y=0.33, label='ab', size=12, hjust='center') +
  annotate('text', x=6, y=0.66, label='c', size=12, hjust='center') +
  annotate('text', x=7, y=0.53, label='bc', size=12, hjust='center')

pushViewport(viewport(layout=grid.layout(1,2)))
print(rich, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(comp, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400x1000


#####10 year or longer-----------------
# 
# #significant lines
# correTotal <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
#   filter(experiment_length>9)%>%
#   select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
#   group_by(variable, trt_type)%>%
#   summarise(total_possible=length(trt_type))%>%
#   ungroup()
# 
# correProp <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
#   filter(experiment_length>9)%>%
#   select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
#   mutate(sig_line=ifelse(N01_shape==0, 0, 1))%>%
#   group_by(variable, trt_type)%>%
#   summarise(count_sig=sum(sig_line))%>%
#   ungroup()%>%
#   left_join(correTotal)%>%
#   mutate(count_nonsig=total_possible-count_sig)%>%
#   filter(trt_type %in% c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp'))%>%
#   mutate(proportion=count_sig/total_possible)
# 
# ###treatment type -- single factor treatments only
# #richness responses
# correPropRich <- correProp%>%filter(variable=='richness')
# prop.test(x=as.matrix(correPropRich[,c('count_sig', 'total_possible')]), alternative='two.sided')
# #compositional responses
# correPropComp <- correProp%>%filter(variable=='mean')
# prop.test(x=as.matrix(correPropComp[,c('count_sig', 'total_possible')]), alternative='two.sided')
# 
# 
# ###treatment type -- single factor treatments only
# #richness responses
# correPropRich <- correProp%>%filter(variable=='richness')%>%select(trt_type,count_sig,count_nonsig)
# trtType <- as.character(correPropRich$trt_type)
# correPropRichMatrix <- as.matrix(correPropRich[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp')), colnames(c('count_sig','count_nonsig'))), byrow=T)
# fisher.test(correPropRichMatrix, workspace = 2e8)
# fisher.multcomp(correPropRichMatrix)
# #compositional responses
# correPropComp <- correProp%>%filter(variable=='mean')%>%select(trt_type,count_sig,count_nonsig)
# trtType <- as.character(correPropComp$trt_type)
# correPropCompMatrix <- as.matrix(correPropComp[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('CO2','drought','herb_removal','irr','mow_clip','N','other_nonresource','other_resource','P','plant_mani','precip_vari','temp')), colnames(c('count_sig','count_nonsig'))), byrow=T)
# fisher.test(correPropCompMatrix, workspace = 2e8)
# fisher.multcomp(correPropCompMatrix)
# 
# 
# 
# ###treatment category
# #significant lines
# correTotalCategory <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
#   filter(experiment_length>9)%>%
#   select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
#   mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('temp','herb_removal','mow_clip','plant_mani','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
#   group_by(variable, trt_category)%>%
#   summarise(total_possible=length(trt_category))%>%
#   ungroup()
# 
# correPropCategory <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
#   filter(experiment_length>9)%>%
#   select(variable, trt_type, site_code, project_name, community_type, treatment, N01_shape)%>%
#   mutate(trt_category=ifelse(trt_type %in% c('threeway','RxRxR','RxR','RxN','NxN'), as.character(trt_type), ifelse(trt_type %in% c('temp','herb_removal','mow_clip','plant_mani','other_nonresource'), 'single_nonresource', 'single_resource')))%>%
#   mutate(sig_line=ifelse(N01_shape==0, 0, 1))%>%
#   group_by(variable, trt_category)%>%
#   summarise(count_sig=sum(sig_line))%>%
#   ungroup()%>%
#   left_join(correTotalCategory)%>%
#   mutate(count_nonsig=total_possible-count_sig)%>%
#   mutate(proportion=count_sig/total_possible)
# 
# #richness responses
# correPropCategoryRich <- correPropCategory%>%filter(variable=='richness')
# prop.test(x=as.matrix(correPropCategoryRich[,c('count_sig', 'total_possible')]), alternative='two.sided')
# #compositional responses
# correPropCategoryComp <- correPropCategory%>%filter(variable=='mean')
# prop.test(x=as.matrix(correPropCategoryComp[,c('count_sig', 'total_possible')]), alternative='two.sided')
# 
# 
# ###treatment type -- single factor treatments only
# #richness responses
# correPropRich <- correPropCategory%>%filter(variable=='richness')%>%select(trt_category,count_sig,count_nonsig)
# trtType <- as.character(correPropRich$trt_category)
# correPropRichMatrix <- as.matrix(correPropRich[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('NxN','RxN','RxR','RxRxR','single_nonresource','single_resource','threeway')), colnames(c('count_sig','count_nonsig'))), byrow=T)
# fisher.test(correPropRichMatrix, workspace = 2e8)
# fisher.multcomp(correPropRichMatrix)
# #compositional responses
# correPropComp <- correPropCategory%>%filter(variable=='mean')%>%select(trt_category,count_sig,count_nonsig)
# trtType <- as.character(correPropComp$trt_category)
# correPropCompMatrix <- as.matrix(correPropComp[,c('count_sig', 'count_nonsig')], ncol=2, dimnames=list(rownames(c('NxN','RxN','RxR','RxRxR','single_nonresource','single_resource','threeway')), colnames(c('count_sig','count_nonsig'))), byrow=T)
# fisher.test(correPropCompMatrix, workspace = 2e8)
# fisher.multcomp(correPropCompMatrix)
