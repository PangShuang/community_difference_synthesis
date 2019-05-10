library(vegan)
library(grid)
library(tidyverse)

#kim's wd
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#ggplot theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#load raw abundance data
trt <- read.csv('ExperimentInformation_March2019.csv')
rawAbundance <- read.csv('SpeciesRelativeAbundance_March2019.csv')%>%
  select(-X)%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  left_join(trt)



###irrigation plots calculations
###calculating bray-curtis and jaccard dissimilarities within and among treatments to get distances between treatment centroids and dispersion among replicate plots within a treatment
#make a new dataframe with just the label
exp_year=rawAbundance%>%
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
  subset=rawAbundance[rawAbundance$exp_year==as.character(exp_year$exp_year[i]),]%>%
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
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances among treatment centroids
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      braycurtis=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="braycurtis" 
  
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
  
  #calculate jaccard dissimilarities
  jd=vegdist(species[,5:ncol(species)], method="jaccard", binary=T)
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(jd, species$treatment, type="centroid")
  
  #getting distances among treatment centroids
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      jaccard=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="jaccard" 
  
  #merging back with labels to get back plot_mani
  centroid2=merge(cent_C_T, labels, by="treatment")
  
  #collecting and labeling distances to centroid from betadisper
  trt_disp=data.frame(data.frame(exp_year=exp_year$exp_year[i], 
                                 plot_id=species$plot_id,
                                 treatment=species$treatment,
                                 dist=disp$distances))%>%
    tbl_df()%>%
    group_by(exp_year, treatment)%>%
    summarize(dispersion=mean(dist))
  
  #merge together change in mean and dispersion data
  distances2<-merge(centroid2, centroid, by=c("exp_year","treatment", "plot_mani"))
  
  
  #pasting dispersions into the dataframe made for this analysis
  for.analysis=rbind(distances2, for.analysis)  
}

# write.csv(for.analysis, 'DiversityMetrics_Nov2017.csv')

rm(list=setdiff(ls(), "for.analysis"))

dist <- for.analysis%>%
  separate(exp_year, into=c('site_code', 'project_name', 'community_type', 'calendar_year'), sep="::")


#figure
irrBray <- ggplot(subset(dist, site_code=='KNZ'&community_type=='u'&treatment=='i'), aes(x=as.integer(calendar_year), y=braycurtis)) +
  geom_point(size=5) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size=2, color='black', se=F) +
  scale_y_continuous(limits=c(0,0.65)) +
  scale_x_continuous(limits=c(1999,2009)) +
  xlab('') +
  ylab('Bray-Curtis Distance') +
  annotate('text', x=1999, y=0.6, label='(a) Konza Irrigation Plots', size=10, hjust='left') +
  theme(legend.position='none')

irrJaccard <- ggplot(subset(dist, site_code=='KNZ'&community_type=='u'&treatment=='i'), aes(x=as.integer(calendar_year), y=jaccard)) +
  geom_point(size=5) +  
  geom_smooth(method = "lm", size=2, color='black', se=F) +
  scale_y_continuous(limits=c(0,0.65)) +
  scale_x_continuous(limits=c(1999,2009)) +
  xlab('') +
  ylab('Jaccard Distance') +
  annotate('text', x=1999, y=0.6, label='(b) Konza Irrigation Plots', size=10, hjust='left') +
  theme(legend.position='none')

e001Bray <- ggplot(subset(dist, site_code=='CDR'&community_type=='D'&treatment=='8'), aes(x=as.integer(calendar_year), y=braycurtis)) +
  geom_point(size=5) +  
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size=2, color='black', se=F) +
  scale_y_continuous(limits=c(0,0.95)) +
  scale_x_continuous(limits=c(1986,2004)) +
  xlab('Year of Experiment') +
  ylab('Bray-Curtis Distance') +
  annotate('text', x=1986, y=0.9, label='(c) Cedar Creek e001', size=10, hjust='left') +
  theme(legend.position='none')

e001Jaccard <- ggplot(subset(dist, site_code=='CDR'&community_type=='D'&treatment=='8'), aes(x=as.integer(calendar_year), y=jaccard)) +
  geom_point(size=5) +  
  geom_smooth(method = "lm", size=2, color='black', se=F) +
  scale_y_continuous(limits=c(0,0.95)) +
  scale_x_continuous(limits=c(1986,2004)) +
  xlab('Year of Experiment') +
  ylab('Jaccard Distance') +
  annotate('text', x=1986, y=0.9, label='(d) Cedar Creek e001', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(irrBray, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(irrJaccard, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(e001Bray, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(e001Jaccard, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export as 1500x1500


