################################################################################
##  time_comparisons.R: Calculating changes and differences between treatment and control plots between first and last years of treatments for only the treatments that exhibited significant curvature.
##
##  Author: Kimberly Komatsu
##  Date created: March 18, 2019
################################################################################

library(grid)
library(tidyverse)
library(codyn)


#kim desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###note: must run utilities file before running updates to codyn functions below
###functions---------------------------
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

#update SERGL function to get lnRR
SERGL2 <- function (df, species.var, abundance.var, abundance.var2) 
{
  out <- c(species.var, "rank", "rank2", abundance.var, abundance.var2)
  out <- unique(df[!(names(df) %in% out)])
  s_t1 <- S(df[[abundance.var]])
  e_t1 <- Evar(as.numeric(df[[abundance.var]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- Evar(as.numeric(df[[abundance.var2]]))
  delta_s <- log(s_t2/s_t1)
  delta_e <- (e_t2 - e_t1)
  df$gain <- ifelse(df[[abundance.var]] == 0, 1, 0)
  df$loss <- ifelse(df[[abundance.var2]] == 0, 1, 0)
  gain <- sum(df$gain)/nrow(df)
  loss <- sum(df$loss)/nrow(df)
  delta_r <- mean(abs(df[["rank"]] - df[["rank2"]]))/nrow(df)
  metrics <- data.frame(richness_change = delta_s, evenness_change = delta_e, 
                        rank_change = delta_r, gains = gain, losses = loss)
  return(cbind(out, metrics))
}

#update RAC_change function to get lnRR
RAC_lnRR_change <- function (df, time.var, species.var, abundance.var, replicate.var = NULL, 
          reference.time = NULL) 
{
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  by <- c(replicate.var)
  allsp <- split_apply_combine(df, by, FUN = fill_species, 
                               species.var, abundance.var)
  by <- c(time.var, replicate.var)
  rankdf <- split_apply_combine(allsp, by, FUN = add_ranks, 
                                abundance.var)
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = "")
  split_by <- c(replicate.var)
  merge_to <- !(names(rankdf) %in% split_by)
  if (is.null(reference.time)) {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = species.var, suffixes = c("", 
                                                          "2"))
      f <- factor(cross[[cross.var]])
      f2 <- factor(cross[[cross.var2]], levels = levels(f))
      idx <- (as.integer(f2) - as.integer(f)) == 1
      cross[idx, ]
    })
  }
  else {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[x[[time.var]] != reference.time, merge_to]
      x <- x[x[[time.var]] == reference.time, ]
      merge(x, y, by = species.var, suffixes = c("", "2"))
    })
  }
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = "")
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  idx <- ranktog[[abundance.var]] != 0 | ranktog[[abundance.var2]] != 
    0
  ranktog <- ranktog[idx, ]
  by <- c(replicate.var, time.var)
  output <- split_apply_combine(ranktog, by, FUN = SERGL2, species.var, 
                                abundance.var, abundance.var2)
  if (any(is.na(output$evenness_change))) 
    warning(paste0("evenness_change values contain NAs because there are plots", 
                   " with only one species"))
  output_order <- c(time.var, paste(time.var, 2, sep = ""), 
                    replicate.var, "richness_change", "evenness_change", 
                    "rank_change", "gains", "losses")
  return(output[intersect(output_order, names(output))])
}

#update SERSp to get lnRR
SERSp2 <- function (df, species.var, abundance.var, abundance.var2) 
{
  out <- c(species.var, "rank", "rank2", abundance.var, abundance.var2)
  out <- unique(df[!(names(df) %in% out)])
  df <- subset(df, df[[abundance.var]] != 0 | df[[abundance.var2]] != 
                 0)
  s_r1 <- S(df[[abundance.var]])
  e_r1 <- Evar(as.numeric(df[[abundance.var]]))
  s_r2 <- S(df[[abundance.var2]])
  e_r2 <- Evar(as.numeric(df[[abundance.var2]]))
  sdiff <- log(s_r2/s_r1)
  ediff <- e_r2 - e_r1
  spdiff <- df[df[[abundance.var]] == 0 | df[[abundance.var2]] == 
                 0, ]
  spdiffc <- nrow(spdiff)/nrow(df)
  rank_diff <- mean(abs(df[["rank"]] - df[["rank2"]]))/nrow(df)
  metrics <- data.frame(richness_diff = sdiff, evenness_diff = ediff, 
                        rank_diff = rank_diff, species_diff = spdiffc)
  return(cbind(out, metrics))
}

#update RAC_difference function to get lnRR
RAC_lnRR_difference <- function (df, time.var = NULL, species.var, abundance.var, replicate.var, 
          treatment.var = NULL, pool = FALSE, block.var = NULL, reference.treatment = NULL) 
{
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  if (pool) {
    rankdf <- pool_replicates(df, time.var, species.var, 
                              abundance.var, replicate.var, treatment.var)
  }
  else {
    by <- c(block.var, time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_species, 
                                 species.var, abundance.var)
    by <- c(block.var, time.var, treatment.var, replicate.var)
    rankdf <- split_apply_combine(allsp, by, FUN = add_ranks, 
                                  abundance.var)
  }
  if (!is.null(block.var)) {
    cross.var <- treatment.var
  }
  else if (pool) {
    cross.var <- treatment.var
  }
  else {
    cross.var <- replicate.var
  }
  to_ordered = is.factor(rankdf[[cross.var]]) & !is.ordered(rankdf[[cross.var]])
  if (to_ordered) {
    class(rankdf[[cross.var]]) <- c("ordered", class(rankdf[[cross.var]]))
  }
  split_by <- c(block.var, time.var)
  merge_to <- !(names(rankdf) %in% split_by)
  cross.var2 <- paste(cross.var, 2, sep = "")
  if (is.null(reference.treatment)) {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = species.var, suffixes = c("", 
                                                          "2"))
      idx <- cross[[cross.var]] < cross[[cross.var2]]
      cross[idx, ]
    })
  }
  else {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[x[[treatment.var]] != reference.treatment, 
             merge_to]
      x <- x[x[[treatment.var]] == reference.treatment, 
             ]
      merge(x, y, by = species.var, suffixes = c("", "2"))
    })
  }
  if (to_ordered) {
    x <- class(ranktog[[cross.var]])
    class(ranktog[[cross.var]]) <- x[x != "ordered"]
    class(ranktog[[cross.var2]]) <- x[x != "ordered"]
  }
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = "")
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  idx <- ranktog[[abundance.var]] != 0 | ranktog[[abundance.var2]] != 
    0
  ranktog <- ranktog[idx, ]
  split_by <- c(block.var, time.var, cross.var, cross.var2)
  output <- split_apply_combine(ranktog, split_by, FUN = SERSp2, 
                                species.var, abundance.var, abundance.var2)
  if (any(is.na(output$evenness_diff))) 
    warning(paste0("evenness_diff values contain NAs because there are plots", 
                   " with only one species"))
  output_order <- c(time.var, block.var, replicate.var, paste(replicate.var, 
                                                              2, sep = ""), treatment.var, paste(treatment.var, 2, 
                                                                                                 sep = ""), "richness_diff", "evenness_diff", "rank_diff", 
                    "species_diff")
  return(output[intersect(output_order, names(output))])
}


###data----------------
#import relative species abundance data
alldata <- read.csv("SpeciesRelativeAbundance_March2019.csv")%>%
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
expinfo <- read.csv("ExperimentInformation_March2019.csv")%>%
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
    select(plot_id, treatment, treatment2)%>%
    unique()
  
  #calculating richness differences
  richnessDifference <- RAC_lnRR_difference(subset, time.var = 'calendar_year', species.var = 'genus_species', abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment2', reference.treatment='TRUECONTROL')%>%
    group_by(calendar_year, treatment2, treatment22)%>%
    summarise(richness_difference=mean(richness_diff))%>%
    ungroup()
  
  #calculating composition difference and abs(dispersion difference)
  difference <- multivariate_difference(subset, time.var = 'calendar_year', species.var = 'genus_species', abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
    left_join(richnessDifference)%>%
    rename(treatment=treatment22)%>%
    select(-treatment2, -trt_greater_disp)%>%
    mutate(site_project_comm=exp_year$site_project_comm[i])
  
  #calculating richness change
  richnessChange <- RAC_lnRR_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', reference.time=min(subset$calendar_year))%>%
    left_join(treatments)%>%
    group_by(calendar_year, calendar_year2, treatment2)%>%
    summarise(richness_difference=mean(richness_change))%>%
    ungroup()
  
  #calculating composition change and abs(dispersion difference)
  change <- multivariate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=min(subset$calendar_year))%>%
    left_join(richnessChange)%>%
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

richnessDifferenceKBSuntilled <- RAC_lnRR_difference(kbsUntilled, time.var = 'calendar_year', species.var = 'genus_species', abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment2', reference.treatment='TRUECONTROL')%>%
  group_by(calendar_year, treatment2, treatment22)%>%
  summarise(richness_difference=mean(richness_diff))%>%
  ungroup()

differenceKBSuntilled <- multivariate_difference(kbsUntilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
  left_join(richnessDifferenceKBSuntilled)%>%
  rename(treatment=treatment22)%>%
  select(-treatment2, -trt_greater_disp)%>%
  mutate(site_project_comm='KBS::T7::0')

richnessChangeKBSuntilled <- RAC_lnRR_change(kbsUntilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', reference.time=1989)%>%
  left_join(treatments)%>%
  group_by(calendar_year, calendar_year2, treatment2)%>%
  summarise(richness_difference=mean(richness_change))%>%
  ungroup()

changeKBSuntilled <- multivariate_change(kbsUntilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=1989)%>%
  left_join(richnessChangeKBSuntilled)%>%
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

richnessDifferenceKBSuntilled <- RAC_lnRR_difference(kbsTilled, time.var = 'calendar_year', species.var = 'genus_species', abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment2', reference.treatment='TRUECONTROL')%>%
  group_by(calendar_year, treatment2, treatment22)%>%
  summarise(richness_difference=mean(richness_diff))%>%
  ungroup()

differenceKBStilled <- multivariate_difference(kbsTilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
  left_join(richnessDifferenceKBSuntilled)%>%
  rename(treatment=treatment22)%>%
  select(-treatment2, -trt_greater_disp)%>%
  mutate(site_project_comm='KBS::T7::0')

richnessChangeKBStilled <- RAC_lnRR_change(kbsTilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', reference.time=1990)%>%
  left_join(treatments)%>%
  group_by(calendar_year, calendar_year2, treatment2)%>%
  summarise(richness_difference=mean(richness_change))%>%
  ungroup()

changeKBStilled <- multivariate_change(kbsTilled, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=1990)%>%
  left_join(richnessChangeKBStilled)%>%
  rename(treatment=treatment2)%>%
  mutate(site_project_comm='KBS::T7::0')

for.analysis.difference=rbind(for.analysis.difference, differenceKBStilled, differenceKBSuntilled) 
# write.csv(for.analysis.difference, 'time_comparisons_difference.csv')
for.analysis.change=rbind(for.analysis.change, changeKBStilled, changeKBSuntilled) 
# write.csv(for.analysis.difference, 'time_comparisons_change.csv')

###final numbers for difference and change
curveShape <- read.csv('stdtimebytrt_shape_classification_04072019.csv')%>%
  select(variable, site_code, project_name, community_type, treatment, N01_shape)%>%
  spread(key=variable, value=N01_shape)%>%
  rename(mean_shape=mean, richness_shape=richness)

#compositional difference
differenceAllMean <- for.analysis.difference%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(mean_shape==7|mean_shape==8)%>%
  group_by(site_project_comm)%>%
  mutate(comparison=ifelse(calendar_year==min(calendar_year), 'trt-ctl_first', 'trt-ctl_last'))%>%
  ungroup()%>%
  rename(value=composition_diff)%>%
  select(comparison, site_project_comm, treatment, calendar_year, value)

changeTrtMean <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  filter(mean_shape==7|mean_shape==8)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

expToUseMean <- changeTrtMean%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(site_project_comm)%>%
  unique()

changeCtlMean <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  full_join(curveShape)%>%
  mutate(mean_shape=ifelse(treatment=='TRUECONTROL', 999, mean_shape))%>%
  filter(mean_shape==999)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(site_project_comm %in% expToUseMean$site_project_comm)

changeAllMean <- rbind(changeTrtMean, changeCtlMean)%>%
  mutate(comparison=ifelse(treatment=='TRUECONTROL', 'ctl_first-last', 'trt_first-last'))%>%
  rename(value=composition_change)%>%
  select(site_project_comm, treatment, comparison, calendar_year, value)

comparisonsMean <- rbind(differenceAllMean, changeAllMean)


#richness difference
differenceAllRich <- for.analysis.difference%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(richness_shape==7|richness_shape==8)%>%
  group_by(site_project_comm)%>%
  mutate(comparison=ifelse(calendar_year==min(calendar_year), 'trt-ctl_first', 'trt-ctl_last'))%>%
  ungroup()%>%
  rename(value=composition_diff)%>%
  select(comparison, site_project_comm, treatment, calendar_year, value)

changeTrtRich <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  right_join(curveShape)%>%
  filter(richness_shape==7|richness_shape==8)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

expToUseRich <- changeTrtRich%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(site_project_comm)%>%
  unique()

changeCtlRich <- for.analysis.change%>%
  separate(site_project_comm, c('site_code', 'project_name', 'community_type'), sep='::')%>%
  full_join(curveShape)%>%
  mutate(richness_shape=ifelse(treatment=='TRUECONTROL', 999, richness_shape))%>%
  filter(richness_shape==999)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  filter(site_project_comm %in% expToUseRich$site_project_comm)

changeAllRich <- rbind(changeTrtRich, changeCtlRich)%>%
  mutate(comparison=ifelse(treatment=='TRUECONTROL', 'ctl_first-last', 'trt_first-last'))%>%
  rename(value=composition_change)%>%
  select(site_project_comm, treatment, comparison, calendar_year, value)

comparisonsRich <- rbind(differenceAllRich, changeAllRich)

###figures
meanComparisonFig <- ggplot(data=barGraphStats(data=comparisonsMean, variable="value", byFactorNames=c("comparison")), aes(x=comparison, y=mean)) +
  geom_bar(stat='identity', color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  ylab('') +
  scale_x_discrete(limits=c('trt-ctl_first', 'ctl_first-last', 'trt_first-last'),
                   labels=c('initial diff', 'ctl change', 'trt change')) +
  theme(axis.title.x=element_blank()) +
  annotate('text', x=0.5, y=0.65, label='(b) Composition Response', size=12, hjust='left') +
  annotate('text', x=1, y=0.29, label='a', size=12, hjust='center') +
  annotate('text', x=2, y=0.52, label='b', size=12, hjust='center') +
  annotate('text', x=3, y=0.57, label='b', size=12, hjust='center')
  

richComparisonFig <- ggplot(data=barGraphStats(data=comparisonsRich, variable="value", byFactorNames=c("comparison")), aes(x=comparison, y=mean)) +
  geom_bar(stat='identity', color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  ylab('Response Magnitude') +
  scale_x_discrete(limits=c('trt-ctl_first', 'ctl_first-last', 'trt_first-last'),
                   labels=c('initial diff', 'ctl change', 'trt change')) +
  theme(axis.title.x=element_blank()) +
  annotate('text', x=0.5, y=0.65, label='(a) Richness Response', size=12, hjust='left') +
  annotate('text', x=1, y=0.27, label='a', size=12, hjust='center') +
  annotate('text', x=2, y=0.60, label='b', size=12, hjust='center') +
  annotate('text', x=3, y=0.55, label='b', size=12, hjust='center')

pushViewport(viewport(layout=grid.layout(1,2)))
print(richComparisonFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanComparisonFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400 x 1200

summary(aov(value~comparison, data=subset(comparisonsMean, comparison!='trt-ctl_last')))
pairwise.t.test(comparisonsMean$value, comparisonsMean$comparison)

summary(aov(value~comparison, data=subset(comparisonsRich, comparison!='trt-ctl_last')))
pairwise.t.test(comparisonsRich$value, comparisonsRich$comparison)

