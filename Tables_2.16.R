#################################################
## READ IN DATA AND DEFINE CLUSTERING FUNCTION ##
#################################################

library(arm)
library(stargazer)
library(pastecs)
library(xtable)
library(sandwich)
library(dplyr)
library(knitr)
library(ggplot2)

#Clustering Function

clx <- function(fm, dfcw, cluster){
  # reweighting the var-cov matrix for the within model
  library(sandwich);library(lmtest)
  M <- length(unique(cluster))   
  N <- length(cluster)           
  K <- fm$rank                        
  dfc <- (M/(M-1))*((N-1)/(N-K))  
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)*dfcw
  return(vcovCL) }


master <- read.csv('/Users/nataliecarlson/Desktop/Microfinance/Data/master_v2_reg.csv', header=TRUE, stringsAsFactors=FALSE)
always_ngo <- filter(master, always_ngo==1)
always_nbfc <- filter(master, always_nbfc==1)
switch <- filter(master, switch==1)

############
## TABLES ##
############

#TABLE 1: PRODUCTIVITY

productivity1 <- lm(productivity_std ~ meaningless_mean, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity1)))])
cov1   <- clx(productivity1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

productivity2 <- lm(productivity_std ~ meaningless_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity2)))])
cov1   <- clx(productivity2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

productivity3 <- lm(productivity_std ~ meaningless_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity3)))])
cov1   <- clx(productivity3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

productivity4 <- lm(productivity_std ~ meaningless_ngo + meaningless_nbfc + meaningless_switch + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity4)))])
cov1   <- clx(productivity4, 1, master.x$mfiuniqueid)
cluster.se4  <- sqrt(diag(cov1))

stargazer(productivity1, productivity2, productivity3, productivity4, title="Productivity Regressions", dep.var.labels = "Productivity measure (standardized)", covariate.labels = c("I sometimes feel my job is meaningless", "I sometimes feel my job is meaningless x NGO", "I sometimes feel my job is meaningless x NBFC", "I sometimes feel my job is meaningless x Switch",  "Always NGO", "Always NBFC"), omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4))



#TABLE 2: WHAT DRIVES MEANING

meaning1 <- lm(meaningless_mean ~ interest_work_branch + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning1)))])
cov1   <- clx(meaning1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

meaning2 <- lm(meaningless_mean ~ community_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning2)))])
cov1   <- clx(meaning2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

meaning3 <- lm(meaningless_mean ~ like_coworkers_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning3)))])
cov1   <- clx(meaning3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

meaning4 <- lm(meaningless_mean ~ dictator_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning4)))])
cov1   <- clx(meaning4, 1, master.x$mfiuniqueid)
cluster.se4   <- sqrt(diag(cov1))

meaning5 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning5)))])
cov1   <- clx(meaning5, 1, master.x$mfiuniqueid)
cluster.se5   <- sqrt(diag(cov1))

meaning6 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning6)))])
cov1   <- clx(meaning6, 1, master.x$mfiuniqueid)
cluster.se6   <- sqrt(diag(cov1))

meaning7 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + bm_total_comp_std + lo_total_comp_std + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning7)))])
cov1   <- clx(meaning7, 1, master.x$mfiuniqueid)
cluster.se7   <- sqrt(diag(cov1))

stargazer(meaning1, meaning2, meaning3, meaning4, meaning5, meaning6, meaning7,  title="Origins of Meaning", dep.var.labels = "I sometimes feel my job is meaningless", covariate.labels = c("Took job because of interest in work","'Our company contributes positively to the community'", "'I like the people I work with'", "Altruism (avg. given in dictator game)", "Productivity (std)","BM total compensation (std)", "LO total compensation (std)",  "Always NGO", "Always NBFC"),omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4, cluster.se5, cluster.se6, cluster.se7))



#Meaning 2: Subsamples

meaning1 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + bm_total_comp_std + lo_total_comp_std + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(meaning1)))])
cov1   <- clx(meaning1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

meaning2 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + bm_total_comp_std + lo_total_comp_std + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, always_ngo, na.action = na.omit)
master.x <- na.omit(always_ngo[ , c("mfiuniqueid", all.vars(formula(meaning2)))])
cov1   <- clx(meaning2, 1, master.x$mfiuniqueid)
cluster.se2  <- sqrt(diag(cov1))

meaning3 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + bm_total_comp_std + lo_total_comp_std + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, always_nbfc, na.action = na.omit)
master.x <- na.omit(always_nbfc[ , c("mfiuniqueid", all.vars(formula(meaning3)))])
cov1   <- clx(meaning3, 1, master.x$mfiuniqueid)
cluster.se3  <- sqrt(diag(cov1))

meaning4 <- lm(meaningless_mean ~ interest_work_branch + community_mean + like_coworkers_mean + dictator_mean + productivity_std + bm_total_comp_std + lo_total_comp_std + bm_branch_age + collect_monthly + fg_n_comp_mfi + stateid + mfiuniqueid, switch, na.action = na.omit)
master.x <- na.omit(switch[ , c("mfiuniqueid", all.vars(formula(meaning4)))])
cov1   <- clx(meaning4, 1, master.x$mfiuniqueid)
cluster.se4  <- sqrt(diag(cov1))

stargazer(meaning1, meaning2, meaning3, meaning4, title="Origins of Meaning: Subsamples", dep.var.labels = "I sometimes feel my job is meaningless", column.labels=c("Full Sample", "Always NGO", "Always NBFC", "Switch"), covariate.labels = c("Took job because of interest in work","'Our company contributes positively to the community'", "'I like the people I work with'", "Altruism (avg. given in dictator game)", "Productivity (std)","BM total compensation (std)", "LO total compensation (std)",  "Always NGO", "Always NBFC"),omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4))



