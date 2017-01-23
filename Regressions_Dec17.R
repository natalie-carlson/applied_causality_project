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


############
## TABLES ##
############

##Meaning --> Productivity##

productivity1 <- lm(productivity_std ~ bm_j_30 + lo_j_26, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity1)))])
cov1   <- clx(productivity1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

productivity2 <- lm(productivity_std ~ bm_j_30 + lo_j_26 + bm_bonus_prop + lo_bonus_prop + bm_trainingdays + lo_totaltraining_days, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity2)))])
cov1   <- clx(productivity2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

productivity3 <- lm(productivity_std ~ bm_j_30 + lo_j_26 + bm_bonus_prop + lo_bonus_prop + bm_trainingdays + lo_totaltraining_days + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity3)))])
cov1   <- clx(productivity3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

productivity4 <- lm(productivity_std ~ bm_j_30 + lo_j_26 + bm_bonus_prop + lo_bonus_prop + bm_trainingdays + lo_totaltraining_days + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity4)))])
cov1   <- clx(productivity4, 1, master.x$mfiuniqueid)
cluster.se4  <- sqrt(diag(cov1))

productivity5 <- lm(productivity_std ~ bm_j_30 + lo_j_26 + bm_bonus_prop + lo_bonus_prop + bm_trainingdays + lo_totaltraining_days +  always_ngo + always_nbfc + default_port_prop + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity5)))])
cov1   <- clx(productivity5, 1, master.x$mfiuniqueid)
cluster.se5  <- sqrt(diag(cov1))

productivity6 <- lm(productivity_std ~ bm_j_30 + lo_j_26 + bm_bonus_prop + lo_bonus_prop + bm_trainingdays + lo_totaltraining_days +  always_ngo + always_nbfc + default_port_prop + avg_loan_size_std +  always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity6)))])
cov1   <- clx(productivity6, 1, master.x$mfiuniqueid)
cluster.se6  <- sqrt(diag(cov1))

stargazer(productivity1, productivity2, productivity3, productivity4, productivity5, productivity6, title="Productivity Regressions", dep.var.labels = "Productivity measure (Standardized)", covariate.labels = c("BM: I sometimes feel my job is meaningless", "LO: I sometimes feel my job is meaningless", "BM bonus as prop. of comp.", "LO bonus as prop. of comp.", "BM training days", "LO training days", "Always NGO", "Always NBFC", "Prop. of clients in default", "Average loan size (std)"), omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "bm_b3_gender","stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls", "Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4, cluster.se5, cluster.se6))


##Meaning --> Hours##

hours1 <- lm(bm_lo_hoursdiff ~ bm_j_30 + bm_total_comp_std + bm_bonus_prop + default_client_prop + bm_j_35 + bm_incentive_loanrepayment + bm_incentive_newclients + bm_lo_informal, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(hours1)))])
cov1   <- clx(hours1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

hours2 <- lm(bm_lo_hoursdiff ~ bm_j_30 + bm_total_comp_std + bm_bonus_prop + default_client_prop + bm_j_35 + bm_incentive_loanrepayment + bm_incentive_newclients + bm_lo_informal + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(hours2)))])
cov1   <- clx(hours2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

hours3 <- lm(bm_lo_hoursdiff ~ bm_j_30 + bm_total_comp_std + bm_bonus_prop + default_client_prop + bm_j_35 + bm_incentive_loanrepayment + bm_incentive_newclients + bm_lo_informal + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid , master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(hours3)))])
cov1   <- clx(hours3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))


stargazer(hours1, hours2, hours3, title="Branch Manager Hours", dep.var.labels = "BM daily hours - LO daily hours", covariate.labels = c("BM: I sometimes feel my job is meaningless", "BM total compensation (std)", "BM bonus as prop. of comp.", "Client default proportion","BM: My job is enjoyable", "BM has incentive for loan repayment", "BM has incentive for new clients","BM/LO informal meetings", "Always NGO", "Always NBFC"),omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "bm_b3_gender","stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls","Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3))


##Hours --> Productivity##

productivity1 <- lm(productivity_std ~ bm_lo_hoursdiff + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity1)))])
cov1   <- clx(productivity1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

productivity2 <- lm(productivity_std ~ bm_lo_hoursdiff + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity2)))])
cov1   <- clx(productivity2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

productivity3 <- lm(productivity_std ~ bm_lo_hoursdiff + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity3)))])
cov1   <- clx(productivity3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

stargazer(productivity1, productivity2, productivity3, title="Productivity Regressions", dep.var.labels = "Productivity measure (Standardized)", covariate.labels = c("BM daily hours - LO daily hours", "BM bonus as prop. of comp.", "LO bonus as prop. of comp.", "Proportion of employees to resign","Prop. of clients in default", "Always NGO", "Always NBFC"), omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "bm_b3_gender","stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls", "Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3))


##Meaning --> Generosity##

generosity1 <- lm(dictator_mean ~ bm_j_30, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(generosity1)))])
cov1   <- clx(generosity1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

generosity2 <- lm(dictator_mean ~ bm_j_30 + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(generosity2)))])
cov1   <- clx(generosity2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

generosity3 <- lm(dictator_mean ~ bm_j_30 + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(generosity3)))])
cov1   <- clx(generosity3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

generosity4 <- lm(dictator_mean ~ bm_j_30 + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(generosity4)))])
cov1   <- clx(generosity4, 1, master.x$mfiuniqueid)
cluster.se4   <- sqrt(diag(cov1))

stargazer(generosity1, generosity2, generosity3, generosity4, title="Generosity Regressions", dep.var.labels = "Average amount given in dictator game (BM/LO)", covariate.labels = c("BM: I sometimes feel my job is meaningless", "BM bonus as prop. of comp.", "LO bonus as prop. of comp.", "Proportion of employees to resign", "Prop. of clients in default", "Always NGO", "Always NBFC"), omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "bm_b3_gender","stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls", "Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4))


## Generosity --> Productivity ##

productivity1 <- lm(productivity_std ~ dictator_mean, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity1)))])
cov1   <- clx(productivity1, 1, master.x$mfiuniqueid)
cluster.se1   <- sqrt(diag(cov1))

productivity2 <- lm(productivity_std ~ dictator_mean + bm_bonus_prop + lo_bonus_prop + emp_prop_resign, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity2)))])
cov1   <- clx(productivity2, 1, master.x$mfiuniqueid)
cluster.se2   <- sqrt(diag(cov1))

productivity3 <- lm(productivity_std ~ dictator_mean + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity3)))])
cov1   <- clx(productivity3, 1, master.x$mfiuniqueid)
cluster.se3   <- sqrt(diag(cov1))

productivity4 <- lm(productivity_std ~ dictator_mean + bm_bonus_prop + lo_bonus_prop + emp_prop_resign + default_client_prop + always_ngo + always_nbfc + bm_branch_age + collect_monthly + fg_n_comp_mfi + bm_b3_gender + stateid + mfiuniqueid, master, na.action = na.omit)
master.x <- na.omit(master[ , c("mfiuniqueid", all.vars(formula(productivity4)))])
cov1   <- clx(productivity4, 1, master.x$mfiuniqueid)
cluster.se4   <- sqrt(diag(cov1))

stargazer(productivity1, productivity2, productivity3, productivity4, title="Productivity Regressions", dep.var.labels = "Productivity measure (Standardized)", covariate.labels = c("Average amount given in dictator game (BM/LO)","BM bonus as prop. of comp.", "LO bonus as prop. of comp.", "Proportion of employees to resign","Prop. of clients in default", "Always NGO", "Always NBFC"), omit = c("bm_branch_age", "collect_monthly", "fg_n_comp_mfi", "bm_b3_gender","stateid", "mfiuniqueid"), omit.labels = c("Controls","Controls","Controls", "Controls","State dummies", "MFI dummies"),no.space = T, digits = 2, df=F, omit.stat=c("ll", "f", "ser", "adj.rsq"), font.size = "footnotesize", se = list(cluster.se1,cluster.se2, cluster.se3, cluster.se4))



