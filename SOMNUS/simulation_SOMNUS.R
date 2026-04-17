#######################################################################################################################
## Goal:                                                                                                             ##
##    1. Function which runs power simulations using simr to determine ideal no. of clinics                          ## 
##       a. Model to simulate data:                                                                                  ##
##          1. DV: per-prescriber average daily dayssply of z drugs                                                  ##
##          2. IVs: def + aj + days + days_ + def*days + aj+days_                                                    ##                                          
##          3. Random intercepts for clinic and provider (nested)                                                    ##
##                                                                                                                   ##
##       b. Function arguments:                                                                                      ##
##           1. c: no. clinics                                                                                       ##
##           2. m: no. prescribers per-clinic; 8 based on NW data                                                    ##                                                              ##
##           3. avg: avg. no. per-prescriber events (dayssply); 36 based on NW data                                  ##
##           4. b0 - b6: fixed effects                                                                               ##
##           5. r1, r2: variances for prescriber and clinic                                                          ##                         
##             i. In absence of pilot data, above effects should be based on strong a priori assumptions.See source. ##                                                 
##           6. no_sim: number of simulations (1000 ideal)                                                           ##     
##                                                                                                                   ##
##    2. Source : Kumle, L., Võ, M.LH. & Draschkow, D. Estimating power in (generalized) linear mixed models:        ##
##          An open introduction and tutorial in R. Behav Res 53, 2528–2543 (2021).                                  ##
########################################################################################################################

#load packages 
library(dplyr)
library(simr)
library(parallel)

#function 
sim <- function(c, m, pr, avg, b0, b1, b2, b3, b4, b5, b6, r1, r2, no_sim) {

#create data
#no. prescribers
p = c*m

#no. clinic and prescribers
sample <- as.data.frame(cbind(1:p, (rep(1:c, each = m))))

#randomize interventions by clinic
set.seed(123)
ru1 <- cbind(rbinom(c, 1, 0.5), 1:c)
set.seed(456)
ru2 <- cbind(rbinom(c, 1, 0.5), 1:c)

#merge randomization id with orginal sample
sample <- merge(sample, ru1, by = "V2") %>% merge(ru2)

#add days; need to center, otherwise model results in errors
#32 visits per day
sample <- as.data.frame(cbind(sample, scale(rep(-546:546, each = p*32))))

#order by clinic and prescriber
sample <- sample[order(sample$V2, sample$V1.x, sample$`scale(rep(-546:546, each = p * 32))`  ),]

#sample size/number of obs 
#1093 = no. of days
#32 = no. visits per day
no <- c*m*1093*32

#zero-inflated Poisson distribution
set.seed(789)
sample <- cbind(sample, ifelse(rbinom(no, 1, p = pr) > 0, 0, rpois(no, avg)))

#rename variables
sample <- rename(sample, pres = V1.x, clinic = V2, aj = V1.y, def = V1, 
                days = `scale(rep(-546:546, each = p * 32))`  , 
                dayssply = `ifelse(rbinom(no, 1, p = pr) > 0, 0, rpois(no, avg))`)


#days truncated
sample$days_ <- ifelse(0 > sample$days, 0, sample$days)

#continuous
samplenz <- subset(sample, dayssply != 0)

#create two dataframes based on zero vs. numeric value  
sample$dayssply <- ifelse(sample$dayssply == 0, 0, 1 )
sample$dayssply <- as.factor(sample$dayssply)

#model formula
formula <- dayssply ~ aj + def + days + days_ + def:days_ +  aj:days_ + (1|clinic/pres)

#fixed and random effects
fixed_effects <- c(b0, b1, b2, b3, b4, b5, b6)
random_variance <- list(r1, r2)

#models
#zero vs. nonzero
artificial_glmer <- makeGlmer(formula, fixef = fixed_effects,
                                  VarCorr = random_variance,
                                 family = "binomial",
                                    data = sample)

#need to use different method with poisson, to avoid "boundary singular fit" warnings/Hessian matrix not positive definite
#control = glmerControl(calc.derivs) to avoid warning 
m1 <- glmer(dayssply ~ aj + def + days + days_ + aj:days_ + def:days_ + 
             (1|clinic/pres), data = samplenz, family = "poisson",
             control = glmerControl(calc.derivs = FALSE))

#change model fixed effect and random variances 
VarCorr(m1)["clinic"] <- r1
VarCorr(m1)["pres:clinic"] <- r2
fixef(m1)["aj:days_"] <- b6
fixef(m1)["def:days_"] <- -b5
fixef(m1)["aj"] <- b4
fixef(m1)["def"] <- b3
fixef(m1)["days"] <- b2
fixef(m1)["days_"] <- b1
fixef(m1)["(Intercept)"] <- b0

#simulations: 1000 ideal, however, takes long time (days)
power_simr_zero <- powerSim(fit = artificial_glmer,
                         test = fixed("aj:days_"), nsim = no_sim)

power_simr_nz <- powerSim(fit = m1,
                          test = fixed("aj:days_"), nsim = no_sim)

#check error and warning log
errors <- lastResult() $error_log
warnings <- lastResult() $warnings
 
#output
return(list(power_simr_nz, errors, warnings))

}

#with parallelization 
#clinic levels 
clinic_no <- c(10, 20, 30, 40)

#number of cores
no_cores <- detectCores() - 4

#parallelization 
clust <- makeCluster(no_cores, type = "FORK")

#results 
res <- parLapply(clust, clinic_no, sim, 8, 36, 0, 0, 0, 0, 0, -.05, -.05, 0.55, 1, 200)})

#shut down cluster
stopCluster(clust)