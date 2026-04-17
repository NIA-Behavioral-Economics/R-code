###################################################################################################
### Power Analysis Polypharmacy Grant 4/15/2026
### Assumptions
#### ICC estimates based on AESOPS data
###### prescriber = 0.09
###### clinic = 0.05
###### Average no. visits = 250 (SD = 40) (May want to check with Ji Young)
###### 90 clinics
###### 300 clinicians
###### Base rate = test 1%-9% by 2%
###### Minimum effect required for 80% power
###### Model glmer: desprescribing (0 vs. 1) = b1(tx) + (1|clinic/clinician)
###################################################################################################

#load libraries 
library(dplyr)
library(lme4)
library(parallel)

#simulation
run_sim <- function(i, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed + i)
  
  #parameters
  nclinic <- 90
  nprov <- 300
  obs_per_prov <- round(rnorm(nprov, mean = 250, sd = 40))
  obs_per_prov[obs_per_prov < 0] <- 0
  obs_per_prov[obs_per_prov > 250] <- 250
  icc_prov_target   <- 0.09
  icc_clinic_target <- 0.05
  #base (control) probability 
  pcont <- 0.05
  #tx probability
  ptx <- 0.068
  #within group variance for logistic
  resid_var <- (pi^2) / 3
  
  #clinic size ('remainder' necessary when number of clinicians is not multiple of clinics, e.g., 300 and 90)
  base_n <- nprov %/% nclinic
  remainder <- nprov %% nclinic
  clinic_sizes <- rep(base_n, nclinic)
  if (remainder > 0) {
    clinic_sizes[1:remainder] <- clinic_sizes[1:remainder] + 1
  }
  
  #id variables
  clinic_id <- rep(1:nclinic, times = clinic_sizes)
  prov_id <- 1:nprov
  df <- data.frame(clinic = clinic_id, provID = prov_id)
  
  #assign clinic tx
  clinic_assign <- data.frame(
    clinic = 1:nclinic,
    tx = rbinom(nclinic, 1, 0.5)
  )
  df <- left_join(df, clinic_assign, by = "clinic")
  
  #expand data by number of rows per clinician 
  data <- df[rep(1:nrow(df), times = obs_per_prov), ]
  
  #random effects
  total_var <- resid_var / (1 - icc_prov_target - icc_clinic_target)
  var_prov   <- icc_prov_target * total_var
  var_clinic <- icc_clinic_target * total_var
  data$prov_re   <- rnorm(nprov, 0, sqrt(var_prov))[data$provID]
  data$clinic_re <- rnorm(nclinic, 0, sqrt(var_clinic))[data$clinic]
  
  #check ICCs
  var_prov_obs <- var(data$prov_re)
  var_clinic_obs <- var(data$clinic_re)
  total_var_obs <- var_prov_obs + var_clinic_obs + resid_var
  icc_prov   <- var_prov_obs / total_var_obs
  icc_clinic <- var_clinic_obs / total_var_obs
  
  #fixed effect-convert to log scale
  #this function is to determine optimal fixed effect intercept to produce 
  #desired control probability (from 0.01 to 0.09)
  re_base <- rnorm(1e5)
  target_fn <- function(beta0, re_var_total, target_p) {
    re <- re_base * sqrt(re_var_total)
    mean(plogis(beta0 + re)) - target_p
  }
  
  beta0_calibrated <- uniroot(
    target_fn,
    interval = c(-10, 0),
    re_var_total = var_prov + var_clinic,
    target_p = pcont
  )$root
  
	#effect on log odds scale
  beta1 <- qlogis(ptx) - qlogis(pcont)
  data$eta <- beta0_calibrated + beta1 * data$tx +
 	data$prov_re + data$clinic_re
  
  #assign rows (0 vs.1) for probabilities converted from eta
  data$deprescribe <- rbinom(nrow(data), 1, plogis(data$eta))
  prob_cont <- mean(data$deprescribe[data$tx == 0])
  prob_tx   <- mean(data$deprescribe[data$tx == 1])
  
  #factors
  data$tx <- as.factor(data$tx)
  data$clinic <- as.factor(data$clinic)
  data$provID <- as.factor(data$provID)
  data$deprescribe <- as.factor(data$deprescribe)
  
  #model 
  m <- glmer(deprescribe ~ tx + (1|clinic/provID),
             data = data,
             family = binomial,   
             control = glmerControl(calc.derivs = FALSE))
  
  #extract estimates
  coefs <- summary(m)$coefficients
  z <- coefs["tx1", "Estimate"] / coefs["tx1", "Std. Error"]
  pval <- 2 * pnorm(-abs(z))
  
  #return as dataframe 
  return(c(
    prob_cont = prob_cont,
    prob_tx = prob_tx,
    icc_prov = icc_prov,
    icc_clinic = icc_clinic,
    Intercept = coefs["(Intercept)", "Estimate"],
    tx1 = coefs["tx1", "Estimate"],
    p_tx1 = pval,
    n_obs = nrow(data)
  ))
}

#this section is for parallel processing 
nsims <- 1000
cl <- makeCluster(detectCores() - 1)

clusterExport(cl, "run_sim")
clusterEvalQ(cl, {
  library(dplyr)
  library(lme4)
})

results <- parLapplyLB(cl, 1:nsims, run_sim)
stopCluster(cl)
results <- do.call(rbind, results)
results <- as.data.frame(results)
results$sig <- ifelse(results$p_tx1 < 0.05, 1, 0)

#mean simulated estimates 
#prob. control
round(mean(results$prob_cont),4)

#prob. tx
round(mean(results$prob_tx),4)

#effect size 
results$efs <- results$prob_tx - results$prob_cont
round(mean(results$efs),4)

#clinician ICC
round(mean(results$icc_prov),4)

#clinic ICC
round(mean(results$icc_clinic),4)

#mean number of obs. 
round(mean(results$n_obs),0)

#power 
round(mean(results$sig),4)
