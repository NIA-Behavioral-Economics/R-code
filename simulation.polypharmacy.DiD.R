
#load libraries 
library(dplyr)
library(lme4)
library(parallel)
library(sqldf)
library(nleqslv)

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
  p00 <- 0.27
  #tx probability
  p01 <- 0.27
  #covariate probability 
  #base group
  p10 <- 0.325
  #cov tx group 
  p11 <- 0.355
  #within group variance for logistic
  resid_var <- (pi^2) / 3
  
  #clinic size ('remainder' necessary when number of clinicians is not a multiple of clinics, e.g., 300 and 90)
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
  data$cov <- rbinom(nrow(data), 1, 0.5)
  
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
  
  #calibrate estimates so that observed probabilities match specified effects 
  data$re_total <- data$prov_re + data$clinic_re
  
  target_fn4 <- function(par, data, targets) {
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    beta3 <- par[4]
    
    eta <- beta0 +
      beta1 * data$tx +
      beta2 * data$cov +
      beta3 * (data$tx * data$cov) +
      data$re_total
    
    p <- plogis(eta)
    
    c(
      mean(p[data$tx==0 & data$cov==0]) - targets[1],
      mean(p[data$tx==0 & data$cov==1]) - targets[2],
      mean(p[data$tx==1 & data$cov==0]) - targets[3],
      mean(p[data$tx==1 & data$cov==1]) - targets[4]
    )
  }
  
  sol <- nleqslv(
    x = c(qlogis(p00), 0, 0, 0),
    fn = target_fn4,
    data = data,
    targets = c(p00, p01, p10, p11)
  )
  
  beta0_calibrated <- sol$x[1]
  beta1 <- sol$x[2]
  beta2 <- sol$x[3]
  beta3 <- sol$x[4]
  
  data$eta <- beta0_calibrated + (beta1 * data$tx) + 
    (beta2 * data$cov) + (beta3 * (data$tx * data$cov)) +
    data$prov_re + data$clinic_re
  
  #assign rows (0 vs.1) for probabilities converted from eta
  data$deprescribe <- rbinom(nrow(data), 1, plogis(data$eta))
  prob_00 <- mean(data$deprescribe[data$tx == 0 & data$cov == 0])
  prob_01   <- mean(data$deprescribe[data$tx == 0 & data$cov == 1])
  prob_10 <- mean(data$deprescribe[data$tx == 1 & data$cov == 0])
  prob_11   <- mean(data$deprescribe[data$tx == 1 & data$cov == 1])
  
  #factors
  data$clinic <- as.factor(data$clinic)
  data$provID <- as.factor(data$provID)
  
  #model 
  m <- glmer(deprescribe ~ tx + cov + tx:cov + (1|clinic/provID),
             data = data,
             family = binomial,   
             control = glmerControl(calc.derivs = FALSE))
  
  #extract estimates
  coefs <- summary(m)$coefficients
  z <- coefs["tx:cov", "Estimate"] / coefs["tx:cov", "Std. Error"]
  pval <- 2 * pnorm(-abs(z))
  
  #return as dataframe 
  return(c(
    prob_00 = prob_00,
    prob_01 = prob_01,
    proc_10 = prob_10,
    prob_11 = prob_11,
    icc_prov = icc_prov,
    icc_clinic = icc_clinic,
    txcov = coefs["tx:cov", "Estimate"],
    p_txcov = pval,
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
  library(nleqslv)
})

results <- parLapply(cl, 1:nsims, run_sim)
stopCluster(cl)
results <- do.call(rbind, results)
results <- as.data.frame(results)
results$sig <- ifelse(results$p_txcov < 0.05, 1, 0)

#mean simulated estimates 
#probabilities in each arm
round(mean(results$prob_00),4)
round(mean(results$prob_01),4)
round(mean(results$proc_10),4)
round(mean(results$prob_11),4)

#tx effect size 
results$p0 <- (results$prob_00+results$prob_01)/2
results$p1 <- (results$prob_11+results$proc_10)/2
results$efs <- results$p1-results$p0
round(mean(results$efs),4)

#Did
results$DiD <- (results$prob_11 - results$proc_10) - (results$prob_01 - results$prob_00)
mean(results$DiD)

#clinician ICC
#round(mean(results$icc_prov),4)

#clinic ICC
#round(mean(results$icc_clinic),4)

#mean number of obs. 
#round(mean(results$n_obs),0)

#power 
round(mean(results$sig),4)