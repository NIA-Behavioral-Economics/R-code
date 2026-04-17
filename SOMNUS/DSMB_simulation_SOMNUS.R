############################################################################################################################################################################
### Power Analysis SOMNUS DSMB 9/6/2024
### Assumptions based on AESOPS R33 data
###### prescriber ICC = 0.09
###### clinic ICC = 0.05
###### Average no. visits = 460 
###### 60 clinics
###### 400 prescribers
###### Base rate = 36 pills
###### Alpha = 0.05
###### 5.2 pill decrease post-intervention in intervention arms 
###### Restricted Model: Y = b1(mnth) + b2(kmnth) + b3(aj) + b4(default) + b5(mnth*aj) + b6(mnth*default) + b7(kmnth*aj) + b8(kmnth*default) + (1|clinic) + (1|clinician)
###########################################################################################################################################################################

#libraries 
library(dplyr)
library(sqldf)
library(broom)
library(lme4)

#create empty datasets to get mean effect sizes and power (# simulations with p value < 0.05)
kdef <- c()
kaj <- c()
def_p <- c()
aj_p <- c()
p_def <- c()
p_aj <- c()
p_aj_or_def <- c()
p_aj_and_def <- c()
mean_p_dist <- c()

for (i in 1:200) {
 
  #no. clinic and prescribers
  p1 <- as.data.frame(cbind(1:100, (rep(1:20, each = 5))))
  p2 <- as.data.frame(cbind(101:240, (rep(21:40, each = 7))))
  p3 <- as.data.frame(cbind(241:400, (rep(41:60, each = 8))))
  sample <- rbind(p1, p2, p3)
 
  #clinic ICC
  f_var_alpha <- function(ICC, var_epsilon){
    var_alpha <- (ICC*var_epsilon)/(1-ICC)
    return(var_alpha)
  }
 
  #clinic variance
  var_clinic <-f_var_alpha(ICC = 0.05, var_epsilon = 1)
  var_clinic
 
  #randomly generate clinic variance
  sample_v2 <- sample %>%
  group_by(V2) %>%
  mutate(clinic_var = rnorm(1, mean = 0, sd = sqrt(var_clinic)))
 
  #prescriber variance
  var_pres <-f_var_alpha(ICC = 0.09, var_epsilon = 1)
  var_pres
 
  #randomly generate prescriber variance
  sample_v2$pres_var <- rnorm(400, mean = 0, sd = sqrt(var_pres))
  sample_v2
 
  #randomly assign intervention
  aj <- as.data.frame(cbind(rbinom(60, 1, 0.5), 1:60))
  def <- as.data.frame(cbind(rbinom(60, 1, 0.5), 1:60))
  int <- merge(aj, def, by = "V2")
  sample_v3 <- merge(sample_v2, int, by = "V2")
 
  #rename variables
  sample_v4 <- rename(sample_v3, clinic = V2, clinician = V1, def = V1.x, aj = V1.y)
  table(sample_v4$def, sample_v4$aj)
 
  #randomly generate per-prescriber number of visits
  no_visits <- round(runif(400, 1, 460), 0)
 
  #add to sample
  sample_v5 <- cbind(sample_v4, no_visits)
 
  #randomly generate months using randomly generated number of per-prescriber visits
  s <- lapply(
    seq_along(sample_v5$no_visits),
    \(i) data.frame(Month = sample(-18:18, sample_v5$no_visits[i], replace = TRUE), clinician = i)
  )
  months <- do.call(rbind.data.frame, s)

  #knotted month and pre vs. post time
  sample_v6 <- merge(sample_v5, months, by = "clinician")
  sample_v6$kmonth <- ifelse(sample_v6$Month <= 0, 0, sample_v6$Month)
  sample_v6$post <- ifelse(sample_v6$kmonth == 0, 0, 1)
 
  #total number of Rxs in each study group pre-to-post
  sample_v6$int <- ifelse(sample_v6$aj == 1 | sample_v6$def == 1, 1, 0)
  cts <- table(sample_v6$int, sample_v6$post)
  precont <- cts[1,1]
  preint <- cts[2,1]
  postcont <- cts[1,2]
  postint <- cts[2,2]
 
  #randomly generate Y for each study arm in pre- and post-period 
  #pre-control
  yprecont <- as.data.frame((rpois(precont, 36))+rnorm(precont, mean = 0, sd = 1))
  yprecont <- rename(yprecont, Y = `(rpois(precont, 36)) + rnorm(precont, mean = 0, sd = 1)`)
 
  #pre-int
  ypreint <- as.data.frame((rpois(preint, 36))+rnorm(preint, mean = 0, sd = 1))
  ypreint <- rename(ypreint, Y = `(rpois(preint, 36)) + rnorm(preint, mean = 0, sd = 1)`)
 
  #post-cont
  ypostcont <- as.data.frame((rpois(postcont, 36))+rnorm(postcont, mean = 0, sd = 1))
  ypostcont <- rename(ypostcont, Y = `(rpois(postcont, 36)) + rnorm(postcont, mean = 0, sd = 1)`)
 
  #post-int
  ypostint <- as.data.frame((rpois(postint, 30.8))+rnorm(postint, mean = 0, sd = 1))
  ypostint <- rename(ypostint, Y = `(rpois(postint, 30.8)) + rnorm(postint, mean = 0, sd = 1)`)
 
  #append y
  y <- rbind(yprecont, ypreint, ypostcont, ypostint)
 
  #order by pre/post and intervention
  sample_v6 <- sqldf("select * from sample_v6 order by post, int")
  sample_v7 <- cbind(sample_v6, y)
  sample_v7$Y <- round((sample_v7$Y+sample_v7$clinic_var+sample_v7$pres_var),0)
 
  #drop unnecessary variables
  drops <- c("clinic_var","no_visits", "int")
  sample_v8 <- sample_v7[ , !(names(sample_v7) %in% drops)]

  #model
  m1 <- glmer(Y ~ Month + kmonth + def + aj +
                  Month:def + Month:aj + kmonth:def +
                  kmonth:aj + (1|clinician),
                  data = sample_v8, family = poisson(link = "log"),  
                  control = glmerControl(calc.derivs = FALSE))
 
  #model results for each simulation
  #estimates
  kdef[i] <- coef(summary(m1))["kmonth:def","Estimate"]
  kaj[i] <- coef(summary(m1))["kmonth:aj","Estimate"]

  #p value
  def_p[i] <- coef(summary(m1))["kmonth:def","Pr(>|z|)"]
  aj_p[i] <- coef(summary(m1))["kmonth:aj","Pr(>|z|)"]

  #number of samples with sig. default
  p_def[i] <- ifelse(def_p[i] < 0.05, 1, 0)

  #number samples with sig. AJ
  p_aj[i] <- ifelse(aj_p[i] < 0.05, 1, 0)

  #number samples with at least one sig. intervention effect
  p_aj_or_def[i] <- ifelse(def_p[i] | aj_p[i] < 0.05, 1, 0)

  #number samples where both interventions are significant
  p_aj_and_def[i] <- ifelse(def_p[i] & aj_p[i] < 0.05, 1, 0)
  mean_p_dist[i] <- mean(p_aj_and_def[i])
}


#mean knotted coefficient for default
mean(kdef)

#95% CI for default
quantile(kdef, 0.05)
quantile(kdef, 0.95)

#mean knotted coefficient for AJ
mean(kaj)

#95% CI for aj
quantile(kaj, 0.05)
quantile(kaj, 0.95)

#power
summary(p_aj_and_def)
power <- mean(p_aj_and_def)
power

#margin of error
me <- 2*sd(p_def<=0.05)/sqrt(1000)
me

#power confidence interval
lcl <- power-me
lcl
ucl <- power+me
ucl