####################################################################################################
### Power Analysis SOMNUS DSMB last update: 2/24/2025
### Assumptions
#### ICC estimates based on AESOPS data 
###### prescriber = 0.09
###### clinic = 0.05
###### 64 clinics
##### 444 prescribers
##### Mean 4.8 Rxs per-month, per-clinician, min = 1, max = 41
##### Base rate = 43 pills 
##### Testing an average 1.58 pill reduction for intervention clinicians during post-period
##### Model: Y = b1(aj) + b2(default) + b3(post) + b4(post*aj) +
##### b7(post*default) + (1|clinic/clinician) 
###################################################################################################


#packages 
library(lme4)
library(dplyr)
library(sqldf)

#create empty datasets for model output
est_def <- c()
est_aj <- c()
def_p <- c()
aj_p <- c()
p_def <- c()
p_aj <- c()
p_aj_or_def <- c()
p_aj_and_def <- c()
mean_p_dist <- c()

for (i in 1:1000) {
  
#no. clinic and prescribers
p1 <- as.data.frame(cbind(1:132, (rep(1:22, each = 6))))
p2 <- as.data.frame(cbind(133:264, (rep(23:44, each = 6))))
p3 <- as.data.frame(cbind(265:444, (rep(45:64, each = 9))))
sample <- rbind(p1, p2, p3)

#Clinic ICC
f_var_alpha <- function(ICC, var_epsilon){
  var_alpha <- (ICC*var_epsilon)/(1-ICC)
  return(var_alpha)
}

#Clinic variance 
var_clinic <-f_var_alpha(ICC = 0.05, var_epsilon = 0.65)
var_clinic

#Randomly generate clinic variance 
sample_v2 <- sample %>% 
  group_by(V2) %>% 
  mutate(clinic_var = rnorm(1, mean = 0, sd = sqrt(var_clinic)))

#Prescriber variance 
var_pres <-f_var_alpha(ICC = 0.09, var_epsilon = 0.65)

#Randomly generate prescriber variance 
sample_v2$pres_var <- rnorm(444, mean = 0, sd = sqrt(var_pres))

#randomly assign intervention 
aj <- as.data.frame(cbind(rbinom(64, 1, 0.5), 1:64))
def <- as.data.frame(cbind(rbinom(64, 1, 0.5), 1:64))

#merge interventions and clinics with sample 
sample_v3 <- sqldf("select t.V1 as clinician, t.V2 as clinic,
                      t.clinic_var, t.pres_var, l.V1 as aj,
                      n.V1 as def from sample_v2 t 
                      left join 
                      aj l 
                      on t.V2 = l.V2 
                      left join 
                      def n
                      on t.V2 = n.V2")

#randomly generate per-prescriber number of Rxs
no_visits <- round(runif(4.8, 1, 41), 0)

#add to sample
sample_v4 <- cbind(sample_v3, no_visits)

#randomly generate months using randomly generated number of per-prescriber visits 
s <- lapply(
  seq_along(sample_v4$no_visits), 
  \(i) data.frame(Month = sample(-18:18, sample_v4$no_visits[i], replace = TRUE), clinician = i)
)
months <- do.call(rbind.data.frame, s)

#knotted month and pre vs. post time
sample_v5 <- merge(sample_v4, months, by = "clinician")
sample_v5$kmonth <- ifelse(sample_v5$Month <= 0, 0, sample_v5$Month)
sample_v5$post <- ifelse(sample_v5$kmonth == 0, 0, 1)

#combined intervention counts
sample_v5$int <- ifelse(sample_v5$aj == 1 | sample_v5$def == 1, 1, 0)
cts <- table(sample_v5$int, sample_v5$post)
precont <- cts[1,1]
preint <- cts[2,1]
postcont <- cts[1,2]
postint <- cts[2,2]

#no. of pills
#pre-control
yprecont <- as.data.frame((rpois(precont, 43))+rnorm(precont, mean = 0, sd = 1))
yprecont <- rename(yprecont, Y = `(rpois(precont, 43)) + rnorm(precont, mean = 0, sd = 1)`)

#pre-int
ypreint <- as.data.frame((rpois(preint, 43))+rnorm(preint, mean = 0, sd = 1))
ypreint <- rename(ypreint, Y = `(rpois(preint, 43)) + rnorm(preint, mean = 0, sd = 1)`)

#post-cont
ypostcont <- as.data.frame((rpois(postcont, 43))+rnorm(postcont, mean = 0, sd = 1))
ypostcont <- rename(ypostcont, Y = `(rpois(postcont, 43)) + rnorm(postcont, mean = 0, sd = 1)`)

#post-int
ypostint <- as.data.frame((rpois(postint, 41.42))+rnorm(postint, mean = 0, sd = 1))
ypostint <- rename(ypostint, Y = `(rpois(postint, 41.42)) + rnorm(postint, mean = 0, sd = 1)`)

#append y
y <- rbind(yprecont, ypreint, ypostcont, ypostint)

#order by pre/post and intervention 
sample_v5 <- sqldf("select t.* from sample_v5 t order by post,int")
sample_v6 <- cbind(sample_v5, y)
sample_v6$Y <- round((sample_v6$Y+sample_v6$clinic_var+sample_v6$pres_var),0)

#drop unnecessary variables
drops <- c("clinic_var","pres_var", "no_visits", "int")
sample_v7 <- sample_v6[ , !(names(sample_v6) %in% drops)]

#diff-in-diff model (originally we used knotted spline)
m1 <- glmer(Y ~ post + def + aj + post:def + post:aj + (1|clinic:clinician), 
            data = sample_v7, family = poisson(link = "log"),
            control = glmerControl(calc.derivs = FALSE))

# Get model results
est_def[i] <- coef(summary(m1))["post:def","Estimate"]
est_aj[i] <- coef(summary(m1))["post:aj","Estimate"]
def_p[i] <- coef(summary(m1))["post:def","Pr(>|z|)"]
aj_p[i] <- coef(summary(m1))["post:aj","Pr(>|z|)"]

#significance default
p_def[i] <- ifelse(def_p[i] < 0.05, 1, 0)
#significant aj
p_aj[i] <- ifelse(aj_p[i] < 0.05, 1, 0)
#significance of either 
p_aj_or_def[i] <- ifelse(def_p[i] | aj_p[i] < 0.05, 1, 0)
#sig. of both 
p_aj_and_def[i] <- ifelse(def_p[i] & aj_p[i] < 0.05, 1, 0)
mean_p_dist[i] <- mean(p_aj_and_def[i])
}

#mean coefficient for default:post
mean(est_def)

#95% CI for default
quantile(est_def, 0.05)
quantile(est_def, 0.95)

#mean knotted coefficient for AJ
mean(est_aj)

#95% CI for aj
quantile(est_aj, 0.05)
quantile(est_aj, 0.95)

#power 
summary(p_aj_and_def)
power <- mean(p_aj_and_def)
power

#margin of error
me <- 2*sd(p_aj_and_def<=0.05)/sqrt(1000)
me

#power confidence interval 
lcl <- power-me
lcl
ucl <- power+me
ucl