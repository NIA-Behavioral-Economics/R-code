#######################################################################################
### SOMNUS Power Analysis: last update 2/14/2025
### Assumptions 
###### prescriber ICC = 0.09
###### clinic ICC = 0.05
###### Average no. visits = 4.8
###### 64 clinics
###### 444 prescribers
###### Base rate = 33 pills
###### Alpha = 0.05
###### 1.5 mean pill decrease pre-to-post intervention for both AJ and default study arms 
##########################################################################################

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
mean(var_pres)

#Randomly generate prescriber variance 
sample_v2$pres_var <- rnorm(444, mean = 0, sd = sqrt(var_pres))

#randomly assign intervention 
aj <- as.data.frame(cbind(rbinom(64, 1, 0.5), 1:64))
def <- as.data.frame(cbind(rbinom(64, 1, 0.5), 1:64))

#merge interventions and clinics with sample 
sample_v3 <- sqldf("select t.V1 as clinician, t.V2 as clinic,
                      t.clinic_var, t.pres_var, l.V1 as aj,
                      n.V1 as def, o.V1 as high_low from sample_v2 t 
                      left join 
                      aj l 
                      on t.V2 = l.V2 
                      left join 
                      def n
                      on t.V2 = n.V2
                      left join 
                      high_low o 
                      on t.V2 = o.V2")

#randomly generate per-prescriber number of Rxs based on Ji Young's aggregate estimates
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
precont
preint <- cts[2,1]
preint
postcont <- cts[1,2]
postcont
postint <- cts[2,2]
postint

#no. of pills
#pre-control
yprecont <- as.data.frame((rpois(precont, 33))+rnorm(precont, mean = 0, sd = 1))
yprecont <- rename(yprecont, Y = `(rpois(precont, 33)) + rnorm(precont, mean = 0, sd = 1)`)

#pre-int
ypreint <- as.data.frame((rpois(preint, 33))+rnorm(preint, mean = 0, sd = 1))
ypreint <- rename(ypreint, Y = `(rpois(preint, 33)) + rnorm(preint, mean = 0, sd = 1)`)

#post-cont
ypostcont <- as.data.frame((rpois(postcont, 33))+rnorm(postcont, mean = 0, sd = 1))
ypostcont <- rename(ypostcont, Y = `(rpois(postcont, 33)) + rnorm(postcont, mean = 0, sd = 1)`)

#post-int
ypostint <- as.data.frame((rpois(postint, 31.5))+rnorm(postint, mean = 0, sd = 1))
ypostint <- rename(ypostint, Y = `(rpois(postint, 31.5)) + rnorm(postint, mean = 0, sd = 1)`)

#append y
y <- rbind(yprecont, ypreint, ypostcont, ypostint)

#order by pre/post and intervention 
sample_v5 <- sqldf("select t.* from sample_v5 t order by post,int")
sample_v6 <- cbind(sample_v5, y)

#add within-group (clinic and clinician) variance to pill count
sample_v6$Y <- round((sample_v6$Y+sample_v6$clinic_var+sample_v6$pres_var),0)

#drop unnecessary variables
drops <- c("clinic_var","pres_var", "no_visits", "int", "mean")
sample_v7 <- sample_v7[ , !(names(sample_v7) %in% drops)]

#diff-in-diff model (originally knotted spline)
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

#95% CI for default coefficient
quantile(est_def, 0.05)
quantile(est_def, 0.95)

#mean coefficient for AJ:post
mean(est_aj)

#95% CI for AJ coefficient 
quantile(est_aj, 0.05)
quantile(est_aj, 0.95)

#power: percentage simulated samples where both AJ + default/post interaction effects are significant 
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