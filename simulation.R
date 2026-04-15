#if not home directory, set to home 
#setwd("/schhome/users/epstewar")

#load libraries
library(sqldf)
library(dplyr)
library(lme4)
library(tidyr)
library(purrr)
library(performance)
library(sjstats)

#number of simulations
nsims <- 1000

#empty df to store glmer results for every simulation (n = 1,000)
estimates <- data.frame(
  Intercept = numeric(nsims),
  assignment1 = numeric(nsims),
  p_assignment1 = numeric(nsims)
)

#store ICC per simulation
ICC_adj <- numeric(nsims)
ICC_unadj <- numeric(nsims)

# data parameters
### 32 classes per subject
### 10 people per class
### base control attendance rate = 0.59
### subject ICC = 0.51, class ICC = 0.05 (specified numbers must be higher to acheive desired model ICC)
##### number of subjects
nsub <- 250
##### number of classes
nclasses <- 25
##### tx probability
ptx <- 0.71
##### subject ICC 
icc_sub <- 0.75
##### class ICC
icc_class <- 0.12

#simulations
for (i in 1:nsims) {
  #set.seed(i)
  #create data (first as counts to add group level variances, will later be elongated for logistic glmer)
  subjects <- as.data.frame(cbind(1:nsub))
  subjects <- subjects %>% rename(subjectID = V1)
  classes <- as.data.frame(cbind((rep(1:nclasses, each = 10))))
  assignment <- as.data.frame(rbinom(nsub, 1, 0.5))
  data <- cbind(subjects, classes, assignment)
  data <- data %>% rename(assignment = `rbinom(nsub, 1, 0.5)`, class = V1)
 
  #randomly assign percent attendance based on 32 'trials'
  data$att <- ifelse(data$assignment == 0, rbinom(nrow(data), 32, 0.59), rbinom(nrow(data), 32, ptx))
 
  #use rnorm to generate between group (random intercept) variance
  ### ICC = random intercept variance/random intercept variance + within group (residual) variance 
  ### random intercept variance = ICC*residual variance/1-ICC
  ### 3.29 is the residual variance for logistic models (pi-squared/3 = 3.29)
  #random intercept variance for subject 
  vars <- (icc_sub*3.29)/(1-icc_sub)
  data <- data %>%
  group_by(subjectID) %>%
  mutate(var_sub = rnorm(1, mean = 0, sd = vars))
  mean(data$var_sub)
 
  #random intercept variance for class
  varc <- (icc_class*3.29)/(1-icc_class)
  data <- data %>%
  group_by(class) %>%
  mutate(var_class = rnorm(1, mean = 0, sd = varc))
 
  #add random intercept variances to y
  data$attend <- round((data$att+data$var_class+data$var_sub),0)
  #cap at 32 and 0 in cases where variance exceeds number of classes
  data$attend <- ifelse(data$attend > 32, 32, data$attend)
  data$attend <- ifelse(data$attend < 0, 0, data$attend)
 
  #number of classes subject didn't attend
  data$notattend <- 32-data$attend
 
  #drop unnecessary variables
  drops <- c("att", "var_sub", "var_class")
  data <- data[ , !(names(data) %in% drops)]
 
  #elongate data to binary form (0 vs. 1) from binomially distributed visit counts
  #attended rows
  attend <- data %>%
    mutate(attend = map(attend, seq_len)) %>%
    unnest(attend)
 
  #change attendance to 1
  attend$attend <- 1
 
  #unattended rows
  notattend <- data %>%
    mutate(notattend = map(notattend, seq_len)) %>%
    unnest(notattend)
 
  #change attendance to 0
  notattend$attend <- 0
 
  #append attended and unattended rows
  final <- rbind(attend, notattend)
  drops <- c("notattend")
  final <- final[ , !(names(final) %in% drops)]
 
  #convert variables to factors
  final$assignment <- as.factor(final$assignment)
  final$class <- as.factor(final$class)
  final$subjectID <- as.factor(final$subjectID)
 
  #glmer model which executes for each simulation (n = 1,000)
  m <- glmer(attend ~ assignment + (1|subjectID/class), data = final, family = binomial,
             control = glmerControl(calc.derivs = FALSE))
  summary(m)
 
  #ICC estimates
  ICC_out <- performance::icc(m)
  ICC_adj[i] <- as.numeric(ICC_out$ICC_adjusted)[1]
  ICC_unadj[i] <- as.numeric(ICC_out$ICC_unadjusted)[1]
 
  #extract p value for assignment 
  coefs <- summary(m)$coefficients
  estimates$p_assignment1[i] <- coefs["assignment1", "Pr(>|z|)"]
}

#check mean model ICC-should approximate 0.51
#mean adjusted ICC
ICC_adj <- as.data.frame(ICC_adj)
mean(ICC_adj$ICC_adj)

#mean unadjusted ICC
ICC_unadj <- as.data.frame(ICC_unadj)
mean(ICC_unadj$ICC_unadj)

#percentage of simulations where p is significant (i.e., power)
estimates$sig <- ifelse(estimates$p_assignment1 < 0.05, 1, 0)
mean(estimates$sig)



