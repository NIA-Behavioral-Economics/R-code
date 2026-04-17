#install packages
#for random number generation
install.packages("random", dependencies = T)
library(random)
#for splitting variable by median
install.packages("sjmisc", dependencies = T)
library(sjmisc)

random <- function (filename="Randomization.csv", 
                    folder="/Users/emilystewart/Desktop/", 
                    v="MME",
                    m="msplit",
                    no=30) 
{
  #import data
  df <- read.csv(file=paste0(folder, filename))
  #split data based on median for block variable (e.g., MME), 0 = <median, 1 = >median
  df[[m]] <- dicho(df[[v]], dich.by="median")
  #order by continous (numeric) block variable and binary variable based on median
  df <- df[order(df[[v]], df[[m]]),]
  #retrieve two lists of binary random numbers from random.org; 0=control, 1=treatment
  random1 <- randomNumbers(n=no, min=0, max=1, col=1)
  random2 <- randomNumbers(n=no, min=0, max=1, col=1)
  #append lists (should remain in same order)
  tx <- as.factor(append(random1, random2))
  #merge with data 
  df <- cbind(df, tx)
  return(df)
}

mydat <- random()