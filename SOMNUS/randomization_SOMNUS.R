#install packages
#for random number generation
#random.org uses atmospheric noise
#with von Neumann randomness
#extractor function
#https://arxiv.org/pdf/2101.02345
install.packages("random", dependencies = T)
install.packages("readr", dependencies =T)
library(random)

#Create Experimental Factors
#nesting justification conditions
#within default conditions
factor1 <- rep(c("No Default","Default"), times=2, each=2, 32)
factor2 <- rep(c("No Justification","Justification"), times=1, each=1, 32)

#create a random sequence of array positions
#for use with each factor; for each of high 
#and low Rx
clinic.assign.low.array.pos <- randomSequence(min=1, max=32, col=1)
clinic.assign.high.array.pos <- randomSequence(min=1, max=32, col=1)

#assign random array positions to full clinic list of 64: 32 low and 32 high Rx rate
#preserving the randomization of eight 2 x 2 factorial cells for each of low and high
clinic.assign.array.pos <- c(clinic.assign.low.array.pos,clinic.assign.high.array.pos)
clinic.assign <- rep("a",64)
#conduct assignment by assigning random array positions
for (i in 1:length(clinic.assign)) {
  
clinic.assign[i] <- 
  paste(factor1[clinic.assign.array.pos[i]], " -- ", factor2[clinic.assign.array.pos[i]])

}
#read in NU file sent from Jeff
library(readr)
SOMNUS.data <- read_csv("Downloads/SOMNUS_Randomization_10222024.csv")
#label columns for final file
col.labels <- c("clinic", "assignment", "rate", "low=0/high=1")
#construct final dataframe
random.assign.df <- data.frame(cbind(SOMNUS.data$deid, clinic.assign, SOMNUS.data$rate, SOMNUS.data$rx_high))
colnames(random.assign.df) <- col.labels
#write results to computer
write.csv(random.assign.df, file="randomizationSOMNUS.csv")