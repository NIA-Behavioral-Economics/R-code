#install and load true random number package
install.packages("random", dependencies = T)
library(random)

#install and load DOE package
install.packages("DoE.base")
library(DoE.base)

## Create a dataframe for a fractional factorial design with 6 factors with 3 levels each
## this is just an example you will have to work with DOE.base to get the right factors 
## and levels
  frac.factorial.list <- oa.design(nfactors=6, nlevels=3)
  
## Turn this list into a data frame
  df2 <- data.frame(sapply(frac.factorial.list,c))
  
## Show Cards
df2

# Use combn() fxn to find all pairs (n choose 2)
# Use t() fxn to transpose matrix so that rows are pairs
paired.comp <-data.frame(t(combn(1:18,2)))
# add names to columns
names(paired.comp)[1] <- "Choice_A"
names(paired.comp)[2] <- "Choice_B"

# add a random order column with which to sort on  
random.order <- randomSequence(min=1, max=length(paired.comp$Choice_A), check=TRUE)
# make random order vector a data.frame
df3 <- data.frame(random.order)
# add random.order to data frame
paired.comp <- cbind(paired.comp, df3)

# Sort pairs by random number
# display paired comparisons
paired.comp <- paired.comp[order(random.order),]
names(paired.comp)[3] <- "Random.Order"
paired.comp

# Display a list of paired comparisons given to each subject 
# that completes the fractional design
# for each divisor of the design

# Function that finds divisors
divisors <- function(x){
  #  Vector of numberes to test against
  y <- seq_len(x)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  y[ x%%y == 0 ]
}

# Put divisors of design into a vector

design.divisors <- divisors(length(paired.comp$Choice_A))

for(i in 1:length(design.divisors)){
  z <- expand.grid(seq(1:(length(paired.comp$Choice_A)/design.divisors[i])),seq(1:design.divisors[i]))
  
  paired.comp <- cbind(paired.comp, z$Var2)
  colnames(paired.comp)[length(paired.comp)] <- paste("SsDesign", i, sep = "")
}
paired.comp