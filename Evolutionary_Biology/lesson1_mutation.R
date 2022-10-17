## Modelling Lurina-Delbrück experiment

x <- matrix(0,1,2)

mut_rate <- 5e-9    # Mutation rate
w <- 1              # Mutant fitness (survival), ranging between 0 and 1.
N_0 <- 1000         # Initial population
gen <- 20           # Number of generations
replicates <- 5000

mutant_matrix <- matrix(0, replicates, 1)   # Retrieve proportion of mutants
for (i in 1:replicates){                    # in each replicate.
  
  x[1,1] <- N_0
  x[1,2] <- 0
  for (g in 1:gen){
    x[1,1] <- x[1,1]*2            # Each generation, both wild-type and mutant
    x[1,2] <- x[1,2]*2*w          # population duplicates.
  
    ## Poisson distribution: discrete distribution with mean lambda.
      
    mutations <- rpois(1, mut_rate*x[1,1])
    x[1,1] <- x[1,1]-mutations
    x[1,2] <- x[1,2]+mutations
    #print(c(g, x))
    mutant_matrix[i, 1] <- x[1, 2]/sum(x)
  }
}

hist(log10(mutant_matrix), breaks = 50, col = "lightblue", ylim = c(0,1000))

## When fitness of the mutants is near to 1, there is a great variety in the
## number of mutant individuals appearing in the last generation between
## the different replicates.
## - Ex: {0, 4, 0, 4, 3, 135, 2, 6, 69, 0, 3}
##
## When fitness is near to 0, the proportion of mutants appearing in the
## last generation in similar to the mutation rate, and its distribution is
## a Poisson distribution.
##
## Thus, if Luria and Delbrück had studied a mutation with much less fitness
## (survival), they would have reached an incorrect conclusion, since the
## number of mutated cells would have been much more homogeneous and the
## statistics would have followed a regular Poisson distribution.

