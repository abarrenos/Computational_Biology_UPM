# This script simulates the evolution of a population of N genomes comprised
# by (n_loci) number of loci.

N <- 1000
n_loci <- 10
u <- 1e-4      # Mutation rate per genome/individual
days <- 50000
t_step <- 200

pop <- matrix(0, n_loci, N)            # Population matrix
S1 <- matrix(0, n_loci, days/t_step)   # single substitutions record
S2 <- matrix(0, n_loci, days/t_step)   # double substitutions record
S3 <- matrix(0, n_loci, days/t_step)   # triple substitutions record
S4 <- matrix(0, n_loci, days/t_step)   # quadruple substitutions record

obs <- 0


# Algorithm for growth and mutation

for (day in (1:days)){
  
  # Reproduction
  
  offspring <- sample(x = 1:N, size = N/2, replace = FALSE)     # Sample 50 random genomes
  pop <- cbind(pop[,offspring], pop[,offspring], deparse.level = 0)   # Duplicate them
  
  # Mutation
  
  mutations <- rpois(1, N*u)                  # Calculate number of mutations. Poisson distribution with mean n of mutations lambda = N*u = 1e-3
  ind <- round(runif(mutations, 1, N))        # Allocate mutations across individuals. Prob of mutation is uniformly distributed through individuals and genes
  gen <- round(runif(mutations, 1, n_loci))   # Allocate mutations across genes.
  coord <- rbind(cbind(c(gen), c(ind)))       # Introduce the mutations in the matrix.
  pop[coord] <- pop[coord]+1
  
  if (day %in% seq(from = 1, to = days, by = t_step)){
    obs <- obs+1
    for (i in (1:n_loci)){
      S1[i, obs] <- sum(pop[i,]==1)
      S2[i, obs] <- sum(pop[i,]==2)
      S3[i, obs] <- sum(pop[i,]==3)
      S4[i, obs] <- sum(pop[i,]==4)
    }
  }
}

#print(obs)
#print(dim(S1))

plot(x = 1:obs, y = S1[1,1:obs], type = "l", ylim = c(1,N))
for (i in (1:n_loci)){
  lines(1:obs, S1[i, 1:obs], col = i, lwd = 1.5)
  lines(1:obs, S2[i, 1:obs], col = i, lwd = 1.5)
  lines(1:obs, S3[i, 1:obs], col = i, lwd = 1.5)
  lines(1:obs, S4[i, 1:obs], col = i, lwd = 1.5)
}

(substitutions <- rbind(S1[,obs], S2[,obs], S3[,obs], S4[,obs]))


