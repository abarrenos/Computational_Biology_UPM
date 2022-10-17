# This script simulates the evolution of a population of N
# genomes comprised by (n_loci) number of loci.

N <- 100
n_loci <- 10
u <- 1e-2      # Mutation rate per genome
days <- 5000

pop <- matrix(0, n_loci, N)         # Population matrix
S1 <- matrix(0, n_loci, days/200)   # single substitutions record
S2 <- matrix(0, n_loci, days/200)   # double substitutions record
S3 <- matrix(0, n_loci, days/200)   # triple substitutions record
S4 <- matrix(0, n_loci, days/200)   # quadruple substitutions record

obs <- 0

# Algorithm for growth and mutation

for (day in 1:days){
  
  # Reproduction
  
  offspring <- sample(x = 1:N, size = N/2, replace = TRUE)
  pop <- cbind(pop[,offspring], pop[,offspring], deparse.level = 0)
  
  # Mutation
  
  mutations <- rpois(1, N*u)                  # Calculate number of mutations
  ind <- round(runif(mutations, 1, N))        # Allocate mutations across individuals
  gen <- round(runif(mutations, 1, n_loci))   # Allocate mutations across genes.
  coord <- rbind(cbind(c(gen), c(ind)))       # Introduce the mutations in the matrix.
  pop[coord] <- pop[coord]+1

  if (day %in% seq(1, 50000, by = 200)){
    obs <- obs+1
    for (i in n_loci){
      S1[i, obs] <- sum(pop[i,]==1)
      S2[i, obs] <- sum(pop[i,]==2)
      S3[i, obs] <- sum(pop[i,]==3)
      S4[i, obs] <- sum(pop[i,]==4)
    }
  }
}

(substitutions <- rbind(S1[,obs], S2[,obs], S3[,obs], S4[,obs]))


plot(x = 1:obs, y = S1[1,1:obs], type = "l", ylim = c(1,N))
for (i in i:n_loci){
  lines(1:obs, S1[i, 1:obs], col = 1, lwd = 1.5)
  lines(1:obs, S2[i, 1:obs], col = 2, lwd = 1.5)
  lines(1:obs, S3[i, 1:obs], col = 3, lwd = 1.5)
  lines(1:obs, S4[i, 1:obs], col = 4, lwd = 1.5)
}



