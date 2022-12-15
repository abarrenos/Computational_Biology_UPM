
# Exercise #1: We have seen in class a computer model that simulates the simple evolution of a
# population of N genomes comprised by neutral genes (see "Neutral evolution predicts the clock").
# Refurbishing this code, estimate empirically how the time to fixation (tfix) and fixation probability (pfix)
# of a single neutral allele changes with:
#   a) population size (N)
#   b) bottleneck size

# Tips:
#   i.  Follow the pseudocode scheme discussed in class. Make sure to seed your population with just
#       one mutation in just one individual, and turn off mutation rate.
#   ii. Suggested range for parameters: days: 10,000; replicates: 100,000; N: 100, 300, 1000 and
#       3000; bottleneck: 1, 0.5, 0.25.
#   iii.To estimate the relationships between tfix and pfix with population and bottleneck size, try fitting
#       a linear or a log-linear model (see ?lm).
#   iv. An acceptable answer will include: the estimated values, the scripts, 2-3 clear plots and a
#       critical discussion of the results (e.g., limitations of the computer code, possible extensions, etc.) .

# We define a function that returns the number of fixation that occur over a certain
# number of replicates from a single random mutation in a population with size N.

library(dplyr)
library(ggplot2)

mutation_fix <- function(N, bottleneck_size, replicates, days){
    
  fixations <- 0                # Counter of fixed genes.
  time_of_fixation <- c()       # Array of times of fixation. 
  
  for (i in 1:replicates){
    
  # Create population and allocate the mutation in a random locus in one random individual.
    
    pop <- matrix(0, 1, N)            # Population matrix rep(x = 0, times = N) 
    
    ind <- round(runif(n = 1, min = 1, max = N))        # Allocate mutation across individuals.
    pop[,ind] <- 1                                      # Introduce the mutation
  
    
  ##### ------------------------ Replication loop ------------------------- #####
  
    for (day in 1:days){
      
  # Sample half of the population genomes and duplicate them to obtain next generation of individuals.
  # A population bottleneck reduces the variation in the gene pool of a population, thus a smaller
  # subset of genetic variants is passed to the future generations. We can introduce this effect by
  # reducing the sampling size by a constant factor (bot_size).
    
      (offspring <- sample(x = (1:N), size = bottleneck_size*N, replace = TRUE))
      
      (pop_off <- pop[,offspring])   # Sample population genomes
  
  # Duplicate the offspring until reaching initial population size.
  
      pop <- rbind(rep(pop_off, 1/bottleneck_size))
  
  # Break the loop if the mutation is fixed or lost. If fixed, annotate the
  # mutation in fixation array.
  
      if (sum(pop == 1) == 0){ #cat("\n")
        break }
      
      if (sum(pop == 1) == N){ #cat("\n")
        fixations <- fixations + 1
        time_of_fixation <- append(x = time_of_fixation, value = day)
        break }
    }
  }
  return(c(fixations, time_of_fixation))
}

### ----------------- Simulate the different conditions ------------------ ###

N_list <- c(100, 300, 1000, 3000)   # Number of individual
bot_size_list <- c(0.25, 0.5, 1)    # Bottleneck size 0 < X < 1
replicates <- 100000                # Number of replicates per condition
days <- 10000                       # Number of days (interpreted as generations)

# Save the data in matrixes
fixations_matrix <- matrix(0, nrow = length(bot_size_list), ncol = length(N_list))
fix_time_avg_matrix <- matrix(0, nrow = length(bot_size_list), ncol = length(N_list))
fix_time_std_matrix <- matrix(0, nrow = length(bot_size_list), ncol = length(N_list))


for (j in 1:length(N_list)){
  for (i in 1:length(bot_size_list)){
    
    (result <- mutation_fix(N_list[j], bot_size_list[i], replicates, days))
    
    (fixations <- result[1])
    (fix_time <- result[-1])
    
    fixations_matrix[i,j] <- fixations
    fix_time_avg_matrix[i,j] <- mean(fix_time)
    fix_time_std_matrix[i,j] <- sd(fix_time)
  }
}

colnames(fixations_matrix) = colnames(fix_time_avg_matrix) = colnames(fix_time_std_matrix) <- c("100", "300", "1000", "3000")
rownames(fixations_matrix) = rownames(fix_time_avg_matrix) = rownames(fix_time_std_matrix) <- c("0.25", "0.5", "1")

print(fixations_matrix)
print(fix_time_avg_matrix)
print(fix_time_std_matrix)

#save(fixations_matrix, fix_time_avg_matrix, fix_time_std_matrix, file = "fixation_data.Rdata")
#load(file = "fixation_data.Rdata", verbose = TRUE)
 
# We observe the number of fixations is correlated with N size, but not with 
# bottleneck size.
summary(fixations_matrix)
summary(t(fixations_matrix))

# The time of fixation is correlated with both N size and bottleneck size.
summary(fix_time_avg_matrix)
summary(t(fix_time_avg_matrix))


### -------------------- Box Plots ------------------------ ###

# Prepare data in matrices for representation
fixation_rate <- fixations_matrix/100000
log_fixation_rate <- log10(fixation_rate)
log_fix_time_avg <- log10(fix_time_avg_matrix)

# Represent the fixation rate as a function of population size.
boxplot(fixation_rate, xlab = "Population size (N)", ylab="Fixation probability")
boxplot(log_fixation_rate, xlab = "Population size (N)", ylab="Fixation probability (Log10)")

# Represent the fixation rate as a function of bottleneck size.
boxplot(fixation_rate, use.cols = FALSE, xlab = "Bottleneck size", ylab="Fixation probability")
boxplot(log_fixation_rate, use.cols = FALSE, xlab = "Bottleneck size", ylab="Fixation probability (Log10)")

# Represent the time to fixation as a function of population size.
boxplot(fix_time_avg_matrix, xlab = "Population size (N)", ylab="Time to fixation")
boxplot(log_fix_time_avg, xlab = "Population size (N)", ylab="Time to fixation (Log10)")

# Represent the time to fixation as a function of bottleneck size.
boxplot(fix_time_avg_matrix, use.cols = FALSE, xlab = "Bottleneck size", ylab="Time to fixation")
boxplot(log_fix_time_avg, use.cols = FALSE, xlab = "Bottleneck size", ylab="Time to fixation (Log10)")


### -------------------- Create Data Frames  ------------------- ###

# I will convert the data into a Data Frame to facilitate model building and representation.

col_N <- c(rep(100,3),rep(300,3),rep(1000,3),rep(3000,3))
col_bot <- rep(c(0.25,0.5,1), 4)
col_fixrate <- c(fixation_rate)
col_total_fix <- c(fixations_matrix)
col_fixtime <- c(fix_time_avg_matrix)
(df <- data.frame(col_N, col_bot, col_fixrate, col_fixtime))
#write.csv(df, file = "fixrate_df.csv", row.names = FALSE)
#df <- read.csv("fixrate_df.csv", header = TRUE)

# Save DataFrame with total number of fixations for presentation in the written report.
write.csv(data.frame(col_N, col_bot, col_total_fix, col_fixtime), file = "total_fix_df.csv", row.names = FALSE)
read.csv(file = "total_fix_df.csv", col.names = c("Population size", "Bottleneck size", "Total Fixations", "Average fixation time"))


### -------------------- Linear regression model ------------------- ###

# Logarithmic regression model N size - log(fixation rate)
lm_N_fixrate <- lm(log10(df$col_fixrate) ~ df$col_N)
a <- c(lm_N_fixrate$coefficients)
(a_r2 <- summary(lm_N_fixrate)$r.squared)
(a_corr <- cor(x = log10(df$col_fixrate), y = df$col_N, method = "spearman"))

# Linear regression model bottleneck size - fixation rate
lm_bot_fixrate <- lm(df$col_fixrate ~ df$col_bot)
b <- c(lm_bot_fixrate$coefficients)
(b_r2 <- summary(lm_bot_fixrate)$r.squared)
(b_corr <- cor(x = df$col_fixrate, y = df$col_bot, method = "spearman"))

# Linear regression model N size - fixation time
lm_N_fixtime <- lm(df$col_fixtime ~ df$col_N)
c <- c(lm_N_fixtime$coefficients)
(c_r2 <- summary(lm_N_fixtime)$r.squared)
(c_corr <- cor(x = df$col_fixtime, y = df$col_N, method = "spearman"))

# Logarithmic regression model bottleneck size - log(fixation time)
lm_bot_fixtime <- lm(log10(df$col_fixtime) ~ df$col_bot)
d <- c(lm_bot_fixtime$coefficients)
(d_r2 <- summary(lm_bot_fixtime)$r.squared)
(d_corr <- cor(x = log10(df$col_fixtime), y = df$col_bot, method = "spearman"))


### ------------------- Scatter Plots -------------------- ###

# Represent the fixation rate as a function of population size.
plot(df$col_N, df$col_fixrate, xlab = "Population size (N)", ylab="Fixation probability")
plot(df$col_N, log10(df$col_fixrate), xlab = "Population size (N)", ylab="Fixation probability (log10)")

# Regression line
abline(coef = a)
legend("topright", legend = paste("R² =", round(a_r2, 2), "\nρ = ", round(a_corr,2)), bty = "n", cex = 0.9)

# Represent the fixation rate as a function of bottleneck size.
plot(df$col_bot, log10(df$col_fixrate), xlab = "Bottleneck size", ylab="Fixation probability (log10)")
plot(df$col_bot, df$col_fixrate, xlab = "Bottleneck size", ylab="Fixation probability")

# Regression line
abline(coef = b)
legend("right", legend = paste("R² =", round(b_r2, 2), "\nρ = ", round(b_corr,2)), bty = "n", cex = 0.9)

# Represent the time to fixation as a function of population size.
plot(df$col_N, log10(df$col_fixtime), xlab = "Population size (N)", ylab="Time to fixation (log10)")
plot(df$col_N, df$col_fixtime, xlab = "Population size (N)", ylab="Time to fixation")

# Regression line 
abline(coef = c)
legend("topleft", legend = paste("R² =", round(c_r2, 2), "\nρ = ", round(c_corr,2)), bty = "n", cex = 0.9)

# Represent the time to fixation as a function of bottleneck size.
plot(df$col_bot, df$col_fixtime, xlab = "Bottleneck size", ylab="Time to fixation")
plot(df$col_bot, log10(df$col_fixtime), xlab = "Bottleneck size", ylab="Time to fixation (log10)")

# Regression line 
abline(coef = d)
legend("bottomright", legend = paste("R² =", round(d_r2, 2), "\nρ = ", round(d_corr,2)), bty = "n", cex = 0.9)


### -------------------- Save the figures ----------------------- ###

jpeg(file="Population_size.jpeg", width = 10, height = 5, units = 'in', res = 500)
graph_param <- par(mfrow=c(1,2))

plot(df$col_N, log10(df$col_fixrate), xlab = "Population size (N)", ylab="Fixation probability (log10)", col = rep(c(1, 2, 4),4))
abline(coef = a)
legend("topright", legend = paste("R² =", round(a_r2, 2), "\nρ = ", round(a_corr,2)), bty = "n", cex = 0.9)
plot(df$col_N, df$col_fixtime, xlab = "Population size (N)", ylab="Time to fixation", col = rep(c(1, 2, 4),4))
abline(coef = c)
legend("topleft", legend = paste("R² =", round(c_r2, 2), "\nρ = ", round(c_corr,2)), bty = "n", cex = 0.9)
mtext("Effect of population size on fixation of neutral mutations", side = 3, line = -2, outer = TRUE)

par(graph_param)
dev.off()

jpeg(file="Bottleneck_size.jpeg", width = 10, height = 5, units = 'in', res = 500)
graph_param <- par(mfrow=c(1,2))

#plot(df$col_bot, df$col_fixrate, xlab = "Bottleneck size", ylab="Fixation probability")
boxplot(fixation_rate, use.cols = FALSE, xlab = "Bottleneck size", ylab="Fixation probability")
legend("topright", legend = paste("\nR² =", round(b_r2, 2), "\nρ = ", round(b_corr,2)), bty = "n", cex = 0.8)
boxplot(log_fix_time_avg, use.cols = FALSE, xlab = "Bottleneck size", ylab="Time to fixation (Log10)")
legend("bottomright", legend = paste("R² =", round(d_r2, 2), "\nρ = ", round(d_corr,2), "\n"), bty = "n", cex = 0.8)

mtext("Effect of bottleneck size on fixation of neutral mutations", side = 3, line = -2, outer = TRUE)

par(graph_param)
dev.off()





