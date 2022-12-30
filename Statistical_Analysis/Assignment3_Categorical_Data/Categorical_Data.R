## Update R
#install.packages("installr")
#library(installr)
#updateR()

# Load data
load("DatosCTC.RData")
df <- DAT4
str(df)

# Variable status can take the levels censored (alive) or event (dead)
# Create 2x2 tables to relate categorical values to death

# Then Pearson's chi-squared test is performed of the null hypothesis that
# the joint distribution of the cell counts in a 2-dimensional contingency
# table is the product of the row and column marginals.

(sex <- table(df$sexo, df$status))      # p-value = 0.5189
chisq.test(sex)
(ctc <- table(df$CTC, df$status))       # p-value = 0.1857
chisq.test(ctc)
(mutant <- table(df$Mutant, df$status)) # p-value = 0.0806
chisq.test(mutant)
(group <- table(relevel(df$group, ref = "R"), df$status))   # p-value = 0.001758
chisq.test(group)

save(sex, ctc, mutant, group, file = "./tables.RData")    # Save contingency table

# Fisher test to calculate Odds Ratio and significance. Odds Ratio can only be calculated
# for 2x2 tables, therefore no OR will be returned for the factor "group".

(fisher_sex <- fisher.test(sex))$p.value
(fisher_ctc <- fisher.test(ctc))
(fisher_mutant <- fisher.test(mutant))
(fisher_group <- fisher.test(group))

fisher  <- cbind(rbind(fisher_sex$estimate, fisher_ctc$estimate, fisher_mutant$estimate), 
          rbind(fisher_sex$conf.int, fisher_ctc$conf.int, fisher_mutant$conf.int, deparse.level = 2),
          rbind(fisher_sex$p.value, fisher_ctc$p.value, fisher_mutant$p.value))

colnames(fisher) <- c("odds raio", "CI low", "CI high", "p-value")

fisher


### --------- Build logistic models for each pair of variables --------- ###

## SEX - STATUS
sex_lm = glm(status ~ sexo, data = df, family = binomial)
summary(sex_lm)
# Convert model coefficients to Odds Ratio
round(exp(cbind(OR = coef(sex_lm), confint(sex_lm))), 3)


## CTC - STATUS
ctc_lm = glm(status ~ CTC, data = df, family = binomial)
summary(ctc_lm)
round(exp(cbind(OR = coef(ctc_lm), confint(ctc_lm))), 3)


## MUTANT - STATUS
mutant_lm = glm(status ~ Mutant, data = df, family = binomial)
summary(mutant_lm)
round(exp(cbind(OR = coef(mutant_lm), confint(mutant_lm))), 3)


## GROUP - STATUS
df$group <- relevel(df$group, ref ="R")   # Relevel "group" factor
group_lm = glm(status ~ group, data = df, family = binomial)
summary(group_lm)
round(exp(cbind(OR = coef(group_lm), confint(group_lm))), 3)

# Interpret results
# https://quantifyinghealth.com/logistic-regression-in-r-with-categorical-variables/

### --------------------- Other Packages ------------------------ ###

## EPITOOLS -> https://jarrettmeyer.com/2019/07/23/odds-ratio-in-r

#install.packages("epitools")
library("epitools")
oddsratio(group)
riskratio(group)


### ----------- Model Selection Using Akaike (AIC) ------------ ###

# AIC value allows us to compare between different models and choose the
# most suitable one, which has the lower AIC.

step(glm(status ~ sexo + CTC + Mutant + group, data = df, family = "binomial"))




### ----------- Model Selection Using Cross Validation ------------ ###

# Load required library
#install.packages("caret")
library(caret)

# Define the model formula
formula <- status ~ sexo + CTC + Mutant + group

# Set up the logistic regression model using the "glm" function
model <- glm(formula, data = df, family = "binomial")

# Set up the cross-validation control using the "trainControl" function
(control <- trainControl(method = "cv", number = 3, savePredictions = TRUE))

# Train the model using the "train" function and the specified cross-validation control
trained_model <- train(formula, data = df, method = "glm", trControl = control, family = "binomial")

# Print the results of the cross-validation
print(trained_model)


library(Rcmdr)











