library(nimble)
library(igraph)
library(coda)
library(R6)
library(ggplot2)
library(cowplot)

# i = student number (= 1, 2, ..., 73)
# Y_i = nr of comprehensive answers on winter exam
# j = quiz number (= 1 (quiz), 2 (nimble))
# z_ij = standardised continuous assesment mark of ith student
#        for quiz j
# x_i = indicator for MT5731 (=1) or MT4531 (=0) of student i

# Assumed model:
# Y_i|B ~ Poisson(lda_i) where B = (B0, B1, B2, B3)
# log(lda_i) = B0 + B1*z_i1 + B2*z_i2 + B3*x_i

####################################################
############ DOWNLOADING AND FIXING DATA ###########
####################################################
load("~/R/Advanced Bayesian Inference/NIMBLE Project (11th)/ExamData0.RData")
data <- ExamData0

# calculating mean of available z_nimb data
z_nimb_nona <- data$z_Nimb[!is.na(data$z_Nimb)]
z_nimb_mean <- mean(z_nimb_nona)

# replacing NA values with mean
data$z_Nimb[is.na(data$z_Nimb)] <- rep(z_nimb_mean, 3)
#-------------------------------------------------------------------------------



####################################################
################## QUESTION 1 CODE #################
####################################################

# // Part 1 //
# Producing a histogram of response Y_i stratified by indicator x_i
Y_5 <- data$y[data$x_module == 1]
Y_4 <- data$y[data$x_module == 0]

# calculating means
mean_Y_5 <- mean(Y_5)
mean_Y_4 <- mean(Y_4)

# see through colors
c1 <- rgb(173, 216, 230, max=255, alpha=90)
c2 <- rgb(216, 238, 192, max=255, alpha=90)

# plotting stratified histogram
hist(Y_5, breaks=8, main='Number of comprehensive answers (Y) in MT5731', 
     col=c1, ylim=c(0, 35), xlim=c(0,10), xlab='Y value', ylab='Count')
hist(Y_4, breaks=8, main='Number of comprehensive answers in MT4531', 
     col=c2, add = TRUE)
legend(4, 30, legend=c('5000 level', '4000 level'), fill=c(c1, c2))


# comparing kernel density estimates of the stratified data
# to the sample's implied poisson distribution
par(mfrow=c(1, 2))

# 5000 level
plot(density(Y_5), ylim=c(0, 0.25), xlim=c(0, 10), col='black',
     xlab='Y value', main='MT5731')
lines(dpois(x=0:10, lambda=mean_Y_5), type='l', col='blue')
legend(4.4, 0.258, legend=c('sample density', 'poisson estimate'), 
       fill=c('black', 'blue'))

# 4000 level
plot(density(Y_4),ylim=c(0, 0.35), xlim=c(0, 10), col='black', 
     xlab='Y value', main = 'MT4531')
lines(dpois(x=0:10, lambda=mean_Y_4), type='l', col='blue')
legend(4.4, 0.361, legend=c('sample density', 'poisson estimate'), 
       fill=c('black', 'blue'))



# // Part 2 //
# Sort data by z_Quiz and z_Nimb separately
z_Q <- data$z_Quiz
z_N <- data$z_Nimb

sort_data_zQ <- data[order(z_Q), ]
sort_data_zN <- data[order(z_N), ]

cut_delta <- length(z_Q)/5 #=14.6 so use 15, 15, 15, 14, 14
# z_Quiz quintiles
Yq1_zQ <- sort_data_zQ$y[1:15] 
Yq2_zQ <- sort_data_zQ$y[16:30] 
Yq3_zQ <- sort_data_zQ$y[31:45] 
Yq4_zQ <- sort_data_zQ$y[46:59] 
Yq5_zQ <- sort_data_zQ$y[60:73] 

# z_Quiz means
mean_q1zQ <- mean(Yq1_zQ)
mean_q2zQ <- mean(Yq2_zQ)
mean_q3zQ <- mean(Yq3_zQ)
mean_q4zQ <- mean(Yq4_zQ)
mean_q5zQ <- mean(Yq5_zQ)

# z_Quiz variances
var_q1zQ <- var(Yq1_zQ)
var_q2zQ <- var(Yq2_zQ)
var_q3zQ <- var(Yq3_zQ)
var_q4zQ <- var(Yq4_zQ)
var_q5zQ <- var(Yq5_zQ)

# z_Nimb quintiles
Yq1_zN <- sort_data_zN$y[1:15] 
Yq2_zN <- sort_data_zN$y[16:30] 
Yq3_zN <- sort_data_zN$y[31:45] 
Yq4_zN <- sort_data_zN$y[46:59] 
Yq5_zN <- sort_data_zN$y[60:73] 

# z_Nimb means
mean_q1zN <- mean(Yq1_zN)
mean_q2zN <- mean(Yq2_zN)
mean_q3zN <- mean(Yq3_zN)
mean_q4zN <- mean(Yq4_zN)
mean_q5zN <- mean(Yq5_zN)

# z_Nimb variances
var_q1zN <- var(Yq1_zN)
var_q2zN <- var(Yq2_zN)
var_q3zN <- var(Yq3_zN)
var_q4zN <- var(Yq4_zN)
var_q5zN <- var(Yq5_zN)


# Table
df_1b <- data.frame(q1=c(mean_q1zQ, var_q1zQ, mean_q1zN, var_q1zN), 
                    q2=c(mean_q2zQ, var_q2zQ, mean_q2zN, var_q2zN), 
                    q3=c(mean_q3zQ, var_q3zQ, mean_q3zN, var_q3zN), 
                    q4=c(mean_q4zQ, var_q4zQ, mean_q4zN, var_q4zN), 
                    q5=c(mean_q5zQ, var_q5zQ, mean_q5zN, var_q5zN))

rownames(df_1b) <- c('Y_sQ mean', 'Y_sQ var', 'Y_sN mean', 'Y_sN var')

df_1b


# // Part 3 //
# z_Quiz log means
logm_q1zQ <- log(mean_q1zQ)
logm_q2zQ <- log(mean_q2zQ)
logm_q3zQ <- log(mean_q3zQ)
logm_q4zQ <- log(mean_q4zQ)
logm_q5zQ <- log(mean_q5zQ)

#z_Nimb log means
logm_q1zN <- log(mean_q1zN)
logm_q2zN <- log(mean_q2zN)
logm_q3zN <- log(mean_q3zN)
logm_q4zN <- log(mean_q4zN)
logm_q5zN <- log(mean_q5zN)

#z_Quiz and z_Nimb log mean quantile plots
logmean_data <- data.frame(quantile=1:5, 
                           zQ = c(logm_q1zQ, logm_q2zQ, logm_q3zQ, 
                                  logm_q4zQ, logm_q5zQ),
                           zN = c(logm_q1zN, logm_q2zN, logm_q3zN, 
                                  logm_q4zN, logm_q5zN))
p1 <- ggplot(logmean_data, aes(quantile, zQ)) + geom_point()
p2 <- ggplot(logmean_data, aes(quantile, zN)) + geom_point()

title <- ggdraw() + draw_label('log mean plots for z_Quiz and z_Nimb quantiles',
                               fontface='bold')
plots <- plot_grid(p1, p2, ncol=2)

plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1))


# Stratifying log means by module
# z_Quiz 5000 level
Y_5_q1zQ <- Yq1_zQ[sort_data_zQ$x_module[1:15]==1]
mean_5_q1zQ <- mean(Y_5_q1zQ)
logm_5_q1zQ <- log(mean_5_q1zQ)

Y_5_q2zQ <- Yq2_zQ[sort_data_zQ$x_module[16:30]==1]
mean_5_q2zQ <- mean(Y_5_q2zQ)
logm_5_q2zQ <- log(mean_5_q2zQ)

Y_5_q3zQ <- Yq3_zQ[sort_data_zQ$x_module[31:45]==1]
mean_5_q3zQ <- mean(Y_5_q3zQ)
logm_5_q3zQ <- log(mean_5_q3zQ)

Y_5_q4zQ <- Yq4_zQ[sort_data_zQ$x_module[46:59]==1]
mean_5_q4zQ <- mean(Y_5_q4zQ)
logm_5_q4zQ <- log(mean_5_q4zQ)

Y_5_q5zQ <- Yq5_zQ[sort_data_zQ$x_module[60:73]==1]
mean_5_q5zQ <- mean(Y_5_q5zQ)
logm_5_q5zQ <- log(mean_5_q5zQ)

# z_Quiz 4000 level
Y_4_q1zQ <- Yq1_zQ[sort_data_zQ$x_module[1:15]==0]
mean_4_q1zQ <- mean(Y_4_q1zQ)
logm_4_q1zQ <- log(mean_4_q1zQ)

Y_4_q2zQ <- Yq2_zQ[sort_data_zQ$x_module[16:30]==0]
mean_4_q2zQ <- mean(Y_4_q2zQ)
logm_4_q2zQ <- log(mean_4_q2zQ)

Y_4_q3zQ <- Yq3_zQ[sort_data_zQ$x_module[31:45]==0]
mean_4_q3zQ <- mean(Y_4_q3zQ)
logm_4_q3zQ <- log(mean_4_q3zQ)

Y_4_q4zQ <- Yq4_zQ[sort_data_zQ$x_module[46:59]==0]
mean_4_q4zQ <- mean(Y_4_q4zQ)
logm_4_q4zQ <- log(mean_4_q4zQ)

Y_4_q5zQ <- Yq5_zQ[sort_data_zQ$x_module[60:73]==0]
mean_4_q5zQ <- mean(Y_4_q5zQ)
logm_4_q5zQ <- log(mean_4_q5zQ)

# z_Nimb 5000 level
Y_5_q1zN <- Yq1_zN[sort_data_zN$x_module[1:15]==1]
mean_5_q1zN <- mean(Y_5_q1zN)
logm_5_q1zN <- log(mean_5_q1zN)

Y_5_q2zN <- Yq2_zN[sort_data_zN$x_module[16:30]==1]
mean_5_q2zN <- mean(Y_5_q2zN)
logm_5_q2zN <- log(mean_5_q2zN)

Y_5_q3zN <- Yq3_zN[sort_data_zN$x_module[31:45]==1]
mean_5_q3zN <- mean(Y_5_q3zN)
logm_5_q3zN <- log(mean_5_q3zN)

Y_5_q4zN <- Yq4_zN[sort_data_zN$x_module[46:59]==1]
mean_5_q4zN <- mean(Y_5_q4zN)
logm_5_q4zN <- log(mean_5_q4zN)

Y_5_q5zN <- Yq5_zQ[sort_data_zN$x_module[60:73]==1]
mean_5_q5zN <- mean(Y_5_q5zN)
logm_5_q5zN <- log(mean_5_q5zN)

# z_Nimb 4000 level
Y_4_q1zN <- Yq1_zN[sort_data_zN$x_module[1:15]==0]
mean_4_q1zN <- mean(Y_4_q1zN)
logm_4_q1zN <- log(mean_4_q1zN)

Y_4_q2zN <- Yq2_zN[sort_data_zN$x_module[16:30]==0]
mean_4_q2zN <- mean(Y_4_q2zN)
logm_4_q2zN <- log(mean_4_q2zN)

Y_4_q3zN <- Yq3_zN[sort_data_zN$x_module[31:45]==0]
mean_4_q3zN <- mean(Y_4_q3zN)
logm_4_q3zN <- log(mean_4_q3zN)

Y_4_q4zN <- Yq4_zN[sort_data_zN$x_module[46:59]==0]
mean_4_q4zN <- mean(Y_4_q4zN)
logm_4_q4zN <- log(mean_4_q4zN)

Y_4_q5zN <- Yq5_zN[sort_data_zN$x_module[60:73]==0]
mean_4_q5zN <- mean(Y_4_q5zN)
logm_4_q5zN <- log(mean_4_q5zN)

logm_data <- data.frame(quantile=1:5, zQ5=c(logm_5_q1zQ, logm_5_q2zQ, 
                                            logm_5_q3zQ, logm_5_q4zQ, 
                                            logm_5_q5zQ), 
                        zQ4=c(logm_4_q1zQ, logm_4_q2zQ, logm_4_q3zQ, 
                              logm_4_q4zQ, logm_4_q5zQ),
                        zN5=c(logm_5_q1zN, logm_5_q2zN, logm_5_q3zN, 
                              logm_5_q4zN, logm_5_q5zN), 
                        zN4=c(logm_4_q1zN, logm_4_q2zN, logm_4_q3zN, 
                              logm_4_q4zN, logm_4_q5zN))


# Plot all four panels 
zQ_5_plot <- ggplot(logm_data, aes(quantile, zQ5)) + geom_point()
zQ_4_plot <- ggplot(logm_data, aes(quantile, zQ4)) + geom_point()
zN_5_plot <- ggplot(logm_data, aes(quantile, zN5)) + geom_point()
zN_4_plot <- ggplot(logm_data, aes(quantile, zN4)) + geom_point()

title <- ggdraw() + draw_label("log mean plots for z_Quiz and z_Nimb quantiles 
                               stratified by module", fontface='bold')

grid <- plot_grid(zQ_5_plot, zQ_4_plot, zN_5_plot, zN_4_plot, ncol=2)

plot_grid(title, grid, ncol=1, rel_heights=c(0.1, 1))
#-------------------------------------------------------------------------------



####################################################
################## QUESTION 2 CODE #################
####################################################


# model from intro
# normal priors: mean = 0, sd = 10 for all B

# Perform Bayesian analysus if the data.
# Store all MCMC samples in object "combinedchains_MQ2" - use mcmc.list


# Part 1: Need to provide details for:
# 1. Number of iterations used
# 2. Convergence assessment (trace plot, )
# 3. Determination of burn-in
# 4. Number of samples collected (post convergence)
# 5. Summary of posterior results in terms of B parameters

# Part 2: Thorough interpretation of results using
# 1. posterior point estimates (of all parameters)
# 2. posterior interval estimates (of all parameters)
# CONSIDER TRANSFORMATIONS OF POSTERIOR DISTRIBUTIONS TO AID INTERPRETABILITY

# Calculate:
# Convergence/Burn-in: trace plot, BGR
# Effective sample size: ACF + fixes (batching/thinning), MC error (batching)
# AutoCorrelation, MC error (connects to sample size)

# // MODEL AND MCMC //
Q2Code <- nimbleCode({
  # Specify the likelihood
  for (i in 1:N) {
    
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0 + B1*z[i, 1] + B2*z[i, 2] + B3*x[i])
    
  }
  
  # Prior specification
  B0 ~ dnorm(0, 1/10^2)
  B1 ~ dnorm(0, 1/10^2)
  B2 ~ dnorm(0, 1/10^2)
  B3 ~ dnorm(0, 1/10^2)
  
})


# Constants of model
Q2Consts <- list(N=73)

# Data of model
Q2Data <- list(x = data$x_module, 
               z = matrix(c(data$z_Quiz, data$z_Nimb), nrow=73, ncol=2), 
               Y = data$y)

# Set of initial values before building model
Q2Inits <- list(B0=1, B1=1, B2=1, B3=1)

# Building the model
Q2 <- nimbleModel(code = Q2Code, name = 'Q2', constants = Q2Consts, 
                  data = Q2Data, inits = Q2Inits)

# Compiling the model
CQ2 <- compileNimble(Q2)

# Set up monitored quantities
Q2Conf <- configureMCMC(Q2, monitors = c('B0', 'B1', 'B2', 'B3'), 
                        print = TRUE)

# Build MCMC algorithm
Q2MCMC <- buildMCMC(Q2Conf)

# Compile MCMC chain
CQ2MCMC <- compileNimble(Q2MCMC, project = Q2)



# // POSTERIOR SAMPLES IN CODA FORMAT //
set.seed(1)
# Initial visualisation of plots
Q2Inits2 <- list(Q2Inits, list(B0=5, B1=5, B2=5, B3=5), 
                 list(B0=3, B1=3, B2=3, B3=3))
posterior <- runMCMC(CQ2MCMC, niter=51000, thin=1, nburnin=0,
                     summary=TRUE, samples=TRUE, nchains=3,
                     samplesAsCodaMCMC=TRUE, inits=Q2Inits2)

combinedchains_MQ2 <- mcmc.list(posterior$samples$chain1, 
                                posterior$samples$chain2, 
                                posterior$samples$chain3)

# Trace and posterior density plots for B0, B1, B2 and B3
plot(combinedchains_MQ2)

# ACF plots
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
autocorr.plot(posterior$samples$chain3)




# Thinning=10 to reduce AC and adding burn-in of 1000
set.seed(1)
posterior <- runMCMC(CQ2MCMC, niter=51000, thin=10, nburnin=1000,
                     summary=TRUE, samples=TRUE, nchains=3,
                     samplesAsCodaMCMC=TRUE, inits=Q2Inits2)

combinedchains_MQ2 <- mcmc.list(posterior$samples$chain1, 
                                posterior$samples$chain2, 
                                posterior$samples$chain3)

# Trace and posterior density plots for B0, B1, B2 and B3
plot(combinedchains_MQ2)

# ACF plots
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
autocorr.plot(posterior$samples$chain3)

# Gelman plots
gelman.diag(combinedchains_MQ2)
gelman.plot(combinedchains_MQ2)


# posterior summary
Q2summary <- posterior$summary$all.chains
Q2summary
summary(combinedchains_MQ2)

effectiveSize(combinedchains_MQ2)



# // Clearer density plots //
par(mfrow = c(2, 2))
# Plot densities and confidence intervals
# BO
#chain1
plot(density(posterior$samples$chain1[ , "B0"]), type = "l", 
     xlab = expression(B0), ylab = "Posterior densities", col='black', 
     xlim=c(-0.1, 1), ylim=c(0, 5), main='B0 density')
#chain2
lines(density(posterior$samples$chain2[ , "B0"]), type = "l", col='blue')
#chain3
lines(density(posterior$samples$chain3[ , "B0"]), type = "l", col='purple')
#CI, mean
lines(x = rep(Q2summary[1,4], 2), y = c(0, 5), col='red')
lines(x = rep(Q2summary[1,5], 2), y = c(0, 5), col='red')
lines(x = rep(Q2summary[1,1], 2), y = c(0, 5), col='orange')
# overlaying implied normal
x <- seq(-0.1, 1, 0.01)
lines(x, dnorm(x, Q2summary[1,1], Q2summary[1,3]), col='green')


#legend
legend(0.7, 5, legend=c('chain 1', 'chain 2', 'chain 3', '95% CI', 'mean'),
       fill=c('black', 'blue', 'purple', 'red', 'orange'))

# B1
plot(density(posterior$samples$chain2[ , "B1"]), type = "l", 
     xlab = expression(B1), ylab = "Posterior densities", col='black', 
     xlim=c(-0.1,1), ylim=c(0,5), main='B1 density')
lines(density(posterior$samples$chain1[ , "B1"]), type = "l", col='blue')
lines(density(posterior$samples$chain3[ , "B1"]), type = "l", col='purple')
lines(x = rep(Q2summary[2,4], 2), y = c(0, 5), col = 'red')
lines(x = rep(Q2summary[2,5], 2), y = c(0, 5), col='red')
lines(x = rep(Q2summary[2,1], 2), y = c(0, 5), col='orange')
x <- seq(-0.1, 1, 0.01)
lines(x, dnorm(x, Q2summary[2,1], Q2summary[2,3]), col='green')

legend(0.7, 5, legend=c('chain 1', 'chain 2', 'chain 3', '95% CI', 'mean'),
       fill=c('black', 'blue', 'purple', 'red', 'orange'))

#B2
plot(density(posterior$samples$chain2[ , "B2"]), type = "l", 
     xlab = expression(B2), ylab = "Posterior densities", col='black', 
     xlim=c(-0.1,1), ylim=c(0,5), main='B2 density')
lines(density(posterior$samples$chain1[ , "B2"]), type = "l", col='blue')
lines(density(posterior$samples$chain3[ , "B2"]), type = "l", col='purple')
lines(x = rep(Q2summary[3,4], 2), y = c(0, 5), col = 'red')
lines(x = rep(Q2summary[3,5], 2), y = c(0, 5), col='red')
lines(x = rep(Q2summary[3,1], 2), y = c(0, 5), col='orange')
x <- seq(-0.1, 1, 0.01)
lines(x, dnorm(x, Q2summary[3,1], Q2summary[3,3]), col='green')

legend(0.7, 5, legend=c('chain 1', 'chain 2', 'chain 3', '95% CI', 'mean'),
       fill=c('black', 'blue', 'purple', 'red', 'orange'))
#B3
plot(density(posterior$samples$chain2[ , "B3"]), type = "l", 
     xlab = expression(B3), ylab = "Posterior densities", col='black', 
     xlim=c(-0.1,1), ylim=c(0, 5), main='B3 density')
lines(density(posterior$samples$chain1[ , "B3"]), type = "l", col='blue')
lines(density(posterior$samples$chain3[ , "B3"]), type = "l", col='purple')
lines(x = rep(Q2summary[4,4], 2), y = c(0, 5), col = 'red')
lines(x = rep(Q2summary[4,5], 2), y = c(0, 5), col='red')
lines(x = rep(Q2summary[4,1], 2), y = c(0, 5), col='orange')
x <- seq(-0.1, 1, 0.01)
lines(x, dnorm(x, Q2summary[4,1], Q2summary[4,3]), col='green')

legend(0.7, 5, legend=c('chain 1', 'chain 2', 'chain 3', '95% CI', 'mean'),
       fill=c('black', 'blue', 'purple', 'red', 'orange'))
#-------------------------------------------------------------------------------



####################################################
################## QUESTION 3 CODE #################
####################################################
# PART 1
# Creating nice dataframe of samples to work with
posterior_c1samples <- data.frame(posterior$samples$chain1)
posterior_c2samples <- data.frame(posterior$samples$chain2)
posterior_c3samples <- data.frame(posterior$samples$chain3)

posterior_samples <- rbind(posterior_c1samples, posterior_c2samples)
posterior_samples <- rbind(posterior_samples, posterior_c3samples)

# All the betas are independent, and so each set of row values for the Betas can
# be used to estimate a lambda_24 sample

# Tranforming Betas into lambda_24 sample
# number of samples
n <- length(posterior_samples[,1])

# initialise empty vector to store values
lambdas <- c()

# getting lambdas
for (i in 1:n) {
  #Beta samples (B3 neglected as x_24 is zero)
  B0 <- posterior_samples[i, 1]
  B1 <- posterior_samples[i, 2]
  B2 <- posterior_samples[i, 3]
  
  # implied lambda_24
  lambda_24 <- exp(B0 + z_Q[24]*B1 + z_N[24]*B2) 
  
  # append
  lambdas <- c(lambdas, lambda_24)
}

# Computing posterior mean and median
# mean
mean_lambda_24 <- mean(lambdas)
mean_lambda_24

# median
med_lambda_24 <- median(lambdas)
med_lambda_24

# Computing the 90% symmetric CI
sort_lambdas <- sort(lambdas)
# want integral -inf -> lda_1 = 0.05
cut1 <- 0.05*n
lda_1 <- sort_lambdas[cut1]
lda_1

# want integral lda_2 -> +inf = 0.05
cut2 <- 0.95*n + 1
lda_2 <- sort_lambdas[cut2]
lda_2


# PART 2
# new MT5731 student: x = 1
# z_{stud,j} = (w_{stud,j} - E[w_{stud,j}])/SD[w_{stud,j}]
#            = 1 for j=1 as z_{stud,1} = E[w_{stud,1}] + SD[w_{stud,1}]
#            = 2 for j=2 as z_{stud,2} = E[w_{stud,2}] + 2*SD[w_{stud,2}]
stud_zQ <- 1 
stud_zN <- 2

# get lambdas then average
new_lambdas <- c()

# getting lambdas
for (i in 1:n) {
  #Beta samples (B3 neglected as x_24 is zero)
  B0 <- posterior_samples[i, 1]
  B1 <- posterior_samples[i, 2]
  B2 <- posterior_samples[i, 3]
  B3 <- posterior_samples[i, 4]
  
  # implied lambda_24
  lambda_24 <- exp(B0 + stud_zQ*B1 + stud_zN*B2 + B3) #x_module = 1
  
  # append
  new_lambdas <- c(new_lambdas, lambda_24)
}

# Getting implied poisson distribution
Y_samples <- c()

for (i in 1:n) {
  # get poisson sample
  pois_sample <- rpois(1, new_lambdas[i])
  
  # append
  Y_samples <- c(Y_samples, pois_sample)
}

# Calculate P(Y > 5) for Y ~ Poisson(stud_lambda) = proportion > 5
prob_g5 <- length(Y_samples[Y_samples > 5]) / n
prob_g5

#-------------------------------------------------------------------------------




####################################################
################## QUESTION 4 CODE #################
####################################################
# // 3 new priors //
Q4Data <- Q2Data # same data

# PRIOR 1: Betas with mean=5 and sd=5
Q4Code1 <- nimbleCode({
  # Specify the likelihood
  for (i in 1:N) {
    
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0 + B1*z[i, 1] + B2*z[i, 2] + B3*x[i])
    
  }
  
  # Prior specification
  B0 ~ dnorm(5, 1/5^2)
  B1 ~ dnorm(5, 1/5^2)
  B2 ~ dnorm(5, 1/5^2)
  B3 ~ dnorm(5, 1/5^2)
  
})


# Constants of model
Q4Consts1 <- list(N=73)

# Set of initial values before building model
Q4Inits1 <- list(B0=1, B1=1, B2=1, B3=1)

# Building the model
Q4_1 <- nimbleModel(code = Q4Code1, name = 'Q4_1', constants = Q4Consts1, 
                    data = Q4Data, inits = Q4Inits1)

# Compiling the model
CQ4_1 <- compileNimble(Q4_1)

# Set up monitored quantities
Q4Conf1 <- configureMCMC(Q4_1, monitors = c('B0', 'B1', 'B2', 'B3'), 
                         print = TRUE)

# Build MCMC algorithm
Q4MCMC_1 <- buildMCMC(Q4Conf1)

# Compile MCMC chain
CQ4MCMC_1 <- compileNimble(Q4MCMC_1, project = Q4_1)



# PRIOR 2: Betas with mean=-5, sd=2
Q4Code2 <- nimbleCode({
  # Specify the likelihood
  for (i in 1:N) {
    
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0 + B1*z[i, 1] + B2*z[i, 2] + B3*x[i])
    
  }
  
  # Prior specification
  B0 ~ dnorm(-5, 1/2^2)
  B1 ~ dnorm(-5, 1/2^2)
  B2 ~ dnorm(-5, 1/2^2)
  B3 ~ dnorm(-5, 1/2^2)
})


# Constants of model
Q4Consts2 <- list(N=73)

# Set of initial values before building model
Q4Inits2 <- list(B0=1, B1=1, B2=1, B3=1)

# Building the model
Q4_2 <- nimbleModel(code = Q4Code2, name = 'Q4_2', constants = Q4Consts2, 
                    data = Q4Data, inits = Q4Inits2)

# Compiling the model
CQ4_2 <- compileNimble(Q4_2)

# Set up monitored quantities
Q4Conf2 <- configureMCMC(Q4_2, monitors = c('B0', 'B1', 'B2', 'B3'), 
                         print = TRUE)

# Build MCMC algorithm
Q4MCMC_2 <- buildMCMC(Q4Conf2)

# Compile MCMC chain
CQ4MCMC_2 <- compileNimble(Q4MCMC_2, project = Q4_2)



# PRIOR 3: New distribution for Betas -> Uniform(-5, 5)
Q4Code3 <- nimbleCode({
  # Specify the likelihood
  for (i in 1:N) {
    
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0 + B1*z[i, 1] + B2*z[i, 2] + B3*x[i])
    
  }
  
  # Prior specification
  B0 ~ dunif(-5, 5)
  B1 ~ dunif(-5, 5)
  B2 ~ dunif(-5, 5)
  B3 ~ dunif(-5, 5)
})


# Constants of model
Q4Consts3 <- list(N=73)

# Set of initial values before building model
Q4Inits3 <- list(B0=1, B1=1, B2=1, B3=1)

# Building the model
Q4_3 <- nimbleModel(code = Q4Code3, name = 'Q4_3', constants = Q4Consts3, 
                    data = Q4Data, inits = Q4Inits3)

# Compiling the model
CQ4_3 <- compileNimble(Q4_3)

# Set up monitored quantities
Q4Conf3 <- configureMCMC(Q4_3, monitors = c('B0', 'B1', 'B2', 'B3'), 
                         print = TRUE)

# Build MCMC algorithm
Q4MCMC_3 <- buildMCMC(Q4Conf3)

# Compile MCMC chain
CQ4MCMC_3 <- compileNimble(Q4MCMC_3, project = Q4_3)




# RUNNING THE MCMC CHAINS
# Using same initial values
Q4Inits_chains <- list(Q4Inits1, list(B0=5, B1=5, B2=5, B3=5), 
                       list(B0=3, B1=3, B2=3, B3=3))

# Results of prior 1
set.seed(1)
posterior_1 <- runMCMC(CQ4MCMC_1, niter=51000, thin=10, nburnin=1000,
                       summary=TRUE, samples=TRUE, nchains=3,
                       samplesAsCodaMCMC=TRUE, inits=Q4Inits_chains)

combinedchains_MQ4_1 <- mcmc.list(posterior_1$samples$chain1, 
                                  posterior_1$samples$chain2, 
                                  posterior_1$samples$chain3)


# Results of prior 2
set.seed(1)
posterior_2 <- runMCMC(CQ4MCMC_2, niter=51000, thin=10, nburnin=1000,
                       summary=TRUE, samples=TRUE, nchains=3,
                       samplesAsCodaMCMC=TRUE, inits=Q4Inits_chains)

combinedchains_MQ4_2 <- mcmc.list(posterior_2$samples$chain1, 
                                  posterior_2$samples$chain2, 
                                  posterior_2$samples$chain3)

# Results of prior 3
set.seed(1)
posterior_3 <- runMCMC(CQ4MCMC_3, niter=51000, thin=10, nburnin=1000,
                       summary=TRUE, samples=TRUE, nchains=3,
                       samplesAsCodaMCMC=TRUE, inits=Q4Inits_chains)

combinedchains_MQ4_3 <- mcmc.list(posterior_3$samples$chain1, 
                                  posterior_3$samples$chain2, 
                                  posterior_3$samples$chain3)



# Posterior summaries
# Prior 1: Normal(mean=5, sd=5)
summary(combinedchains_MQ4_1)
Q4summary1 <- posterior_1$summary$all.chains # with error
Q4summary1

# Prior 2: Normal(mean=-5, sd=2)
summary(combinedchains_MQ4_2)
Q4summary2 <- posterior_2$summary$all.chains # with error
Q4summary2

# Prior 3: Unif(-1, 1)
summary(combinedchains_MQ4_3)
Q4summary3 <- posterior_3$summary$all.chains # with error
Q4summary3

Q2summary
summary(combinedchains_MQ2)
#-------------------------------------------------------------------------------



####################################################
################## QUESTION 5 CODE #################
####################################################

# Supposedly: E[lambda_i] (x_i=0) = [2, 5] with 95% certainty IF
# > w_i1 = E[w_i1] -> z_i1 = 0 
# > w_i2 = E[w_i2] -> z_i2 = 0

# Model: log(lambda_i) = B0  -> as x_i = z_i1 = z_i2 = 0
#        lambda_i = exp(B0)

# PART 2
# Simulate from the prior predictive distribution of p(y_i*) for
# students with average scores in the quiz and nimble.
N <- 5

# B0 prior: N(1.1513, 0.2338^2)
# lambda_i* = exp(B0)
# Y*|B0 ~ Poisson(lambda_i*)

# From now on, i refers to i*.

# Model 
Q5Code <- nimbleCode({
  # likelihood
  for (i in 1:N) {
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0)
  }
  
  # prior
  B0 ~ dnorm(1.1513, 1/(0.2338)^2)
})

# Constants
Q5Consts <- list(N=5)

# Initial values
Q5Inits <- list(B0 = 1)

# Build model
Q5Model <- nimbleModel(code = Q5Code, 
                       constants = Q5Consts,
                       inits = Q5Inits)

# Simulate it = 1000 sets of N=5 values from the prior predictive distribution
set.seed(5)
# dependencies to simulate
simnodes <- Q5Model$getDependencies('B0', self=FALSE, downstream=TRUE)
simnodes

# iterations
it <- 10000

# storing results
Y_sim <- list()

for (i in 1:it){
  # Simulate from prior predictive distribution
  Q5Model$simulate(simnodes)
  sim_vals <- c(Q5Model$Y)
  
  #append
  Y_sim <- c(Y_sim, list(sim_vals))
}


# Distribution of g(y_i*) 
means <- c()

for (i in 1:it) {
  # extract simulated values of iteration
  it_vals <- Y_sim[[i]]
  
  # mean value
  Y_bar <- mean(it_vals)
  
  # append
  means <- c(means, Y_bar)
}

# density plot of mean values
par(mfrow=c(1,1))
plot(density(means), xlab='Mean value')



# PART 3
Y_data <- data$y[(abs(data$z_Quiz) < 0.25) & (abs(data$z_Nimb) < 0.25)]
Y_data

mean_data <- mean(Y_data)
mean_data

mean_sim <- mean(means)
mean_sim

sort_means <- sort(means)
lowcut <- it*0.025
low95ci <- sort_means[lowcut]
low95ci

# Proposed prior is NOT compatible with the data
# mean Y of z_abs<0.25 data "<" mean Y implied by prior for z_vals = 0
# the data values of z_vals = 0 should be included in z_abs<0.25, yet
# the mean from the data is still smaller -> implies prior is wrong, and
# assumes there are too many students with average scores in both.
#-------------------------------------------------------------------------------



####################################################
################## QUESTION 6 CODE #################
####################################################
# Original data
data_na <- ExamData0
pyt <- data_na$Pyt
z_Q <- data_na$z_Quiz # should be the same

# z_i2 = y0 + y1*pyt_i + e_i   -> for i = 1, 2, ..., 73
# y0 ~ N(0, 10^2)
# y1 ~ N(0, 10^2)
# e_i ~ N(0, sig^2)
# sig ~ U(0, 5)

Q6Code <- nimbleCode({
  # Specify the likelihood
  for (i in 1:N) {
    Y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(B0 + B1*zQ[i] + B2*zN[i] + B3*x[i])
    
    zN[i] <- y0 + y1*pyt[i] + e[i] # 44, 57, 65 are NA
    e[i] ~ dnorm(0, 1/(sig*sig))
    
  }
  
  # sample new values
  z_N_44.new <- y0 + y1*pyt[44] + e[44] 
  z_N_57.new <- y0 + y1*pyt[57] + e[57] 
  z_N_65.new <- y0 + y1*pyt[65] + e[65] 
  
  # Prior specification
  B0 ~ dnorm(0, 1/10^2)
  B1 ~ dnorm(0, 1/10^2)
  B2 ~ dnorm(0, 1/10^2)
  B3 ~ dnorm(0, 1/10^2)
  
  # linear regression parameters
  y0 ~ dnorm(0, 1/10^2)
  y1 ~ dnorm(0, 1/10^2)
  sig ~ dunif(0, 5)
})



# Constants of model
Q6Consts <- list(N=73)

# Data of model
Q6Data <- list(x = data_na$x_module, zQ = data_na$z_Quiz, zN = data_na$z_Nimb,
               pyt = data_na$Pyt, Y = data_na$y)

# Set of initial values before building model
Q6Inits <- list(B0=1, B1=1, B2=1, B3=1, y0=1, y1=1, sig=1)

# Building the model
Q6 <- nimbleModel(code = Q6Code, name = 'Q6', constants = Q6Consts, 
                  data = Q6Data, inits = Q6Inits)

# Compiling the model
CQ6 <- compileNimble(Q6)

# Set up monitored quantities (also tracking e samples)
Q6Conf <- configureMCMC(Q6, monitors=c('B0', 'B1', 'B2', 'B3', 'sig', 'y0', 'y1',
                                       'z_N_44.new', 'z_N_57.new', 'z_N_65.new'), 
                        print=TRUE) 

# Build MCMC algorithm
Q6MCMC <- buildMCMC(Q6Conf)

# Compile MCMC chain
CQ6MCMC <- compileNimble(Q6MCMC, project = Q6)


# TESTING
set.seed(1)
Q6Inits2 <- list(Q6Inits, list(B0=4, B1=3, B2=2, B3=1, y0=2, y1=2, sig=3), 
                 list(B0=-1, B1=-1, B2=-1, B3=-1, y0=-1, y1=-1, sig=2))
posterior1 <- runMCMC(CQ6MCMC, niter=50000, thin=1, nburnin=0,
                      summary=TRUE, samples=TRUE, nchains=3,
                      samplesAsCodaMCMC=TRUE, inits=Q6Inits2)

combinedchains_MQ6_test <- mcmc.list(posterior1$samples$chain1, 
                                     posterior1$samples$chain2, 
                                     posterior1$samples$chain3)

plot(combinedchains_MQ6_test)
effectiveSize(combinedchains_MQ6)

# ACF plots
autocorr.plot(posterior1$samples$chain1)
autocorr.plot(posterior1$samples$chain2)
autocorr.plot(posterior1$samples$chain3)



# ADDING BURNIN AND THINNING
set.seed(1)
# Initial visualisation of plots
posterior2 <- runMCMC(CQ6MCMC, niter=2015000, thin=2000, nburnin=15000,
                     summary=TRUE, samples=TRUE, nchains=3,
                     samplesAsCodaMCMC=TRUE, inits=Q6Inits2)

combinedchains_MQ6 <- mcmc.list(posterior2$samples$chain1, 
                                posterior2$samples$chain2, 
                                posterior2$samples$chain3)



plot(combinedchains_MQ6)
effectiveSize(combinedchains_MQ6)


# ACF plots
autocorr.plot(posterior2$samples$chain1)
autocorr.plot(posterior2$samples$chain2)
autocorr.plot(posterior2$samples$chain3)

# Gelman plots
gelman.diag(combinedchains_MQ6)
gelman.plot(combinedchains_MQ6)

# RESULTS
# posterior summary
summary(combinedchains_MQ6)
Q6summary <- posterior2$summary$all.chains
Q6summary






# // PART 2: Posterior predictive distribution of model vs Q2 for i=57 //
# f(Y_57) ~ Poisson(lambda_57)

# Q6 MODEL: 
# Get posterior data
posterior_chain1 <- data.frame(combinedchains_MQ6[[1]])
posterior_chain2 <- data.frame(combinedchains_MQ6[[2]])
posterior_chain3 <- data.frame(combinedchains_MQ6[[3]])

posterior_chains <- rbind(posterior_chain1, posterior_chain2)
posterior_chains <- rbind(posterior_chains, posterior_chain3)


# get posterior lambda_57 and sample a Y_57 value from it 
n <- length(posterior_chains$B0)
Y_57s_6 <- c()

for (i in 1:n) {
  # get z_N_57 value
  z_N_57_val <- posterior_chains$z_N_57.new[i]
  
  # get lambda
  lambda_57 <- exp(posterior_chains$B0[i] + posterior_chains$B1[i]*z_Q[57] +
                   posterior_chains$B2[i]*z_N_57_val)
  
  # sample from poisson
  Y_57_val <- rpois(1, lambda_57)
  
  # append
  Y_57s_6 <- c(Y_57s_6, Y_57_val)
}



# Q2 MODEL:
# Get posterior data
posterior_chain1 <- data.frame(combinedchains_MQ2[[1]])
posterior_chain2 <- data.frame(combinedchains_MQ2[[2]])
posterior_chain3 <- data.frame(combinedchains_MQ2[[3]])

posterior_chains <- rbind(posterior_chain1, posterior_chain2)
posterior_chains <- rbind(posterior_chains, posterior_chain3)


# get posterior lambda_57 and sample a Y_57 value from it 
n <- length(posterior_chains$B0)
Y_57s_2 <- c()

for (i in 1:n) {
  # get lambda
  lambda_57 <- exp(posterior_chains$B0[i] + posterior_chains$B1[i]*z_Q[57] +
                     posterior_chains$B2[i]*z_N[57]) #z_N from data in Q2
  
  # sample from poisson
  Y_57_val <- rpois(1, lambda_57)
  
  # append
  Y_57s_2 <- c(Y_57s_2, Y_57_val)
}


# Plotting density 
par(mfrow=c(1,2))
hist(Y_57s_2, col=c1, xlim=c(-0.3,18), probability=TRUE, main='Q2', xlab='Y value')
# Plotting density 
hist(Y_57s_6, col=c1, xlim=c(-0.3,18), probability=TRUE, main='Q6', xlab='Y value')



