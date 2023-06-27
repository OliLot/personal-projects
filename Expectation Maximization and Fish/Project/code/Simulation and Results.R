# ---------------------- Simulation to test EM Algorithm -----------------------

# // Simulate data from Gaussian mixture and demonstrate that EM algorithm
# // can recover the correct parameter values

# nice looking parameters for mix
lambda1 = 0.2
mean1 = 15
sd1 = 6

lambda2 = 0.5
mean2 = 40
sd2 = 9

lambda3 = 0.3
mean3 = 70
sd3 = 10



# mix probability function
dmix_sim_original <- function(x) {
  
  return(lambda1 * dnorm(x, mean1, sd1) + lambda2 * dnorm(x, mean2, sd2) + 
           lambda3 * dnorm(x, mean3, sd3))
} 

# plot to visualise data_sim spread of fish lengths
curve(dmix_sim_original, xlim = c(0, 100), ylim = c(0, 0.03))



# choose sample size
sampsize <- 1000

# fish from age group with (with associated probability)
Age <- sample(c(1, 2, 3), sampsize, replace = TRUE, prob = c(lambda1, lambda2, 
                                                             lambda3))

# creating simulated data frame
data_sim <- data.frame("FishID" = seq(1, sampsize), "Length" = NaN, "Age" = Age)

# for each group, sample from the simulated normal distribution
for (i in 1:sampsize) {
  if (data_sim$Age[i] == 1) {
    data_sim$Length[i] <- rnorm(1, mean1, sd1)
  } else if (data_sim$Age[i] == 2) {
    data_sim$Length[i] <- rnorm(1, mean2, sd2)
  } else if (data_sim$Age[i] == 3) {
    data_sim$Length[i] <- rnorm(1, mean3, sd3)
  }
}


# Plot showing typical fit to simulated data
# extract parameter values
parameters <- teamEM(data_sim, 1e-10)$estimates
parameters

# probability density function with estimate parameters
dmix_sim_EM <- function(x) {
  
  return(parameters[1, 1] * dnorm(x, parameters[1, 2], parameters[1, 3]) + 
           parameters[2, 1] * dnorm(x, parameters[2, 2], parameters[2, 3]) + 
           parameters[3, 1] * dnorm(x, parameters[3, 2], parameters[3, 3]))
}

# histogram for fish lengths overlayed with mix probability density function
hist(data_sim$Length, prob = TRUE, main = "Histogram and Probability Density", 
     xlab = "Fish Length", ylab = "Probability Density", xlim = c(0, 100), 
     ylim = c(0, 0.025), breaks = 25)
# estimate from EM algorithm as blue
curve(dmix_sim_EM, lwd = 2, add = TRUE, col = "blue")
# initial parameters as red
curve(dmix_sim_original, lwd = 2, add = TRUE, col = "red")




# ------------------------ Results of the assignment ---------------------------

# Results for the assignment
# Table of estimates:
results <- teamEM(data_real, epsilon = 1e-12)$estimates
results

# Plot of model fits data:
# probability density function with estimate parameters
dmix_real <- function(x) {
  
  return(results[1, 1] * dnorm(x, results[1, 2], results[1, 3]) + 
           results[2, 1] * dnorm(x, results[2, 2], results[2, 3]) + 
           results[3, 1] * dnorm(x, results[3, 2], results[3, 3]))
}

# histogram for fish lengths overlayed with mix probability density function
hist(data_real$Length, prob = TRUE, main = "Histogram and Probability Density", 
     xlab = "Fish Length", ylab = "Probability Density", xlim = c(0, 100), 
     ylim = c(0, 0.04))
curve(dmix_real, lwd = 2, add = TRUE)

