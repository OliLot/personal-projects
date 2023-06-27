install.packages("dslabs") # to get dataset
install.packages("sn") # to calculate skewed normal distributions
install.packages('plyr') # to count frequencies of occurrences in data

library("sn")
library("plyr")
library("dslabs")
help(package = "dslabs")

set.seed(190005300)
# ---------------------- Research Question/Preparing Data ----------------------
# RQ: Are self-reported male heights normally distributed?
# i.e. test for normality

# H_0: the self-reported heights are normally distributed
# H_1: the self-reported heights are not normally distributed


data <- heights
head(data)

# // Clean up the data //
# delete any data rows with NA values
data <- na.omit(data)
data


# ----------------------------- Test Statistics --------------------------------
# Parametric: Jarque-Bera (JB) test
# JB statistic tests whether data has a normal distribution (using its skewness
# S and kurtosis K)
# + asymptotically approaches chi squared when sample data assumed to be normal, 
# which will be taken to be the case here (parametric).

# Non-parametric: Lilliefors (LF) test
# Lilliefors statistic tests whether the data (which can have any distribution)
# fits the specified distribution based on estimate parameters from the data.
# In this case, the data distribution will be compared to the normal
# distribution with estimated sample mean and sd of the data. 
# It works by comparing the supremum of the absolute difference between the
# empriical cumulative distribution of the data, and that of the specified 
# distribution that is being tested for.



# --------------------- Initial look at data/analysis --------------------------
# // Plots of data //
hist(data$height[data$sex == "Male"], breaks = 20, probability = TRUE)


# Looks like a normal distribution would be a reasonable approximation here 
# (with the exception of a few outliers)
mean_heights = mean()
mean_heights
sd_heights = sd(data$height[data$sex == "Male"])
sd_heights 

# A standard case of a normal distribution with mean 70 and standard deviation
# of 4 seems apt (heights are measured in inches).
x <- seq(0, 100)
lines(x, dnorm(x, mean_heights, sd_heights))


# ------------------------ Statistical test functions --------------------------
# Functions to calculate JB and LF test statistics

# // Jarque-Bera test statistic function //
# Function calculating JB statistic for set of data
# INPUT:
#   - vector of simulated test score data
# OUTPUT:
#   - JB statistic value
JB_stat <- function(data) {
  
  # prerequisites
  n <- length(data)
  data_mean <- mean(data)
  
  S_sum_num <- 0
  K_sum_num <- 0
  sum_denom <- 0
  
  # calculating from formula
  for (i in 1:n) {
    
    S_sum_num <- S_sum_num + (data[i] - data_mean)^3 
    K_sum_num <- K_sum_num + (data[i] - data_mean)^4
    
    sum_denom <- sum_denom + (data[i] - data_mean)^2
  }
  
  S <- ((1/n) * S_sum_num) / ((1/n) * sum_denom)^(3/2)
  K <- ((1/n) * K_sum_num) / ((1/n) * sum_denom)^2
  
  JB <- (n/6) * (S^2 + (K - 3)^2 / 4)
  return(JB)
}

# JB critical value
# chi squared distribution (2 df) for finding critical value of JB statistic
# The null hypothesis is rejected if the JB statistic is greater that or equal
# to this critical value
alpha = 0.05

JB_critval <- qchisq(1 - alpha, 2)


# // Lilliefors test statistic functions //
# Function calculating empirical distribution values for sorted data
emp_cdf <- function(data) {
  
  # prerequisites (getting data as a frequency table, and its length)
  freq_data <- count(data)
  n <- length(freq_data$freq)
  
  # calculating the empirical distribution function
  v <- c()
  freq <- 0
  for (i in 1:n) {
    freq <- freq + freq_data$freq[i]
    v <- c(v, rep(freq/length(data), freq_data$freq[i]))
  }
  
  return(v)
}

# Function calculating the theoretical cumulative distribution values for the
# sorted data
distr_cdf <- function(data) {
  
  # prerequisites
  data_mean <- mean(data)
  data_sd <- sd(data)
  n <- length(data)
  
  # calculating theoretical cumulative distribution values of the data
  v <- c()
  for (i in 1:n) {
    
    v <- c(v, pnorm(data[i], mean = data_mean, sd = data_sd))
    
  }
  return(v)
}


# Function computing LF test statistic value given input data
# INPUTS:
#   - data = vector of simulated data for test scores
# OUTPUTS:
#   - LF test statistic
#   - Plot of empirical cumulative distribution and theoretical cumulative 
#     distribution for comparison
LF_stat <- function(data) {
  
  # sorting the data from smallest to largest
  data <- sort(data)
  
  # Empirical distribution
  data_cdf <- emp_cdf(data)
  
  # Cumulative distribution
  theo_cdf <- distr_cdf(data)
  
  # Plot the CDFs against each other (if run once: to visualise comparison)
  # x <- data
  # plot(x, data_cdf, type="l", col="blue", lwd=2)
  # lines(x, theo_cdf, col="red", lwd=2)
  
  D_n <- max(abs(data_cdf - theo_cdf))
  return(D_n)
}



# LF critical value (5% level will be used)
# The null hypothesis is rejected if the LF statistic is greater than these
# critical values for given % levels - which depends on the number of trials

# alpha = 0.01 (1% level) - 1% chance that the LF value is this or greater
LF_crit_1 <- function(n) {
  return(1.035/((0.83 + n)/sqrt(n) - 0.01))
}

# alpha = 0.05 (5% level) - 5% chance that the LF value is this or greater
LF_crit_5 <- function(n) {
  return(0.895/((0.83 + n)/sqrt(n) - 0.01))
}

# alpha = 0.1 (10% level) - 10% chance that the LF value is this or greater
LF_crit_10 <- function(n) {
  return(0.819/((0.83 + n)/sqrt(n) - 0.01))
}


# --------------------------- Simulations/Scenarios ----------------------------
# Will now examine the effect of changing sample size, spread of the tails, and
# skewness on the size and power of the JB and LF tests

# H_O: The distribution of the data is the normal distribution
# H_1: The distribution of the data is not the normal distribution


# // Scenario 1: Changing sample size //
# Will change the sample size looking at small, medium and large samples.
# To calculate size: Will simulate data from the a normal distribution with mean
#                    height 70 and a standard deviation of 4. Any mean and sd 
#                    that could reasonably represent the data for male heights
#                    can be used here.
# To calculate power: The uniform distribution (far from normal) in addition to  
#                     the cauchy distribution (changing scale and location 
#                     parameters - similar to the normal)



# Function for simulating n heights from a particular distribution
# INPUTS:
#   - n = number of samples to take
#   - distr = "function type" (e.g. "norm", "cauchy", "unif", "sn")
#   - ... = extra arguments to be passed into the distribution function
# OUTPUTS;
#   - Vector of n values samples from the input distribution
n_data <- function(n, distr, ...) {
  
  f <- get(paste("r", distr, sep = ""), mode = "function")
  
  return(f(n, ...))
}

# ---------------------- SCENARIO 1: CHANGING SAMPLE SIZE ----------------------
# Calculating size/power for the JB/LF statistic varying the sample size n
# (calculated at the 5% confidence level i.e. alpha = 0.05)
# INPUTS:
#   - n = number of samples
#   - test = input as either "JB" or "LF" to specify statistic to test
#   - iterations = number of iterations in accepting/rejecting the hypothesis
#                  (default as 1000)
#   - distr = name of distribution type (e.g. "norm", "cauchy", "unif"). If this
#             distribution is "norm" -> testing for size, else for power.
#   - ... = extra arguments required to specify the distr input (e.g. mean, sd 
#           for distr = "norm").
# OUTPUTS:
#   - size/power of the JB/LF statistic for a sample size of n at the 5% level
sp_samples <- function(n, test, iterations = 1000, distr, ...) {
  
  if (test == "JB") {
    
    # critical JB value at alpha = 0.05 level
    crit <- qchisq(0.95, 2)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      # get sample data
      data <- n_data(n, distr, ...)
      
      # get statistic value
      stat <- JB_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
    
  } else if (test == "LF") {
    
    # critical LF value at alpha = 0.05 level
    crit <- LF_crit_5(n)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      # get sample data - must also be sorted
      data <- sort(n_data(n, distr, ...))
      
      # get statistic value
      stat <- LF_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
  }
  
  return(sum(v)/length(v))
}


# Plot and table of the size of JB/LF statistic for varying sample sizes 
# n = 25, 50, 100, 500, 1000, 2000, 5000, 10000, 20000.
# INPUT:
#   - test = "JB" or "LF" to specify test being conducted
#   - iterations = number of times to conduct the test (accept or reject normal)
#                  (default as 1000)
#   - distr = distribution class as string, e.g. "norm", "cauchy", "unif".
#             If the choice of distribution is "norm", it is a test of size. If
#             it is any other distribution it is a test of power.
#   - ... = any further arguments/parameters to be passed to the distr function
# OUTPUT:
#   - prints a plot and table of the sample sizes with corresponding JB/LF sizes
s_calc_samples <- function(test, iterations = 1000, distr, ...) {
  
  # sample sizes to test
  small <- c(25, 50, 100)
  medium <- c(500, 1000, 2000)
  large <- c(5000, 10000, 20000)
  
  # creating a vector of sample sizes
  x <- c(small, medium, large)
  
  # creating a vector of corresponding test size/power values
  y <- c() 
  for (i in 1:length(x)) {
    
    val <- sp_samples(x[i], test, iterations, distr, ...)
    y <- c(y, val)
  
  }
  
  # plot the results
  plot(x,y)
  
  # table of values
  size <- data.frame(x, y, row.names = NULL)
  colnames(size) <- c("sample", "proportion")
  size
}

# Same as above, but calculates power: The range of sample sizes is chosen to
# better analyse the outcome.
p_calc_samples <- function(test, iterations = 1000, distr, ...) {
  
  # sample sizes to test
  small <- c(5, 10, 15, 20, 25, 30)
  medium <- c(100, 200, 500, 1000, 2000)
  large <- c(4000, 10000)
  
  # creating a vector of sample sizes
  x <- c(small, medium, large)
  
  # creating a vector of corresponding test size/power values
  y <- c() 
  for (i in 1:length(x)) {
    
    val <- sp_samples(x[i], test, iterations, distr, ...)
    y <- c(y, val)
  }
  
  # plot the results
  plot(x,y)
  
  # table of the data
  size <- data.frame(x, y, row.names = NULL)
  colnames(size) <- c("sample", "proportion")
  size
}



# ////////////// Calculating size for both statistical tests //////////////
# Will use normal distribution mentioned above (mean = 70, sd = 4), though any
# reasonable normal distribution to represent male heights could be used here


# plot of the distribution that is being used to sample data from
x <- seq(40, 100)
y <- dnorm(x, mean = 70, sd = 4)
plot(x, y)


# - SIZE for JB statistic -
# with data sampled from normal distribution
s_calc_samples(test = "JB", iterations = 1000, 
               distr = "norm", mean = 70, sd = 4)

# ambiguous below 2000 sample size as expected, since it is only asymptotically
# approximated by a chisq (2df). Large samples above 2000: the size is steady 
# at about 5-6%.


# - SIZE for LF statistic -
# with data sampled from normal distribution
s_calc_samples(test = "LF", iterations = 100, 
               distr = "norm", mean = 70, sd = 4)


# small samples size fluctuates quite heavily - which should
# be good, however, making accurate predictions for a distribution based on 
# merely a few samples is clearly not accurate. The size of test is varies
# between 4-6% for larger samples.

# It would of course be better to use more iterations, but it is too
# intense computationally, and requires a too-long runtime. Table values for 
# iterations with 1000, but leave code with 100. Even better as well to simulate
# from extremely large sample sizes (like 100000), but again too computationally
# intense. 


# For even larger sample size of  50'000
# JB size
sp_samples(50000, test = "JB", iterations = 1000, 
           distr = "norm", mean = 70, sd = 4)
# result: 0.052, 0.038, 

# LF size
sp_samples(30000, test = "LF", iterations = 500, 
           distr = "norm", mean = 70, sd = 4)
# result: 0.04


# both results support the aforementioned trend





# ////////////// Calculating power for both statistical tests //////////////
# 1. Simulating data from a cauchy distribution (similar to normal distribution)
# with mean 70 and sd 4. The distribution that most looks like the normal is the
# Cauchy distribution with location = 70 and scale = 3.1.

# plot to visualise similarity of the two distributions
x <- seq(40, 100)
y_norm <- dnorm(x, mean = 70, sd = 4)
plot(x, y_norm)
y_cauchy <- dcauchy(x, location = 70, scale = 3.1)
lines(x, y_cauchy)



# - POWER of JB statistic -
p_calc_samples(test = "JB", iterations = 1000, 
               distr = "cauchy", location = 70, scale = 3.1)

# asymptotically has a power of 1 as the sample size is increased, for very 
# small samples - where it is harder to predict what the distribution is - the 
# power is quite small. It does, however, pick up very quickly with a 70% power
# already achieved around 15 samples, and about 90% at 25.


# - POWER of LF statistic - 
p_calc_samples(test = "LF", iterations = 100, 
               distr = "cauchy", location = 70, scale = 3.1)

# similar power to the JB statistic, asymptotically it is 1, and lower for small
# sample sizes. The power for the test seems to be better for very
# small samples (sample of 5, 10 are better), but then is slightly worse
# for the greater samples, also seems to approach 1 later (0.8 for 25)




# 2. Simulating data from a uniform distribution (far from the normal)
# will take the uniform distribution to stretch from 57.5 - 82.5
x <- seq(40, 100)
y_norm <- dnorm(x, mean = 70, sd = 4)
plot(x, y_norm)
y_unif <- dunif(x, 57.5, 82.5)
lines(x, y_unif)


# - POWER of JB statistic -
p_calc_samples(test = "JB", iterations = 1000, distr = "unif", 57.5, 82.5)


# - POWER of LF statistic -
p_calc_samples(test = "LF", iterations = 100, distr = "unif", 57.5, 82.5)





# ------------------ Scenario 2: Changing tail length/size ---------------------
# Will change the tail length/size from short/small to long/large
# To calculate size: Will simulate data from normal distribution with varying sd,
#                    mean chosen as the same of 70 (as standard reference).
# To calculate power: Cauchy distribution with varying scale parameter
#                     (location chosen as 70). Also uniform distribution across
#                     differing intervals.


# Calculating size/power for the JB/LF statistic varying the tail length/size
# (calculated at the 5% confidence level i.e. alpha = 0.05)
# INPUTS:
#   - distr = name of distribution type (either "norm", "cauchy", "unif"). If 
#             distribution is "norm" -> testing for size, else testing for power.
#   - spread = input value of parameter that causes spreading (If distr = "norm"
#              the spread value must be a single integer representing sd. If
#              distr = "cauchy" the spread value must be a single integer for
#              scale. If distr = "unif" the spread must be a vector of the 
#              start and end value.
#   - test = input as either "JB" or "LF" to specify statistic to test
#   - iterations = number of iterations in accepting/rejecting the hypothesis
#                  (default as 1000)
#   - n = number of samples to simulate (default as 3000)
# OUTPUTS:
#   - size/power of the JB/LF statistic for given spread parameter at the 5% 
#     level
sp_tails <- function(distr, spread, test, iterations = 1000, n = 3000) {
  
  if (test == "JB") {
    
    # critical JB value at alpha = 0.05 level
    crit <- qchisq(0.95, 2)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      
      # get sample data
      if (distr == "norm") {
        
        data <- n_data(n, distr, mean = 70, sd = spread)
        
      } else if (distr == "cauchy") {
        
        data <- n_data(n, distr, location = 70, scale = spread)
        
      } else if (distr == "unif") {
        
        data <- n_data(n, distr, spread[1], spread[2])
        
      }
      
      # get statistic value
      stat <- JB_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
    
  } else if (test == "LF") {
    
    # critical LF value at alpha = 0.05 level
    crit <- LF_crit_5(n)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      
      # get sample data
      if (distr == "norm") {
        
        data <- sort(n_data(n, distr, mean = 70, sd = spread))
        
      } else if (distr == "cauchy") {
        
        data <- sort(n_data(n, distr, location = 70, scale = spread))
        
      } else if (distr == "unif") {
        
        data <- sort(n_data(n, distr, spread[1], spread[2]))
        
      }
      
      # get statistic value
      stat <- LF_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
  }
  
  return(sum(v)/length(v))
}




# Plot and table of the size/power of JB/LF statistic for varying tail spreads
# If simulation distribution is "norm":
#   - Have values ranging from sd = 1, 2, 3, 4, 5, 10, 15, 20, 50
#     (clearly the large sd's are unreasonable, but they still give insight)
# If simulation distribution is "cauchy":
#   - Have values ranging from scale = 1, 2, 3, 4, 5, 10, 15, 20, 50
# If simulation distribution is "unif":
#   - Have intervals ranging from between (60, 80), (50, 90), (40, 100), 
#     (30, 110), (20, 120), (10, 130)
# INPUT:
#   - test = statistic to test (either "JB" or "LF")
#   - iterations = number of times to calculate the test statistic, from which
#                  the proportions are calculated to get size/power 
#                  (default is 1000)
#   - distr = distribution class as string (either "norm", "cauchy", "unif").
#             If the choice of distribution is "norm", it is a test of size. If
#             it is any other distribution it is a test of power.
#   - n = number of samples to base calculation on (default as 3000)
# OUTPUT:
#   - prints a plot and table of the spread parameter with corresponding test
#     size/power
sp_calc_spread <- function(test, iterations = 1000, distr, n = 3000) {
  
  if (distr == "norm" | distr == "cauchy") {
    
    # spread values to test
    spread <- c(1, 2, 3, 4, 5, 10, 15, 20, 50)
    
    # creating a vector of corresponding test size/power values
    y <- c()
    for (i in 1:length(spread)) {
      
      val <- sp_tails(distr, spread[i], test, iterations, n)
      y <- c(y, val)
    }
    
    # plot the values
    plot(spread,y)
    
    # table for size/power values
    size <- data.frame(spread, y, row.names = NULL)
    colnames(size) <- c("spread", "proportion")
    size
    
  } else if (distr == "unif") {
    
    # spread values to test
    spread <- array(c(c(60, 80), c(50, 90), c(40, 100), c(30, 110), c(20, 120)), 
                    dim = c(2, 5))
    spread
    
    # creating a vecto of corresponding test size/power values
    y <- c()
    for (i in 1:length(spread[1, ])) {
      
      val <- sp_tails(distr, spread[, i], test, iterations, n)
      y <- c(y, val)
      
    }
    
    # table of the values
    size <- data.frame(y, row.names = c("(60, 80)", "(50, 90)", "(40, 100)", 
                                        "(30, 110)", "(20, 120)"))
    colnames(size) <- c("proportion")
    size
  }
}



# ///// Calculating sizes of JB/LF statistic with varying tail length/size /////
# - SIZE of JB - 
sp_calc_spread(test = "JB", iterations = 1000, distr = "norm", n = 3000)

# for very small spread
sp_tails(distr = "norm", spread = 0.0001, test = "JB", 
         iterations = 1000, n = 3000)
# similar result

# for very large spread
sp_tails(distr = "norm", spread = 10000, test = "JB", 
         iterations = 1000, n = 3000)
# similar result


# - SIZE of LF -
sp_calc_spread(test = "LF", iterations = 500, distr = "norm", n = 3000)

# for very small spread
sp_tails(distr = "norm", spread = 0.0001, test = "LF", 
         iterations = 1000, n = 3000)
# similar result

# for very large spread
sp_tails(distr = "norm", spread = 10000, test = "LF", 
         iterations = 500, n = 3000)
# similar result


# CONCLUSION: size seems to be invariant when changing the spread of the normal 
#             distribution that the data is simulated from. Size values are 
#             constant between around 4-6%.




# //// Calculating powers of JB/LF statistic with varying tail length/size ////
# 1. For a distribution similar to the normal, e.g. the cauchy distribution,
#    also chosen to have location 70 and varying scale, choose n = 15 to be able
#    to discern effect of varying tail spread.
# - POWER of JB -
sp_calc_spread(test = "JB", iterations = 5000, distr = "cauchy", n = 15)

# for very small spread
sp_tails(distr = "cauchy", spread = 0.0001, test = "JB", 
         iterations = 5000, n = 200)
# similar result

# for very large spread
sp_tails(distr = "cauchy", spread = 10000, test = "JB", 
         iterations = 5000, n = 15)
# similar result



# - POWER of LF -
sp_calc_spread(test = "LF", iterations = 5000, distr = "cauchy", n = 15)

# for very small spread
sp_tails(distr = "cauchy", spread = 0.0001, test = "LF", 
         iterations = 5000, n = 15)
# similar result

# for very large spread
sp_tails(distr = "cauchy", spread = 10000, test = "LF", 
         iterations = 5000, n = 15)
# similar result


# CONCLUSION: The similar cauchy distribution has no variation in its power for
#             different spreads




# 2. For a distribution far from the gaussian e.g. the uniform distribution. 
# Spread varied by varying its intervals as previously specified.
# - POWER of JB - 
sp_calc_spread(test = "JB", iterations = 5000, distr = "unif", n = 100)

# for very small spread
sp_tails(distr = "unif", spread = c(69, 71), test = "JB", 
         iterations = 5000, n = 100)
# similar result

# for very large spread
sp_tails(distr = "unif", spread = c(0, 200), test = "JB", 
         iterations = 5000, n = 100)
# similar result



# - POWER of LF -
sp_calc_spread(test = "LF", iterations = 5000, distr = "unif", n = 100)

# for very small spread
sp_tails(distr = "unif", spread = c(69, 71), test = "LF", 
         iterations = 5000, n = 100)
# similar result

# for very large spread
sp_tails(distr = "unif", spread = c(0, 200), test = "LF", 
         iterations = 5000, n = 100)
# similar result


# CONCLUSION: The far-off uniform distribution has no variation in its power for
#             different spreads





# -------------------------- Scenario 3: Skewness ------------------------------
# Since a normal distribution cannot be skewed, this is merely an exploration of
# the power of the tests
# Power varying skewness: Use skewed normal distributions starting from the
#                         normal with mean = 70 and sd = 4, increasing skewness
#                         progressively



# Calculating power for the JB/LF statistic varying the skewness of the normal
# distribution (calculated at the 5% confidence level i.e. alpha = 0.05)
# INPUTS:
#   - skew = skewness parameter, the greater it is, the greater the skewness
#            (can be any real number)
#   - test = input as either "JB" or "LF" to specify test
#   - iterations = number of iterations in accepting/rejecting the hypothesis
#                  (default as 1000)
#   - n = number of samples to simulate (default as 3000)
# OUTPUTS:
#   - power of the JB/LF statistic for specified skewness at the 5% level
p_skew <- function(skew, test, iterations = 1000, n = 3000) {
  
  if (test == "JB") {
    
    # critical JB value at alpha = 0.05 level
    crit <- qchisq(0.95, 2)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      
      # get sample data
      data <- n_data(n, distr = "sn", xi = 70, omega = 4, alpha = skew)
      
      # get statistic value
      stat <- JB_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
    
  } else if (test == "LF") {
    
    # critical LF value at alpha = 0.05 level
    crit <- LF_crit_5(n)
    
    # repeat testing appending 0 or 1 if normality is accepted or rejected
    v <- c()
    for (i in 1:iterations) {
      
      # get sample data
      data <- sort(n_data(n, distr = "sn", xi = 70, omega = 4, alpha = skew))
      
      # get statistic value
      stat <- LF_stat(data)
      
      # accept normality (0) or reject it (1)
      if (stat <= crit) {
        
        # accept that it is normal
        v <- c(v, 0)
        
      } else if ( stat > crit) {
        
        # reject that it is normal
        v <- c(v, 1)
        
      }
    }
  }
  
  return(sum(v)/length(v))
}




# Plot and table of the size of JB/LF statistic for varying skew parameter alpha
# (called skew as the input, skew = 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25,
#  2, 3, 5, 10, 20)
#   - test = specifying whether to conduct the "JB" or "LF" test
#   - iterations = number of times to conduct the test (reject or accept) 
#                  from which the power is computed (default as 1000)
#   - n = sample size to use for the test (default as 3000)
# OUTPUT:
#   - prints a plot and table of the skew parameter with corresponding power
p_calc_skew <- function(test, iterations = 1000, n = 3000) {
  
  # skew = alpha values to test
  x <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 2, 3, 5, 10, 20)
  
  # negative values simply mirror these in the other direction, so will have the
  # same results
  
  # creating a vector of corresponding test size/power values
  y <- c() 
  for (i in 1:length(x)) {
    
    val <- p_skew(x[i], test, iterations, n)
    y <- c(y, val)
  }
  
  # plot the results
  plot(x,y)
  
  # table of the data
  size <- data.frame(x, y, row.names = NULL)
  colnames(size) <- c("skew", "proportion")
  size
}

# - Calculating the power for varying skewness
# Low skew number means closer to regular normal, large means greatly skewed
# POWER of the JB test
p_calc_skew(test = "JB", iterations = 1000, n = 3000)

# POWER of the LF test
p_calc_skew(test = "LF", iterations = 500, n = 3000)




# -------------------- ADDITIONAL CALCULATIONS THAT WERE MADE ------------------

# size measurement for JB test for 100000 sample size
sp_samples(100000, test = "JB", iterations = 1000, 
           distr = "norm", mean = 70, sd =4)
# results: within rounding of 4-6 range 


# size measurement for LF test for 100000 sample size
sp_samples(100000, test = "JB", iterations = 100, 
           distr = "norm", mean = 70, sd =4)
# results: within rounding of 4-6 range 

x <- seq(40, 100)
y <- dnorm(x, mean = 70, sd = 4)
plot(x, y, ylim = c(0, 0.15))
lines(x, dsn(x, xi = 70, omega = 4, alpha = 1), col = "blue")
lines(x, dsn(x, xi = 70, omega = 4, alpha = 2), col = "red")
lines(x, dsn(x, xi = 70, omega = 4, alpha = 0.5), col = "orange")
