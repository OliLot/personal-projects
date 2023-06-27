# I confirm that the attached is my own work, except where clearly indicated in the text.

# ------------------------ FUNCTION 1: my_rnorm() ----------------------------
# my_rnorm() function returning a vector of pseudo-random values from a normal
# distribution using the Box-Muller algorithm.
# Inputs:
#   - n, number of values to return (the length of the vector)
#   - mean, the mean of the normal distribution
#   - sd, the standard deviation of the normal distribution
# Outputs:
#   - a vector of pseudo-random values from normal distribution N(mean, sd^2),
#     of length n
my_rnorm <- function(n, mean = 0, sd = 1) {

  # Testing arguments n, mean, and sd are valid: stops execution if not
  if (
    # testing n
    length(n) != 1 || n <= 0 || is.numeric(n) == FALSE ||
    abs(n - round(n)) > 1e-10 || is.finite(n) == FALSE ||
    
    # testing mean
    length(mean) != 1 || is.numeric(mean) == FALSE || 
    is.finite(mean) == FALSE ||
    
    
    # testing sd
    length(sd) != 1 || is.numeric(sd) == FALSE || is.finite(sd) == FALSE || 
    sd < 0 
    ) {
    stop("invalid arguments")
  } else {
    # Defining function that runs Box-Muller algorithm m times, with
    # all pairs of results as elements of a vector
    # Input:
    #   - m, number of times to run algorithm
    # Output:
    #   - vector of all pairs of results
    bm_algorithm <- function(m) {
      
      # define initially empty vector
      v_bm <- c()
    
      # run bm algorithm m times
      for (i in 1:m) {
        # extract two values from uniform distribution with min = 0, max = 1
        a <- runif(1, 0, 1)
        b <- runif(1, 0, 1)
      
        # pair of pseudo-random values from normal distribution (with: mean, sd)
        x1 <- sd*(sin(2*pi*a)*sqrt(-2*log(b, base = exp(1)))) + mean
        x2 <- sd*(cos(2*pi*a)*sqrt(-2*log(b, base = exp(1)))) + mean
      
        # append x1, x2 pair of values to v_bm
        v_bm <- c(v_bm, x1, x2)
      }
      # return vector of sample pairs
      return(v_bm)
    }
  
    # For even n: the bm algorithm can be run n/2 times and a vector with n
    # samples from the specified normal distribution will be obtained.
    if (n %% 2 == 0) {
      v_norm <- bm_algorithm(n/2)
    
      # return the vector
      return(v_norm)
    } else {
      # For odd n: the bm algorithm can be run (n+1)/2 times, and a vector with
      # n+1 samples from the specified normal distribution will be obtained
      v_norm <- bm_algorithm((n+1)/2)
    
      # last entry of the vector is deleted s.t. it only contains n samples
      v_norm <- head(v_norm, -1)
    
      # return the vector
      return(v_norm)
    }
  }
}



#-------------------------- FUNCTION 2: my_rchisq() ---------------------------

# my_rchisq() function returning a vector of pseudo-random chi squared
# distributed deviates
# Inputs:
#   - n, number of values to return (the length of the vector output)
#   - df, degrees of freedom for chi squared distribution
# Outputs:
#   - vector of pseudo-random values from chi squared distribution with df
#     degrees of freedom, of length n
my_rchisq <- function(n, df = 1) {
  
  # Testing arguments n, df are valid: stops execution if not
  if (
    # testing n
    length(n) != 1 || n <= 0 || is.numeric(n) == FALSE ||
    abs(n - round(n)) > 1e-10 || is.finite(n) == FALSE ||
    
    # testing df
    length(df) != 1 || is.numeric(df) == FALSE || is.finite(df) == FALSE ||
    df < 0  || abs(df - round(df)) > 1e-10
  ) {
    stop("invalid arguments")
  } else {
    # Function composing one sample of the chi squared distribution with df 
    # degrees of freedom (out of standard normal values)
    # Input:
    #   - df, degrees of freedom
    # Output:
    #   - one pseudo-random deviate from the chi squared distribution with df
    #     degrees of freedom
    chisq_sample <- function(df) {
      
      # define initially "empty"/zero value
      chi_deviate <- 0
      
      # avoid error message for n = df = 0 in my_rnorm(), as df = 0 is possible
      # for the chi squared distribution, but n = 0 is not possible for my_rnorm
      if (df != 0) {
        # sample df standard normal values
        snorm_values <- my_rnorm(df)
        
        # compute chi squared deviate
        for (i in 1:df) {
          elmt_sq <- (snorm_values[i])**2
          
          chi_deviate <- chi_deviate + elmt_sq
        }
        # return chi squared distribution sample (df != 0)
        return(chi_deviate)
      } else {
        # return chi squared distribution sample (df = 0, all deviates = 0)
        return(chi_deviate)
      }
    }
    
    # Define initially empty vector
    v_chisq <- c()
    
    # fill vector with n samples from the specified chi squared distribution
    for (i in 1:n) {
      v_chisq <- c(v_chisq, chisq_sample(df))
    }
    
    # return the vector
    return(v_chisq)
  }
}



#---------------------------- FUNCTION 3: my_rt() ----------------------------

# my_rt() function returning a vector of pseudo-random t-distributed deviates
# Inputs:
#   - n, number of values to return (the length of the vector)
#   - df, degrees of freedom of the t-distribution 
# Outputs:
#   - vector of pseudo-random values from t-distribution with df degrees of
#     freedom, of length n
my_rt <- function(n, df = 1) {
  
  # Testing arguments n, df are valid: stops execution if not
  if (
    # testing n
    length(n) != 1 || n <= 0 || is.numeric(n) == FALSE ||
    abs(n - round(n)) > 1e-10 || is.finite(n) == FALSE ||
    
    # testing df
    length(df) != 1 || is.numeric(df) == FALSE || is.finite(df) == FALSE ||
    df <= 0 || abs(df - round(df)) > 1e-10
  ) {
    stop("invalid arguments")
  } else {
    # Function extracting a single sample value from the t-distribution
    # with df degrees of freedom
    # Input:
    #   - df, degrees of freedom
    # Output:
    #   - sample value from the specified t-distribution
    t_sample <- function(df) {
      
      # compute t sample value
      t_val <- my_rnorm(1)/(sqrt(my_rchisq(1, df)/df))
      
      # return the sample
      return(t_val)
    }
    
    # Define initially empty vector
    v_tdist <- c()
    
    # fill vector with n samples from the specified t-distribution
    for (i in 1:n) {
      v_tdist <- c(v_tdist, t_sample(df))
    }
    
    # return the vector
    return(v_tdist)
  }
}



#------------ TESTING OUTPUTS OF FUNCTIONS/VALIDITY OF DISTR. ------------------

# Testing output type and similarity of my_distribution functions to their
# R counter parts. The Z-test is used in comparing whether the two distributions
# are similar. Note that this is only accurate for large sample sizes. Also,
# you can choose how many times to run the tests. If a large sample size is
# chosen, then use a smaller number of testing iterations so the computer can 
# handle it. The sample size and test iterations must be chosen with each other
# in mind.



# 1. Function testing the output of my_rnorm() with given n, mean, and sd
# Input:
#   - n, mean, sd as for my_rnorm()
#   - testrep, number of times to repeat testing
# Output:
#   - TEST PASSED or TEST FAILED, see inside function for tests
test_myrnorm <- function(n, mean = 0, sd = 1, testrep = 100) {
  
  # Define empty matrix (used for determining success of tests)
  m <- matrix(0, nrow = testrep, ncol = 4)
  
  # Testing my_rnorm() function test rep times
  for (i in 1:testrep) {
    # Get output
    x <- my_rnorm(n, mean, sd)
    
    # TEST 1: output is a vector of length n
    pass.test <- length(x) == n
    if (pass.test == TRUE) {
      # pass
      m[i, 1] <- 1
    } else {
      # fail
      m[i, 1] <- 0
    }
    
    # TEST 2: output is numeric
    pass.test <- is.numeric(x)
    if (pass.test == TRUE) {
      # pass
      m[i, 2] <- 1
    } else {
      # fail
      m[i, 2] <- 0
    }
    
    # TEST 3: output is finite
    pass.test <- is.finite(x)
    if (setequal(pass.test, rep(TRUE, times = n)) == TRUE) {
      # pass
      m[i, 3] <- 1
    } else {
      # fail
      m[i, 3] <- 0
    }
    
    if (mean != 0 & sd != 0) {
      # TEST 4: Z-statistic of my_rnorm relative to rnorm (<2 implies similar)
      # defining values (Note: n+1 voids error with n = 1; sd(...) = NA )
      my_mean <- mean(my_rnorm(n + 1, mean, sd))
      my_sd <- sd(my_rnorm(n + 1, mean, sd))
      r_mean <- mean(rnorm(n + 1, mean, sd))
      r_sd <- sd(rnorm(n + 1, mean, sd))
      
      # calculating z-statistic
      z_stat <- (my_mean - r_mean)/sqrt(my_sd**2 + r_sd**2)
      
      # running test
      if (z_stat < 2) {
        # pass
        m[i, 4] <- 1
      } else {
        # fail
        m[i, 4] <- 0
      }
    } else {
      # when mean = 0 and sd = 0 the z-statistic = 0, so it passes the test
      m[i, 4] <- 1
    }
  }
  
  # Test 1 result
  if (sum(m[, 1]) == testrep) {
    print("TEST 1 PASSED")
  } else {
    print("TEST 1 FAILED: Not vector of length n")
  }
  
  # Test 2 result
  if (sum(m[, 2]) == testrep) {
    print("TEST 2 PASSED")
  } else {
    print("TEST 2 FAILED: Not numeric")
  }
  
  # Test 3 result
  if (sum(m[, 3]) == testrep) {
    print("TEST 3 PASSED")
  } else {
    print("TEST 3 FAILED: Not finite")
  }
  
  # Test 4 result
  if (sum(m[, 4]) == testrep) {
    print("TEST 4 PASSED")
  } else {
    print("TEST 4 FAILED: Z-statistic is not less than two")
  }
}




# 2. Function testing output of my_rchisq() with given n and df
# Input:
#   - n and df, as for my_rchisq()
#   - testrep, number of times to repeat testing
# Output:
#   - TEST PASSED or TEST FAILED, see inside function for tests
test_myrchisq <- function(n, df = 1, testrep = 100) {
  
  # Define empty matrix (used for determining success of tests)
  m <- matrix(0, nrow = testrep, ncol = 4)
  
  # Testing my_rchisq() function testrep times
  for (i in 1:testrep) {
    # Get output
    x <- my_rchisq(n, df)
    
    # TEST 1: output is a vector of length n
    pass.test <- length(x) == n
    if (pass.test == TRUE) {
      # pass
      m[i, 1] <- 1
    } else {
      # fail
      m[i, 1] <- 0
    }
    
    # TEST 2: output is numeric
    pass.test <- is.numeric(x)
    if (pass.test == TRUE) {
      # pass
      m[i, 2] <- 1
    } else {
      # fail
      m[i, 2] <- 0
    }
    
    # TEST 3: output is finite
    pass.test <- is.finite(x)
    if (setequal(pass.test, rep(TRUE, times = n)) == TRUE) {
      # pass
      m[i, 3] <- 1
    } else {
      # fail
      m[i, 3] <- 0
    }
    
    if (df != 0) {
      # TEST 4: Z-statistic of my_rchisq relative to rchisq (<2 implies similar)
      # defining values (Note: n+1 avoids error with n = 1; sd(...) = NA )
      my_mean <- mean(my_rchisq(n + 1, df))
      my_sd <- sd(my_rchisq(n + 1, df))
      r_mean <- mean(rchisq(n + 1, df))
      r_sd <- sd(rchisq(n + 1, df))
  
      # calculating z-statistic
      z_stat <- (my_mean - r_mean)/sqrt(my_sd**2 + r_sd**2)
      
      # running test
      if (z_stat < 2) {
        # pass
        m[i, 4] <- 1
      } else {
        # fail
        m[i, 4] <- 0
      }
    } else {
      # when df = 0 the z-statistic = 0, so it passes the test
      m[i, 4] <- 1
    }
  }
  
  # Test 1 result
  if (sum(m[, 1]) == testrep) {
    print("TEST 1 PASSED")
  } else {
    print("TEST 1 FAILED: Not vector of length n")
  }
  
  # Test 2 result
  if (sum(m[, 2]) == testrep) {
    print("TEST 2 PASSED")
  } else {
    print("TEST 2 FAILED: Not numeric")
  }
  
  # Test 3 result
  if (sum(m[, 3]) == testrep) {
    print("TEST 3 PASSED")
  } else {
    print("TEST 3 FAILED: Not finite")
  }
  
  # Test 4 result
  if (sum(m[, 4]) == testrep) {
    print("TEST 4 PASSED")
  } else {
    print("TEST 4 FAILED: Z-statistic is not less than two")
  }
}




# 3. Function testing output of my_rt() given n and df
# Input:
#   - n and df, as for my_rt()
#   - testrep, number of times to repeat testing
# Output:
#   - TEST PASSED or TEST FAILED, see inside function for tests
test_myrt <- function(n, df = 1, testrep = 100) {
  
  # Define empty matrix (used for determining success of tests)
  m <- matrix(0, nrow = testrep, ncol = 4)
  
  # Testing my_rt() function testrep times
  for (i in 1:testrep) {
    # Get output
    x <- my_rt(n, df)
    
    # TEST 1: output is a vector of length n
    pass.test <- length(x) == n
    if (pass.test == TRUE) {
      # pass
      m[i, 1] <- 1
    } else {
      # fail
      m[i, 1] <- 0
    }
    
    # TEST 2: output is numeric
    pass.test <- is.numeric(x)
    if (pass.test == TRUE) {
      # pass
      m[i, 2] <- 1
    } else {
      # fail
      m[i, 2] <- 0
    }
    
    # TEST 3: output is finite
    pass.test <- is.finite(x)
    if (setequal(pass.test, rep(TRUE, times = n)) == TRUE) {
      # pass
      m[i, 3] <- 1
    } else {
      # fail
      m[i, 3] <- 0
    }
    
    # TEST 4: Z-statistic of my_rt relative to rt (<2 implies similar)
    # defining values (Note: n+1 avoids error with n = 1; sd(...) = NA)
    my_mean <- mean(my_rt(n + 1, df))
    my_sd <- sd(my_rt(n + 1, df))
    r_mean <- mean(rt(n + 1, df))
    r_sd <- sd(rt(n + 1, df))
      
    # calculating z statistic
    z_stat <- (my_mean - r_mean)/sqrt(my_sd**2 + r_sd**2)
      
    # running test
    if (z_stat < 2) {
      # pass
      m[i, 4] <- 1
    } else {
      # fail
      m[i, 4] <- 0
    }
  }
  
  # Test 1 result
  if (sum(m[, 1]) == testrep) {
    print("TEST 1 PASSED")
  } else {
    print("TEST 1 FAILED: Not vector of length n")
  }
  
  # Test 2 result
  if (sum(m[, 2]) == testrep) {
    print("TEST 2 PASSED")
  } else {
    print("TEST 2 FAILED: Not numeric")
  }
  
  # Test 3 result
  if (sum(m[, 3]) == testrep) {
    print("TEST 3 PASSED")
  } else {
    print("TEST 3 FAILED: Not finite")
  }
  
  # Test 4 result
  if (sum(m[, 4]) == testrep) {
    print("TEST 4 PASSED")
  } else {
    print("TEST 4 FAILED: Z-statistic is not less than two")
  }
}
