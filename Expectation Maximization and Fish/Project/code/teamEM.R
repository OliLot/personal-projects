# --------------------------  EM algorithm function ---------------------------

# Function that runs the EM algorithm to estimate parameters
# INPUTS:
#   - data = a dataframe of observations for fish (FishID, Length, Age)
#   - epsilon = tolerance value (default = 1e-8)
#   - maxit = maximum number of iterations (default = 1000)
# OUTPUTS:
#   - estimates: dataframe of estimates with rows (Age1, Age2, Age3), and
#                columns (lambda, mu, sigma)
#   - inits: dataframe of initial values with rows (Age1, Age2, Age3) and
#            columns (lambda, mu, sigma)
#   - converged: TRUE/FALSE boolean
#   - posterior: N x k dataframe of posterior probabilities (col = group)
#   - likelihood: vector of log-likelihood values (length is nr of iterations)
teamEM <- function(data, epsilon = 1e-8, maxit = 1000) {
  
  # // Testing inputs: //
  if (is.null(data$FishID) == TRUE || is.null(data$Length) == TRUE || 
      is.null(data$Age) == TRUE || 
      length(data$FishID[is.na(data$FishID) == TRUE]) != 0 ||
      length(data$Length[is.na(data$Length) == TRUE]) != 0 ||
      length(data$Age[is.na(data$Age) == TRUE]) != 0 ||
      is.numeric(epsilon) == FALSE || is.finite(epsilon) == FALSE ||
      is.numeric(maxit) == FALSE || is.finite(maxit) == FALSE || maxit <= 0) {
    
    # invalid inputs, stop iteration
    stop("Invalid arguments")
  }
  
  # // Initialise: //
  # Given rough age group assignments, compute initial values of mean_k and 
  # sd_k for each age group k = 1, 2, 3.
  mean_init <- aggregate(data$Length, list(data$Age), FUN = mean)
  colnames(mean_init) <- c("Group", "Mean guess")
  
  sd_init <- aggregate(data$Length, list(data$Age), FUN = sd)
  colnames(sd_init) <- c("Group", "Sd guess")
  
  # Compute initial values of lamdda_k = probability that fish is in group k
  fishes_per_group <- aggregate(data$Age, by = list(data$Age), FUN = length)
  colnames(fishes_per_group) <- c("Group", "Fish")
  
  lambda_init <- c()
  for (i in 1:3) {
    lambda_init[i] <- fishes_per_group[i, 2]/length(data$FishID)
  }
  lambda_init <- data.frame("Group" = c(1, 2, 3), "lambda" = lambda_init)
  
  # Create tables containing lambda, mean, and sd values for each group
  # (new estimates will be added as new columns)
  for (i in 1:3) {
    initial <- c(lambda_init[i, 2], mean_init[i, 2], sd_init[i, 2])
    assign(paste0("group", i), 
           data.frame(initial, row.names = c("lambda", "mean", "sd")))
  }
  
  # Log likelihood for initial estimates (will append new values from iteration)
  loglikelihoods <- c(loglikelihood(data, lambda_init[, 2], mean_init[, 2], 
                                    sd_init[, 2]))
  
  # // Iterate: //
  for (i in 1:maxit) {
    
    # // Expectation step: get posterior probabilities for each group //
    # Function that calculates the posterior probability for all fish
    # Inputs:
    #   - data = the data frame including columns for FishID, Length, Age
    #   - lambda_k = estimate lambda value for the kth group
    #   - mean_k = estimate mean fish length for the kth group
    #   - sd:k = estimate sd value for the kth group
    # OUTPUTS:
    #   - vector of posterior probabilities for fish being in the kth group, 
    #     where the ith element represents the probability of the ith fish
    expectation <- function(data, lambda_k, mean_k, sd_k) {
      
      # Testing inputs
      if (lambda_k < 0 || lambda_k > 1 || mean_k <= 0 || sd_k < 0 || 
          is.numeric(c(lambda_k, mean_k, sd_k)) == FALSE ||
          all(is.finite(c(lambda_k, mean_k, sd_k))) == FALSE) {
        
        stop("Invalid arguments")
      }
      
      
      p_expectation <- (lambda_k * dnorm(data$Length, mean_k, sd_k)) / (
        (group1[1, i] * dnorm(data$Length, group1[2, i], group1[3, i])) +
          (group2[1, i] * dnorm(data$Length, group2[2, i], group2[3, i])) +  
          (group3[1, i] * dnorm(data$Length, group3[2, i], group3[3, i])) 
      )
      
      return(p_expectation)
    }
    
    # the posterior probabilities
    group1_exp <- expectation(data, group1[1, i], group1[2, i], group1[3, i])
    group2_exp <- expectation(data, group2[1, i], group2[2, i], group2[3, i])
    group3_exp <- expectation(data, group3[1, i], group3[2, i], group3[3, i])
    
    
    # // Maximization step: get new estimates for parameters //
    # new mean estimates
    new_mean1 <- calc_weighted_mean(data, group1_exp)
    new_mean2 <- calc_weighted_mean(data, group2_exp)
    new_mean3 <- calc_weighted_mean(data, group3_exp)
    
    # new sd estimates
    new_sd1 <- calc_weighted_sd(data, group1_exp, new_mean1)
    new_sd2 <- calc_weighted_sd(data, group2_exp, new_mean2)
    new_sd3 <- calc_weighted_sd(data, group3_exp, new_mean3)
    
    # new lambda estimates
    new_lambda1 <- calc_weight(data, group1_exp)
    new_lambda2 <- calc_weight(data, group2_exp)
    new_lambda3 <- calc_weight(data, group3_exp)
    
    
    # // Test convergence: //
    # vectors of new estimates
    new_lambdas <- c(new_lambda1, new_lambda2, new_lambda3)
    new_means <- c(new_mean1, new_mean2, new_mean3)
    new_sds <- c(new_sd1, new_sd2, new_sd3)
    
    # if converged: retrieve values, else: continue iteration with new estimates
    loglikelihoods <- c(loglikelihoods, loglikelihood(data, new_lambdas, 
                                                      new_means, new_sds))
    
    if (loglikelihoods[i + 1] - loglikelihoods[i] < epsilon) {
      
      # It has converged
      converged <- TRUE
      
      # Stop iteration
      break
      
    } else if (i == maxit) {
      
      # Not converged after maxit iterations
      converged <- FALSE
      
    } else {
      
      # Continue iterating
      # add new estimates for each age group as columns in the group data frame
      new_estimates1 <- c(new_lambda1, new_mean1, new_sd1)
      group1 <- cbind(group1, new_estimates1)
      
      new_estimates2 <- c(new_lambda2, new_mean2, new_sd2)
      group2 <- cbind(group2, new_estimates2)
      
      new_estimates3 <- c(new_lambda3, new_mean3, new_sd3)
      group3 <- cbind(group3, new_estimates3)
    }
  }
  
  # // Get required outputs as a single list //
  estimates <- data.frame(new_lambdas, new_means, new_sds)
  colnames(estimates) <- c("lambda", "mu", "sigma")
  rownames(estimates) <- c("Age1", "Age2", "Age3")
  
  inits <- data.frame(lambda_init[, 2], mean_init[, 2], sd_init[, 2])
  colnames(inits) <- c("lambda", "mu", "sigma")
  rownames(inits) <- c("Age1", "Age2", "Age3")
  
  posterior <- cbind(group1_exp, group2_exp, group3_exp)
  
  likelihood <- loglikelihoods[-1]
  
  output <- list("estimates" = estimates, "inits" = inits,
                 "converged" = converged, "posterior" = posterior, 
                 "likelihood" = likelihood)
  
  return(output)
}
