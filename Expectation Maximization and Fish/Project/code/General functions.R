# --------------------- General functions to be used ---------------------------

# Function calculating the weighted lambda value for the kth group
# INPUTS:
#   - data = the data frame including columns for FishID, Length, Age
#   - expectations_k = posterior probability for all fish to be in kth group
# OUTPUTS:
#   - weighted lambda value for kth group
calc_weight <- function(data, expectations_k) {
  
  # Testing inputs
  if (is.numeric(expectations_k) == FALSE || 
      length(expectations_k) != length(data$FishID) ||
      all(expectations_k >= 0) == FALSE || 
      all(expectations_k <= 1) == FALSE) {
    
    stop("Invalid arguments")
  }
  
  # Calculating using formula
  sum(expectations_k)/length(data$FishID)
}



# Function calculating the weighted mean length of fish for the kth group
# INPUTS:
#   - data = the data frame including columns for FishID, Length, Age
#   - expectations_k = posterior probability for all fish to be in kth group 
# OUTPUTS:
#   - weighted mean of fish length for kth group
calc_weighted_mean <- function(data, expectations_k) {
  
  # Testing inputs
  if (is.numeric(expectations_k) == FALSE || 
      length(expectations_k) != length(data$FishID) ||
      all(expectations_k >= 0) == FALSE || 
      all(expectations_k <= 1) == FALSE) {
    
    stop("Invalid arguments")
  }
  
  # Calculating using formula
  return(sum(expectations_k * data$Length) / sum(expectations_k))
}



# Function calculating the weighted standard deviation for the kth group
# INPUTS:
#   - data = the data frame including columns for FishID, Length, Age
#   - expectations_k = posterior probability for all fish to be in kth group
#   - weighted_mean_k = the weighted mean fish length of the kth group
# OUTPUTS:
#   - weighted standard deviation of fish lengths for the kth group
calc_weighted_sd <- function(data, expectations_k, weighted_mean_k) {
  
  # Testing inputs
  if (is.numeric(expectations_k) == FALSE || 
      length(expectations_k) != length(data$FishID) ||
      all(expectations_k >= 0) == FALSE || 
      all(expectations_k <= 1) == FALSE ||
      is.numeric(weighted_mean_k) == FALSE || 
      is.finite(weighted_mean_k) == FALSE || weighted_mean_k <= 0) {
    
    stop("Invalid arguments")
  }
  
  # Calculating using formula
  return(sqrt(sum(expectations_k * (data$Length - weighted_mean_k)^2 / 
                    sum(expectations_k))))
}



# Function calculating the log likelihood of data given the new estimates
# INPUTS:
#   - data = data set being used (with columns FishID, Length, Age)
#   - lambda_ks = vector of new estimates to lambda_k of each group
#   - mean_ks = vector of new estimates to mean_k of each group
#   - variance_ks = vector of new estimates to variance_k of each group
# OUTPUTS:
#   - Log likelihood of data given the input estimates
loglikelihood <- function(data, lambda_ks = new_lambdas, mean_ks = new_means, 
                          sd_ks = new_sds) {
  
  # Testing Inputs
  if (length(lambda_ks) != 3 || length(mean_ks) != 3 || length(sd_ks) != 3 ||
      is.numeric(c(lambda_ks, mean_ks, sd_ks)) == FALSE || 
      all(is.finite(c(lambda_ks, mean_ks, sd_ks))) == FALSE || 
      all(mean_ks > 0) == FALSE ||  all(sd_ks >= 0) == FALSE || 
      1 - sum(lambda_ks) > 1e-8 || all(lambda_ks >= 0) == FALSE || 
      all(lambda_ks <= 1) == FALSE) {
    
    print(lambda_ks)
    stop("Invalid arguments")
  }
  
  
  # Defining f(y) as shown in the assignment script
  f <- function(y_i) {
    
    sum <- 0
    for (k in 1:3) {
      sum <- sum + lambda_ks[k] * dnorm(y_i, mean_ks[k], sd_ks[k])
    }
    
    return(sum)
  }
  
  # Calculating log likelihood
  log_likelihood <- 0
  for (i in 1:length(data$FishID)) {
    
    log_likelihood <- log_likelihood + log(f(data$Length[i]))
  }
  
  # after iterating over all fish, return the result
  return(log_likelihood)
}

