# Group 14

################################################################################
# SETUP
################################################################################

# Set working directory where data is located
setwd("~/Documents/Studie/MSc Econometrics/Multivariate Statistics/Assignments/Week 1")

# Load the data
load("supermarket1996.RData")

# Overview of the variables included in the dataframe
str(supermarket1996)
summary(supermarket1996)

# Select predictor variables age9 to shpindx
x = as.matrix(supermarket1996[ ,c(6:50)])
# Define dependant variable
y = as.vector(supermarket1996$GROCCOUP_sum)

# PART A
################################################################################
# Majorization Function for Elastic Net
################################################################################
# DESCRIPTION: Majorization algorithm for the elastic net
# INPUTS
# x: A general matrix containing the predictor variables (not a df)
# y: The dependant variable
# lambda: Strength of the penalty term, default is set equal to 1
# alpha: mixing paramater of Ridge and LASSO, default is set equal to 1 (LASSO)
# verbose: control output of updating during algorithm, default is TRUE
# eps: stopping criterium, if the relative imrpovement falls below this threshold
# the algorithm is stopped. Default is set equal to 10^-9 
# OUTPUTS
# b: vector with regression weights

elasticnet <- function(x, y, lambda = 1, alpha = 1, verbose=T, eps = 1e-09) {
  
  # Define amount of variables in x
  n = ncol(x)
  
  # Initialize some b vector, we use a vector with ones
  b <- as.matrix(rep(1, n))
  
  # Compute L_en(b)
  c = (2*n^-1)*(t(y)%*%y) + lambda*(0.5*alpha*sum(abs(b)))
  d_jj = 1/(2*max(abs(b), eps))
  D = diag(d_jj, n)
  
  L_en_0 = t(b) %*% ( (2*n)^-1 * t(x)%*%x + (lambda*0.5*(1-alpha))*diag(n) + 
                      lambda*alpha*D ) %*% b - (2*n)^-1 * (t(b)%*%t(x)) %*% y + c
  
  # Update estimates
  t = 0
  while (t == 0 | L_en[ ,t] - (L_en[ ,t] / L_en[ ,t+1]) > eps ) {
    t = t + 1
    d_jj = 1/(2*max(abs(b), eps))
    D = diag(d_jj, n)
    A = (2*n)^-1 * (t(x)%*%x) + 0.5*(lambda*(1-alpha))*diag(n) + lambda*alpha*D
    b = solve(A, t(x)%*%y*(2*n)^-1)
    
    L_en = cbind(L_en, t(b) %*% ( (2*n)^-1 * t(x)%*%x + (lambda*0.5*(1-alpha))*diag(n) + 
                                    lambda*alpha*D ) %*% b - 
                   (2*n)^-1 * (t(b)%*%t(x)) %*% y + c)
    
    # Print steps if verbose is set to TRUE
    if (verbose == T) {
      print(paste0("t=:", t, "   L_en=", L_en[ ,t]))
    }
  }
  return(b)
}

################################################################################
# TEST FUNCTION
################################################################################

elasticnet(x, y)

################################################################################
# USE GLMNET PACKAGE FOR COMPARISON
################################################################################

# Load package
library(glmnet)
# Fit elastic net (LASSO, alpha=1)
check = cv.glmnet(x, y, alpha=1)
# Plot lambda's after CV
plot(check)
# Coefficient estimates after shrinking
coef(check, s = "lambda.min")
