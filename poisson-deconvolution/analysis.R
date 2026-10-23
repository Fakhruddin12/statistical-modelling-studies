# This assignment was completed collaboratively by 
#   Mansi Neb: model matrices, penalty matrix, likelihood, and gradient
#   Riccardo Monforte: fitting and lambda grid
#   Fakhruddin Hakim Hussain: bootstrap, plotting and write-up
# Each member contributed approximately one-third (33%) to the total effort.

# ------------------------------------------------------------------------------#

# GitHub Repo Link: 
# https://github.com/fakhruddin-12/Extended_Statistical_programming_Assignment_1/tree/main

# ------------------------------------------------------------------------------#

# This code reconstructs the trajectory of new infections each day that eventually
# resulted in deaths in English hospitals from Covid-19 during 2020. It uses a 
# Poisson deconvolution model with B-splines to represent infections, applies a 
# smoothness penalty to avoid overfitting, and selects the best smoothing parameter
# via BIC. The final output includes fitted death curves, estimated infections, 
# and 95% confidence intervals.

# ------------------------------------------------------------------------------#

## Setup and Load Data
library(splines) # for splineDesign

# Reading the raw data
# setwd("add path")
a <- read.table("engcov.txt", header = TRUE)

y <- a$nhs      # observed deaths per day (vector length n)
t <- a$julian   # day-of-year values
n <- length(y)  # number of observations

# Evaluate the model probability of each fatal disease duration
# log(d) ~ N(3.151, 0.469^2)
d <- 1:80; edur <- 3.151; sdur <- 0.469
pd <- dlnorm(d, edur, sdur)   # probability of dying after j days after infection
pd <- pd / sum(pd)    

## Build model matrices X_tilde, X, and S

build_matrices <- function(t, K, pd) {
  # t: vector of julian days (observed death days)
  # K: number of B-spline basis functions
  # pd: probability distribution for infection-to-death duration
  # Returns a list of matrices:
  #   X_tilde: (n+30 x K) B-spline basis evaluated at infection times
  #   X: (n x K) model matrix for deaths
  #   S: (K x K) penalty matrix
  
  n <- length(t)
  t1 <- min(t); tn <- max(t)  # first and last observation day
  
  # Time grid for f(t): from (t1 - 30) to tn
  # We need f(t) for all infection times that could lead to observed deaths
  time_f <- seq(t1 - 30, tn); n_f <- length(time_f)
  
  # Define K+4 knots: middle K-2 span the infection period, 2 on each side for boundaries
  tmin <- min(time_f); tmax <- max(time_f); trange <- tmax - tmin
  pad <- 0.01 * trange  # add small padding on ends
  
  knots <- seq(tmin - pad, tmax + pad, length.out = K + 4)
  
  # Build X_tilde: row is a time point, column is a basis function
  X_tilde <- splineDesign(knots = knots, x = time_f, ord = 4, outer.ok = TRUE)
  
  # Build X matrix by convolving X_tilde with pd
  # For death on day ti, sum over all possible infection times:
  # mu_i = sum{j=1}^{min(29+i,80)} f(ti - j) * pd[j]
  X <- matrix(0, n, K); max_dur <- length(pd)
  
  for (i in 1:n) {
    jmax <- min(29 + i, max_dur)
    
    # compute all infection times
    t_inf <- t[i] - 1:jmax
    
    # convert to X_tilde row indices
    idx <- t_inf - time_f[1] + 1
    
    # find rows of X by calculating weighted sum of basis rows
    X[i, ] <- colSums(X_tilde[idx, , drop=FALSE] * pd[1:jmax])
  }
  
  # Build penalty matrix S for second differences
  S <- crossprod(diff(diag(K), differences = 2))
  
  list(X_tilde = X_tilde, X = X, S = S, time_f = time_f)
}

## Build matrices
K <- 80  # number of B-spline basis functions

mat <- build_matrices(t = t, K = K, pd = pd)
X_tilde <- mat$X_tilde
X <- mat$X
S <- mat$S

## Build Likelihood and Gradient Functions

# Penalized negative log-likelihood
pnll <- function(gamma, y, X, S, lambda , w = NULL) {
  # gamma: log-transformed parameters (length K)
  # y: observed deaths (length n)
  # X: model matrix (n x K)
  # S: penalty matrix (K x K)
  # lambda: smoothing parameter; w: weight (needed while bootstrapping)
  ## Returns a scalar penalized nll value
  
  if(is.null(w)) w <- rep(1 , length(y))
  beta <- exp(gamma)   # to ensure positive f(t)
  mu <- X %*% beta
  
  mu[mu <= 0] <- 1e-9   # add numerical stability check
  
  # Poisson nll (dropped constant factorial term)
  nll <- -sum(w*(y * log(mu) - mu))
  
  # Penalty term
  penalty <- (lambda / 2) * sum(beta *(S %*% beta))
  
  nll + penalty
}

# Gradient of penalized negative log-likelihood w.r.t. gamma
grad_pnll <- function(gamma, y, X, S, lambda , w = NULL) {
  # Inputs: gamma0 (params), y (death data), X (model matrix), S (penalty matrix), 
  #         lambda (smoothing param), w (weight)
  ## Returns a gradient vector (length K)
  
  if(is.null(w)) w <- rep(1 , length(y))
  beta <- exp(gamma)
  mu <- X %*% beta
  mu[mu <= 0] <- 1e-9
  
  # Gradient of nll w.r.t. gamma
  # d(nll)/d(gamma) = diag(beta) * d(nll)/d(beta)
  # where d(nll)/d(beta) = -X^T * (y/mu - 1)
  resid <- w*(y / mu - 1)
  grad_nll <- -beta * as.vector(crossprod(X, resid))
  
  # Gradient of penalty w.r.t. gamma
  # d(P)/d(gamma) = diag(beta) * d(P)/d(beta) = diag(beta) * lambda * S * beta
  grad_penalty <- lambda * beta * (S %*% beta)
  
  grad_nll + grad_penalty
}


## Gradient Test Function
test_gradient <- function(gamma0, y, X, S, lambda, eps = 1e-7) {
  # Purpose: Verifies the analytical gradient (grad_pnll) using finite differences.
  # Inputs: gamma0 (params), y (deathdata), X (model matrix), S (penalty matrix), 
  #         lambda (smoothing param), eps (step size)
  ## Returns: Prints the test result (PASSED/FAILED) and 
  #           max difference between analytical and FD gradient
  
  # Compute analytical gradient
  grad_a <- grad_pnll(gamma0, y, X, S, lambda)
  
  # Compute finite difference gradient
  nll0 <- pnll(gamma0, y = y, X = X, S = S, lambda = lambda)
  grad_fd <- numeric(length(gamma0))
  
  for (i in 1:length(gamma0)) {  #loop over parameters
    gamma1 <- gamma0
    gamma1[i] <- gamma1[i] + eps  # increase gamma1[i] by eps
    nll1 <- pnll(gamma1, y = y, X = X, S = S, lambda = lambda) # compute resulting nll
    grad_fd[i] <- (nll1 - nll0) / eps  # Approximate gradient (slope)
  }
  
  # find max difference between analytical and FD gradient
  max_diff <- max(abs(grad_a - grad_fd))
    cat("Max difference between analytic and FD gradient:", max_diff, "\n")
  
  # Threshold check
  if (max_diff < 1e-3) cat("Gradient test PASSED\n")
   else cat("Gradient test FAILED\n")
}

## Run gradient test

# Set lambda and starting values for testing
lambda0 <- 5e-5
# Initial guess for beta say average daily deaths, mean(y)
beta0 <- rep(mean(y), K) 
gamma0 <- log(beta0)  #as beta = exp(gamma) --> gamma = log(beta)
test_result <- test_gradient(gamma0, y, X, S, lambda0)

## Sanity Check Fit
# Fit the model to estimate parameters gamma that minimizes pnll with Lambda = 5e-5
# Now we fit the model to estimate the parameters gamma that minimizes the 
# penalized negative log likelihood with Lambda = 5e-5
# Then, plot the observed deaths, fitted deaths, and estimated infections
# Where observed deaths are y, fitted deaths are mi = X * beta_hat, and 
# estimated infections are f(t) = X_tilde * beta_hat

#Fit the model using optim function
fit0 <- optim(gamma0,pnll,gr=grad_pnll, y = y, X=X,S=S,lambda=lambda0,
              method="BFGS",hessian=FALSE)  #used initial guess of gamma,lambda from gradient testing

# Extract estimated parameters
gamma0.hat <- fit0$par
beta0.hat <- exp(gamma0.hat)

# Expected deaths per day
mhat <- drop((X %*% beta0.hat))

# Expected infections per day
fhat <- drop((X_tilde %*% beta0.hat))

# ---  Sanity Check Plot ---
# Plotting actual vs fitted deaths along with fitted daily infection curve f

# Setting up the time range of infections, time grid of f(t)
t.inf <- mat$time_f

# Setup Plot Limits
ylim <- range(c(y, fhat))
xlim <- range(t.inf)

# 1. Base Plot: Observed Deaths
plot(t, y, 
     type = "p", pch = 19, col = "grey", cex = 0.8,
     xlim = xlim, 
     ylim = ylim,
     xlab = "Day of Year (Julian)", 
     ylab = "Count",
     main = "Preliminary Sanity Check: Deaths vs Est. Infections")

grid() # Add grid for readability

# 2. Add Curves
lines(t.inf, fhat, col = "blue", lwd = 1.5) # Estimated Infections
lines(t, mhat, col = "black", lwd = 1.5)         # Fitted Deaths

# 3. Add Legend
legend("topright", 
       legend = c("Observed Deaths", "Fitted Deaths", "Est. Infections"),
       col = c("grey", "black", "blue"),
       pch = c(19, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       bty = "n", cex = 0.8)


# --- BIC Grid Search ---
# Implement a grid search over log(lambda) to find the optimal smoothing parameter
# BIC = 2*nll + log(n)*edf
# Where nll is the Poisson negative log-likelihood (without penalty term)
# edf is the effective degrees of freedom defined as: edf = trace(H_lambda_inv * H_0)
# H_lambda is defined as: X^T * W * X + lambda * S
# H_0 is defined as: X^T * W * X
# W is a diagonal matrix with entries W_ii = y_i / mhat_i^2

log_lambda_grid <- seq(-13,-7,length=50) #Given log(lambda) grid

# Function to compute optimal lambda using BIC
optim_lambda <- function(log_lambda_grid, gamma0, y, X, S, grad_pnll, pnll) {
  
  # Create vector to store BIC values for each lambda
  bic_vector <- rep(0, length(log_lambda_grid))
  
  # Loop over each log(lambda) in the grid
  for (i in 1:length(log_lambda_grid)) {
    
    log_lambda <- log_lambda_grid[i]   #Extract current log(lambda)
    lambda <- exp(log_lambda)    #Convert to lambda scale
    
    # Fit the model using optim
    #(gamma.hat --> beta.hat --> mhat --> W --> H_lambda, H_lambda_inv, H_0 --> edf --> BIC)
    fit <- optim(gamma0, pnll, grad_pnll, y, X, S, lambda = lambda,
                 method = "BFGS", hessian = FALSE)
    
    #Extract estimated parameters
    gamma.hat <- fit$par
    beta.hat <- exp(gamma.hat)
    
    # Compute fitted values
    mhat <- drop((X %*% beta.hat))
    
    # Compute W diagonal entries and weighted X matrix
    # Weight matrix W = diag(y / mhat^2)
    W <- y / mhat^2
    wsqrt <- sqrt(W)
    X_w <- X * wsqrt
    
    # Compute H_0 = X^T * W * X
    H_0 <- crossprod(X_w)
    
    # Compute H_lambda = H_0 + lambda * S
    H_lambda <- H_0 + lambda * S
    
    # Compute H_lambda_inv using Cholesky
    H_lambda_inv <- chol2inv(chol(H_lambda))
    
    # Instead of sum(diag(H_lambda_inv %*% H_0)) which is O(N^3),
    # we use sum(H_lambda_inv * H_0) which is O(N^2) and mathematically equivalent
    edf <- sum(H_lambda_inv * H_0)
    
    # Compute mu and nll(without penalty)
    mu <- X %*% beta.hat
    mu[mu <= 0] <- 1e-9
    nll <- -sum(y * log(mu) - mu)
    
    # Compute BIC and store in bic_vector for the current lambda
    n <- length(y)
    bic_vector[i] <- 2*nll + log(n) * edf
  }  
  
  #Find the lambda that minimizes BIC 
  best_idx <- which.min(bic_vector)
  best_log_lambda <- log_lambda_grid[best_idx]
  best_lambda <- exp(best_log_lambda)
  
  list(best_lambda = best_lambda, bic_vector = bic_vector, log_lambda_grid = log_lambda_grid)
}

# Run the optim_lambda function to find the best lambda
best_lambda <- optim_lambda(log_lambda_grid, gamma0.hat, y, X, S, 
                            grad_pnll, pnll)$best_lambda

## Refit model with best lambda
fit_best <- optim(gamma0.hat, pnll, grad_pnll, y = y, X = X, S = S, 
                  lambda = best_lambda, method = "BFGS")
f_best <- X_tilde %*% exp(fit_best$par)


## Bootstrapping to estimate uncertainty (95% CI)

## Finding the uncertainty in our estimate of the function 'f'.
# (f(t) signifies new infections at time t)

# We need to estimate the uncertainty in our estimate as it is build using
# noisy observed death count data.
# This means that there would be statistical uncertainty associated with our
# Beta estimates of the B-splines and hence in the fitted curve.

# This can be done using bootstrapping.
# Here, we will sample(with replacement) from our original data and will fit
# the curve using this data. We will do this multiple times (lets say B times).
# Now for each time t, we will have B values of the curve. With these B values,
# we can find the 95% confidence interval for that particular time t.

bootstrapping_f <- function(samples , gamma0 , pnll , grad_pnll , y ,X_tilde, X , S ,
                             lambda , time_boots)
{
  # Defining an empty matrix to store the bootstrapping results
  f_boot_matrix <- matrix(NA , nrow = time_boots , ncol = samples)
  
  # Iterating
  for (i in 1:samples){
    
    # Defining weights which is just the number of times each row is sampled
    weights = tabulate(sample(length(y) , replace = TRUE) , length(y))
    
    # fitting the model to using the weights.
    fit_b <- optim(gamma0, pnll, grad_pnll, y, X, S, lambda = lambda, w = weights,
                   method = "BFGS", hessian = FALSE)
    # Best gamma
    gamma.hat <- fit_b$par
    # Best beta corresponding to the best gamma
    beta.hat <- exp(gamma.hat)
    
    # The expected infection function corresponding to the best beta
    fhat <- drop((X_tilde %*% beta.hat))
    
    # Storing the expected infection function
    f_boot_matrix[, i] <- fhat
  }
  f_boot_matrix
}

f_mat <- bootstrapping_f(200 , fit_best$par , pnll , grad_pnll , y ,X_tilde, X , S , 
                          best_lambda , length(t.inf))

# Calculate 2.5% and 97.5% quantiles
f_ci <- apply(f_mat , 1 , quantile , probs = c(0.025 , 0.975) , na.rm = TRUE)

## Final Plotting
# Plot observed deaths, fitted deaths, and estimated infections with 95% CI

#Set plot limit
ylim1 <- range(c(y, f_ci))   #xlim same as perviously plotted graphs: t.inf

plot(t, y, 
     type = "p", pch = 19, col = "grey", cex = 0.8,
     xlim = xlim, 
     ylim = ylim1,
     main = "Covid-19: Deaths vs Estimated Fatal Infections (2020)",
     xlab = "Day of Year (Julian)", 
     ylab = "Count")

grid() 

# Infection Confidence Interval (Shaded)
polygon(c(t.inf, rev(t.inf)), 
        c(f_ci[1, ], rev(f_ci[2, ])),
        col = rgb(0, 0, 1, 0.3), # blue with transparency
        border = NA)

# Estimated Infection Rate f(t) (Blue line)
lines(t.inf, f_best, col = "blue", lwd = 1.5)

# Fitted Death Rate mu(t) (Black line)
mu_best <- X %*% exp(fit_best$par)
lines(t, mu_best, col = "black", lwd = 1.5)

legend("topright", 
       legend = c("Observed Deaths", "Fitted Deaths", "Est. Infections", "95% CI (Infections)"),
       col = c("grey", "black", "blue", NA),
       pch = c(19, NA, NA, NA),
       lty = c(NA, 1, 1, NA),
       lwd = c(NA, 2, 2, NA),
       fill = c(NA, NA, NA, rgb(0, 0, 1, 0.2)),
       border = NA,
       bty = "n", cex = 0.8)

# find the specific days where the curves peak
peak_inf <- t.inf[which.max(f_best)]  # day of max infections
peak_dth <- t[which.max(mu_best)]          # day of max fitted deaths

cat("Estimated time lag:", peak_dth - peak_inf, "days\n")

# Results:
#- The fitted death curve closely matches the observed daily deaths,
#  indicating the model fits well.
#- The estimated infection curve peaks ~20 days before the death peak,
#  consistent with known infection-to-death delays for COVID-19.
#- The infection peak is noticeably higher than the death peak, as
#  expected since only a fraction of infections lead to death.
#- Bootstrap confidence bands are tight around the epidemic peak and
#  wider in early and late periods, reflecting uncertainty where fewer
#  data points are informative.