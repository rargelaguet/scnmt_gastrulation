#' MLE estimate of Beta-Binomial model
#'
#' @param x An n x 2 data.table, where 1st column keeps total number of trials
#'  and 2nd column number of successes, n is the total number of samples.
#' @param w Vector with initial values of \alpha and \beta, if NULL the method
#'  of moments is used to initialize them.
#'
bb_mle <- function(x, w = NULL, n_starts = 10, lower_thresh = 1e-3){
  require(VGAM)
  require(MASS)
  require(data.table)
  # Ensure that x is a data.table object with 2 columns
  assertthat::assert_that(NCOL(x) == 2)
  if (!is.data.table(x)) { x <- as.data.table(x) }

  is_binomial <- FALSE  # Assume that data follow a Beta-Binomial
  is_conv <- FALSE      # Assume we have no convergence
  best_ll <- -Inf       # Best model is no model...

  # Create a grid of possible initial values around MM estimates
  grid <- c(0.5,-0.5,1,-1)
  grid <- matrix(rep(grid, 2), ncol = 2)

  rho = a = b <- 0  # Set alpha, beta, overdispersion parameters to 0
  bb_ll = AIC_bb = BIC_bb <- 0  # Set all values of the Beta Binomial to zero
  mu <- mean(x[[2]] / x[[1]]) # Compute mean proportion

  # Check if we have data with proportion of success (mu) is ~0 or ~1
  # then assume data coming from a Binomial distribution
  if (mu < lower_thresh) {
    mu <- lower_thresh
    is_binomial <- TRUE
  }else if (mu > 1 - lower_thresh) {
    mu <- 1 - lower_thresh
    is_binomial <- TRUE
  }else{
    # If we have no initial values for \alpha and \beta ...
    if (is.null(w)) {
      # ... compute them using method of moments
      w <- bb_mm(n = x[[1]], k = x[[2]])
      # Check if any parameter is ~0 or even negative, then probably we have
      # un- or under-dispersion or the MM estimate has issues (e.g. different n_i).
      # This might cause problems and might not converge to ML estimate.
      # TODO: Again, assume data coming from a Binomial distribution???
      if (any(w < lower_thresh)) { 
        w[w < lower_thresh] = lower_thresh 
      } else {
        mu <- w[1] / sum(w) # Compute mean from MM
      } 
    }
  }

  # If we assume a Beta-Binomial distribution
  if (!is_binomial) {
    # Initially try the vglm fit function of VGA package
    fit_try <- try(fit_vglm <- vglm(cbind(x[[2]], x[[1]]-x[[2]]) ~ 1, betabinomial, trace = FALSE))
    if (!inherits(fit_try, "try-error")) {
      # Results of the vglm optimization are for \mu and \rho
      mu <- Coef(fit_vglm)[1]
      rho <- Coef(fit_vglm)[2]
      # TODO: Should we really need to to that???
      # Convert them to \alpha and \beta
      tmp <- (1 / rho) - 1
      a <- mu * tmp
      # If mu and rho are close enought to 0, then \alpha will be NaN rather than 1
      if (is.na(a)) { a <- 1 }
      b <- tmp - a
      is_conv <- TRUE   # LABEL that the algorithm has converged
    } else {
      warning("Error in VGLM method. Trying with different initial values.")
      # Create grid around MM estimates of possible initial values
      w_mat <- t(t(grid) + w + rnorm(1, 0, 0.1))
      # Be careful to not introduce negative initial values
      w_mat[w_mat < 0.001] <- 0.001
      w_mat <- rbind(w, w_mat, c(0.1,0.1), c(0.5,0.5), c(1,1), c(3,3), c(4,4))
      # Maximum number of possible initial values
      iter <- min(NROW(w_mat), n_starts)
      for (k in 1:iter) {
        # Fit a Beta-Binomial distribution using Newton's method
        fit_try <- try(fit <- bb_newton(x, w_mat[k, ]), silent = TRUE)
        if (inherits(fit_try, "try-error")) { warning("Error in Newton's method. Trying with different initial values.")}
        # Newton's method is diverging?
        if (fit$conv == 10) { warning("Newton's method not converging. Trying again with different initial values.") }
        # Did not converge yet. Start again using the last value
        else if (fit$conv == 1) {
          if (k < iter) { w_mat[k + 1, ] <- fit$w }
        }else if (any(fit$w < lower_thresh)) { next } # If valuer <=0 try with other initial values
        else{
          # Compute log-likelihood under these parameters
          est_ll <- bb_lik(fit$w, x)
          if (est_ll >= best_ll) {
            # TODO: Should we check if MLE mu change much from MM mu?
            mu <- fit$w[1] / sum(fit$w)  # MLE estimate of mu
            rho <- 1 / (sum(fit$w) + 1)  # rho = 1 / (\alpha + \beta + 1)
            a <- fit$w[1]     # MLE estimates of \alpha and \beta
            b <- fit$w[2]
            is_conv <- TRUE   # LABEL that the algorithm has converged
            best_ll <- est_ll # Update best LL estimate
          }
        }
      }
    }

    # Function did not converge?
    if (!is_conv) {
      warning("Could not find MLE estimates, using method of moments estimates.")
      # Compute parameters using method of moments
      w <- bb_mm(n = x[[1]], k = x[[2]])
      # Check if any parameter is ~0
      if (any(w < lower_thresh)) {
        w[w < lower_thresh] = lower_thresh
        is_binomial <- TRUE
      }else{
        mu <- w[1] / sum(w)
        rho <- 1 / (sum(fit$w) + 1)
        a <- fit$w[1]
        b <- fit$w[2]
      }
    }
  }
  # Compute log likelihood of Binomial with MLE estiamte of probability of success
  bin_ll <- sum(dbinom(x[[2]], x[[1]], prob = mean(x[[2]] / x[[1]]), log = TRUE))
  # AIC_bin <- -2 * bin_ll + 2            # Compute AIC
  BIC_bin <- -2 * bin_ll + log(NROW(x))   # Compute BIC
  if (!is_binomial) {
    bb_ll <- bb_lik(c(a, b), x)             # log likelihood of Beta Binomial
    # AIC_bb <- -2 * bb_ll + 2 * length(w)  # AIC : -2 * logLL + 2 * K
    BIC_bb <- -2 * bb_ll + log(NROW(x)) * length(w)  # BIC : -2 * logLL + log(n) * K
  }
  # Compute likelihood ratio statistic
  lrt <- -2 * (bin_ll - bb_ll)
  # If LRT < 0, the Binomial is better fit, so rho should be better set to 0.
  if (lrt < 0) { rho <- 0 }
  # Compute chi^2 test with 1 degree of freedom
  chi2_test <- pchisq(q = lrt, df = 1, lower.tail = FALSE)

  # Compute Tarone's Z statistic for overdispersion
  p_hat = sum(x[[2]]) / sum(x[[1]])
  S = sum( (x[[2]] - x[[1]] * p_hat)^2 / (p_hat * (1 - p_hat)) )
  Z_score = (S - sum(x[[1]])) / sqrt(2 * sum(x[[1]] * (x[[1]] - 1)))
  # Compute p-values under standard normal for the NULL hypothesis
  z_test <- 2 * pnorm(-abs(Z_score))

  obj <- structure(list(rho = rho, mu = mu, alpha = a, beta = b,
                        is_conv = is_conv, lrt = lrt, chi2_test = chi2_test,
                        Z_score = Z_score, z_test = z_test, bb_ll = bb_ll,
                        BIC_bb = BIC_bb, bin_ll = bin_ll, BIC_bin = BIC_bin),
                   class = "bb_mle")
  return(obj)
}


#' Fit Beta-Binomial model parameters using method of moments approach
#'
#' TODO: We do not take into accout different N for each sample, but we take
#' just the average N of all samples..
#'
#' @param n Vector of total number of trials per sample
#' @param k Vector of successes per subject
#'
bb_mm <- function(n, k){
  # Total number of samples (i.e. cells)
  size <- length(n)
  # Sample moment 1
  m1 <- sum(k) / size
  # Sample moment 2
  m2 <- sum(k^2) / size
  # Number of possible outputs of Binomial distribution.
  # TODO: How to compute this N, mean(n) or max(n) or how??
  # TODO: Why MM sometimes returns negative values, under-dispersion??
  N <- mean(n)
  a <- (N * m1 - m2) / (N * (m2 / m1 - m1 - 1) + m1)
  b <- (N - m1) * (N - m2 / m1) / (N * (m2 / m1 - m1 - 1) + m1)
  return(c(a, b))
}


#' Likelihood of Beta-Binomial wrt to \code{\alpha} and \code{\beta}
#'
#' @param w Vector with the values of \code{\alpha} and \code{\beta} to compute
#'  the gradient.
#' @param x An n x 2 matrix, where 1st column keeps total number of trials
#'  and 2nd column number of successes.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
#'
bb_lik <- function(w, x, is_NLL = FALSE) {
  # Check if both parameters \alpha and \beta are positive
  w <- as.vector(w)
  assertthat::assert_that(w[1] > 0)
  assertthat::assert_that(w[2] > 0)
  if (!is.data.table(x)) { x <- as.data.table(x) }
  # Compute NLL
  f <- sum(dbetabinom.ab(x[[2]], x[[1]], shape1 = w[1], shape2 = w[2], log = TRUE))
  # f <- sum(log(choose(N, K)) + lbeta(K + w[1], N -K + w[2]) - lbeta(w[1], w[2]) )
  # If we required the Negative Log Likelihood
  if (is_NLL) { f <- (-1) * f }
  return(f)
}


#' Gradient of Beta-Binomial wrt to \alpha and \beta
#'
#' TODO: Check that this is right?
#'
#' @param w Vector with the values of \alpha and \beta to compute the gradient.
#' @param x An n x 2 matrix, where 1st column keeps total number of trials
#'   and 2nd column number of successes.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
#'
bb_grad <- function(w, x, is_NLL = FALSE) {
  if (!is.data.table(x)) { x <- as.data.table(x) }
  n <- nrow(x)       # sample size
  N <- x[[1]]        # Total number of trials
  K <- x[[2]]        # Number of successes
  ab_sum <- sum(w)   # Precompute sum of \alpha and \beta
  # Derivative wrt \alpha
  da <- n * (digamma(ab_sum) - digamma(w[1])) + sum(digamma(w[1] + K)     - digamma(ab_sum + N))
  # Derivative wrt \beta
  db <- n * (digamma(ab_sum) - digamma(w[2])) + sum(digamma(w[2] + N - K) - digamma(ab_sum + N))
  gr <- c(da, db)
  # If we required the Negative Log Likelihood
  if (is_NLL) { gr <- (-1) * gr }
  return(gr)
}

#' Compute Hessian matrix of Beta-Binomial wrt to \alpha and \beta
#'
#' TODO: Check that this is right?
#'
#' @param w Vector with the values of \alpha and \beta to compute the gradient.
#' @param x An n x 2 matrix, where 1st column keeps total number of trials
#'   and 2nd column number of successes.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
bb_hessian <- function(w, x, is_NLL = FALSE) {
  if (!is.data.table(x)) { x <- as.data.table(x) }
  n <- nrow(x)       # sample size
  N <- x[[1]]        # Total number of trials
  K <- x[[2]]        # Number of successes
  ab_sum <- sum(w)   # Precompute sum of \alpha and \beta
  # Derivative wrt \alpha
  dada <- n * (trigamma(ab_sum) - trigamma(w[1])) + sum(trigamma(w[1] + K)     - trigamma(ab_sum + N))
  dbda <- n * trigamma(ab_sum)  - sum(trigamma(ab_sum + N))
  # Derivative wrt \beta
  dbdb <- n * (trigamma(ab_sum) - trigamma(w[2])) + sum(trigamma(w[2] + N - K) - trigamma(ab_sum + N))
  h <- matrix(c(dada, dbda, dbda, dbdb), ncol = 2, byrow = T)
  # If we required the Negative Log Likelihood
  if (is_NLL) { h <- (-1) * h }
  return(h)
}


#' Newton method for maximizing Beta Binomial model parameters
#'
#' @param x An n x 2 data.table, where 1st column keeps total number of trials
#'   and 2nd column number of successes, n is the total number of samples.
#' @param w Vector with initial values of \alpha and \beta.
#' @param gamma Step size parameter for Newton's method. Default to 1, other
#'   values gamma in (0, 1) is referred to relaxed Newton's method.
#' @param max_iter Maximum number of iterations for the Newton's method.
#' @param epsilon Numeric denoting the difference of parameters between
#'   consecutive iterations.
bb_newton <- function(x, w = c(0.1, 0.1), gamma = 1, max_iter = 100, epsilon = 1e-5) {
  conv <- 1
  w_prev <- w
  for (i in 1:max_iter) {
    w_cur <- w_prev - gamma * ginv(bb_hessian(w_prev, x)) %*% bb_grad(w_prev, x)
    if (any(w_cur > 2e7)) {
      conv <- 10
      break
    }
    if ((sum(w_cur - w_prev))^2 < epsilon) {
      conv <- 0
      break
    }
    w_prev <- w_cur
  }
  return(list(w = as.vector(w_cur), conv = conv))
}
