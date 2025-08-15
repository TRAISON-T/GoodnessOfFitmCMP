#' @title Simulate Random Samples from the mCMP Distribution
#'
#' @description Generates random samples from a mean-reparametrized Conway-Maxwell Poisson (mCMP) distribution
#'   using inverse transform sampling based on its probability mass function.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param mu A numeric value for the location parameter.
#' @param phi A numeric value for the dispersion parameter.
#' @param max_terms An integer for the maximum number of terms to sum for the normalizing constant in `dcmp`.
#' @param tol A numeric tolerance for the normalizing constant sum in `dcmp`.
#'
#' @return A numeric vector of random samples.
#' @examples
#' # Simulate 10 samples from an equi-dispersed mCMP distribution (similar to Poisson)
#' simulate_mcmp_sample(n = 10, mu = 5, phi = 0)
#'
#' # Simulate 15 samples from an under-dispersed mCMP distribution
#' simulate_mcmp_sample(n = 15, mu = 3, phi = 1)
#'
#' # Simulate 20 samples from an over-dispersed mCMP distribution
#' simulate_mcmp_sample(n = 20, mu = 7, phi = -0.5)
simulate_mcmp_sample <- function(n, mu, phi, max_terms = 1000, tol = 1e-12) {
  # Step 1: Compute probabilities (PMF) using the dcmp function within the package
  x_vals <- 0:max_terms
  probs <- dcmp(x_vals, mu, phi, max_terms = max_terms, tol = tol)

  # Replace NA/NaN/negative probabilities with 0 to ensure valid processing
  probs[!is.finite(probs) | probs < 0] <- 0

  # Normalize probabilities so they sum exactly to 1, as required for a valid PMF
  prob_sum <- sum(probs)
  if (prob_sum == 0) stop("All probabilities are zero. Check mu, phi values.")
  probs <- probs / prob_sum

  # Step 2: Compute Cumulative Distribution Function (CDF)
  cdf_vals <- cumsum(probs)

  # Step 3: Inverse transform sampling (vectorized for efficiency)
  # Generate 'n' random numbers from a uniform distribution [0, 1]
  u <- runif(n)
  # Map each uniform random number to an integer from x_vals based on where it falls
  # within the CDF intervals. '+ 1' is needed because findInterval returns 0 for values
  # less than the first interval boundary.
  samples <- x_vals[findInterval(u, cdf_vals) + 1]

  return(samples)
}
