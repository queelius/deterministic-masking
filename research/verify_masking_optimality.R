#!/usr/bin/env Rscript
# verify_masking_optimality.R
#
# Verify that uniform masking maximizes Fisher information among all masking
# mechanisms satisfying C1 (failed component in candidate set) with fixed
# cardinality w = m-1, for exponential series systems.
#
# General FIM for any masking mechanism:
#   I_{jk} = (1/Lambda) * sum_{c: j in c, k in c} p_{c|j} * p_{c|k} / sigma_c
# where sigma_c = sum_{k in c} p_{c|k} * lambda_k
#
# For w=m-1, candidate sets are c_j = {1,...,m}\{j}, j=1,...,m.
# Masking probabilities: p_{c_j|k} for k != j, with sum_{j!=k} p_{c_j|k} = 1.

library(MASS)  # for ginv if needed

# Compute FIM for w=m-1 given lambda and masking matrix P
# P[j,k] = Pr(C = c_j | K = k) for j != k (0 on diagonal)
# c_j = {1,...,m}\{j}
compute_fim_wm1 <- function(lambda, P) {
  m <- length(lambda)
  Lambda <- sum(lambda)

  # sigma_{c_j} = sum_{k != j} P[j,k] * lambda[k]
  sigma <- numeric(m)
  for (j in 1:m) {
    sigma[j] <- sum(P[j, -j] * lambda[-j])
  }

  # FIM: I_{ell, ell'} = (1/Lambda) * sum_{j: ell in c_j, ell' in c_j} P[j,ell]*P[j,ell'] / sigma[j]
  # c_j contains ell iff j != ell. So the sum is over j != ell AND j != ell'.
  I <- matrix(0, m, m)
  for (ell in 1:m) {
    for (ellp in 1:m) {
      s <- 0
      for (j in 1:m) {
        if (j != ell && j != ellp) {
          s <- s + P[j, ell] * P[j, ellp] / sigma[j]
        }
      }
      I[ell, ellp] <- s / Lambda
    }
  }
  return(I)
}

# Build uniform masking matrix for w=m-1
uniform_masking <- function(m) {
  P <- matrix(1/(m-1), m, m)
  diag(P) <- 0
  return(P)
}

# Generate random masking matrix satisfying constraints
random_masking <- function(m) {
  P <- matrix(0, m, m)
  for (k in 1:m) {
    # For each failed component k, distribute probability over m-1 candidate sets
    probs <- runif(m - 1)
    probs <- probs / sum(probs)  # normalize to sum to 1
    idx <- (1:m)[-k]
    P[idx, k] <- probs
  }
  return(P)
}

# Compare masking mechanisms
compare_mechanisms <- function(lambda, n_random = 1000) {
  m <- length(lambda)

  # Uniform masking FIM
  P_unif <- uniform_masking(m)
  I_unif <- compute_fim_wm1(lambda, P_unif)
  det_unif <- det(I_unif)
  tr_inv_unif <- sum(diag(solve(I_unif)))  # trace of inverse = sum of variances
  max_var_unif <- max(diag(solve(I_unif)))

  cat(sprintf("m=%d, lambda=(%s), Lambda=%.1f\n", m, paste(lambda, collapse=","), sum(lambda)))
  cat(sprintf("Uniform masking:  det(I)=%.6e, tr(I^-1)=%.4f, max_var=%.4f\n",
              det_unif, tr_inv_unif, max_var_unif))

  # Random masking mechanisms
  det_better <- 0
  tr_better <- 0
  max_better <- 0

  worst_det <- det_unif
  worst_tr <- tr_inv_unif
  best_det <- det_unif
  best_tr <- tr_inv_unif

  for (i in 1:n_random) {
    P_rand <- random_masking(m)
    I_rand <- compute_fim_wm1(lambda, P_rand)
    det_rand <- det(I_rand)
    tryCatch({
      inv_rand <- solve(I_rand)
      tr_rand <- sum(diag(inv_rand))
      max_rand <- max(diag(inv_rand))

      if (det_rand > det_unif * 1.001) det_better <- det_better + 1
      if (tr_rand < tr_inv_unif * 0.999) tr_better <- tr_better + 1
      if (max_rand < max_var_unif * 0.999) max_better <- max_better + 1

      if (det_rand < worst_det) worst_det <- det_rand
      if (tr_rand > worst_tr) worst_tr <- tr_rand
      if (det_rand > best_det) best_det <- det_rand
      if (tr_rand < best_tr) best_tr <- tr_rand
    }, error = function(e) NULL)
  }

  cat(sprintf("Random mechanisms beating uniform (out of %d):\n", n_random))
  cat(sprintf("  det(I) larger:    %d\n", det_better))
  cat(sprintf("  tr(I^-1) smaller: %d\n", tr_better))
  cat(sprintf("  max_var smaller:  %d\n", max_better))
  cat(sprintf("Range of det(I):    [%.6e, %.6e] (uniform=%.6e)\n",
              worst_det, best_det, det_unif))
  cat(sprintf("Range of tr(I^-1): [%.4f, %.4f] (uniform=%.4f)\n",
              best_tr, worst_tr, tr_inv_unif))
  cat("\n")

  return(list(det_better=det_better, tr_better=tr_better,
              max_better=max_better, det_unif=det_unif,
              tr_inv_unif=tr_inv_unif))
}

cat("=== Verifying uniform masking optimality for w=m-1 ===\n\n")

# Test 1: m=3, equal rates
cat("--- Test 1: m=3, equal rates ---\n")
compare_mechanisms(c(1, 1, 1))

# Test 2: m=3, unequal rates (paper's example)
cat("--- Test 2: m=3, unequal rates ---\n")
compare_mechanisms(c(1, 3, 5))

# Test 3: m=3, highly unequal rates
cat("--- Test 3: m=3, highly unequal rates ---\n")
compare_mechanisms(c(1, 1, 100))

# Test 4: m=4, equal rates
cat("--- Test 4: m=4, equal rates ---\n")
compare_mechanisms(c(1, 1, 1, 1))

# Test 5: m=4, unequal rates
cat("--- Test 5: m=4, unequal rates ---\n")
compare_mechanisms(c(1, 2, 3, 4))

# Test 6: m=5, unequal rates
cat("--- Test 6: m=5, unequal rates ---\n")
compare_mechanisms(c(1, 2, 3, 4, 5))

# Test 7: m=5, highly unequal rates
cat("--- Test 7: m=5, highly unequal rates ---\n")
compare_mechanisms(c(1, 1, 1, 1, 100))

# Test 8: m=6
cat("--- Test 8: m=6, moderate heterogeneity ---\n")
compare_mechanisms(c(1, 2, 3, 5, 8, 13))

cat("\n=== Also verify: w=1 (no masking) always beats w=m-1 uniform ===\n\n")
for (lam_set in list(c(1,3,5), c(1,1,1), c(1,1,100), c(1,2,3,4))) {
  m <- length(lam_set)
  Lambda <- sum(lam_set)

  # w=1: exact identification. FIM = diag(p_j / lambda_j^2) / Lambda... no wait
  # For w=1, each obs identifies the failed component exactly.
  # Score: d/d lambda_j [log lambda_k - Lambda*t] = 1_{k=j}/lambda_j - t
  # FIM_{jk} = E[1_{K=j}/lambda_j^2] for j=k, and cross terms
  # Actually: d^2/d lambda_j d lambda_k [log lambda_K - Lambda t] = -1_{K=j}/lambda_j^2 * 1_{j=k}
  # So FIM = diag(p_j / lambda_j^2) = diag(1/(lambda_j * Lambda))
  # Var(lambda_j) = lambda_j * Lambda / n
  I_w1 <- diag(1 / (lam_set * Lambda))
  vars_w1 <- lam_set * Lambda

  # w=m-1 uniform:
  P_unif <- uniform_masking(m)
  I_wm1 <- compute_fim_wm1(lam_set, P_unif)
  vars_wm1 <- diag(solve(I_wm1))

  cat(sprintf("lambda=(%s), Lambda=%.0f\n", paste(lam_set, collapse=","), Lambda))
  cat(sprintf("  w=1 vars (×n):     %s\n", paste(sprintf("%.2f", vars_w1), collapse=", ")))
  cat(sprintf("  w=m-1 unif vars:   %s\n", paste(sprintf("%.2f", vars_wm1), collapse=", ")))
  cat(sprintf("  Ratio (wm1/w1):    %s\n", paste(sprintf("%.2f", vars_wm1/vars_w1), collapse=", ")))
  cat("\n")
}

cat("=== Summary ===\n")
cat("If uniform masking is optimal, no random mechanism should beat it.\n")
cat("If w=1 always beats w=m-1, all ratios should be > 1.\n")
