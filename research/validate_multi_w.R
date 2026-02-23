#!/usr/bin/env Rscript
# validate_multi_w.R
#
# Verify Theorem 3.2 for m=5, lambda=(1,2,3,5,7) at multiple cardinalities
# w=2,3,4. For each w, construct a separating-deterministic injection phi
# and confirm the MLE variances match w=1 (no masking).
#
# NOT included in the paper — just a sanity check.

set.seed(42)

lambda <- c(1, 2, 3, 5, 7)
m <- length(lambda)
Lambda <- sum(lambda)
probs <- lambda / Lambda

R <- 50000   # replications
n <- 2000    # sample size

cat(sprintf("m=%d, lambda=(%s), Lambda=%.0f, R=%s, n=%d\n\n",
            m, paste(lambda, collapse=","), Lambda,
            format(R, big.mark=","), n))

# Theoretical w=1 variances (scaled by n): Lambda * lambda_j
theo_w1 <- Lambda * lambda
cat(sprintf("Theoretical n*Var (w=1): %s\n\n",
            paste(sprintf("%.1f", theo_w1), collapse=", ")))

# ============================================================================
# Build separating-deterministic injections for each w
# ============================================================================
# We need phi: {1,...,m} -> C(m,w) such that k in phi(k) and phi is injective.
# Strategy: for w < m-1, assign each k a unique w-subset containing k.

build_injection <- function(m, w) {
  # For each k, build phi(k) = a w-subset of {1,...,m} containing k.
  # Ensure all phi(k) are distinct.
  #
  # Simple approach for small m: for k=1,...,m, start with {k}, then add
  # (w-1) elements from {1,...,m}\{k} using a deterministic pattern that
  # ensures distinctness.

  phi <- vector("list", m)

  if (w == 1) {
    for (k in 1:m) phi[[k]] <- k
    return(phi)
  }

  if (w == m - 1) {
    # Derangement: phi(k) = {1,...,m}\{sigma(k)} where sigma is cyclic shift
    sigma <- c(2:m, 1L)
    for (k in 1:m) phi[[k]] <- setdiff(1:m, sigma[k])
    return(phi)
  }

  # General w: use a systematic construction.
  # For each k, include k and the next (w-1) elements cyclically (skipping k).
  others_from <- function(k, m) {
    # Elements {1,...,m}\{k} in cyclic order starting after k
    idx <- ((k:(k + m - 2)) %% m) + 1
    idx[idx != k]  # remove k if present
  }

  for (k in 1:m) {
    pool <- others_from(k, m)
    phi[[k]] <- sort(c(k, pool[1:(w - 1)]))
  }

  # Verify containment and injectivity
  for (k in 1:m) stopifnot(k %in% phi[[k]], length(phi[[k]]) == w)
  phi_strings <- sapply(phi, paste, collapse = ",")
  if (length(unique(phi_strings)) < m) {
    stop(sprintf("Injection failed at w=%d: phi values not all distinct", w))
  }

  return(phi)
}

# ============================================================================
# Simulation for a given w
# ============================================================================

simulate_w <- function(lambda, n, R, w) {
  m <- length(lambda)
  Lambda <- sum(lambda)
  probs <- lambda / Lambda

  phi <- build_injection(m, w)

  # Build lookup: for each candidate set (as string), what is k*?
  phi_strings <- sapply(phi, paste, collapse = ",")
  lookup <- setNames(1:m, phi_strings)

  cat(sprintf("  w=%d: phi mappings:\n", w))
  for (k in 1:m) cat(sprintf("    phi(%d) = {%s}\n", k, phi_strings[k]))

  # Generate data and compute MLEs
  mle_mat <- matrix(0, R, m)

  for (r in 1:R) {
    # Generate n observations of (S, K)
    t_vec <- rexp(n, rate = Lambda)
    K_vec <- sample.int(m, n, replace = TRUE, prob = probs)

    # Apply deterministic masking: C = phi(K)
    # Recover K via phi^{-1}
    recovered_K <- integer(n)
    for (i in 1:n) {
      c_str <- paste(phi[[K_vec[i]]], collapse = ",")
      recovered_K[i] <- lookup[c_str]
    }

    # MLE: lambda_j = count(recovered_K == j) / (n * tbar)
    tbar <- mean(t_vec)
    for (j in 1:m) {
      mle_mat[r, j] <- sum(recovered_K == j) / (n * tbar)
    }
  }

  # Empirical variances scaled by n
  emp_var <- apply(mle_mat, 2, var) * n
  return(emp_var)
}

# ============================================================================
# Run for w = 1, 2, 3, 4
# ============================================================================

cat("=" , rep("=", 59), "\n", sep = "")

for (w in 1:(m - 1)) {
  cat(sprintf("\n--- w = %d ---\n", w))
  emp <- simulate_w(lambda, n, R, w)
  cat(sprintf("  n*Var empirical:    %s\n", paste(sprintf("%7.2f", emp), collapse=", ")))
  cat(sprintf("  n*Var theory (w=1): %s\n", paste(sprintf("%7.2f", theo_w1), collapse=", ")))
  cat(sprintf("  Ratio (emp/theo):   %s\n",
              paste(sprintf("%7.4f", emp / theo_w1), collapse=", ")))
}

cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("If Theorem 3.2 holds, all ratios should be ~1.00.\n")
cat("Done.\n")
