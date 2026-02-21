#!/usr/bin/env Rscript
# validate_deterministic_masking.R
#
# Monte Carlo verification of Proposition 4.7: deterministic (derangement)
# masking at w=m-1 recovers the same Fisher information as w=1 (no masking).
#
# Three masking mechanisms are compared on *identical* (S, K) data:
#   1. w=1:            Observe K directly. MLE: lambda_j = count(K=j) / (n*tbar)
#   2. Deterministic:  C = phi(K) via derangement. Since phi is injective,
#                      phi^{-1}(C) recovers K exactly. Same MLE as w=1.
#   3. Uniform w=m-1:  Random size-(m-1) set containing K. MLE via Theorem 5.1.
#
# Prediction: mechanisms (1) and (2) produce identical empirical covariance
# (both match I^{-1}/n = diag(Lambda * lambda_j)/n), while (3) is strictly
# worse (I^{-1}/n = Lambda^2 * [structured matrix] / n).
#
# Self-contained: requires only base R.

set.seed(2026)

# ============================================================================
# Theoretical covariance matrices
# ============================================================================

# w=1 (or deterministic): I^{-1}/n = diag(Lambda * lambda_j) / n
cov_w1 <- function(lambda, n) {
  Lambda <- sum(lambda)
  diag(Lambda * lambda) / n
}

# Uniform w=m-1: I^{-1}/n
# For m=3: Lambda * matrix(Lambda, -l3, -l2; -l3, Lambda, -l1; -l2, -l1, Lambda) / n
# General m: closed-form from Prop 4.6 + matrix inversion
cov_uniform_wm1 <- function(lambda, n) {
  m <- length(lambda)
  Lambda <- sum(lambda)

  # Build FIM using Prop 4.6 for w=m-1
  # Candidate sets: c_j = {1,...,m}\{j}
  # I_{ell,ell'} = [1 / (binom(m-1,m-2) * Lambda)] *
  #   sum_{c: ell in c, ell' in c} (sum_{p in c} lambda_p)^{-1}
  # binom(m-1, m-2) = m-1
  # c_j contains ell iff j != ell
  # sum_{p in c_j} lambda_p = Lambda - lambda_j

  I_mat <- matrix(0, m, m)
  for (ell in 1:m) {
    for (ellp in 1:m) {
      s <- 0
      for (j in 1:m) {
        if (j != ell && j != ellp) {
          s <- s + 1 / (Lambda - lambda[j])
        }
      }
      I_mat[ell, ellp] <- s / ((m - 1) * Lambda)
    }
  }
  solve(I_mat) / n
}

# ============================================================================
# Derangement construction
# ============================================================================

# Simple cyclic derangement: sigma(k) = (k mod m) + 1
# i.e., 1->2, 2->3, ..., m->1
make_derangement <- function(m) {
  sigma <- c(2:m, 1L)
  # Verify: no fixed points
  stopifnot(all(sigma != 1:m))
  sigma
}

# ============================================================================
# Main simulation
# ============================================================================

BATCH_SIZE <- 10000L

run_simulation <- function(lambda, sample_sizes, R) {
  m <- length(lambda)
  Lambda <- sum(lambda)
  probs <- lambda / Lambda

  sigma <- make_derangement(m)
  # phi(k) = {1,...,m} \ {sigma(k)}
  # phi_inv: candidate set c_j (excluding j) -> phi^{-1}(c_j)
  # If c excludes j, then sigma(K) = j, so K = sigma^{-1}(j)
  sigma_inv <- integer(m)
  sigma_inv[sigma] <- 1:m  # sigma_inv[sigma[k]] = k

  cat(sprintf("m=%d, lambda=(%s), Lambda=%.1f\n",
              m, paste(lambda, collapse=", "), Lambda))
  cat(sprintf("Derangement sigma: (%s)\n", paste(sigma, collapse=", ")))
  cat(sprintf("R=%s per sample size\n\n", format(R, big.mark = ",")))

  results <- list()

  for (si in seq_along(sample_sizes)) {
    n <- sample_sizes[si]
    cat(sprintf("  n = %d ... ", n))

    # Accumulate MLEs for all three mechanisms
    mle_w1_chunks <- list()
    mle_det_chunks <- list()
    mle_unif_chunks <- list()
    chunk_idx <- 0L
    R_remaining <- R
    n_rejected_unif <- 0L

    while (R_remaining > 0L) {
      B <- min(BATCH_SIZE, R_remaining)
      R_remaining <- R_remaining - B
      chunk_idx <- chunk_idx + 1L

      # --- Generate (S, K) data: B samples of n observations ---
      # System lifetimes
      t_mat <- matrix(rexp(B * n, rate = Lambda), nrow = B, ncol = n)
      t_bars <- rowMeans(t_mat)

      # Failed components
      K_mat <- matrix(sample.int(m, B * n, replace = TRUE, prob = probs),
                      nrow = B, ncol = n)

      # --- Mechanism 1: w=1 (no masking) ---
      # Count how many times each component failed in each sample
      counts_w1 <- matrix(0L, B, m)
      for (j in 1:m) {
        counts_w1[, j] <- rowSums(K_mat == j)
      }
      mle_w1_chunks[[chunk_idx]] <- counts_w1 / (n * t_bars)

      # --- Mechanism 2: Deterministic masking ---
      # C = {1,...,m}\{sigma(K)}, so the excluded component is sigma(K)
      # To recover K: K = sigma_inv[excluded component]
      # But excluded = sigma(K), so sigma_inv[sigma(K)] = K.
      # The "identified" component is sigma_inv applied to the excluded index.
      # Since this is a bijection, count(C=phi(j)) = count(K=j).
      # So the MLE is identical to w=1!
      #
      # But let's compute it explicitly through the candidate set to verify:
      excluded_mat <- matrix(sigma[K_mat], nrow = B, ncol = n)
      recovered_K <- matrix(sigma_inv[excluded_mat], nrow = B, ncol = n)

      counts_det <- matrix(0L, B, m)
      for (j in 1:m) {
        counts_det[, j] <- rowSums(recovered_K == j)
      }
      mle_det_chunks[[chunk_idx]] <- counts_det / (n * t_bars)

      # --- Mechanism 3: Uniform masking at w=m-1 ---
      # For each observation, candidate set is a random (m-1)-subset containing K
      # Equivalently, exclude one of the m-1 non-failed components uniformly
      non_K <- matrix(0L, B * n, m - 1)
      K_vec <- as.vector(K_mat)
      for (i in seq_along(K_vec)) {
        non_K[i, ] <- (1:m)[-K_vec[i]]
      }
      # Pick which non-failed component to exclude
      excl_idx <- sample.int(m - 1, B * n, replace = TRUE)
      excluded_unif <- non_K[cbind(seq_along(excl_idx), excl_idx)]
      excluded_unif_mat <- matrix(excluded_unif, nrow = B, ncol = n)

      # B_j = count of candidate sets excluding component j
      B_counts <- matrix(0L, B, m)
      for (j in 1:m) {
        B_counts[, j] <- rowSums(excluded_unif_mat == j)
      }

      # Theorem 5.1: lambda_j = (A_j - (m-2)*B_j) / (n*tbar)
      # A_j = n - B_j (every obs either contains j or excludes it)
      A_counts <- n - B_counts
      mle_unif_raw <- (A_counts - (m - 2) * B_counts) / (n * t_bars)

      # Filter out non-identifiable samples (any lambda_j <= 0)
      valid <- apply(mle_unif_raw, 1, function(x) all(x > 0))
      n_rejected_unif <- n_rejected_unif + sum(!valid)

      mle_unif_chunks[[chunk_idx]] <- mle_unif_raw[valid, , drop = FALSE]
    }

    mle_w1_mat <- do.call(rbind, mle_w1_chunks)
    mle_det_mat <- do.call(rbind, mle_det_chunks)
    mle_unif_mat <- do.call(rbind, mle_unif_chunks)

    cov_emp_w1 <- cov(mle_w1_mat)
    cov_emp_det <- cov(mle_det_mat)
    cov_emp_unif <- cov(mle_unif_mat)

    cov_theo_w1 <- cov_w1(lambda, n)
    cov_theo_unif <- cov_uniform_wm1(lambda, n)

    # Relative Frobenius errors
    frob <- function(A) sqrt(sum(A^2))
    err_w1 <- frob(cov_emp_w1 - cov_theo_w1) / frob(cov_theo_w1)
    err_det <- frob(cov_emp_det - cov_theo_w1) / frob(cov_theo_w1)
    err_unif <- frob(cov_emp_unif - cov_theo_unif) / frob(cov_theo_unif)

    # Also: how close is deterministic to w=1 empirically?
    err_det_vs_w1 <- frob(cov_emp_det - cov_emp_w1) / frob(cov_emp_w1)

    reject_pct <- n_rejected_unif / R * 100

    cat(sprintf("w1=%.3f%%, det=%.3f%%, unif=%.3f%%, det_vs_w1=%.4f%%, reject=%.1f%%\n",
                err_w1 * 100, err_det * 100, err_unif * 100,
                err_det_vs_w1 * 100, reject_pct))

    results[[si]] <- list(
      n = n,
      var_w1 = diag(cov_emp_w1),
      var_det = diag(cov_emp_det),
      var_unif = diag(cov_emp_unif),
      var_theo_w1 = diag(cov_theo_w1),
      var_theo_unif = diag(cov_theo_unif),
      err_w1 = err_w1,
      err_det = err_det,
      err_unif = err_unif,
      err_det_vs_w1 = err_det_vs_w1,
      reject_pct = reject_pct
    )
  }

  results
}

# ============================================================================
# Run simulations
# ============================================================================

cat("=" , rep("=", 69), "\n", sep = "")
cat("Proposition 4.7 Verification: Deterministic Masking Recovers w=1 FIM\n")
cat("=", rep("=", 69), "\n\n", sep = "")

sample_sizes <- c(50, 100, 200, 500, 1000, 2000)
R <- 100000

# --- Configuration 1: m=3 asymmetric ---
cat("--- Configuration 1: m=3, lambda=(1,3,5) ---\n")
res_m3 <- run_simulation(c(1, 3, 5), sample_sizes, R)

# --- Configuration 2: m=5 ---
cat("\n--- Configuration 2: m=5, lambda=(1,2,3,5,7) ---\n")
res_m5 <- run_simulation(c(1, 2, 3, 5, 7), sample_sizes, R)

# ============================================================================
# Summary tables
# ============================================================================

print_summary <- function(results, lambda) {
  m <- length(lambda)
  Lambda <- sum(lambda)

  cat("\n  Empirical variances (×n) at largest n:\n")
  r <- results[[length(results)]]
  n <- r$n
  cat(sprintf("  %-14s", "Component:"))
  for (j in 1:m) cat(sprintf("  j=%d    ", j))
  cat("\n")

  cat(sprintf("  %-14s", "Theo w=1:"))
  for (j in 1:m) cat(sprintf("  %7.2f", r$var_theo_w1[j] * n))
  cat("\n")

  cat(sprintf("  %-14s", "Emp w=1:"))
  for (j in 1:m) cat(sprintf("  %7.2f", r$var_w1[j] * n))
  cat("\n")

  cat(sprintf("  %-14s", "Emp derange:"))
  for (j in 1:m) cat(sprintf("  %7.2f", r$var_det[j] * n))
  cat("\n")

  cat(sprintf("  %-14s", "Theo uniform:"))
  for (j in 1:m) cat(sprintf("  %7.2f", r$var_theo_unif[j] * n))
  cat("\n")

  cat(sprintf("  %-14s", "Emp uniform:"))
  for (j in 1:m) cat(sprintf("  %7.2f", r$var_unif[j] * n))
  cat("\n")

  cat("\n  Convergence (relative Frobenius error vs theory):\n")
  cat(sprintf("  %6s  %10s  %10s  %10s  %12s  %8s\n",
              "n", "w=1", "derange", "uniform", "det vs w1", "reject%"))
  for (r in results) {
    cat(sprintf("  %6d  %9.3f%%  %9.3f%%  %9.3f%%  %11.4f%%  %7.1f%%\n",
                r$n, r$err_w1 * 100, r$err_det * 100, r$err_unif * 100,
                r$err_det_vs_w1 * 100, r$reject_pct))
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY: m=3, lambda=(1,3,5)\n")
cat(rep("=", 70), "\n", sep = "")
print_summary(res_m3, c(1, 3, 5))

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY: m=5, lambda=(1,2,3,5,7)\n")
cat(rep("=", 70), "\n", sep = "")
print_summary(res_m5, c(1, 2, 3, 5, 7))

# ============================================================================
# Generate figure
# ============================================================================

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
figpath <- file.path(script_dir, "..", "paper", "fig_deterministic_masking.pdf")
if (!dir.exists(dirname(figpath))) figpath <- "fig_deterministic_masking.pdf"

cat(sprintf("\nWriting figure to %s\n", normalizePath(figpath, mustWork = FALSE)))

pdf(figpath, width = 8, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 5, 2.5, 1), cex.lab = 1.05, cex.axis = 0.85)

# --- Panel (a): m=3 variance of lambda_1 (the rarest component) ---
lambda_m3 <- c(1, 3, 5)
Lambda_m3 <- sum(lambda_m3)
j_plot <- 1  # component 1 (rarest, most interesting)

ns <- sapply(res_m3, `[[`, "n")
var_w1 <- sapply(res_m3, function(r) r$var_w1[j_plot])
var_det <- sapply(res_m3, function(r) r$var_det[j_plot])
var_unif <- sapply(res_m3, function(r) r$var_unif[j_plot])
theo_w1 <- sapply(res_m3, function(r) r$var_theo_w1[j_plot])
theo_unif <- sapply(res_m3, function(r) r$var_theo_unif[j_plot])

ylim <- range(c(var_w1, var_det, var_unif, theo_w1, theo_unif)) * c(0.7, 1.4)

plot(ns, var_w1, type = "b", pch = 16, col = "#2166AC", lwd = 2,
     log = "xy", xlab = "Sample size (n)",
     ylab = expression(Var(hat(lambda)[1])),
     ylim = ylim, main = expression("(a)" ~ m == 3 ~ ", " ~ lambda == "(1,3,5)"),
     axes = FALSE)
axis(1, at = ns, labels = c("50", "100", "200", "500", "1K", "2K"))
axis(2, las = 1)
box()

lines(ns, var_det, type = "b", pch = 17, col = "#E66101", lwd = 2, lty = 2)
lines(ns, var_unif, type = "b", pch = 15, col = "#B2182B", lwd = 2)
lines(ns, theo_w1, col = "#2166AC", lwd = 1, lty = 3)
lines(ns, theo_unif, col = "#B2182B", lwd = 1, lty = 3)

legend("topright",
       legend = c("w=1 (no masking)", "Deterministic w=m-1",
                  "Uniform w=m-1", "Theory"),
       col = c("#2166AC", "#E66101", "#B2182B", "gray40"),
       pch = c(16, 17, 15, NA), lwd = c(2, 2, 2, 1),
       lty = c(1, 2, 1, 3), cex = 0.7, bg = "white")

# --- Panel (b): m=5 variance of lambda_1 ---
lambda_m5 <- c(1, 2, 3, 5, 7)
Lambda_m5 <- sum(lambda_m5)

var_w1_5 <- sapply(res_m5, function(r) r$var_w1[j_plot])
var_det_5 <- sapply(res_m5, function(r) r$var_det[j_plot])
var_unif_5 <- sapply(res_m5, function(r) r$var_unif[j_plot])
theo_w1_5 <- sapply(res_m5, function(r) r$var_theo_w1[j_plot])
theo_unif_5 <- sapply(res_m5, function(r) r$var_theo_unif[j_plot])

ylim5 <- range(c(var_w1_5, var_det_5, var_unif_5, theo_w1_5, theo_unif_5)) * c(0.7, 1.4)

plot(ns, var_w1_5, type = "b", pch = 16, col = "#2166AC", lwd = 2,
     log = "xy", xlab = "Sample size (n)",
     ylab = expression(Var(hat(lambda)[1])),
     ylim = ylim5, main = expression("(b)" ~ m == 5 ~ ", " ~ lambda == "(1,2,3,5,7)"),
     axes = FALSE)
axis(1, at = ns, labels = c("50", "100", "200", "500", "1K", "2K"))
axis(2, las = 1)
box()

lines(ns, var_det_5, type = "b", pch = 17, col = "#E66101", lwd = 2, lty = 2)
lines(ns, var_unif_5, type = "b", pch = 15, col = "#B2182B", lwd = 2)
lines(ns, theo_w1_5, col = "#2166AC", lwd = 1, lty = 3)
lines(ns, theo_unif_5, col = "#B2182B", lwd = 1, lty = 3)

legend("topright",
       legend = c("w=1 (no masking)", "Deterministic w=m-1",
                  "Uniform w=m-1", "Theory"),
       col = c("#2166AC", "#E66101", "#B2182B", "gray40"),
       pch = c(16, 17, 15, NA), lwd = c(2, 2, 2, 1),
       lty = c(1, 2, 1, 3), cex = 0.7, bg = "white")

invisible(dev.off())

cat("\n", rep("=", 70), "\n", sep = "")
cat("CONCLUSION\n")
cat(rep("=", 70), "\n", sep = "")
cat("If Prop 4.7 is correct:\n")
cat("  - 'w=1' and 'derange' columns should match (same theory curve)\n")
cat("  - 'det vs w1' should be near zero (sampling noise only)\n")
cat("  - 'uniform' should converge to a DIFFERENT (larger) theory curve\n")
cat("Done.\n")
