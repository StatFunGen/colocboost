library(testthat)

# ---- Shared test data ----

generate_xref_test_data <- function(n = 200, p = 30, L = 2, n_ref = 15, seed = 42) {
  set.seed(seed)
  sigma <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) sigma[i, j] <- 0.9^abs(i - j)
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  rownames(X) <- paste0("sample", 1:n)

  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5
  true_beta[5, 2] <- 0.4

  Y <- matrix(0, n, L)
  for (l in 1:L) Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)

  LD <- cor(X)

  # Reference panel (subset of samples or independent)
  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  colnames(X_ref) <- paste0("SNP", 1:p)

  list(X = X, Y = Y, LD = LD, X_ref = X_ref, true_beta = true_beta)
}

make_sumstat <- function(X, Y) {
  p <- ncol(X); L <- ncol(Y)
  ss <- list()
  for (i in 1:L) {
    z <- se <- b <- rep(0, p)
    for (j in 1:p) {
      fit <- summary(lm(Y[, i] ~ X[, j]))$coef
      if (nrow(fit) == 2) { b[j] <- fit[2, 1]; se[j] <- fit[2, 2]; z[j] <- b[j] / se[j] }
    }
    ss[[i]] <- data.frame(beta = b, sebeta = se, z = z, n = nrow(X), variant = colnames(X))
  }
  ss
}

# ---- Test: X_ref with N_ref < P produces valid results ----

test_that("X_ref with N_ref < P produces valid colocboost results", {
  # Use p=20 with n_ref=15 for a realistic but not too extreme ratio
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 15)
  ss <- make_sumstat(td$X, td$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = td$X_ref, M = 50, output_level = 1)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
  expect_equal(length(result$data_info$variables), 20)
})

# ---- Test: X_ref with N_ref >= P auto-converts to LD ----

test_that("X_ref with N_ref >= P auto-converts to LD path", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 25)
  ss <- make_sumstat(td$X, td$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = td$X_ref, M = 20, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: X_ref and LD both provided -> error ----

test_that("providing both LD and X_ref returns error", {
  td <- generate_xref_test_data()
  ss <- make_sumstat(td$X, td$Y)

  result <- suppressWarnings(colocboost(sumstat = ss, LD = td$LD, X_ref = td$X_ref, M = 10))
  expect_null(result)
})

# ---- Test: Neither LD nor X_ref -> LD-free mode ----

test_that("neither LD nor X_ref triggers LD-free mode", {
  td <- generate_xref_test_data()
  ss <- make_sumstat(td$X, td$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, M = 1, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
})

# ---- Test: X_ref numerical equivalence with LD (N_ref >= P case) ----

test_that("X_ref produces identical results to precomputed LD when N_ref >= P", {
  set.seed(123)
  n <- 200; p <- 20
  sigma <- diag(p)
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n * 2), n, 2)
  Y[, 1] <- Y[, 1] + X[, 5] * 0.5
  Y[, 2] <- Y[, 2] + X[, 5] * 0.4

  # Use X itself as reference (N_ref = 200 >= P = 20)
  X_ref <- X
  LD <- get_cormat(X)
  rownames(LD) <- colnames(LD) <- colnames(X)

  ss <- make_sumstat(X, Y)

  suppressWarnings(suppressMessages({
    res_ld <- colocboost(sumstat = ss, LD = LD, M = 50, output_level = 2)
    res_xref <- colocboost(sumstat = ss, X_ref = X_ref, M = 50, output_level = 2)
  }))

  # Both should have same structure
  expect_equal(res_ld$data_info$n_outcomes, res_xref$data_info$n_outcomes)
  expect_equal(length(res_ld$data_info$variables), length(res_xref$data_info$variables))

  # Model info should be identical (since N_ref >= P auto-converts to same LD)
  expect_equal(res_ld$model_info, res_xref$model_info)
})

# ---- Test: X_ref with list of matrices ----

test_that("X_ref works as list of matrices", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)
  X_ref_list <- list(td$X_ref, td$X_ref)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = X_ref_list, M = 20, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: X_ref with dictionary mapping ----

test_that("X_ref works with dict_sumstatLD", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)
  X_ref_list <- list(td$X_ref, td$X_ref)
  dict <- cbind(1:2, c(1, 2))

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = X_ref_list, dict_sumstatLD = dict, M = 20, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
})

# ---- Test: X_ref with focal outcome ----

test_that("X_ref works with focal outcome", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = td$X_ref, focal_outcome_idx = 1, M = 20, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$outcome_info$is_focal[1], TRUE)
})

# ---- Test: X_ref with missing variants ----

test_that("X_ref handles partially overlapping variants", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss1 <- make_sumstat(td$X, td$Y[, 1, drop = FALSE])[[1]]
  ss2 <- make_sumstat(td$X, td$Y[, 2, drop = FALSE])[[1]]
  ss2 <- ss2[1:15, ]  # Remove 5 variants from second sumstat

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = list(ss1, ss2), X_ref = td$X_ref, M = 20, output_level = 1)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 2)
})

# ---- Test: X_ref with 3 outcomes ----

test_that("X_ref works with 3 outcomes", {
  set.seed(99)
  n <- 150; p <- 20; n_ref <- 10
  sigma <- diag(p)
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n * 3), n, 3)
  Y[, 1] <- Y[, 1] + X[, 5] * 0.5
  Y[, 2] <- Y[, 2] + X[, 5] * 0.4
  Y[, 3] <- Y[, 3] + X[, 10] * 0.3

  X_ref <- MASS::mvrnorm(n_ref, rep(0, p), sigma)
  colnames(X_ref) <- paste0("SNP", 1:p)

  ss <- make_sumstat(X, Y)

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = X_ref, M = 20, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
  expect_equal(result$data_info$n_outcomes, 3)
})

# ---- Test: X_ref without sample size ----

test_that("X_ref works without sample size in sumstat", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)
  for (i in seq_along(ss)) ss[[i]]$n <- NULL

  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = td$X_ref, M = 10, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
})

# ---- Test: X_ref with HyPrColoc format ----

test_that("X_ref works with HyPrColoc format", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  X <- td$X; Y <- td$Y

  beta <- se <- matrix(0, ncol(X), ncol(Y))
  for (i in 1:ncol(Y)) for (j in 1:ncol(X)) {
    fit <- summary(lm(Y[, i] ~ X[, j]))$coef
    if (nrow(fit) == 2) { beta[j, i] <- fit[2, 1]; se[j, i] <- fit[2, 2] }
  }
  rownames(beta) <- rownames(se) <- colnames(X)
  colnames(beta) <- colnames(se) <- paste0("Y", 1:ncol(Y))

  suppressWarnings(suppressMessages({
    result <- colocboost(
      effect_est = beta, effect_se = se, effect_n = rep(nrow(X), ncol(Y)),
      X_ref = td$X_ref, M = 20, output_level = 2
    )
  }))

  expect_s3_class(result, "colocboost")
})

# ---- Test: Single matrix X_ref (not in list) ----

test_that("X_ref accepts single matrix (auto-wraps in list)", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)

  # Pass as matrix directly, not as list
  suppressWarnings(suppressMessages({
    result <- colocboost(sumstat = ss, X_ref = td$X_ref, M = 10, output_level = 2)
  }))

  expect_s3_class(result, "colocboost")
})

# ---- Test: X_ref colnames validation ----

test_that("X_ref without colnames returns error", {
  td <- generate_xref_test_data(n = 200, p = 20, n_ref = 10)
  ss <- make_sumstat(td$X, td$Y)
  X_ref_nonames <- td$X_ref
  colnames(X_ref_nonames) <- NULL

  result <- suppressWarnings(colocboost(sumstat = ss, X_ref = X_ref_nonames, M = 10))
  expect_null(result)
})
