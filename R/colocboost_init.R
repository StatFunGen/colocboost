#' @title Set of Internal functions for initial colocboost objects
#'
#' @description
#' The `colocboost_inits` function provides an interface for initializing different colocboost objects, including data, model,
#'      and model parameters. This documentation serves as a summary for all related initialization functions.
#'
#' @usage
#' colocboost_init_data(X, Y, dict_YX, Z, LD, N_sumstat, dict_sumstatLD, Var_y, SeBhat)
#' colocboost_init_model(cb_data)
#' colocboost_init_para(cb_data, cb_model)
#'
#' @details
#' The following functions are included in this set:
#' `colocboost_init_data` initial colocboost data object.
#' `colocboost_init_model` initial colocboost model object.
#' `colocboost_init_para` initial colocboost model global parameters object.
#'
#' These functions are not exported individually and are accessed via `colocboost_inits`.
#'
#' @keywords cb_objects
#' @rdname colocboost_objects
#' @noRd
colocboost_inits <- function() {
  message("This function initializes colocboost objects. See details for more information.")
}


#' @noRd
#' @keywords cb_objects
colocboost_init_data <- function(X, Y, dict_YX,
                                 Z, LD, X_ref = NULL, N_sumstat, dict_sumstatLD,
                                 Var_y, SeBhat,
                                 keep_variables,
                                 focal_outcome_idx = NULL,
                                 focal_outcome_variables = TRUE,
                                 overlap_variables = FALSE,
                                 intercept = TRUE, standardize = TRUE,
                                 residual_correlation = NULL) {
  #################  initialization #######################################
  cb_data <- list("data" = list())
  class(cb_data) <- "colocboost"
  # if (!is.null(dict_YX) & !is.null(dict_sumstatLD)) {
  #   dict <- c(dict_YX, max(dict_YX) + dict_sumstatLD)
  #   n_ind_variable <- max(dict_YX)
  # } else if (!is.null(dict_YX) & is.null(dict_sumstatLD)) {
  #   dict <- dict_YX
  # } else if (is.null(dict_YX) & !is.null(dict_sumstatLD)) {
  #   dict <- dict_sumstatLD
  #   n_ind_variable <- 0
  # }
  
  # note here and remove later: dict is designed by X and LD, but we need additional dictionary for keep variables.
  # keep variables for individual level data is based on X - there is no issue.
  # keep variables for sumstat data is based on sumstat (not LD) - there is a issue to index the focal outcome based on dict later.
  if (!is.null(dict_YX) & !is.null(dict_sumstatLD)) {
    dict <- c(dict_YX, max(dict_YX) + dict_sumstatLD)
    n_ind_variable <- max(dict_YX)
    dict_keep_variables <- c(dict_YX, 1:length(dict_sumstatLD) + n_ind_variable)
  } else if (!is.null(dict_YX) & is.null(dict_sumstatLD)) {
    dict <- dict_YX
    dict_keep_variables <- dict_YX
  } else if (is.null(dict_YX) & !is.null(dict_sumstatLD)) {
    dict <- dict_sumstatLD
    n_ind_variable <- 0
    dict_keep_variables <- 1:length(dict_sumstatLD)
  }
  
  if (focal_outcome_variables & !is.null(focal_outcome_idx)) {
    if (focal_outcome_idx > length(dict)) {
      stop("Target outcome index is over the total number of outcomes! please check!")
    }
    # keep_variable_names <- keep_variables[[dict[focal_outcome_idx]]]
    keep_variable_names <- keep_variables[[dict_keep_variables[focal_outcome_idx]]]
    if (overlap_variables) {
      # keep_tmp <- lapply(keep_variables[-dict[focal_outcome_idx]], function(tm) {
      #   intersect(keep_variable_names, tm)
      # })
      keep_tmp <- lapply(keep_variables[-dict_keep_variables[focal_outcome_idx]], function(tm) {
        intersect(keep_variable_names, tm)
      })
      keep_variable_names <- Reduce(union, keep_tmp)
    }
  } else {
    if (overlap_variables) {
      keep_variable_names <- Reduce(intersect, keep_variables)
    } else {
      keep_variable_names <- get_merge_ordered_with_indices(keep_variables)
    }
  }
  cb_data$variable.names <- keep_variable_names
  flag <- 1
  # if individual: X, Y
  if (!is.null(X) & !is.null(Y)) {
    
    ind_formated <- process_individual_data(
      X, Y, dict_YX, target_variants = keep_variable_names,
      intercept = intercept, standardize = standardize
    )
    for (i in 1:length(Y)) {
      cb_data$data[[flag]] <- ind_formated$result[[i]]
      names(cb_data$data)[flag] <- paste0("ind_outcome_", i)
      flag <- flag + 1
    }
    cb_data$dict <- c(ind_formated$new_dict)
  } 
  n_ind <- flag - 1
  # if summary: XtX XtY, YtY
  if (!is.null(Z) & (!is.null(LD) | !is.null(X_ref))) {
    ####################### need to consider more #########################
    # ------ only code up one sumstat
    variant_lists <- keep_variables[c((n_ind_variable+1):length(keep_variables))]
    sumstat_formated <- process_sumstat(
      Z, N_sumstat, Var_y, SeBhat, LD,
      variant_lists, dict_sumstatLD,
      keep_variable_names,
      X_ref = X_ref
    )
    for (i in 1:length(Z)) {
      cb_data$data[[flag]] <- sumstat_formated$results[[i]]
      names(cb_data$data)[flag] <- paste0("sumstat_outcome_", i)
      flag <- flag + 1
    }
    cb_data$dict <- c(cb_data$dict, sumstat_formated$unified_dict + n_ind)
  }
  # - if residual correlation matrix is not NULL, we need to adjust the study dependence
  if (is.null(residual_correlation)) {
    for (i in 1:length(cb_data$data)) {
      cb_data$data[[i]]$dependency <- 1
    }
  } else {
    pseudo_inverse <- function(residual_correlation) {
      eigen_Sigma <- eigen(residual_correlation)
      L <- which(cumsum(eigen_Sigma$values) / sum(eigen_Sigma$values) > 0.999)[1]
      return(eigen_Sigma$vectors[, 1:L] %*% diag(1 / eigen_Sigma$values[1:L]) %*% t(eigen_Sigma$vectors[, 1:L]))
    }
    Theta <- pseudo_inverse(residual_correlation)
    for (i in 1:length(cb_data$data)) {
      cb_data$data[[i]]$dependency <- Theta[i, i]
    }
  }


  return(cb_data)
}

#' @noRd
#' @keywords cb_objects
colocboost_init_model <- function(cb_data,
                                  learning_rate_init = 0.01,
                                  stop_thresh = 1e-06,
                                  func_multi_test = "lfdr",
                                  stop_null = 0.05,
                                  multi_test_max = 1,
                                  ash_prior = "normal",
                                  p.adjust.methods = "fdr",
                                  focal_outcome_idx = NULL) {
  #################  initialization #######################################
  cb_model <- list()
  class(cb_model) <- "colocboost"

  P <- if (!is.null(cb_data$data[[1]]$X)) ncol(cb_data$data[[1]]$X) else length(cb_data$data[[1]]$XtY)
  L <- length(cb_data$data)
  for (i in 1:length(cb_data$data)) {
    tmp <- list(
      "res" = NULL,
      "beta" = rep(0, P),
      # "weights_path" = c(),
      "weights_path" = list(),
      "profile_loglike_each" = NULL,
      "obj_path" = 999999,
      "obj_single" = 999999,
      "change_loglike" = NULL,
      "correlation" = NULL,
      "z" = NULL,
      "learning_rate_init" = learning_rate_init,
      "stop_thresh" = stop_thresh,
      "ld_jk" = c(),
      "jk" = c(),
      "scaling_factor" = if (!is.null(cb_data$data[[i]]$N)) (cb_data$data[[i]]$N - 1) else 1,
      "beta_scaling" = if (!is.null(cb_data$data[[i]]$N)) 1 else 100
    )

    data_each <- cb_data$data[[i]]
    X_dict <- cb_data$dict[i]
    # - calculate change of loglikelihood for data
    tmp$change_loglike <- estimate_change_profile(
      X = cb_data$data[[X_dict]]$X, Y = data_each[["Y"]], N = data_each$N,
      YtY = data_each$YtY, XtY = data_each$XtY
    )
    # - initial profile loglikelihood
    tmp$profile_loglike_each <- estimate_profile_loglike(Y = data_each[["Y"]], N = data_each$N, YtY = data_each$YtY) * data_each$dependency
    # - initial residual: res for ind; xtr for sumstat
    tmp$res <- inital_residual(Y = data_each[["Y"]], XtY = data_each$XtY)
    # - initial correlation between X and residual
    tmp$correlation <- get_correlation(
      X = cb_data$data[[X_dict]]$X, res = tmp$res, XtY = data_each$XtY,
      N = data_each$N, YtY = data_each$YtY,
      XtX = cb_data$data[[X_dict]]$XtX,
      beta_k = tmp$beta,
      miss_idx = data_each$variable_miss
    )
    # - initial z-score between X and residual based on correlation
    tmp$z <- get_z(tmp$correlation, n = data_each$N, tmp$res)
    tmp$z_univariate <- tmp$z
    tmp$beta_hat_univariate <- get_beta(z = tmp$z, n = data_each$N)
    # - add-hoc stop method
    if (P <= 1) {
      multiple_testing_correction <- 1
    } else {
      multiple_testing_correction <- get_multiple_testing_correction(
        z = tmp$z, miss_idx = data_each$variable_miss,
        func_multi_test = func_multi_test,
        ash_prior = ash_prior,
        p.adjust.methods = p.adjust.methods
      )
    }
    tmp$multi_correction <- multiple_testing_correction
    tmp$multi_correction_univariate <- multiple_testing_correction
    if (all(multiple_testing_correction == 1)) {
      tmp$stop_null <- 1
    } else if (min(multiple_testing_correction) > multi_test_max) {
      tmp$stop_null <- min(multiple_testing_correction)
    } else {
      if (min(multiple_testing_correction) < multi_test_max) {
        tmp$stop_null <- multi_test_max
      } else {
        tmp$stop_null <- min(min(multiple_testing_correction) + 0.1, multi_test_max)
      }
    }

    cb_model[[i]] <- tmp
  }
  names(cb_model) <- names(cb_data$data)
  return(cb_model)
}

#' @noRd
#' @keywords cb_objects
#' @importFrom utils tail
colocboost_init_para <- function(cb_data, cb_model, tau = 0.01,
                                 func_simplex = "z2z",
                                 lambda = 0.5, 
                                 lambda_focal_outcome = 1,
                                 learning_rate_decay = 1,
                                 multi_test_thresh = 1,
                                 multi_test_max = 1,
                                 func_multi_test = "lfdr",
                                 LD_free = FALSE,
                                 outcome_names = NULL,
                                 focal_outcome_idx = NULL,
                                 dynamic_learning_rate = TRUE,
                                 prioritize_jkstar = TRUE,
                                 jk_equiv_corr = 0.8,
                                 jk_equiv_loglik = 1,
                                 func_compare = "min_max",
                                 coloc_thresh =  0.1) {
  #################  initialization #######################################
  # - sample size
  N <- sapply(cb_data$data, function(dt) dt$N)
  # - number of variables
  P <- if (!is.null(cb_data$data[[1]]$X)) ncol(cb_data$data[[1]]$X) else length(cb_data$data[[1]]$XtY)
  # - number of outcomes
  L <- length(cb_data$data)
  # - initial profile loglikelihood
  profile_loglike <- sum(sapply(1:length(cb_model), function(i) tail(cb_model[[i]]$profile_loglike_each, n = 1)))
  # - check initial update outcome
  stop_null <- sapply(cb_model, function(tmp) min(tmp$multi_correction_univariate))
  pos_stop <- which(stop_null >= multi_test_thresh)
  update_y <- rep(1, L)
  if (length(pos_stop) != 0) {
    update_y[pos_stop] <- 0
  } else {
    pos_stop <- NULL
  }

  if (!is.null(outcome_names)) {
    if (length(outcome_names) != L) {
      warning(paste(
        "There are", L, "outcomes inputed, but only", length(outcome_names), "provided.",
        "Using default outcome_names as Y1,...,YL."
      ))
      outcome_names <- paste0("Y", 1:L)
    } else {
      outcome_names <- outcome_names
    }
  } else {
    outcome_names <- paste0("Y", 1:L)
  }
  # - initial coloc_thresh
  coloc_thresh <- (1 - coloc_thresh) * max(sapply(1:length(cb_model), function(i) max(cb_model[[i]]$change_loglike)))
  # - initial num_updates_outcome and coveraged outcome
  num_updates_outcome <- lapply(1:L, function(i) 1)
  coveraged_outcome <- lapply(1:L, function(i) TRUE)
  names(num_updates_outcome) <- names(coveraged_outcome) <- names(cb_model)
  
  cb_model_para <- list(
    "L" = L,
    "P" = P,
    "N" = N,
    "tau" = tau,
    "func_simplex" = func_simplex,
    "lambda" = lambda,
    "lambda_focal_outcome" = lambda_focal_outcome,
    "learning_rate_decay" = learning_rate_decay,
    "profile_loglike" = profile_loglike,
    "dynamic_learning_rate" = dynamic_learning_rate,
    "prioritize_jkstar" = prioritize_jkstar,
    "jk_equiv_corr" = jk_equiv_corr,
    "jk_equiv_loglik" = jk_equiv_loglik,
    "func_compare" = func_compare,
    "coloc_thresh" = coloc_thresh,
    "update_status" = c(),
    "jk" = c(),
    "update_y" = update_y,
    "true_stop" = pos_stop,
    "LD_free" = LD_free,
    "real_update_jk" = c(),
    "outcome_names" = outcome_names,
    "variables" = cb_data$variable.names,
    "focal_outcome_idx" = focal_outcome_idx,
    "coveraged" = TRUE,
    "num_updates" = 0,
    "coveraged_outcome" = coveraged_outcome,
    "num_updates_outcome" = num_updates_outcome,
    "func_multi_test" = func_multi_test,
    "multi_test_thresh" = multi_test_thresh,
    "multi_test_max" = multi_test_max,
    "model_used" = "original"
  )
  class(cb_model_para) <- "colocboost"

  return(cb_model_para)
}


estimate_profile_loglike <- function(Y = NULL, N = NULL, YtY = NULL) {
  if (!is.null(Y)) {
    return(mean(Y^2) * N / (N - 1))
  } else if (!is.null(YtY)) {
    if (is.null(N)) {
      return(YtY)
    } else {
      return(YtY / (N - 1))
    }
  }
}


estimate_change_profile <- function(X = NULL, Y = NULL, N = NULL,
                                    YtY = NULL, XtY = NULL) {
  if (!is.null(X)) {
    yty <- sum(Y^2) / (N - 1)
    xty <- t(X) %*% Y / (N - 1)
  } else if (!is.null(XtY)) {
    if (is.null(N)) {
      yty <- YtY
      xty <- XtY
    } else {
      yty <- YtY / (N - 1)
      xty <- XtY / (N - 1)
    }
  }
  numerator <- xty^2 / (2 * yty)
  denominator <- 0.5 * log(2 * pi * yty) + 0.5
  change_loglike <- numerator / denominator
  change_loglike <- (change_loglike - min(change_loglike)) / (max(change_loglike) - min(change_loglike))
  return(change_loglike)
}

#' Get genotype matrix from a data entry for LD/correlation computation.
#'
#' Each data entry is either individual-level (has \code{$X}) or summary statistics
#' (has \code{$XtX} and/or \code{$X_ref}). These are mutually exclusive:
#' \code{$X} and \code{$X_ref} never coexist on the same entry.
#'
#' This helper is used where a genotype matrix is needed to compute correlations
#' (e.g., purity checks via \code{get_cormat}). For summary stats entries that have
#' \code{$XtX} (precomputed LD), callers typically use \code{Xcorr = $XtX} directly
#' and this function returns NULL — the Xcorr path takes priority in \code{get_purity}.
#' For summary stats entries with \code{$X_ref} (reference panel, no precomputed LD),
#' this returns the reference panel so correlations can be computed on the fly.
#'
#' @param data_entry A single element of \code{cb_data$data}
#' @return A genotype matrix (individual-level X or reference panel X_ref), or NULL
#' @noRd
get_genotype_matrix <- function(data_entry) {
  if (!is.null(data_entry$X)) return(data_entry$X)
  if (!is.null(data_entry$X_ref)) return(data_entry$X_ref)
  return(NULL)
}

#' Compute the product XtX \%*\% v from whichever LD source is available.
#'
#' Dispatches to the appropriate computation:
#' \itemize{
#'   \item \code{XtX} matrix (P x P): direct matrix-vector multiply, O(P^2)
#'   \item \code{X_ref} reference panel (N_ref x P): two matrix-vector multiplies, O(N_ref * P)
#'   \item Neither (LD-free / scalar case): returns \code{v} unchanged (identity)
#' }
#' @param v Numeric vector to multiply
#' @param XtX LD matrix (P x P), or scalar 1 for LD-free, or NULL
#' @param X_ref Standardized reference panel matrix (N_ref x P), or NULL
#' @return Numeric vector of same length as v
#' @noRd
compute_xtx_product <- function(v, XtX = NULL, X_ref = NULL) {
  if (!is.null(XtX) && length(XtX) > 1) {
    as.vector(XtX %*% v)
  } else if (!is.null(X_ref)) {
    as.vector(crossprod(X_ref, X_ref %*% v) / (nrow(X_ref) - 1))
  } else {
    v
  }
}

inital_residual <- function(Y = NULL, XtY = NULL) {
  if (!is.null(Y)) {
    return(Y)
  } else if (!is.null(XtY)) {
    return(XtY)
  }
}


# - Calculate correlation between X and res
get_correlation <- function(X = NULL, res = NULL, XtY = NULL, N = NULL,
                            YtY = NULL, XtX = NULL, X_ref = NULL, beta_k = NULL, miss_idx = NULL,
                            XtX_beta_cache = NULL) {
  if (!is.null(X)) {
    corr <- suppressWarnings({
      correls(res, X)[, "correlation"]
    })
    corr[which(is.na(corr))] <- 0
    return(corr)
  } else if (!is.null(XtY)) {
    corr <- rep(0, length(XtY))
    scaling_factor <- if (!is.null(N)) (N - 1) else 1
    beta_scaling <- if (!is.null(N)) 1 else 100
    beta_k <- beta_k / beta_scaling
    YtY <- YtY / scaling_factor
    XtX <- XtX
    if (length(miss_idx) != 0) {
      XtY <- XtY[-miss_idx] / scaling_factor
      Xtr <- res[-miss_idx] / scaling_factor
      beta_k <- beta_k[-miss_idx]
    } else {
      Xtr <- res / scaling_factor
      XtY <- XtY / scaling_factor
    }
    if (!is.null(XtX_beta_cache)) {
      XtX_beta_k <- XtX_beta_cache
    } else {
      XtX_beta_k <- compute_xtx_product(beta_k, XtX = XtX, X_ref = X_ref)
    }
    var_r <- YtY - 2 * sum(beta_k * XtY) + sum(XtX_beta_k * beta_k)
    if (var_r > 1e-6) {
      corr_nomiss <- Xtr / sqrt(var_r)
      if (length(miss_idx) != 0) {
        corr[-miss_idx] <- corr_nomiss
      } else {
        corr <- corr_nomiss
      }
      if (max(abs(corr)) > 1 & !is.null(N)) {
        corr <- rep(0, length(corr))
      }
    }
    return(corr)
  }
}

# - calculate z-score from correlation
get_z <- function(cor_vals, n = NULL, res = NULL) {
  z <- if (!is.null(n)) cor_vals * sqrt(n - 2) / sqrt(1 - cor_vals^2) else res
  return(z)
}

get_beta <- function(z, n = NULL) {
  if (is.null(n)) {
    return(z)
  } else {
    return(z / sqrt(n - 2 + z^2))
  }
}

get_lfsr <- function(z, miss_idx = NULL, ash_prior = "normal") {
  P <- length(z)
  if (length(miss_idx) != 0) {
    z <- z[-miss_idx]
    lfsr_nomissing <- ashr::ash(drop(z), rep(1, P - length(miss_idx)), mixcompdist = ash_prior)$result$lfsr
    lfsr <- rep(1, P)
    lfsr[-miss_idx] <- lfsr_nomissing
  } else {
    lfsr <- ashr::ash(drop(z), rep(1, P), mixcompdist = ash_prior)$result$lfsr
  }
  return(lfsr)
}

#' @importFrom stats pchisq
get_lfdr <- function(z, miss_idx = NULL) {
  P <- length(z)
  lambda_max <- 0.95
  if (length(miss_idx) != 0) {
    z <- z[-miss_idx]
    try_run <- 1
    while (try_run == 1 && lambda_max >= 0.05) {
      # result <- try(
      #   {
      #     lfdr_nomissing <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), lambda = seq(0.05, lambda_max, 0.05))$lfdr
      #   },
      #   silent = TRUE
      # )
      # if (inherits(result, "try-error")) {
      #   lambda_max <- lambda_max - 0.05 # Decrement lambda_max if error occurs
      # } else {
      #   try_run <- 0
      # }
      result <- tryCatch(
        {
          lfdr_nomissing <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), 
                                   lambda = seq(0.05, lambda_max, 0.05))$lfdr
          list(success = TRUE, value = lfdr_nomissing)
        },
        warning = function(w) {
          list(success = FALSE, message = w$message)
        },
        error = function(e) {
          list(success = FALSE, message = e$message)
        }
      )
      
      if (!result$success) {
        lambda_max <- lambda_max - 0.05
      } else {
        lfdr_nomissing <- result$value
        try_run <- 0
      }
    }
    if (try_run == 1) {
      lfdr_nomissing <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), lambda = 0)$lfdr
    }
    lfdr <- rep(1, P)
    lfdr[-miss_idx] <- lfdr_nomissing
  } else {
    try_run <- 1
    while (try_run == 1 && lambda_max >= 0.05) {
      result <- try(
        {
          lfdr <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), lambda = seq(0.05, lambda_max, 0.05))$lfdr
        },
        silent = TRUE
      )
      if (inherits(result, "try-error")) {
        lambda_max <- lambda_max - 0.05 # Decrement lambda_max if error occurs
      } else {
        try_run <- 0
      }
    }
    if (try_run == 1) {
      lfdr <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), lambda = 0)$lfdr
    }
  }
  return(lfdr)
}

#' @importFrom stats pchisq
get_padj <- function(z, miss_idx = NULL, p.adjust.methods = "fdr") { # test also p.adjust.methods = "BY"
  P <- length(z)
  if (length(miss_idx) != 0) {
    z <- z[-miss_idx]
    pv <- pchisq(drop(z^2), 1, lower.tail = FALSE)
    fdr_nomissing <- stats::p.adjust(pv, method = p.adjust.methods, n = length(pv))
    fdr <- rep(1, P)
    fdr[-miss_idx] <- fdr_nomissing
  } else {
    pv <- pchisq(drop(z^2), 1, lower.tail = FALSE)
    fdr <- stats::p.adjust(pv, method = p.adjust.methods, n = length(pv))
  }
  return(fdr)
}



get_multiple_testing_correction <- function(z, miss_idx = NULL, func_multi_test = "lfdr",
                                            ash_prior = "normal",
                                            p.adjust.methods = "fdr") {
  if (func_multi_test == "lfsr") {
    lfsr <- get_lfsr(z = z, miss_idx = miss_idx, ash_prior = ash_prior)
    return(lfsr)
  } else if (func_multi_test == "lfdr") {
    lfdr <- get_lfdr(z = z, miss_idx = miss_idx)
    return(lfdr)
  } else if (func_multi_test == "padj") {
    fdr <- get_padj(z = z, miss_idx = miss_idx, p.adjust.methods = p.adjust.methods)
    return(fdr)
  } else {
    stop("Invalid option for func_multi_test")
  }
}


#' Process LD matrices and variant lists with optimized storage
#'
#' @param ld_matrices List of LD matrices (e.g., list(PP=PP_matrix, QQ=QQ_matrix))
#' @param variant_lists List of variant lists (e.g., list(P1=P1_variants, P2=P2_variants, ...))
#' @param dict Vector mapping variant lists to LD matrices (e.g., c(1,1,2,2,2))
#' @param target_variants Vector of variants to be considered (M variants)
#'
#' @return List containing processed data with optimized LD submatrix storage
#' @noRd
process_sumstat <- function(Z, N, Var_y, SeBhat, ld_matrices, variant_lists, dict, target_variants, X_ref = NULL) {
  
  
  # Step 1: Identify unique combinations of (variant list, LD matrix)
  unified_dict <- integer(length(variant_lists))

  # First item is always assigned its own position
  unified_dict[1] <- 1
  
  # Process remaining items
  if (length(variant_lists) > 1) {
    for (i in 2:length(variant_lists)) {
      # Check if current combination is duplicate of any previous one
      is_duplicate <- FALSE
      for (j in 1:(i - 1)) {
        if (identical(variant_lists[[i]], variant_lists[[j]]) && dict[i] == dict[j]) {
          unified_dict[i] <- unified_dict[j]
          is_duplicate <- TRUE
          break
        }
      }

      if (!is_duplicate) {
        # If not a duplicate, assign its exact index
        unified_dict[i] <- i
      }
    }
  }

  # Step 2: Process each variant list
  results <- list()

  for (i in 1:length(variant_lists)) {
    tmp <- list(
      "XtX" = NULL,
      "XtY" = NULL,
      "YtY" = NULL,
      "N" = N[[i]],
      "variable_miss" = NULL
    )

    # Get current status
    current_variants <- variant_lists[[i]]
    current_z <- Z[[i]]
    current_n <- N[[i]]

    # Get corresponding LD/X_ref from original dictionary mapping
    ld_index <- dict[i]
    use_xref <- !is.null(X_ref) && is.null(ld_matrices)
    if (use_xref) {
      current_xref <- X_ref[[ld_index]]
      current_ld_matrix <- NULL
    } else {
      current_ld_matrix <- ld_matrices[[ld_index]]
      current_xref <- NULL
    }

    # Find common variants between current list and target variants
    common_variants <- intersect(current_variants, target_variants)

    # Find variants in target but not in current list
    missing_variants <- setdiff(target_variants, current_variants)
    tmp$variable_miss <- which(target_variants %in% missing_variants)

    # - creat extend Z by setting 0 to missing variants
    Z_extend <- rep(0, length(target_variants))
    pos_z <- match(common_variants, current_variants)
    pos_target <- match(common_variants, target_variants)
    Z_extend[pos_target] <- current_z[pos_z]

    # Calculate submatrix for each unique entry (not duplicates)
    ld_submatrix <- NULL
    xref_submatrix <- NULL
    if (use_xref) {
      # X_ref path: extract column submatrix
      if (length(common_variants) > 0 && i == unified_dict[i]) {
        matched_indices <- match(common_variants, colnames(current_xref))
        xref_submatrix <- current_xref[, matched_indices, drop = FALSE]
        colnames(xref_submatrix) <- common_variants
      }
    } else if (length(current_ld_matrix) == 1) {
      ld_submatrix <- current_ld_matrix
    } else {
      if (length(common_variants) > 0) {
        # Only include the submatrix if this entry is unique or is the first occurrence
        if (i == unified_dict[i]) {
          # Check if common_variants and rownames have identical order
          if (identical(common_variants, rownames(current_ld_matrix))) {
            # If order is identical, use the matrix directly without reordering
            ld_submatrix <- current_ld_matrix
          } else {
            # If order is different, reorder using matched indices
            matched_indices <- match(common_variants, rownames(current_ld_matrix))
            ld_submatrix <- current_ld_matrix[matched_indices, matched_indices, drop = FALSE]
            rownames(ld_submatrix) <- common_variants
            colnames(ld_submatrix) <- common_variants
          }
        }
      }
    }

    # Organize data
    if (is.null(current_n)) {
      if (!is.null(xref_submatrix)) {
        tmp$X_ref <- xref_submatrix
      } else {
        tmp$XtX <- ld_submatrix
      }
      tmp$XtY <- Z_extend
      tmp$YtY <- 1
    } else {
      if (!is.null(SeBhat[[i]]) & !is.null(Var_y[[i]])) {
        # var_y, shat (and bhat) are provided, so the effects are on the
        # *original scale*.
        adj <- 1 / (Z_extend^2 + current_n - 2)
        if (!is.null(xref_submatrix)) {
          # X_ref path with scaling: compute LD from X_ref, then scale
          ld_from_xref <- get_cormat(xref_submatrix)
          rownames(ld_from_xref) <- colnames(ld_from_xref) <- common_variants
          XtXdiag <- Var_y[[i]] * adj / (SeBhat[[i]]^2)
          xtx <- t(ld_from_xref * sqrt(XtXdiag)) * sqrt(XtXdiag)
          tmp$XtX <- (xtx + t(xtx)) / 2
        } else if (!is.null(ld_submatrix)) {
          XtXdiag <- Var_y[[i]] * adj / (SeBhat[[i]]^2)
          xtx <- t(ld_submatrix * sqrt(XtXdiag)) * sqrt(XtXdiag)
          tmp$XtX <- (xtx + t(xtx)) / 2
        }
        tmp$YtY <- (current_n - 1) * Var_y[[i]]
        tmp$XtY <- Z_extend * sqrt(adj) * Var_y[[i]] / SeBhat[[i]]
      } else {
        if (!is.null(xref_submatrix)) {
          tmp$X_ref <- xref_submatrix
        } else if (!is.null(ld_submatrix)) {
          tmp$XtX <- ld_submatrix
        }
        tmp$YtY <- (current_n - 1)
        tmp$XtY <- sqrt(current_n - 1) * Z_extend
      }
    }

    # Store results for current list
    results[[i]] <- tmp
  }


  # Return results with the unified dictionary
  return(list(
    results = results,
    unified_dict = unified_dict,
    original_dict = dict
  ))
}

#' process individual level input format
#' @noRd
process_individual_data <- function(X, Y, dict_YX, target_variants,
                                    intercept = TRUE,
                                    standardize = TRUE) {
  
  # Step 0: Check if sample IDs match between Y and corresponding X
  for (i in 1:length(Y)) {
    current_matrix_type <- dict_YX[i]
    # If row counts match, we assume samples are in the same order
    if (nrow(Y[[i]]) != nrow(X[[current_matrix_type]])) {
      # Row counts don't match, so check rownames
      ind_id_Y <- rownames(Y[[i]])
      ind_id_X <- rownames(X[[current_matrix_type]])
      if (is.null(ind_id_X) || is.null(ind_id_Y)) {
        stop("Please provide the sample index of X and Y, since they do not have the same samples!")
      }
      # Find matching samples
      pos <- match(ind_id_Y, ind_id_X)
      if (sum(!is.na(pos)) == 0) {
        stop("No samples in Y match any samples in the corresponding X matrix!")
      }
    }
  }
  
  # Step 1: Update dictionary to handle duplicates samples
  if (length(Y) == 1){
    new_dict = 1
  } else {
    sample_lists <- lapply(Y, rownames)
    new_dict <- 1:length(dict_YX)
    # For each pair of Y indices
    for (i in 1:(length(dict_YX)-1)) {
      for (j in (i+1):length(dict_YX)) {
        # Check if they map to the same X matrix AND have the same samples
        if (dict_YX[i] == dict_YX[j] && identical(sample_lists[[i]], sample_lists[[j]])) {
          # If same matrix and same samples, use the smaller index
          new_dict[j] <- new_dict[i]
        }
      }
    }
  }
  
  # Step 2: Create result list
  result <- list()
  
  for (i in 1:length(Y)) {
    tmp <- list(
      "X" = NULL,
      "Y" = NULL,
      "N" = NULL,
      "variable_miss" = NULL
    )
    matrix_type <- dict_YX[i]
    # Get the appropriate matrix from X list
    current_X <- X[[matrix_type]]
    current_Y <- Y[[i]]
    # Check if we need to match samples or if we can use as-is
    if (nrow(current_Y) == nrow(current_X)) {
      # Same number of rows, assume same order
      matched_X <- current_X
      matched_Y <- scale(current_Y, center = intercept, scale = standardize)
    } else {
      # Different number of rows, find matching samples
      overlap_samples <- intersect(rownames(current_Y), rownames(current_X))
      matched_X <- current_X[match(overlap_samples, rownames(current_X)), , drop = FALSE]
      matched_Y <- current_Y[match(overlap_samples, rownames(current_Y)), , drop = FALSE]
      matched_Y <- scale(matched_Y, center = intercept, scale = standardize)
    }
    tmp$Y <- matched_Y
    tmp$N <- length(matched_Y)
    
    # - if missing X
    variable.name <- colnames(matched_X)
    if (length(variable.name) != length(target_variants)) {
      x_extend <- matrix(0, nrow = nrow(matched_X), ncol = length(target_variants),
                         dimnames = list(rownames(matched_X), target_variants) )
      variable.tmp <- intersect(target_variants, variable.name)
      pos <- match(variable.tmp, target_variants)
      tmp$variable_miss <- setdiff(1:length(target_variants), pos)
      poss <- match(variable.tmp, variable.name)
      x_extend[, pos] <- matched_X[, poss]
      matched_X <- x_extend
    }
    if (new_dict[i] == i) {
      if (!intercept & !standardize) {
        x_stand <- matched_X
      } else {
        x_stand <- standardise(matched_X, center = intercept, scale = standardize)
      }
      x_stand[which(is.na(x_stand))] <- 0
      tmp$X <- x_stand
    }
    
    # Create components for each list
    result[[i]] <- tmp
  }
  # Return results with the unified dictionary
  return(list(
    result = result,
    new_dict = new_dict,
    original_dict = dict_YX
  ))
}
