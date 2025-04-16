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
                                 Z, LD, N_sumstat, dict_sumstatLD,
                                 Var_y, SeBhat,
                                 keep_variables = NULL,
                                 focal_outcome_idx = NULL,
                                 focal_outcome_variables = TRUE,
                                 overlap_variables = FALSE,
                                 intercept = TRUE, standardize = TRUE,
                                 residual_correlation = NULL) {
  #################  initialization #######################################
  cb_data <- list("data" = list())
  class(cb_data) <- "colocboost"
  if (!is.null(dict_YX) & !is.null(dict_sumstatLD)) {
    dict <- c(dict_YX, max(dict_YX) + dict_sumstatLD)
  } else if (!is.null(dict_YX) & is.null(dict_sumstatLD)) {
    dict <- dict_YX
  } else if (is.null(dict_YX) & !is.null(dict_sumstatLD)) {
    dict <- dict_sumstatLD
  }
  if (focal_outcome_variables & !is.null(focal_outcome_idx)) {
    if (focal_outcome_idx > length(dict)) {
      stop("Target outcome index is over the total number of outcomes! please check!")
    }
    keep_variable_names <- keep_variables[[dict[focal_outcome_idx]]]
    if (overlap_variables) {
      keep_tmp <- lapply(keep_variables[-dict[focal_outcome_idx]], function(tm) {
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
    drop_lowfreq <- c()
    dict_YX_final <- dict_YX
    for (ij in 1:length(X)) {
      index <- which(dict_YX == ij)
      dict_YX_final[index] <- index[1]
    }
    for (i in 1:length(Y)) {
      tmp <- list(
        "X" = NULL,
        "Y" = scale(Y[[i]], center = intercept, scale = standardize),
        "N" = length(Y[[i]]),
        "variable_miss" = NULL
      )
      x_tmp <- X[[dict_YX[i]]]
      change_x <- if (dict_YX_final[i] == i) TRUE else FALSE
      # - if sample different
      if (nrow(x_tmp) != length(Y[[i]])) {
        change_x <- TRUE
        ind_id_Y <- rownames(Y[[i]])
        ind_id_X <- rownames(x_tmp)
        if (is.null(ind_id_X) | is.null(ind_id_Y)) {
          stop("Please provide the sample index of X and Y, since they do not have the same samples!")
        } else {
          pos <- match(ind_id_Y, ind_id_X)
          x_tmp <- x_tmp[pos, , drop = FALSE]
          if (sum(is.na(pos)) != 0) {
            tmp$Y <- tmp$Y[-which(is.na(pos))]
          }
        }
      }
      # - if missing X
      variable.name <- keep_variables[[dict_YX[i]]]
      if (length(variable.name) != length(keep_variable_names)) {
        x_extend <- matrix(0,
          nrow = nrow(x_tmp), ncol = length(keep_variable_names),
          dimnames = list(rownames(x_tmp), keep_variable_names)
        )
        variable.tmp <- intersect(keep_variable_names, variable.name)
        pos <- match(variable.tmp, keep_variable_names)
        tmp$variable_miss <- setdiff(1:length(keep_variable_names), pos)
        poss <- match(variable.tmp, variable.name)
        x_extend[, pos] <- x_tmp[, poss]
        x_tmp <- x_extend
      }
      if (change_x) {
        dict_YX_final[i] == i
        if (!intercept & !standardize) {
          x_stand <- x_tmp
        } else {
          x_stand <- Rfast::standardise(x_tmp, center = intercept, scale = standardize)
        }
        x_stand[which(is.na(x_stand))] <- 0
        tmp$X <- x_stand
      }
      cb_data$data[[flag]] <- tmp
      names(cb_data$data)[flag] <- paste0("ind_outcome_", i)
      flag <- flag + 1
    }
    cb_data$dict <- c(dict_YX_final)
    ind_idx <- max(dict_YX)
  } else {
    ind_idx <- 0
  }
  n_ind <- flag - 1
  # if summary: XtX XtY, YtY
  if (!is.null(Z) & !is.null(LD)) {
    ####################### need to consider more #########################
    # ------ only code up one sumstat
    variant_lists <- keep_variables[c(flag:length(keep_variables))]
    sumstat_formated <- process_sumstat(Z, N_sumstat, Var_y, SeBhat, LD, 
                                        variant_lists, dict_sumstatLD, 
                                        keep_variable_names)
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
      "weights_path" = c(),
      "profile_loglike_each" = NULL,
      "obj_path" = 999999,
      "obj_single" = 999999,
      "change_loglike" = NULL,
      "correlation" = NULL,
      "z" = NULL,
      "learning_rate_init" = learning_rate_init,
      "stop_thresh" = stop_thresh,
      "ld_jk" = c(),
      "jk" = c()
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
    if (length(multiple_testing_correction) == 1) print(class(multiple_testing_correction))
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
                                 lambda = 0.5, lambda_focal_outcome = 1,
                                 multi_test_thresh = 1,
                                 func_multi_test = "lfdr",
                                 LD_free = FALSE,
                                 outcome_names = NULL,
                                 focal_outcome_idx = NULL) {
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

  cb_model_para <- list(
    "L" = L,
    "P" = P,
    "N" = N,
    "tau" = tau,
    "func_simplex" = func_simplex,
    "lambda" = lambda,
    "lambda_focal_outcome" = lambda_focal_outcome,
    "profile_loglike" = profile_loglike,
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
    "num_updates" = 1
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

inital_residual <- function(Y = NULL, XtY = NULL) {
  if (!is.null(Y)) {
    return(Y)
  } else if (!is.null(XtY)) {
    return(XtY)
  }
}


# - Calculate correlation between X and res
get_correlation <- function(X = NULL, res = NULL, XtY = NULL, N = NULL,
                            YtY = NULL, XtX = NULL, beta_k = NULL, miss_idx = NULL) {
  if (!is.null(X)) {
    corr <- suppressWarnings({
      Rfast::correls(res, X)[, "correlation"]
    })
    corr[which(is.na(corr))] <- 0
    return(corr)
  } else if (!is.null(XtY)) {
    corr <- rep(0, length(XtY))
    scaling_factor <- if (!is.null(N)) (N - 1) else 1
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
    var_r <- YtY - 2 * sum(beta_k * XtY) + sum((XtX %*% as.matrix(beta_k)) * beta_k)
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
      result <- try(
        {
          lfdr_nomissing <- qvalue(pchisq(drop(z^2), 1, lower.tail = FALSE), lambda = seq(0.05, lambda_max, 0.05))$lfdr
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
#' ld_matrices = LD
#' variant_lists = variant_lists
#' dict = dict_sumstatLD
#' target_variants = keep_variable_names
#' N = N_sumstat
process_sumstat <- function(Z, N, Var_y, SeBhat, ld_matrices, variant_lists, dict, target_variants) {
  
  # Step 1: Identify unique combinations of (variant list, LD matrix)
  unified_dict <- integer(length(variant_lists))
  
  # First item is always assigned its own position
  unified_dict[1] <- 1
  
  # Process remaining items
  if (length(variant_lists) > 1){
    for (i in 2:length(variant_lists)) {
      # Check if current combination is duplicate of any previous one
      is_duplicate <- FALSE
      for (j in 1:(i-1)) {
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
    
    # Get corresponding LD matrix from original dictionary mapping
    ld_index <- dict[i]
    current_ld_matrix <- ld_matrices[[ld_index]]
    
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
    
    # Organize data
    if (is.null(current_n)) {
      tmp$XtX <- ld_submatrix
      tmp$XtY <- Z_extend
      tmp$YtY <- 1
    } else {
      if (!is.null(SeBhat[[i]]) & !is.null(Var_y[[i]])) {
        # var_y, shat (and bhat) are provided, so the effects are on the
        # *original scale*.
        adj <- 1 / (Z_extend^2 + current_n - 2)
        if (!is.null(LD_tmp)) {
          XtXdiag <- Var_y[[i]] * adj / (SeBhat[[i]]^2)
          xtx <- t(ld_submatrix * sqrt(XtXdiag)) * sqrt(XtXdiag)
          tmp$XtX <- (xtx + t(xtx)) / 2
        }
        tmp$YtY <- (current_n - 1) * Var_y[[i]]
        tmp$XtY <- Z_extend * sqrt(adj) * Var_y[[i]] / SeBhat[[i]]
      } else {
        if (!is.null(ld_submatrix)) {
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
