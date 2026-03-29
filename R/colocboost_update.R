#' Joint boosting algorithm in ColocBoost
#'
#' @details
#' The gradient boosting algorithm for multiple outcomes
#'
#' @importFrom utils head tail
#' @return colocboost object after gradient boosting update
#' @noRd
colocboost_update <- function(cb_model, cb_model_para, cb_data) {
  
  # - clear which outcome need to be updated at which jk
  pos.update <- which(cb_model_para$update_temp$update_status != 0)
  focal_outcome_idx <- cb_model_para$focal_outcome_idx
  tau = cb_model_para$tau

  for (i in pos.update) {
    update_jk <- cb_model_para$update_temp$real_update_jk[i]
    X_dict <- cb_data$dict[i]
    # adj_dependence
    adj_dep <- cb_data$data[[i]]$dependency

    ########## BEGIN: MAIN CALCULATION ###################

    # - calucalate LD between update_jk and other variables
    if (update_jk %in% unlist(cb_model[[i]]$jk)) {
      pos <- which(unlist(cb_model[[i]]$jk) == update_jk)
      ld_jk <- cb_model[[i]]$ld_jk[pos, ]
    } else {
      cb_model[[i]]$jk <- c(cb_model[[i]]$jk, update_jk)
      ld_jk <- get_LD_jk(update_jk,
        X = get_genotype_matrix(cb_data$data[[X_dict]]),
        XtX = cb_data$data[[X_dict]]$XtX,
        N = cb_data$data[[i]]$N,
        remain_idx = setdiff(1:cb_model_para$P, cb_data$data[[i]]$variable_miss),
        P = cb_model_para$P
      )
      cb_model[[i]]$ld_jk <- rbind(cb_model[[i]]$ld_jk, ld_jk)
    }
    ld_feature <- sqrt(abs(ld_jk))

    # - calculate delta
    if (is.null(focal_outcome_idx)) {
      lambda_outcome <- cb_model_para$lambda
    } else {
      lambda_outcome <- ifelse(i == focal_outcome_idx, 
                               cb_model_para$lambda_focal_outcome, 
                               cb_model_para$lambda)
    }
    delta <- boost_KL_delta(
      z = cb_model[[i]]$z,
      ld_feature = ld_feature, adj_dep = adj_dep,
      func_simplex = cb_model_para$func_simplex, 
      lambda = lambda_outcome
    )

    x_tmp <- cb_data$data[[X_dict]]$X
    scaling_factor <- cb_model[[i]]$scaling_factor
    cov_Xtr <- if (!is.null(x_tmp)) {
      t(x_tmp) %*% as.matrix(cb_model[[i]]$res) / scaling_factor
    } else {
      cb_model[[i]]$res / scaling_factor
    }
    obj_ld <- if (cb_model_para$LD_free) ld_feature else rep(1, length(ld_feature))
    if (length(cb_data$data[[i]]$variable_miss) != 0) {
      obj_ld[cb_data$data[[i]]$variable_miss] <- 0
    }
    exp_term <- adj_dep * obj_ld * (abs(cov_Xtr))
    # - calculate individual objective function
    cb_model[[i]]$obj_path <- c(cb_model[[i]]$obj_path, tau * matrixStats::logSumExp(exp_term / tau + log(delta)))
    cb_model[[i]]$obj_single <- c(cb_model[[i]]$obj_single, abs(cov_Xtr[update_jk]))

    exp_term <- exp_term - max(exp_term)
    exp_abs_cor <- delta * exp(exp_term / tau)
    # weights <- adj_dep * ld_feature / scaling_factor * exp_abs_cor / sum(exp_abs_cor)
    weights <- adj_dep * obj_ld * exp_abs_cor / sum(exp_abs_cor)
    weights <- weights / sum(weights)
    # cb_model[[i]]$weights_path <- rbind(cb_model[[i]]$weights_path, as.vector(weights))
    cb_model[[i]]$weights_path <- c(cb_model[[i]]$weights_path, list(as.vector(weights)))
    ########## END: MAIN CALCULATION ###################


    ########## BEGIN: MAIN UPDATE ######################
    # - Gradient ascent on beta
    beta_grad <- weights * sign(cb_model[[i]]$correlation)
    if (cb_model_para$dynamic_learning_rate) {
      if (tail(cb_model[[i]]$obj_path, n = 1) > 0.5) {
        step1 <- max(0.5 * (1 / (1 + cb_model_para$learning_rate_decay * (length(cb_model[[i]]$obj_path) - 1))), cb_model[[i]]$learning_rate_init)
      } else {
        step1 <- cb_model[[i]]$learning_rate_init
      }
    } else {
      step1 <- cb_model[[i]]$learning_rate_init
    }
    cb_model[[i]]$beta <- cb_model[[i]]$beta + step1 * beta_grad

    ## Gradient descent on residuals
    if (!is.null(cb_data$data[[X_dict]]$X)) {
      # - individual level data
      prediction_beta <- cb_data$data[[X_dict]]$X %*% (beta_grad)
      cb_model[[i]]$res <- cb_model[[i]]$res - step1 * prediction_beta
      # - profile-loglikelihood
      x <- cb_data$data[[X_dict]]$X
      y <- cb_data$data[[i]]$Y
      beta <- cb_model[[i]]$beta
      profile_log <- mean((y - x %*% beta)^2) * adj_dep
    } else if (!is.null(cb_data$data[[X_dict]]$XtX) || !is.null(cb_data$data[[X_dict]]$X_ref)) {
      beta_scaling <- cb_model[[i]]$beta_scaling
      # - summary statistics
      xtx <- cb_data$data[[X_dict]]$XtX
      xref <- cb_data$data[[X_dict]]$X_ref
      cb_model[[i]]$res <- rep(0, cb_model_para$P)
      if (length(cb_data$data[[i]]$variable_miss) != 0) {
        beta <- cb_model[[i]]$beta[-cb_data$data[[i]]$variable_miss]  / beta_scaling
        xty <- cb_data$data[[i]]$XtY[-cb_data$data[[i]]$variable_miss]
        XtX_beta <- compute_xtx_product(beta, XtX = xtx, X_ref = xref)
        cb_model[[i]]$res[-cb_data$data[[i]]$variable_miss] <- xty - scaling_factor * XtX_beta
      } else {
        beta <- cb_model[[i]]$beta / beta_scaling
        xty <- cb_data$data[[i]]$XtY
        XtX_beta <- compute_xtx_product(beta, XtX = xtx, X_ref = xref)
        cb_model[[i]]$res <- xty - scaling_factor * XtX_beta
      }
      # - cache XtX %*% beta for reuse in get_correlation (avoids redundant O(P^2) computation)
      cb_model[[i]]$XtX_beta_cache <- XtX_beta
      # - profile-loglikelihood (reuses cached XtX_beta)
      yty <- cb_data$data[[i]]$YtY / scaling_factor
      xty <- xty / scaling_factor
      if (!is.null(xtx) && length(xtx) == 1){
        profile_log <- (yty - 2 * sum(beta * xty) + sum(beta^2)) * adj_dep
      } else {
        profile_log <- (yty - 2 * sum(beta * xty) + sum(XtX_beta * beta)) * adj_dep
      }
    }
    cb_model[[i]]$profile_loglike_each <- c(cb_model[[i]]$profile_loglike_each, profile_log)
  }

  return(cb_model)
}


# - calculate LD for update_jk
get_LD_jk <- function(jk1, X = NULL, XtX = NULL, N = NULL, remain_idx = NULL, P = NULL) {
  if (!is.null(X)) {
    corr <- suppressWarnings({
      correls(X[, jk1], X)[, "correlation"]
    })
    corr[which(is.na(corr))] <- 0
  } else if (!is.null(XtX)) {
    jk1.remain <- which(remain_idx == jk1)
    corr <- rep(0, P)
    if (length(XtX) == 1 | length(jk1.remain)==0){
      corr[remain_idx] <- 1
    } else {
      corr[remain_idx] <- XtX[, jk1.remain]
    }
  }
  corr
}


boost_KL_delta <- function(z, ld_feature, adj_dep,
                           func_simplex = "LD_z2z",
                           lambda = 0.5) {
  # if (!is.null(n)){ z <- z * sqrt( (n-1)/(z^2+n-2) ) }


  if (func_simplex == "Range_Z") {
    z2z <- lambda * 0.5 * z^2 + (1 - lambda) * abs(z)
    z2z <- (z2z - min(z2z)) / (max(z2z) - min(z2z)) * 5
    z2z <- adj_dep * ld_feature * z2z
    delta <- exp(z2z - max(z2z))
  } else if (func_simplex == "LD_z2z") {
    z2z <- lambda * 0.5 * z^2 + (1 - lambda) * abs(z)
    z2z <- adj_dep * ld_feature * z2z
    delta <- exp(z2z - max(z2z))
  } else if (func_simplex == "only_z2z") {
    z2z <- lambda * 0.5 * z^2 + (1 - lambda) * abs(z)
    z2z <- adj_dep * z2z
    delta <- exp(z2z - max(z2z))
  } else if (func_simplex == "entropy") {
    delta <- rep(1, length(z))
  } else {
    stop("please specific the data-driven prior probability style for calculating delta in each iteration!")
  }
  delta <- delta / sum(delta)
  return(delta)
}


boost_check_stop <- function(cb_model, cb_model_para, pos_stop, stop_no_coverage) {
  # - check the iteration for the stop outcome (pos_stop has the same jk with original data)
  iter_each <- sapply(pos_stop, function(i) {
    length(cb_model[[i]]$obj_path) - 1
  })
  lfsr_each <- sapply(pos_stop, function(i) cb_model[[i]]$stop_null < cb_model_para$multi_test_max)
  pos_need_more <- which(iter_each <= 10 & lfsr_each)

  # pos_need_more <- which(iter_each <= 10)
  if (length(pos_need_more) == 0) {
    # --- stop all pos_stop outcomes
    cb_model_para$update_y[pos_stop] <- 0
    cb_model_para$true_stop <- pos_stop
    cb_model_para$no_coverage_stop <- which(stop_no_coverage==TRUE)
    cb_model_para$need_more <- NULL
  } else {
    # keep boosting for outcome with <= 10 iterations
    update_need_more <- pos_stop[pos_need_more]
    cb_model_para$need_more <- update_need_more
    for (i in update_need_more) {
      cb_model[[i]]$stop_thresh <- cb_model[[i]]$stop_thresh * 0.5
      if (cb_model_para$func_multi_test == "Z") {
        cb_model[[i]]$stop_null <- cb_model[[i]]$stop_null - 0.1
      } else if (cb_model_para$func_multi_test == "lfsr" | cb_model_para$func_multi_test == "lfdr") {
        cb_model[[i]]$stop_null <- cb_model[[i]]$stop_null + 0.05
      }
      cb_model[[i]]$learning_rate_init <- cb_model[[i]]$learning_rate_init * 0.5
    }

    # stop update for outcome with > 10 iterations
    if (length(pos_need_more) != length(pos_stop)) {
      pos_true_stop <- pos_stop[-pos_need_more]
      cb_model_para$true_stop <- pos_true_stop
      cb_model_para$update_y[pos_true_stop] <- 0
    } else {
      cb_model_para$true_stop <- NULL
    }
  }
  return(list("cb_model" = cb_model, "cb_model_para" = cb_model_para))
}


boost_obj_last <- function(cb_data, cb_model, cb_model_para) {
  
  pos.stop <- cb_model_para$true_stop
  focal_outcome_idx <- cb_model_para$focal_outcome_idx
  tau <- cb_model_para$tau

  for (i in pos.stop) {
    # - check which jk update
    correlation <- cb_model[[i]]$correlation
    if (sum(correlation) != 0) {
      abs_cor <- abs(correlation)
      jk <- which(abs_cor == max(abs_cor))
      jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
      # adj_dependence
      adj_dep <- cb_data$data[[i]]$dependency

      ########## MAIN CALCULATION ###################
      X_dict <- cb_data$dict[i]
      ld_jk <- get_LD_jk(jk,
        X = get_genotype_matrix(cb_data$data[[X_dict]]),
        XtX = cb_data$data[[X_dict]]$XtX,
        N = cb_data$data[[i]]$N,
        remain_idx = setdiff(1:cb_model_para$P, cb_data$data[[i]]$variable_miss),
        P = cb_model_para$P
      )
      ld_feature <- sqrt(abs(ld_jk))

      # - calculate delta
      if (is.null(focal_outcome_idx)) {
        lambda_outcome <- cb_model_para$lambda
      } else {
        lambda_outcome <- ifelse(i == focal_outcome_idx, 
                                 cb_model_para$lambda_focal_outcome, 
                                 cb_model_para$lambda)
      }
      delta <- boost_KL_delta(
        z = cb_model[[i]]$z,
        ld_feature = ld_feature, adj_dep = adj_dep,
        func_simplex = cb_model_para$func_simplex, 
        lambda = lambda_outcome
      )

      x_tmp <- cb_data$data[[X_dict]]$X
      scaling_factor <- cb_model[[i]]$scaling_factor
      cov_Xtr <- if (!is.null(x_tmp)) {
        t(x_tmp) %*% as.matrix(cb_model[[i]]$res) / scaling_factor
      } else {
        cb_model[[i]]$res / scaling_factor
      }

      obj_ld <- if (cb_model_para$LD_free) ld_feature else rep(1, length(ld_feature))
      if (length(cb_data$data[[i]]$variable_miss) != 0) {
        obj_ld[cb_data$data[[i]]$variable_miss] <- 0
      }
      exp_term <- adj_dep * obj_ld * (abs(cov_Xtr))
      # - calculate individual objective function
      cb_model[[i]]$obj_path <- c(cb_model[[i]]$obj_path, tau * matrixStats::logSumExp(exp_term / tau + log(delta)))
      cb_model[[i]]$obj_single <- c(cb_model[[i]]$obj_single, abs(cov_Xtr[jk]))
    }
  }
  return(cb_model)
}
