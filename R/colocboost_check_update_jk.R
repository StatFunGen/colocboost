#' @title Model selection scheme in ColocBoost
#'
#' @details
#' Model selection scheme in proximity smoothed gradient boosting algorithm to select update traits and update variants at each iteration.
#'
#' @return update_status and real_update_jk for each trait
#' @noRd
colocboost_check_update_jk <- function(cb_model, cb_model_para, cb_data,
                                       prioritize_jkstar = TRUE,
                                       jk_equiv_corr = 0.8, ##### more than 2 traits
                                       jk_equiv_loglik = 1, ## more than 2 traits
                                       func_compare = "min_max", ##### more than 3 traits
                                       coloc_thresh = 0.1) {
  pos.update <- which(cb_model_para$update_y == 1)
  focal_outcome_idx <- cb_model_para$focal_outcome_idx
  if (is.null(focal_outcome_idx)) {
    cb_model_para <- boost_check_update_jk_nofocal(cb_model, cb_model_para, cb_data,
      prioritize_jkstar = prioritize_jkstar,
      jk_equiv_corr = jk_equiv_corr,
      jk_equiv_loglik = jk_equiv_loglik,
      func_compare = func_compare,
      coloc_thresh = coloc_thresh
    )
  } else {
    if (focal_outcome_idx %in% pos.update) {
      cb_model_para <- boost_check_update_jk_focal(cb_model, cb_model_para, cb_data,
        prioritize_jkstar = prioritize_jkstar,
        jk_equiv_corr = jk_equiv_corr,
        jk_equiv_loglik = jk_equiv_loglik,
        func_compare = func_compare,
        coloc_thresh = coloc_thresh,
        focal_outcome_idx = focal_outcome_idx
      )
    } else {
      cb_model_para <- boost_check_update_jk_nofocal(cb_model, cb_model_para, cb_data,
        prioritize_jkstar = prioritize_jkstar,
        jk_equiv_corr = jk_equiv_corr,
        jk_equiv_loglik = jk_equiv_loglik,
        func_compare = func_compare,
        coloc_thresh = coloc_thresh
      )
    }
  }
  return(cb_model_para)
}



#' @importFrom stats median
boost_check_update_jk_nofocal <- function(cb_model, cb_model_para, cb_data,
                                           prioritize_jkstar = TRUE,
                                           jk_equiv_corr = 0.8, ##### more than 2 traits
                                           jk_equiv_loglik = 1, ## more than 2 traits
                                           func_compare = "min_max", ##### more than 3 traits
                                           coloc_thresh = 0.1) {
  ############# Output #################
  # we will obtain and update the cb_model_para$update_status and cb_model_para$real_update_jk
  ######################################

  # - initial update_status and jk
  update_status <- rep(0, cb_model_para$L)
  update_jk <- rep(NA, cb_model_para$L + 1)
  real_update_jk <- rep(NA, cb_model_para$L)
  if (is.null(cb_model_para$coloc_thresh)) {
    cb_model_para$coloc_thresh <- (1 - coloc_thresh) * max(sapply(1:length(cb_model), function(i) max(cb_model[[i]]$change_loglike)))
  }

  # - update only Ys which is not stop
  pos.update <- which(cb_model_para$update_y == 1)

  if (length(pos.update) == 1) {
    # -- update status
    update_status[pos.update] <- -1
    correlation <- abs(cb_model[[pos.update]]$correlation)
    jk <- which(correlation == max(correlation))
    jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
    update_jk[c(1, pos.update + 1)] <- c(jk, jk)
    real_update_jk[pos.update] <- jk
  } else {
    # -- extract data and model based on pos.update
    model_update <- cb_model[pos.update]
    # data_update = cb_data$data[pos.update]
    X_dict <- cb_data$dict[pos.update]
    # adj_dependence
    adj_dep <- sapply(pos.update, function(ii) cb_data$data[[ii]]$dependency)

    # -- define jk and jk_each
    cor_vals_each <- do.call(cbind, lapply(model_update, function(cc) cc$correlation))
    abs_cor_vals_each <- sweep(abs(cor_vals_each), 2, adj_dep, `*`)
    abs_cor_vals <- rowSums(abs_cor_vals_each)

    jk <- which(abs_cor_vals == max(abs_cor_vals))
    jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
    jk_each <- apply(abs_cor_vals_each, 2, function(temp) {
      jk_temp <- which(temp == max(temp))
      return(ifelse(length(jk_temp) == 1, jk_temp, sample(jk_temp, 1)))
    })
    update_jk[c(1, pos.update + 1)] <- c(jk, jk_each)
    judge_each <- check_jk_jkeach(jk, jk_each,
      pos.update,
      model_update = model_update,
      cb_data = cb_data,
      X_dict = X_dict,
      jk_equiv_corr = jk_equiv_corr,
      jk_equiv_loglik = jk_equiv_loglik
    )

    if (all(judge_each)) {
      # -- colocalization region
      # -- update all Y at jk
      update_status[pos.update] <- 1
      if (prioritize_jkstar) {
        real_update_jk[pos.update] <- rep(jk, length(pos.update))
      } else {
        real_update_jk[pos.update] <- jk_each
      }
    } else if (all(!judge_each)) {
      # -- uncolocalization region
      # -- update all Y based on jk_each, respectively
      # -- Note: do not update the pair of Y_i and Y_j if the jk_each[i] ~ jk_each[j]

      # ----- relative change of profile loglikelihood of each pair of Y
      change_each_pair <- check_pair_jkeach(jk_each, pos.update,
        model_update = model_update,
        cb_data = cb_data,
        X_dict = X_dict,
        jk_equiv_corr = jk_equiv_corr,
        jk_equiv_loglik = jk_equiv_loglik
      )
      pos <- which(colSums(change_each_pair) != 0)
      if (length(pos) == length(pos.update)) {
        # -- we only need to update paired Y_i at jk_i with higher change of loglikelihood
        data_update <- cb_data$data[pos.update]
        change_res_each <- lapply(1:length(jk_each), function(ii) {
          estimate_change_profile_res(jk_each[ii],
            X = cb_data$data[[X_dict[ii]]]$X,
            res = model_update[[ii]]$res,
            N = data_update[[ii]]$N,
            XtX = cb_data$data[[X_dict[ii]]]$XtX,
            YtY = data_update[[ii]]$YtY,
            XtY = data_update[[ii]]$XtY,
            beta_k = model_update[[ii]]$beta,
            miss_idx = data_update[[ii]]$variable_miss
          )
        })
        change_res_each <- as.numeric(unlist(change_res_each))

        # define category with same jk
        temp <- sapply(1:nrow(change_each_pair), function(x) {
          tt <- c(x, which(change_each_pair[x, ] != 0))
          return(paste0(sort(tt), collapse = ";"))
        })
        temp <- merge_sets(temp)
        cate <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
        max_change <- sapply(1:length(cate), function(ii) {
          poss <- as.numeric(unlist(cate[[ii]]))
          return(max(change_res_each[poss]))
        })
        true_update <- as.numeric(unlist(cate[[which.max(max_change)]]))
        update_status[pos.update[true_update]] <- 1
        # - re-define jk_star
        abs_cor_vals_redefine <- rowSums(abs_cor_vals_each[, true_update, drop = FALSE])
        jk_redefine <- which(abs_cor_vals_redefine == max(abs_cor_vals_redefine))
        jk_redefine <- ifelse(length(jk_redefine) == 1, jk_redefine, sample(jk_redefine, 1))
        real_update_jk[pos.update[true_update]] <- jk_redefine
      } else {
        # -- we only need to update unpaired Y_i at jk_i
        true_update <- setdiff(1:length(judge_each), pos)
        update_status[pos.update[true_update]] <- -1
        real_update_jk[pos.update[true_update]] <- jk_each[true_update]
      }
    } else {
      # -- define two sets of jk_each: one for jk = jk_each, one for jk != jk_each
      pos <- which(judge_each) # for jk = jk_each
      pos_not <- which(!judge_each) # for jk != jk_each
      check_coloc <- sapply(pos_not, function(i) {
        check_temp <- model_update[[i]]$change_loglike[jk] > cb_model_para$coloc_thresh
        return(check_temp)
      })

      # - jk is colocalized for which trait
      pos_not_coloc <- which(check_coloc)
      if (length(pos_not_coloc) == 1) {
        # - jk is the colocalized variant only for one trait
        # - if jk is colocalized variant for the trait with jk_each != jk
        # - update this trait at jk_each
        true_update <- pos_not[pos_not_coloc]
        # -- update status
        update_status[pos.update[true_update]] <- -1
        update_jk[c(1, pos.update[true_update] + 1)] <- c(jk_each[true_update], jk_each[true_update])
        real_update_jk[pos.update[true_update]] <- jk_each[true_update]
      } else if (length(pos_not_coloc) >= 2) {
        # - jk is the colocalized variant for at least two traits
        # - update the one with higher relative change of loglikelihood of each Y at each jk_each
        # - still we need to check is there is paired jk_each

        # - relative change of profile loglikelihood of each pair of Y in pos_not_coloc set
        pos_index <- pos_not[pos_not_coloc]
        change_each_pair <- check_pair_jkeach(jk_each[pos_index], pos.update[pos_index],
          model_update = model_update[pos_index],
          cb_data = cb_data,
          X_dict = X_dict[pos_index],
          jk_equiv_corr = jk_equiv_corr,
          jk_equiv_loglik = jk_equiv_loglik
        )
        pos_similar <- pos_not[which(colSums(change_each_pair) != 0)]
        # -- we only need to update paired Y_i at jk_i with higher change of loglikelihood
        data_update <- cb_data$data[pos.update]
        change_res_each <- lapply(1:length(jk_each), function(ii) {
          estimate_change_profile_res(jk_each[ii],
            X = cb_data$data[[X_dict[ii]]]$X,
            res = model_update[[ii]]$res,
            N = data_update[[ii]]$N,
            XtX = cb_data$data[[X_dict[ii]]]$XtX,
            YtY = data_update[[ii]]$YtY,
            XtY = data_update[[ii]]$XtY,
            beta_k = model_update[[ii]]$beta,
            miss_idx = data_update[[ii]]$variable_miss
          )
        })
        change_res_each <- as.numeric(unlist(change_res_each))

        if (length(pos_similar) == length(pos_not_coloc)) {
          # define category with same jk
          temp <- sapply(1:nrow(change_each_pair), function(x) {
            tt <- c(x, which(change_each_pair[x, ] != 0))
            return(paste0(sort(tt), collapse = ";"))
          })
          temp <- merge_sets(temp)
          cate <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
          max_change <- sapply(1:length(cate), function(ii) {
            poss <- as.numeric(unlist(cate[[ii]]))
            true_up <- pos_not[pos_not_coloc[poss]]
            return(max(change_res_each[true_up]))
          })
          true_update <- pos_not[pos_not_coloc[as.numeric(unlist(cate[[which.max(max_change)]]))]]
          update_status[pos.update[true_update]] <- 1
          # - re-define jk_star
          abs_cor_vals_redefine <- rowSums(abs_cor_vals_each[, true_update, drop = FALSE])
          jk_redefine <- which(abs_cor_vals_redefine == max(abs_cor_vals_redefine))
          jk_redefine <- ifelse(length(jk_redefine) == 1, jk_redefine, sample(jk_redefine, 1))
          real_update_jk[pos.update[true_update]] <- jk_redefine
        } else {
          # -- we only need to update unpaired Y_i at jk_i
          true_update <- pos_not[setdiff(pos_not_coloc, pos_similar)]
          update_status[pos.update[true_update]] <- -1
          real_update_jk[pos.update[true_update]] <- jk_each[true_update]
        }
      } else {
        # - update the ones with higher relative change of loglikelihood
        # - calculate change of loglikelihood of each Y at each jk_each
        data_update <- cb_data$data[pos.update]
        change_res_each <- lapply(1:length(jk_each), function(ii) {
          estimate_change_profile_res(jk_each[ii],
            X = cb_data$data[[X_dict[ii]]]$X,
            res = model_update[[ii]]$res,
            N = data_update[[ii]]$N,
            XtX = cb_data$data[[X_dict[ii]]]$XtX,
            YtY = data_update[[ii]]$YtY,
            XtY = data_update[[ii]]$XtY,
            beta_k = model_update[[ii]]$beta,
            miss_idx = data_update[[ii]]$variable_miss
          )
        })
        change_res_each <- as.numeric(unlist(change_res_each))
        change_res_each_jk <- change_res_each[pos] # R*
        change_res_each_jk_not <- change_res_each[-pos] # R!*

        # - relative change of profile loglikelihood of each pair of Y in R!* set
        change_each_pair <- check_pair_jkeach(jk_each[pos_not], pos.update[pos_not],
          model_update = model_update[pos_not],
          cb_data = cb_data,
          X_dict = X_dict[pos_not],
          jk_equiv_corr = jk_equiv_corr,
          jk_equiv_loglik = jk_equiv_loglik
        )
        pos_similar <- pos_not[which(colSums(change_each_pair) != 0)]

        # - compare relative change of loglikelihood
        if (func_compare == "median") {
          check_update <- median(change_res_each_jk) >= median(change_res_each_jk_not)
        } else if (func_compare == "max_min") {
          check_update <- max(change_res_each_jk) >= min(change_res_each_jk_not)
        } else if (func_compare == "max_max") {
          check_update <- max(change_res_each_jk) >= max(change_res_each_jk_not)
        } else if (func_compare == "min_max") {
          check_update <- min(change_res_each_jk) >= max(change_res_each_jk_not)
        } else if (func_compare == "min_min") {
          check_update <- min(change_res_each_jk) >= min(change_res_each_jk_not)
        }
        # thoughts for default=min_max: put jk at very end to ensure all traits are affected by jk can be updated together.
        # prioritize to update R!* but still check relative change of loglikelihood.

        if (check_update | (!check_update && length(pos_similar) != 0 && all(pos_not %in% pos_similar))) {
          # - if R* set has the higher relative change loglikelihood, we need to update traits in R* set at jk
          true_update <- pos
          update_status[pos.update[true_update]] <- 1
          if (prioritize_jkstar) {
            real_update_jk[pos.update[true_update]] <- rep(jk, length(true_update))
          } else {
            real_update_jk[pos.update[true_update]] <- jk_each[true_update]
          }
        } else {
          # - if R!* set has the higher relative change loglikelihood, we need to update traits in R!* set carefully at jk_each
          # !!!!!!! - prefer: update only one trait with highest relative change loglikelihood and not in any pairs
          candidate <- setdiff(pos_not, pos_similar)
          i_max <- candidate[order(change_res_each[candidate], decreasing = TRUE)[1]]
          true_update <- i_max
          update_status[pos.update[true_update]] <- -1
          real_update_jk[pos.update[true_update]] <- jk_each[true_update]
        }
      }
    }
  }

  # - update cb_model and report results
  cb_model_para$jk <- rbind(cb_model_para$jk, update_jk)
  cb_model_para$update_status <- cbind(cb_model_para$update_status, as.matrix(update_status))
  cb_model_para$real_update_jk <- rbind(cb_model_para$real_update_jk, real_update_jk)

  update_temp <- list(
    "update_status" = update_status,
    "real_update_jk" = real_update_jk
  )
  cb_model_para$update_temp <- update_temp
  return(cb_model_para)
}


boost_check_update_jk_focal <- function(cb_model, cb_model_para, cb_data,
                                         prioritize_jkstar = TRUE,
                                         jk_equiv_corr = 0.8, ##### more than 2 traits
                                         jk_equiv_loglik = 1, ## more than 2 traits
                                         func_compare = "min_max", ##### more than 3 traits
                                         coloc_thresh = 0.1,
                                         focal_outcome_idx = 1) {
  ############# Output #################
  # we will obtain and update the cb_model_para$update_status and cb_model_para$real_update_jk
  ######################################

  # - initial update_status and jk
  update_status <- rep(0, cb_model_para$L)
  update_jk <- rep(NA, cb_model_para$L + 1)
  real_update_jk <- rep(NA, cb_model_para$L)
  if (is.null(cb_model_para$coloc_thresh)) {
    cb_model_para$coloc_thresh <- (1 - coloc_thresh) * max(sapply(1:length(cb_model), function(i) max(cb_model[[i]]$change_loglike)))
  }

  # - update only Ys which is not stop
  pos.update <- which(cb_model_para$update_y == 1)

  if (length(pos.update) == 1) {
    # -- update status
    update_status[pos.update] <- -1
    correlation <- abs(cb_model[[pos.update]]$correlation)
    jk <- which(correlation == max(correlation))
    jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
    update_jk[c(1, pos.update + 1)] <- c(jk, jk)
    real_update_jk[pos.update] <- jk
  } else {
    # -- extract data and model based on pos.update
    model_update <- cb_model[pos.update]
    # data_update = cb_data$data[pos.update]
    X_dict <- cb_data$dict[pos.update]
    # adj_dependence
    adj_dep <- sapply(pos.update, function(ii) cb_data$data[[ii]]$dependency)

    # -- define jk and jk_each
    cor_vals_each <- do.call(cbind, lapply(model_update, function(cc) as.vector(cc$correlation)))
    abs_cor_vals_each <- sweep(abs(cor_vals_each), 2, adj_dep, `*`)
    jk_each <- apply(abs_cor_vals_each, 2, function(temp) {
      jk_temp <- which(temp == max(temp))
      return(ifelse(length(jk_temp) == 1, jk_temp, sample(jk_temp, 1)))
    })
    pp_focal <- which(pos.update == focal_outcome_idx)
    jk_focal <- jk_each[pp_focal]

    # - if jk_focal missing in all other traits
    data_update <- cb_data$data[pos.update]
    variable_missing <- Reduce("intersect", lapply(data_update[-pp_focal], function(d) d$variable_miss))
    if (jk_focal %in% variable_missing) {
      # ---- first, check LD between jk_focal and jk_each based on focal LD
      ld <- sapply(jk_each[-pp_focal], function(jki) {
        get_LD_jk1_jk2(jk_focal, jki,
          X = cb_data$data[[X_dict[pp_focal]]]$X,
          XtX = cb_data$data[[X_dict[pp_focal]]]$XtX,
          N = cb_data$data[[X_dict[pp_focal]]]$N,
          remain_jk = 1:cb_model_para$P
        )
      })
      # ----- second, if within the same LD buddies, select the following variants
      if (max(ld) > jk_equiv_corr) {
        cor_focal <- abs_cor_vals_each[, pp_focal]
        order_cor <- order(cor_focal, decreasing = TRUE)
        jk_focal_tmp <- setdiff(order_cor, variable_missing)[1]
        # ----- third, if picked variant within the same LD buddies
        ld_tmp <- get_LD_jk1_jk2(jk_focal, jk_focal_tmp,
          XtX = cb_data$data[[X_dict[pp_focal]]]$XtX,
          remain_jk = 1:cb_model_para$P
        )
        if (ld_tmp > jk_equiv_corr) {
          jk_focal <- jk_focal_tmp
          jk_each[pp_focal] <- jk_focal
        }
      }
    }


    # -- check jk_focal and jk_each
    judge_each <- check_jk_jkeach(jk_focal, jk_each,
      pos.update,
      model_update = model_update,
      cb_data = cb_data,
      X_dict = X_dict,
      jk_equiv_corr = jk_equiv_corr,
      jk_equiv_loglik = jk_equiv_loglik
    )

    judge_nofocal <- judge_each[-pp_focal]
    update_jk[c(1, pos.update + 1)] <- c(jk_focal, jk_each)

    if (all(judge_nofocal)) {
      # -- all strongest signals for non-focal are the same as strongest signal for focal trait
      # -- update all Y at jk
      update_status[pos.update] <- 1
      abs_cor_vals <- rowSums(abs_cor_vals_each)
      jk <- which(abs_cor_vals == max(abs_cor_vals))
      jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
      update_jk[1] <- jk
      if (prioritize_jkstar) {
        real_update_jk[pos.update] <- rep(jk, length(pos.update))
      } else {
        real_update_jk[pos.update] <- jk_each
      }
    } else if (all(!judge_nofocal)) {
      # -- all strongest signals for non-focal are different from the strongest signal for focal trait
      check_coloc <- sapply(pos.update, function(i) {
        check_temp <- cb_model[[i]]$change_loglike[jk_focal] > cb_model_para$coloc_thresh
        return(check_temp)
      })
      check_coloc <- check_coloc[-pp_focal]

      # -- jk_focal is colocalized for which trait
      pos_not_coloc <- which(check_coloc)
      if (length(pos_not_coloc) == 0) {
        # - update focal trait at jk_focal only
        update_status[pos.update[pp_focal]] <- -1
        real_update_jk[pos.update[pp_focal]] <- jk_focal
      } else {
        # - if jk_focal is colocalized variant for the trait with jk_each != jk
        # - update this trait at jk_each
        pp_tmp <- pos.update[-pp_focal]
        update_status[pp_tmp[pos_not_coloc]] <- -1
        real_update_jk[pp_tmp[pos_not_coloc]] <- jk_each[-pp_focal][pos_not_coloc]
      }
    } else {
      # -- define two sets of traits: one for jk_focal = jk_each, one for jk_focal != jk_each
      pos <- which(judge_each) # for jk_focal = jk_each
      pos_not <- which(!judge_each) # for jk_focal != jk_each
      check_coloc <- sapply(pos_not, function(i) {
        check_temp <- model_update[[i]]$change_loglike[jk_focal] > cb_model_para$coloc_thresh
        return(check_temp)
      })

      # - jk_focal is colocalized for which trait
      pos_not_coloc <- which(check_coloc)
      if (length(pos_not_coloc) != 0) {
        # - jk_focal is the colocalized variant for at least one trait
        # - update those traits at jk_each
        true_update <- pos_not[pos_not_coloc]
        update_status[pos.update[true_update]] <- -1
        real_update_jk[pos.update[true_update]] <- jk_each[true_update]
      } else {
        # - jk_focal is not the colocalized variant for other traits
        true_update <- pos.update[pos]
        update_status[true_update] <- 1
        abs_cor_vals <- rowSums(abs_cor_vals_each[, pos, drop = FALSE])
        jk <- which(abs_cor_vals == max(abs_cor_vals))
        jk <- ifelse(length(jk) == 1, jk, sample(jk, 1))
        # - if jk is missing in focal trait, we will prioritize jk_focal
        focal_trait_missing <- cb_data$data[[pos.update[pp_focal]]]$variable_miss
        if (jk %in% focal_trait_missing) {
          jk <- jk_focal
        }
        update_jk[1] <- jk
        if (prioritize_jkstar) {
          real_update_jk[true_update] <- rep(jk, length(true_update))
        } else {
          real_update_jk[true_update] <- jk_each[pos]
        }
      }
    }
  }

  # - update cb_model and report results
  cb_model_para$jk <- rbind(cb_model_para$jk, update_jk)
  cb_model_para$update_status <- cbind(cb_model_para$update_status, as.matrix(update_status))
  cb_model_para$real_update_jk <- rbind(cb_model_para$real_update_jk, real_update_jk)

  update_temp <- list(
    "update_status" = update_status,
    "real_update_jk" = real_update_jk
  )
  cb_model_para$update_temp <- update_temp
  return(cb_model_para)
}

#' @importFrom stats cor
get_LD_jk1_jk2 <- function(jk1, jk2,
                           X = NULL, XtX = NULL, N = NULL,
                           remain_jk = NULL) {
  if (!is.null(X)) {
    LD_temp <- suppressWarnings({
      cor(X[, c(jk1, jk2)])
    })
    LD_temp[which(is.na(LD_temp))] <- 0
    LD_temp <- LD_temp[1, 2]
  } else if (!is.null(XtX)) {
    jk1.remain <- which(remain_jk == jk1)
    jk2.remain <- which(remain_jk == jk2)
    # scaling <- if (!is.null(N)) N-1 else 1
    LD_temp <- XtX[jk1.remain, jk2.remain] # / scaling
  }
  return(LD_temp)
}

check_jk_jkeach <- function(jk, jk_each,
                            pos.update,
                            model_update,
                            cb_data, X_dict,
                            jk_equiv_corr = 0.8,
                            jk_equiv_loglik = 0.001) {
  data_update <- cb_data$data[pos.update]
  # -- check if jk ~ jk_each
  judge <- c()
  for (i in 1:length(jk_each)) {
    if (!(jk %in% data_update[[i]]$variable_miss) & !(jk_each[i] %in% data_update[[i]]$variable_miss)) {
      change_log_jk <- model_update[[i]]$change_loglike[jk]
      change_log_jkeach <- model_update[[i]]$change_loglike[jk_each[i]]
      change_each <- abs(change_log_jk - change_log_jkeach)
      LD_temp <- get_LD_jk1_jk2(jk, jk_each[i],
        X = cb_data$data[[X_dict[i]]]$X,
        XtX = cb_data$data[[X_dict[i]]]$XtX,
        N = data_update[[i]]$N,
        remain_jk = setdiff(1:length(model_update[[i]]$res), data_update[[i]]$variable_miss)
      )
      judge[i] <- (change_each <= jk_equiv_loglik) & (abs(LD_temp) >= jk_equiv_corr)
    } else {
      judge[i] <- FALSE
    }
  }
  return(judge)
}


check_pair_jkeach <- function(jk_each,
                              pos.update,
                              model_update,
                              cb_data, X_dict,
                              jk_equiv_corr = 0.8,
                              jk_equiv_loglik = 0.001) {
  data_update <- cb_data$data[pos.update]
  # -- check if jk_i ~ jk_j
  change_each_pair <- matrix(NA, nrow = length(jk_each), ncol = length(jk_each))
  for (i in 1:length(jk_each)) {
    jk_i <- jk_each[i]
    change_log_jk_i <- model_update[[i]]$change_loglike[jk_i]
    for (j in 1:length(jk_each)) {
      if (j != i) {
        jk_j <- jk_each[j]

        if (!(jk_i %in% data_update[[i]]$variable_miss) & !(jk_j %in% data_update[[i]]$variable_miss)) {
          change_log_jk_j <- model_update[[i]]$change_loglike[jk_j]
          change_each <- abs(change_log_jk_i - change_log_jk_j)
          LD_temp <- get_LD_jk1_jk2(jk_i, jk_j,
            X = cb_data$data[[X_dict[i]]]$X,
            XtX = cb_data$data[[X_dict[i]]]$XtX,
            N = data_update[[i]]$N,
            remain_jk = setdiff(1:length(model_update[[i]]$res), data_update[[i]]$variable_miss)
          )
          change_each_pair[i, j] <- (change_each <= jk_equiv_loglik) & (abs(LD_temp) >= jk_equiv_corr)
        } else {
          change_each_pair[i, j] <- FALSE
        }
      } else {
        change_each_pair[i, j] <- FALSE
      }
    }
  }
  change_each_pair <- change_each_pair + t(change_each_pair)
  return(change_each_pair)
}


estimate_change_profile_res <- function(jk,
                                        X = NULL, res = NULL, N = NULL,
                                        XtX = NULL, YtY = NULL, XtY = NULL,
                                        beta_k = NULL,
                                        miss_idx = NULL) {
  if (!is.null(X)) {
    rtr <- sum(res^2) / (N - 1)
    xtr <- t(X[, jk]) %*% res / (N - 1)
  } else if (!is.null(XtY)) {
    scaling_factor <- if (!is.null(N)) (N - 1) else 1
    yty <- YtY / scaling_factor
    xtr <- res[jk] / scaling_factor
    xtx <- XtX # / scaling_factor
    if (length(miss_idx) != 0) {
      xty <- XtY[-miss_idx] / scaling_factor
      beta_k <- beta_k[-miss_idx]
    } else {
      xty <- XtY / scaling_factor
    }
    rtr <- yty - 2 * sum(beta_k * xty) + sum((xtx %*% as.matrix(beta_k)) * beta_k)
  }
  numerator <- xtr^2 / (2 * rtr)
  denominator <- 0.5 * log(2 * pi * rtr) + 0.5
  change_loglike <- numerator / denominator

  return(change_loglike)
}
