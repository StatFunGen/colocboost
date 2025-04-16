#' @importFrom stats as.dist cutree hclust
colocboost_assemble_cos <- function(cb_obj,
                                    coverage = 0.95,
                                    weight_fudge_factor = 1.5,
                                    check_null = 0.1,
                                    check_null_method = "profile",
                                    dedup = TRUE,
                                    overlap = TRUE,
                                    n_purity = 100,
                                    min_abs_corr = 0.5,
                                    sec_coverage_thresh = 0.8,
                                    median_abs_corr = NULL,
                                    min_cluster_corr = 0.8,
                                    median_cos_abs_corr = 0.8,
                                    tol = 1e-9) {
  if (!inherits(cb_obj, "colocboost")) {
    stop("Input must from colocboost function!")
  }

  cb_model <- cb_obj$cb_model
  cb_model_para <- cb_obj$cb_model_para
  cb_data <- cb_obj$cb_data

  # define the confident sets for colocalization
  update <- cb_model_para$update_status
  pos.coloc <- which(colSums(update) > 1)
  if (length(pos.coloc) == 0) {
    ll <- list(
      "cos" = NULL,
      "evidence_strength" = NULL,
      "requested_coverage" = coverage
    )
  } else if (length(pos.coloc) == 1) {
    # try pos.coloc = pos.coloc[15]
    coloc_outcomes <- which(update[, pos.coloc] == 1)
    avWeight <- get_avWeigth(cb_model, coloc_outcomes, update, pos.coloc, name_weight = T)

    # - check purity for each outcome
    check_purity <- c()
    for (iiii in 1:length(coloc_outcomes)) {
      X_dict <- cb_data$dict[coloc_outcomes[iiii]]
      tmp <- w_purity(avWeight[, iiii, drop = FALSE],
        X = cb_data$data[[X_dict]]$X, Xcorr = cb_data$data[[X_dict]]$XtX,
        N = cb_data$data[[coloc_outcomes[iiii]]]$N, n = n_purity, coverage = sec_coverage_thresh,
        min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr,
        miss_idx = cb_data$data[[coloc_outcomes[iiii]]]$variable_miss
      )
      check_purity[iiii] <- length(tmp) == 1
    }
    if (sum(check_purity) < 2) {
      is_non_null <- NULL
    } else {
      pos_purity <- which(check_purity)
      avWeight <- avWeight[, pos_purity, drop = FALSE]
      coloc_outcomes <- coloc_outcomes[pos_purity]
      weights <- get_integrated_weight(avWeight, weight_fudge_factor = weight_fudge_factor)
      coloc_cos <- get_in_cos(weights, coverage = coverage)
      evidence_strength <- sum(weights[coloc_cos[[1]]])

      # ----- null filtering
      res_temp <- check_null_post(cb_obj, coloc_cos, coloc_outcomes,
        check_null = check_null,
        check_null_method = check_null_method
      )
      cs_change <- res_temp$cs_change
      is_non_null <- res_temp$is_non_null
    }

    if (length(is_non_null) == 0) {
      ll <- list(
        "cos" = NULL,
        "evidence_strength" = NULL,
        "requested_coverage" = coverage
      )
    } else {
      purity <- c()
      pos <- as.numeric(unlist(coloc_cos))
      for (i in coloc_outcomes) {
        X_dict <- cb_data$dict[i]
        if (!is.null(cb_data$data[[X_dict]]$XtX)) {
          pos <- match(pos, setdiff(1:cb_model_para$P, cb_data$data[[i]]$variable_miss))
        }
        p_tmp <- matrix(get_purity(pos,
          X = cb_data$data[[X_dict]]$X,
          Xcorr = cb_data$data[[X_dict]]$XtX,
          N = cb_data$data[[i]]$N, n = n_purity
        ), 1, 3)
        purity <- c(purity, list(p_tmp))
      }
      purity_all <- Reduce(pmax, purity)
      purity_all <- as.data.frame(purity_all)
      colnames(purity_all) <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
      if (is.null(median_abs_corr)) {
        is_pure <- (purity_all[, 1] >= min_abs_corr)
      } else {
        is_pure <- (purity_all[, 1] >= min_abs_corr | purity_all[, 3] >= median_abs_corr)
      }
      if (!is_pure) {
        ll <- list(
          "cos" = NULL,
          "evidence_strength" = NULL,
          "requested_coverage" = coverage
        )
      } else {
        row_names <- "cos1"
        names(coloc_cos) <- row_names
        rownames(purity_all) <- row_names
        coloc_outcomes <- list(coloc_outcomes)

        ll <- list(
          "cos" = coloc_cos,
          "purity" = purity_all,
          "evidence_strength" = evidence_strength,
          "requested_coverage" = coverage,
          "cs_change" = cs_change,
          "avWeight" = list(avWeight),
          "coloc_outcomes" = coloc_outcomes
        )
      }
    }
  } else {
    coloc_candidate <- update[, pos.coloc]
    coloc_candidate <- apply(coloc_candidate, 2, paste0, collapse = ",")
    coloc_temp <- table(coloc_candidate)
    # iterations for each colocalization sets
    pos_coloc_sets <- lapply(1:length(coloc_temp), function(x) {
      which(coloc_candidate == names(coloc_temp)[x])
    })
    names(pos_coloc_sets) <- names(coloc_temp)
    # - define coloc_sets
    coloc_sets <- avWeight_coloc_sets <-
      total_change_Loglik_coloc <- evidence_strength_coloc <-
      cs_change_coloc <- coloc_outcomes_sets <- list()
    flag <- 0
    for (i in 1:length(coloc_temp)) {
      pos_temp_coloc_each <- pos_coloc_sets[[i]]
      coloc_outcomes <- which(unlist(strsplit(names(coloc_temp)[i], ",")) == 1)

      # - if only one iteration for this coloc_set
      if (length(pos_temp_coloc_each) == 1) {
        avWeight <- get_avWeigth(cb_model, coloc_outcomes, update, pos.coloc[pos_temp_coloc_each], name_weight = T)

        # - check purity for each outcome
        check_purity <- c()
        for (iiii in 1:length(coloc_outcomes)) {
          X_dict <- cb_data$dict[coloc_outcomes[iiii]]
          tmp <- w_purity(avWeight[, iiii, drop = FALSE],
            X = cb_data$data[[X_dict]]$X, Xcorr = cb_data$data[[X_dict]]$XtX,
            N = cb_data$data[[coloc_outcomes[iiii]]]$N, n = n_purity, coverage = sec_coverage_thresh,
            min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr,
            miss_idx = cb_data$data[[coloc_outcomes[iiii]]]$variable_miss
          )
          check_purity[iiii] <- length(tmp) == 1
        }
        if (sum(check_purity) < 2) {
          is_non_null <- NULL
        } else {
          pos_purity <- which(check_purity)
          avWeight <- avWeight[, pos_purity, drop = FALSE]
          coloc_outcomes <- coloc_outcomes[pos_purity]
          weights <- get_integrated_weight(avWeight, weight_fudge_factor = weight_fudge_factor)
          coloc_cos <- get_in_cos(weights, coverage = coverage)
          evidence_strength <- sum(weights[coloc_cos[[1]]])

          # ----- null filtering
          res_temp <- check_null_post(cb_obj, coloc_cos, coloc_outcomes,
            check_null = check_null,
            check_null_method = check_null_method
          )

          cs_change <- res_temp$cs_change
          is_non_null <- res_temp$is_non_null
        }

        if (length(is_non_null) > 0) {
          flag <- flag + 1
          coloc_sets[[flag]] <- as.numeric(unlist(coloc_cos))
          total_change_Loglik_coloc[[flag]] <- "ONE"
          avWeight_coloc_sets[[flag]] <- avWeight
          evidence_strength_coloc[[flag]] <- evidence_strength
          cs_change_coloc[[flag]] <- cs_change
          coloc_outcomes_sets[[flag]] <- coloc_outcomes
        }
      } else {
        av <- lapply(coloc_outcomes, function(ii) {
          get_avWeigth(cb_model, ii, update, pos.coloc[pos_temp_coloc_each])
        })
        weight_coloc <- do.call(cbind, av)

        # Hierachical Clustering iteration based on sequenced weights
        cormat <- get_cormat(t(weight_coloc))
        hc <- hclust(as.dist(1 - cormat))
        n_cluster <- get_n_cluster(hc, cormat, min_cluster_corr = min_cluster_corr)$n_cluster
        index <- cutree(hc, n_cluster)
        B <- sapply(1:n_cluster, function(t) as.numeric(index == t))
        B <- as.matrix(B)

        # --- define the confident sets
        LogLik_change <- abs(diff(cb_model_para$profile_loglik)[pos.coloc[pos_temp_coloc_each]])
        avWeight <- apply(B, 2, function(t) {
          idx <- which(t == 1)
          w <- LogLik_change[idx]
          check_purity <- list()
          for (iiii in 1:length(coloc_outcomes)) {
            weight_cluster <- t(av[[iiii]][idx, , drop = FALSE])
            X_dict <- cb_data$dict[coloc_outcomes[iiii]]
            tmp <- w_purity(weight_cluster,
              X = cb_data$data[[X_dict]]$X, Xcorr = cb_data$data[[X_dict]]$XtX,
              N = cb_data$data[[coloc_outcomes[iiii]]]$N, n = n_purity, coverage = sec_coverage_thresh,
              min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr,
              miss_idx = cb_data$data[[coloc_outcomes[iiii]]]$variable_miss
            )
            check_purity[[iiii]] <- tmp
          }
          # check_purity <- unique(as.numeric(check_purity))
          check_purity <- Reduce("intersect", check_purity)
          if (length(check_purity) != 0) {
            w <- w[check_purity]
            w_tmp <- t(weight_coloc[idx, , drop = FALSE])
            weight_cluster <- w_tmp[, check_purity]
            avW <- as.matrix(weight_cluster) %*% as.matrix(w) / sum(w)
          } else {
            avW <- rep(0, cb_model_para$P * length(coloc_outcomes))
          }
          return(avW)
        })
        avWeight <- as.matrix(avWeight)
        avWeight_each <- lapply(1:length(coloc_outcomes), function(x) {
          start <- (x - 1) * cb_model_para$P + 1
          end <- x * cb_model_para$P
          return(as.matrix(avWeight[c(start:end), ]))
        })
        avWeight_coloc <- lapply(1:ncol(B), function(x) {
          a <- lapply(avWeight_each, function(xx) {
            xx[, x]
          })
          a <- do.call(cbind, a)
          colnames(a) <- paste0("outcome", coloc_outcomes)
          return(a)
        })
        total_change <- t(B) %*% as.matrix(LogLik_change) / colSums(B)
        coloc_cos <- evidence_strength <- avWeight_csets <- list()
        total_change_csets <- c()
        fl <- 1
        for (i.w in 1:length(avWeight_coloc)) {
          w <- avWeight_coloc[[i.w]]
          if (sum(w) != 0) {
            weights <- get_integrated_weight(w, weight_fudge_factor = weight_fudge_factor)
            csets <- get_in_cos(weights, coverage = coverage)
            coloc_cos[[fl]] <- unlist(csets)
            evidence_strength[[fl]] <- sum(weights[coloc_cos[[fl]]])
            avWeight_csets[[fl]] <- w
            total_change_csets <- c(total_change_csets, total_change[i.w])
            fl <- fl + 1
          }
        }
        if (length(coloc_cos) != 0) {
          avWeight_coloc <- avWeight_csets
          total_change <- total_change_csets
          # ----- null filtering
          res_temp <- check_null_post(cb_obj, coloc_cos, coloc_outcomes,
            check_null = check_null,
            check_null_method = check_null_method
          )
          cs_change <- res_temp$cs_change
          is_non_null <- res_temp$is_non_null
          avWeight_coloc <- avWeight_csets

          if (length(is_non_null) > 0) {
            coloc_cos <- coloc_cos[is_non_null]
            evidence_strength <- evidence_strength[is_non_null]
            cs_change <- cs_change[is_non_null, ]
            avWeight_coloc <- avWeight_coloc[is_non_null]
            total_change <- total_change[is_non_null]
            for (i.b in 1:length(is_non_null)) {
              flag <- flag + 1
              coloc_sets[[flag]] <- as.numeric(coloc_cos[[i.b]])
              total_change_Loglik_coloc[[flag]] <- total_change[i.b]
              avWeight_coloc_sets[[flag]] <- avWeight_coloc[[i.b]]
              evidence_strength_coloc[[flag]] <- evidence_strength[i.b]
              cs_change_coloc[[flag]] <- cs_change[i.b, ]
              coloc_outcomes_sets[[flag]] <- coloc_outcomes
            }
          }
        }
      }
    }

    # --- filter 2*: remove overlap confidence sets
    if (length(coloc_sets) >= 2) {
      if (overlap) {
        is_overlap <- c()
        ncsets <- length(coloc_sets)
        for (i2 in 1:(ncsets - 1)) {
          for (j in (i2 + 1):ncsets) {
            cset1 <- coloc_sets[[i2]]
            cset2 <- coloc_sets[[j]]
            if (all(cset1 %in% cset2) | all(cset2 %in% cset1)) {
              if (cb_model_para$L == 2) {
                is_overlap <- c(
                  is_overlap,
                  check_two_overlap_sets(total_change_Loglik_coloc, i2, j)
                )
              } else {
                outcome1 <- sort(colnames(avWeight_coloc_sets[[i2]]))
                outcome2 <- sort(colnames(avWeight_coloc_sets[[j]]))
                if (identical(outcome1, outcome2)) {
                  is_overlap <- c(
                    is_overlap,
                    check_two_overlap_sets(total_change_Loglik_coloc, i2, j)
                  )
                } else {
                  if (all(outcome1 %in% outcome2)) {
                    is_overlap <- c(is_overlap, i2)
                  } else if (all(outcome2 %in% outcome1)) {
                    is_overlap <- c(is_overlap, j)
                  }
                }
              }
            }
          }
        }
        is_overlap <- unique(is_overlap)
        if (length(is_overlap) > 0) {
          coloc_sets <- coloc_sets[-is_overlap]
          evidence_strength_coloc <- evidence_strength_coloc[-is_overlap]
          cs_change_coloc <- cs_change_coloc[-is_overlap]
          avWeight_coloc_sets <- avWeight_coloc_sets[-is_overlap]
          total_change_Loglik_coloc <- total_change_Loglik_coloc[-is_overlap]
          coloc_outcomes_sets <- coloc_outcomes_sets[-is_overlap]
        }
      }
    }



    # --- filter 2*: remove overlap confidence sets based on median_cos_abs_corr
    if (length(coloc_sets) >= 2) {
      if (overlap) {
        # calculate between purity
        ncsets <- length(coloc_sets)
        min_between <- max_between <- ave_between <- matrix(0, nrow = ncsets, ncol = ncsets)
        for (i.between in 1:(ncsets - 1)) {
          for (j.between in (i.between + 1):ncsets) {
            cset1 <- coloc_sets[[i.between]]
            cset2 <- coloc_sets[[j.between]]
            res <- list()
            for (i in 1:cb_model_para$L) {
              X_dict <- cb_data$dict[i]
              res[[i]] <- get_between_purity(cset1, cset2,
                X = cb_data$data[[X_dict]]$X,
                Xcorr = cb_data$data[[X_dict]]$XtX,
                N = cb_data$data[[i]]$N,
                miss_idx = cb_data$data[[i]]$variable_miss,
                P = cb_model_para$P
              )
            }
            res <- Reduce(pmax, res)
            min_between[i.between, j.between] <- min_between[j.between, i.between] <- res[1]
            max_between[i.between, j.between] <- max_between[j.between, i.between] <- res[2]
            ave_between[i.between, j.between] <- ave_between[j.between, i.between] <- res[3]
          }
        }
        is.between <- (min_between > min_abs_corr) * (abs(max_between - 1) < tol) * (ave_between > median_cos_abs_corr)
        if (sum(is.between) != 0) {
          temp <- sapply(1:nrow(is.between), function(x) {
            tt <- c(x, which(is.between[x, ] != 0))
            return(paste0(sort(tt), collapse = ";"))
          })
          temp <- merge_sets(temp)
          potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
          potential_merged <- potential_merged[which(sapply(potential_merged, length) >= 2)]
          coloc_sets_merged <- avWeight_merged <-
            cs_change_merged <- evidence_strength_merged <-
            coloc_outcomes_merged <- list()
          is_merged <- c()
          for (i.m in 1:length(potential_merged)) {
            temp_set <- as.numeric(potential_merged[[i.m]])
            is_merged <- c(is_merged, temp_set)
            # define merged set
            coloc_sets_merged <- c(coloc_sets_merged, list(unique(unlist(coloc_sets[temp_set]))))
            # refine avWeight
            merged <- do.call(cbind, avWeight_coloc_sets[temp_set])
            unique_coloc_outcomes <- unique(colnames(merged))
            coloc_p_merged <- as.numeric(gsub("[^0-9.]+", "", unique_coloc_outcomes))
            coloc_outcomes_merged <- c(
              coloc_outcomes_merged,
              list(coloc_p_merged[order(coloc_p_merged)])
            )
            temp <- lapply(unique_coloc_outcomes, function(tt) {
              pos <- which(colnames(merged) == tt)
              tmp <- as.matrix(rowMeans(as.matrix(merged[, pos])))
              colnames(tmp) <- tt
              return(tmp)
            })
            temp <- do.call(cbind, temp)
            temp <- temp[, order(unique_coloc_outcomes)]
            # colnames(temp) <- unique_coloc_outcomes
            avWeight_merged <- c(avWeight_merged, list(temp))
            # refine cs_change
            cs_change_merged <- c(
              cs_change_merged,
              list(do.call(pmax, cs_change_coloc[temp_set]))
            )
            # refine evidence strength
            evidence_strength_merged <- c(
              evidence_strength_merged,
              list(min(unlist(evidence_strength_coloc[temp_set])))
            )
          }
          coloc_sets_single <- coloc_sets[-is_merged]
          avWeight_sets_single <- avWeight_coloc_sets[-is_merged]
          cs_change_single <- cs_change_coloc[-is_merged]
          evidence_strength_single <- evidence_strength_coloc[-is_merged]
          coloc_outcomes_single <- coloc_outcomes_sets[-is_merged]

          # --- combine merged and single
          coloc_sets <- c(coloc_sets_single, coloc_sets_merged)
          evidence_strength_coloc <- c(evidence_strength_single, evidence_strength_merged)
          cs_change_coloc <- c(cs_change_single, cs_change_merged)
          avWeight_coloc_sets <- c(avWeight_sets_single, avWeight_merged)
          coloc_outcomes_sets <- c(coloc_outcomes_single, coloc_outcomes_merged)
        }
      }
    }

    # --- filter 3: purity
    if (length(coloc_sets) != 0) {
      purity <- vector(mode = "list", length = length(coloc_sets))
      for (ee in 1:length(coloc_sets)) {
        coloc_t <- coloc_outcomes_sets[[ee]]
        p_tmp <- c()
        for (i3 in coloc_t) {
          pos <- coloc_sets[[ee]]
          X_dict <- cb_data$dict[i3]
          if (!is.null(cb_data$data[[X_dict]]$XtX)) {
            pos <- match(pos, setdiff(1:cb_model_para$P, cb_data$data[[i3]]$variable_miss))
          }
          tmp <- matrix(get_purity(pos,
            X = cb_data$data[[X_dict]]$X,
            Xcorr = cb_data$data[[X_dict]]$XtX,
            N = cb_data$data[[i3]]$N, n = n_purity
          ), 1, 3)
          p_tmp <- rbind(p_tmp, tmp)
        }
        purity[[ee]] <- matrix(apply(p_tmp, 2, max), 1, 3)
      }
      purity_all <- do.call(rbind, purity)
      purity_all <- as.data.frame(purity_all)
      colnames(purity_all) <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
      if (is.null(median_abs_corr)) {
        is_pure <- which(purity_all[, 1] >= min_abs_corr)
      } else {
        is_pure <- which(purity_all[, 1] >= min_abs_corr | purity_all[, 3] >= median_abs_corr)
      }

      if (length(is_pure) > 0) {
        coloc_sets <- coloc_sets[is_pure]
        evidence_strength_coloc <- unlist(evidence_strength_coloc[is_pure])
        cs_change_coloc <- do.call(rbind, cs_change_coloc[is_pure])
        avWeight_coloc_sets <- avWeight_coloc_sets[is_pure]
        coloc_outcomes_sets <- coloc_outcomes_sets[is_pure]
        purity_all <- purity_all[is_pure, ]
        row_names <- paste0("cos", is_pure)
        names(coloc_sets) <- row_names
        rownames(purity_all) <- row_names
        rownames(cs_change_coloc) <- row_names


        if (length(coloc_sets) == 1) {
          ll <- list(
            "cos" = coloc_sets,
            "purity" = purity_all,
            "evidence_strength" = evidence_strength_coloc,
            "requested_coverage" = coverage,
            "cs_change" = as.matrix(cs_change_coloc),
            "avWeight" = avWeight_coloc_sets,
            "coloc_outcomes" = coloc_outcomes_sets
          )
        } else {
          # Re-order CS list and purity rows based on purity.
          ordering <- order(purity_all[, 1], decreasing = TRUE)
          ll <- list(
            "cos" = coloc_sets[ordering],
            "purity" = purity_all[ordering, ],
            "evidence_strength" = evidence_strength_coloc[ordering],
            "requested_coverage" = coverage,
            "cs_change" = cs_change_coloc[ordering, ],
            "avWeight" = avWeight_coloc_sets[ordering],
            "coloc_outcomes" = coloc_outcomes_sets[ordering]
          )
        }
      } else {
        ll <- list(
          "cos" = NULL,
          "evidence_strength" = NULL,
          "requested_coverage" = coverage
        )
      }
    } else {
      ll <- list(
        "cos" = NULL,
        "evidence_strength" = NULL,
        "requested_coverage" = coverage
      )
    }
  }
  out <- list("cos" = ll)
  return(out)
}
