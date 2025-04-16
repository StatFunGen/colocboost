#' @importFrom stats as.dist cutree hclust
colocboost_assemble_ucos <- function(cb_obj_single,
                                     coverage = 0.95,
                                     check_null = 0.1,
                                     check_null_method = "profile",
                                     dedup = TRUE,
                                     overlap = TRUE,
                                     squared = FALSE,
                                     n_purity = 100,
                                     min_abs_corr = 0.5,
                                     median_abs_corr = NULL,
                                     min_cluster_corr = 0.8,
                                     median_cos_abs_corr = 0.5,
                                     weaker_effect = TRUE,
                                     tol = 1e-9) {
  if (!inherits(cb_obj_single, "colocboost")) {
    stop("Input must from colocboost function!")
  }

  cb_data <- cb_obj_single$cb_data
  cb_model <- cb_obj_single$cb_model
  cb_model_para <- cb_obj_single$cb_model_para

  weights <- cb_model[[1]]$weights_path
  temp <- as.matrix(weights)
  if (sum(temp) == 0) {
    ll <- list(
      "ucos" = NULL,
      "evidence_strength" = NULL,
      "requested_coverage" = coverage
    )
    pip_av <- rep(0, cc$P)
  } else if (ncol(temp) == 1 | nrow(temp) == 1) {
    weights <- as.vector(weights)
    avWeight <- weights
    temp <- order(weights, decreasing = T)
    confidence_sets <- temp[1:min(which(cumsum(weights[temp]) > coverage))]
    evidence_strength <- sum(weights[confidence_sets])
    confidence_sets <- list(unique(confidence_sets))

    # ----- null filtering
    res_temp <- check_null_post(cb_obj_single, confidence_sets,
      coloc_outcomes = 1,
      check_null = check_null,
      check_null_method = check_null_method,
      weaker_effect = weaker_effect
    )
    cs_change <- res_temp$cs_change
    is_non_null <- res_temp$is_non_null

    if (length(is_non_null) > 0) {
      confidence_sets <- confidence_sets
      evidence_strength <- evidence_strength
      cs_change <- cs_change
    } else {
      confidence_sets <- NULL
    }

    # --- filter 3: purity
    if (length(confidence_sets) != 0) {
      pos <- confidence_sets[[1]]
      if (!is.null(cb_data$data[[1]]$XtX)) {
        pos <- match(pos, setdiff(1:cb_model_para$P, cb_data$data[[1]]$variable_miss))
      }
      purity <- matrix(get_purity(pos,
        X = cb_data$data[[1]]$X, Xcorr = cb_data$data[[1]]$XtX,
        N = cb_data$data[[1]]$N, n = n_purity
      ), 1, 3)
      purity <- as.data.frame(purity)
      colnames(purity) <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
      if (is.null(median_abs_corr)) {
        is_pure <- (purity[, 1] >= min_abs_corr)
      } else {
        is_pure <- (purity[, 1] >= min_abs_corr | purity[, 3] >= median_abs_corr)
      }
      if (is_pure) {
        confidence_sets <- confidence_sets
        evidence_strength <- evidence_strength
        cs_change <- cs_change
        purity <- purity
        row_names <- "ucos1"
        names(confidence_sets) <- row_names
        rownames(purity) <- row_names

        # --- report pip
        pip_av <- weights

        ll <- list(
          "ucos" = confidence_sets,
          "purity" = purity,
          "evidence_strength" = evidence_strength,
          "requested_coverage" = coverage,
          "cs_change" = cs_change,
          "avWeight" = as.matrix(avWeight)
        )
      } else {
        ll <- list(
          "ucos" = NULL,
          "evidence_strength" = NULL,
          "requested_coverage" = coverage
        )
        pip_av <- rep(0, cb_model_para$P)
      }
    } else {
      ll <- list(
        "ucos" = NULL,
        "evidence_strength" = NULL,
        "requested_coverage" = coverage
      )
      pip_av <- rep(0, cb_model_para$P)
    }
  } else {
    LogLik_change <- abs(diff(cb_model[[1]]$profile_loglike_each))

    # Hierachical Clustering iteration based on weights
    cormat <- get_cormat(t(weights))
    hc <- hclust(as.dist(1 - cormat))
    n_cluster <- get_n_cluster(hc, cormat, min_cluster_corr = min_cluster_corr)$n_cluster
    # n_cluster = 6
    index <- cutree(hc, n_cluster)
    B <- sapply(1:n_cluster, function(t) as.numeric(index == t))
    B <- as.matrix(B)

    #### calculate the weighted average of weights in each cluster
    total_change <- t(B) %*% as.matrix(LogLik_change)
    avWeight <- apply(B, 2, function(t) {
      idx <- which(t == 1)
      w <- LogLik_change[idx]
      weight_cluster <- t(weights[idx, , drop = FALSE])
      check_purity <- w_purity(weight_cluster,
        X = cb_data$data[[1]]$X, Xcorr = cb_data$data[[1]]$XtX,
        N = cb_data$data[[1]]$N, n = n_purity, coverage = coverage,
        min_abs_corr = min_abs_corr, median_abs_corr = median_abs_corr,
        miss_idx = cb_data$data[[1]]$variable_miss
      )
      if (length(check_purity) != 0) {
        w <- w[check_purity]
        weight_cluster <- weight_cluster[, check_purity]
        av <- as.matrix(weight_cluster) %*% as.matrix(w) / sum(w)
      } else {
        av <- rep(0, cb_model_para$P)
      }
      return(av)
    })
    avWeight <- as.matrix(avWeight)

    # Initial define the confidence sets
    confidence_sets <- avWeight_csets <- list()
    fl <- 1
    for (mm in 1:ncol(B)) {
      av <- avWeight[, mm]
      if (sum(av) != 0) {
        temp <- order(av, decreasing = T)
        confidence_sets[[fl]] <- temp[1:min(which(cumsum(avWeight[, mm][temp]) > coverage))]
        avWeight_csets[[fl]] <- av
        fl <- fl + 1
      }
    }

    if (length(confidence_sets) == 0) {
      confidence_sets <- NULL
    } else {
      avWeight <- do.call(cbind, avWeight_csets)
      # calculate evidence_strength
      evidence_strength <- unlist(lapply(1:length(confidence_sets), function(t) {
        cs <- confidence_sets[[t]]
        w <- avWeight[, t]
        return(sum(w[cs]) / sum(w))
      }))

      # ----- null filtering
      res_temp <- check_null_post(cb_obj_single, confidence_sets,
        coloc_outcomes = 1,
        check_null = check_null,
        check_null_method = check_null_method,
        weaker_effect = weaker_effect
      )
      cs_change <- res_temp$cs_change
      is_non_null <- res_temp$is_non_null
      if (length(is_non_null) > 0) {
        confidence_sets <- confidence_sets[is_non_null]
        evidence_strength <- evidence_strength[is_non_null]
        cs_change <- cs_change[is_non_null, , drop = FALSE]
        avWeight <- as.matrix(avWeight[, is_non_null])
      } else {
        confidence_sets <- NULL
      }
    }
    # --- filter 2: remove duplicate confidence sets
    if (length(confidence_sets) != 0) {
      if (dedup) {
        is_dup <- which(duplicated(confidence_sets) == TRUE)
        if (length(is_dup) > 0) {
          confidence_sets <- confidence_sets[-is_dup]
          evidence_strength <- evidence_strength[-is_dup]
          cs_change <- cs_change[-is_dup, , drop = FALSE]
          avWeight <- as.matrix(avWeight[, -is_dup])
        }
      }
    }
    # --- filter 2*: remove overlap confidence sets
    if (length(confidence_sets) >= 2) {
      if (overlap) {
        is_overlap <- c()
        ncsets <- length(confidence_sets)
        for (i in 1:(ncsets - 1)) {
          for (j in (i + 1):ncsets) {
            cset1 <- confidence_sets[[i]]
            cset2 <- confidence_sets[[j]]
            if (all(cset1 %in% cset2) | all(cset2 %in% cset1)) {
              total1 <- total_change[i]
              total2 <- total_change[j]
              is_overlap <- c(is_overlap, ifelse(total1 > total2, j, i))
            }
          }
        }
        if (length(is_overlap) > 0) {
          confidence_sets <- confidence_sets[-is_overlap]
          evidence_strength <- evidence_strength[-is_overlap]
          cs_change <- cs_change[-is_overlap, , drop = FALSE]
          avWeight <- as.matrix(avWeight[, -is_overlap])
        }
      }
    }

    # --- filter 2*: remove overlap confidence sets based on median_cos_abs_corr
    if (length(confidence_sets) >= 2) {
      if (overlap) {
        # calculate between purity
        ncsets <- length(confidence_sets)
        min_between <- max_between <- ave_between <- matrix(0, nrow = ncsets, ncol = ncsets)
        for (i.between in 1:(ncsets - 1)) {
          for (j.between in (i.between + 1):ncsets) {
            cset1 <- confidence_sets[[i.between]]
            cset2 <- confidence_sets[[j.between]]
            res <- get_between_purity(cset1, cset2,
              X = cb_data$data[[1]]$X,
              Xcorr = cb_data$data[[1]]$XtX,
              N = cb_data$data[[1]]$N,
              miss_idx = cb_data$data[[1]]$variable_miss,
              P = cb_model_para$P
            )
            min_between[i.between, j.between] <- min_between[j.between, i.between] <- res[1]
            max_between[i.between, j.between] <- max_between[j.between, i.between] <- res[2]
            ave_between[i.between, j.between] <- ave_between[j.between, i.between] <- res[3]
          }
        }
        is.between <- (min_between > min_abs_corr) * (abs(max_between - 1) < tol) * (ave_between > median_cos_abs_corr)
        if (sum(is.between) != 0) {
          # potential_merged <- find_merged_csets(is.between)
          temp <- sapply(1:nrow(is.between), function(x) {
            tt <- c(x, which(is.between[x, ] != 0))
            return(paste0(sort(tt), collapse = ";"))
          })
          potential_merged <- sapply(names(table(temp)), strsplit, split = ";")
          confidence_sets_merged <- list()
          avWeight_merged <- cs_change_merged <- evidence_strength_merged <- c()
          is_merged <- c()
          for (i.m in 1:length(potential_merged)) {
            temp_set <- as.numeric(potential_merged[[i.m]])
            is_merged <- c(is_merged, temp_set)
            # define merged set
            confidence_sets_merged <- c(
              confidence_sets_merged,
              list(unique(unlist(confidence_sets[temp_set])))
            )
            # refine avWeight
            avWeight_merged <- cbind(avWeight_merged, rowMeans(as.matrix(avWeight[, temp_set])))
            # refine cs_change
            cs_change_merged <- c(cs_change_merged, max(cs_change[temp_set, ]))
            # refine evidence strength
            evidence_strength_merged <- c(
              evidence_strength_merged,
              min(evidence_strength[temp_set])
            )
          }
          confidence_sets_single <- confidence_sets[-is_merged]
          avWeight_sets_single <- avWeight[, -is_merged]
          cs_change_single <- cs_change[-is_merged, , drop = FALSE]
          evidence_strength_single <- evidence_strength[-is_merged]

          # --- combine merged and single
          confidence_sets <- c(confidence_sets_single, confidence_sets_merged)
          evidence_strength <- c(evidence_strength_single, evidence_strength_merged)
          cs_change <- rbind(cs_change_single, as.matrix(cs_change_merged))
          avWeight <- cbind(avWeight_sets_single, avWeight_merged)
        }
      }
    }


    # --- filter 3: purity
    if (length(confidence_sets) != 0) {
      purity <- NULL
      for (ee in 1:length(confidence_sets)) {
        pos <- confidence_sets[[ee]]
        if (!is.null(cb_data$data[[1]]$XtX)) {
          pos <- match(pos, setdiff(1:cb_model_para$P, cb_data$data[[1]]$variable_miss))
        }
        purity <-
          rbind(
            purity,
            matrix(get_purity(pos,
              X = cb_data$data[[1]]$X, Xcorr = cb_data$data[[1]]$XtX,
              N = cb_data$data[[1]]$N, n = n_purity
            ), 1, 3)
          )
      }
      purity <- as.data.frame(purity)
      colnames(purity) <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
      if (is.null(median_abs_corr)) {
        is_pure <- which(purity[, 1] >= min_abs_corr)
      } else {
        is_pure <- which(purity[, 1] >= min_abs_corr | purity[, 3] >= median_abs_corr)
      }

      if (length(is_pure) > 0) {
        confidence_sets <- confidence_sets[is_pure]
        evidence_strength <- evidence_strength[is_pure]
        cs_change <- cs_change[is_pure, , drop = FALSE]
        avWeight <- as.matrix(avWeight[, is_pure])
        purity <- purity[is_pure, ]
        row_names <- paste0("ucos", is_pure)
        names(confidence_sets) <- row_names
        rownames(purity) <- row_names

        # --- report pip
        pip_av <- as.vector(1 - apply(1 - avWeight, 1, prod))


        # Re-order CS list and purity rows based on purity.
        ordering <- order(purity[, 1], decreasing = TRUE)
        ll <- list(
          "ucos" = confidence_sets[ordering],
          "purity" = purity[ordering, ],
          "evidence_strength" = evidence_strength[ordering],
          "requested_coverage" = coverage,
          "cs_change" = cs_change[ordering, , drop = FALSE],
          "avWeight" = as.matrix(avWeight[, ordering])
        )
      } else {
        ll <- list(
          "ucos" = NULL,
          "evidence_strength" = NULL,
          "requested_coverage" = coverage
        )
        pip_av <- rep(0, cb_model_para$P)
      }
    } else {
      ll <- list(
        "ucos" = NULL,
        "evidence_strength" = NULL,
        "requested_coverage" = coverage
      )
      pip_av <- rep(0, cb_model_para$P)
    }
  }

  out <- list("ucos" = ll, "pip_av" = pip_av)
  return(out)
}
