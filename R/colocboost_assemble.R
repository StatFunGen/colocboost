#' Main function for colocboost post aggregating analysis
#'
#' @details
#' The following functions are included in the post-hoc analysis:
#'
#' Colocalization signal - `colocboost_assemble_cos` - identify the colocalized confidence sets and the corresponding causal configurations.
#'
#' Un-colocalization signal - `colocboost_assemble_ucos` - identify the causal confidence sets for each outcome only.
#'
#' Add-hoc merge_cos functions including
#'
#' \itemize{
#'  \item{merge_coloc_single}{merge the colocalized sets and the single causal set if pass the \code{median_cos_abs_corr}}
#'  \item{merge_single}{merge the single causal sets for different outcomes if pass the \code{median_cos_abs_corr}}
#' }
#'
#' Refine of the colocalization sets (TO-DO-LIST)
#'
#' Summary of the colocboost results and get the output of colocboost (TO-DO-LIST)
#'
#' @noRd
colocboost_assemble <- function(cb_obj,
                                coverage = 0.95,
                                weight_fudge_factor = 1.5,
                                check_null = 0.1,
                                check_null_method = "profile",
                                check_null_max = 2e-5,
                                dedup = TRUE,
                                overlap = TRUE,
                                n_purity = 100,
                                min_abs_corr = 0.5,
                                sec_coverage_thresh = 0.8,
                                median_abs_corr = NULL,
                                min_cluster_corr = 0.8,
                                median_cos_abs_corr = 0.8,
                                weaker_effect = TRUE,
                                merge_cos = TRUE,
                                tol = 1e-9,
                                output_level = 1) {
  if (!inherits(cb_obj, "colocboost")) {
    stop("Input must from colocboost object!")
  }

  # - data information
  data_info <- get_data_info(cb_obj)
  model_info <- get_model_info(cb_obj, outcome_names = data_info$outcome_info$outcome_names)
  if (data_info$n_outcomes == 1 & output_level == 1) {
    output_level <- 2
  }
  if (cb_obj$cb_model_para$num_updates == 1) {
    cb_output <- list(
      "cos_summary" = NULL,
      "vcp" = NULL,
      "cos_details" = NULL,
      "data_info" = data_info,
      "model_info" = model_info
    )
    # - save model and all coloc and single information for diagnostic
    if (output_level != 1) {
      tmp <- get_full_output(cb_obj = cb_obj, past_out = NULL, variables = NULL)
      if (output_level == 2) {
        cb_output$ucos_details <- tmp$ucos_detials
        cb_output <- cb_output[c("cos_summary", "vcp", "cos_details", "data_info", "model_info", "ucos_details")]
      } else {
        cb_output$ucos_details <- tmp$ucos_detials
        cb_output$diagnostic_details <- tmp[-1]
        cb_output <- cb_output[c("cos_summary", "vcp", "cos_details", "data_info", "model_info", "ucos_details", "diagnostic_details")]
      }
      if (data_info$n_outcome == 1) {
        cb_output <- list("ucos_summary", "vpa" = NULL, "ucos_details" = NULL, "data_info" = data_info)
      }
    }
  } else {
    if (cb_obj$cb_model_para$M == 1) {
      check_null_method <- "obj"
      check_null <- check_null_max <- 0.1
    }
    cb_obj <- get_max_profile(cb_obj, check_null_max = check_null_max, check_null_method = check_null_method)
    # --------- about colocalized confidence sets ---------------------------------
    out_cos <- colocboost_assemble_cos(cb_obj,
      coverage = coverage,
      weight_fudge_factor = weight_fudge_factor,
      check_null = check_null,
      check_null_method = check_null_method,
      dedup = dedup,
      overlap = overlap,
      n_purity = n_purity,
      min_abs_corr = min_abs_corr,
      sec_coverage_thresh = sec_coverage_thresh,
      median_abs_corr = median_abs_corr,
      min_cluster_corr = min_cluster_corr,
      median_cos_abs_corr = median_cos_abs_corr,
      tol = tol
    )

    # --------- about non-colocalized confidence sets ---------------------------------
    L <- cb_obj$cb_model_para$L
    if (L == 1) {
      weaker_effect <- FALSE
    }
    update <- cb_obj$cb_model_para$update_status
    ucos_each <- list()
    change_obj_each <- purity_each <- vector(mode = "list", length = L)
    avWeight_ucos_each <- ucos_outcome <- c()
    for (i in 1:L) {
      pos.each.single <- which((colSums(update) == 1) * (update[i, ] == 1) + (update[i, ] == -1) == 1)
      pos.each.all <- which(update[i, ] != 0)
      # - if length(pos.each.all) = 1, idea check lfsr
      if (length(pos.each.single) >= 1) {
        pos_temp <- match(pos.each.single, pos.each.all)
        cb_obj_single <- list(
          "cb_data" = list(),
          "cb_model" = cb_obj$cb_model[i],
          "cb_model_para" = cb_obj$cb_model_para
        )
        cb_obj_single$cb_data <- list("data" = cb_obj$cb_data$data[i], "dict" = 1)
        # change R number
        cb_obj_single$cb_model_para$L <- 1
        # for single outcome, if missing we need to extract from other
        if (!is.null(cb_obj_single$cb_data$data[[1]][["Y"]])) {
          if (is.null(cb_obj_single$cb_data$data[[1]]$X)) {
            X_dict <- cb_obj$cb_data$dict[i]
            cb_obj_single$cb_data$data[[1]]$X <- cb_obj$cb_data$data[[X_dict]]$X
          }
        }
        if (!is.null(cb_obj_single$cb_data$data[[1]][["XtY"]])) {
          if (is.null(cb_obj_single$cb_data$data[[1]]$XtX)) {
            X_dict <- cb_obj$cb_data$dict[i]
            cb_obj_single$cb_data$data[[1]]$X <- cb_obj$cb_data$data[[X_dict]]$XtX
          }
        }
        class(cb_obj_single) <- "colocboost"
        out_ucos_each <- colocboost_assemble_ucos(cb_obj_single,
          coverage = coverage,
          check_null = check_null,
          check_null_method = check_null_method,
          dedup = dedup,
          overlap = overlap,
          n_purity = n_purity,
          min_abs_corr = min_abs_corr,
          median_abs_corr = median_abs_corr,
          min_cluster_corr = min_cluster_corr,
          median_cos_abs_corr = median_cos_abs_corr,
          weaker_effect = weaker_effect,
          tol = tol
        )
        aaa <- out_ucos_each$ucos$ucos
        if (length(aaa) != 0) {
          ucos_outcome <- c(ucos_outcome, rep(i, length(aaa)))
          names(aaa) <- paste0("sets:", "Y", i, ":", names(aaa))
          ucos_each <- c(ucos_each, aaa)
          bbb <- out_ucos_each$ucos$avWeight
          colnames(bbb) <- names(aaa)
          avWeight_ucos_each <- cbind(avWeight_ucos_each, bbb)
          temp <- matrix(0, ncol = L, nrow = length(aaa))
          temp[, i] <- as.numeric(unlist(out_ucos_each$ucos$cs_change))
          rownames(temp) <- names(aaa)
          change_obj_each[[i]] <- temp
          purity_tmp <- out_ucos_each$ucos$purity
          rownames(purity_tmp) <- names(aaa)
          purity_each[[i]] <- purity_tmp
        }
      }
    }
    out_ucos <- list(
      "ucos_each" = ucos_each,
      "avW_ucos_each" = avWeight_ucos_each,
      "change_obj_each" = as.data.frame(do.call(rbind, change_obj_each)),
      "purity_each" = as.data.frame(do.call(rbind, purity_each)),
      "ucos_outcome" = ucos_outcome
    )

    # -------------------------- Post Processing --------------------------
    # ----- Remove high correlated confident sets between cos and ucos
    if (length(out_cos$cos$cos) != 0 & length(out_ucos$ucos_each) != 0) {
      past_out <- merge_cos_ucos(cb_obj, out_cos, out_ucos,
        coverage = coverage,
        min_abs_corr = min_abs_corr, tol = tol,
        median_cos_abs_corr = median_cos_abs_corr
      )
    } else if (length(out_cos$cos$cos) != 0 & length(out_ucos$ucos_each) == 0) {
      past_out <- list("ucos" = NULL, "cos" = out_cos)
    } else if (length(out_cos$cos$cos) == 0 & length(out_ucos$ucos_each) != 0) {
      past_out <- list("ucos" = out_ucos, "cos" = NULL)
    } else {
      past_out <- list("ucos" = NULL, "cos" = NULL)
    }
    # ----- Merge two ucos sets
    if (length(past_out$ucos$ucos_each) > 1) {
      if (merge_cos) {
        past_out <- merge_ucos(cb_obj, past_out,
          min_abs_corr = min_abs_corr,
          median_abs_corr = median_abs_corr,
          n_purity = n_purity,
          median_cos_abs_corr = median_cos_abs_corr,
          tol = tol
        )
      }
    }

    ############# - extract colocboost output - ####################
    # - colocalization results
    cb_obj$cb_model_para$weight_fudge_factor <- weight_fudge_factor
    cb_obj$cb_model_para$coverage <- coverage
    if (cb_obj$cb_model_para$M == 1) {
      # fixme for single iteration model
      cb_obj <- get_max_profile(cb_obj, check_null_max = 0.01, check_null_method = "profile")
    }
    cos_results <- get_cos_details(cb_obj, coloc_out = past_out$cos$cos, data_info = data_info)
    cb_output <- list(
      "vcp" = cos_results$vcp,
      "cos_details" = cos_results$cos_results,
      "data_info" = data_info,
      "model_info" = model_info
    )
    class(cb_output) <- "colocboost"

    ### - extract summary table
    focal_outcome_idx <- cb_obj$cb_model_para$focal_outcome_idx
    summary_table <- get_cos_summary(cb_output)
    cb_output <- c(cb_output, list(cos_summary = summary_table))
    cb_output <- cb_output[c("cos_summary", "vcp", "cos_details", "data_info", "model_info")]

    # - save model and all coloc and single information for diagnostic
    if (output_level != 1) {
      tmp <- get_full_output(
        cb_obj = cb_obj, past_out = past_out, variables = data_info$variables,
        cb_output = cb_output
      )
      if (output_level == 2) {
        cb_output <- c(cb_output, list("ucos_details" = tmp$ucos_details))
        cb_output <- cb_output[c("cos_summary", "vcp", "cos_details", "data_info", "model_info", "ucos_details")]
      } else {
        cb_output <- c(cb_output, list("ucos_details" = tmp$ucos_details))
        cb_output$diagnostic_details <- tmp[-1]
        cb_output <- cb_output[c("cos_summary", "vcp", "cos_details", "data_info", "model_info", "ucos_details", "diagnostic_details")]
      }
      # - if fine-boost, the summary table will be the summary of finemapping
      if (data_info$n_outcomes == 1) {
        cb_output <- cb_output[-match(c("cos_summary", "vcp", "cos_details"), names(cb_output))]
        remain_obj <- names(cb_output)
        if (!is.null(cb_output$ucos_details$ucos)) {
          cb_output$vpa <- apply(do.call(cbind, cb_output$ucos_details$ucos_weight), 1, function(w0) 1 - prod(1 - w0))
          names(cb_output$vpa) <- data_info$variables
          cb_output$ucos_summary <- get_ucos_summary(cb_output)
        } else {
          tmp <- list("vpa" = NULL, "ucos_summary" = NULL)
          cb_output <- c(cb_output, tmp)
        }
        cb_output <- cb_output[c("ucos_summary", "vpa", remain_obj)]
      }
    }
  }
  return(cb_output)
}
