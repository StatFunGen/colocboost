#' @rdname get_colocboost_summary
#'
#' @title Get summary tables from a ColocBoost output.
#'
#' @description `get_colocboost_summary` get colocalization and trait-specific summary table 
#' with or without the outcomes of interest.
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param summary_level When \code{summary_level = 1}, return basic summary table for colocalization results. See details in `get_ucos_summary` function when \code{summary_level = 2}.
#' @param outcome_names Optional vector of names of outcomes, which has the same order as Y in the original analysis.
#' @param interest_outcome Optional vector specifying a subset of outcomes from \code{outcome_names} to focus on. When provided, only colocalization events that include at least one of these outcomes will be returned.
#' @param region_name Optional character string. When provided, adds a column with this gene name to the output table for easier filtering in downstream analyses.
#' @param min_abs_corr_between_ucos Minimum absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.5.
#' @param median_abs_corr_between_ucos Median absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.8.
#' 
#' @return A list containing results from the ColocBoost analysis:
#' \itemize{
#'   \item When \code{summary_level = 1} (default):
#'   \itemize{
#'     \item \code{cos_summary}: A summary table for colocalization events with the following columns:
#'     \itemize{
#'       \item \code{focal_outcome}: The focal outcome being analyzed if exists. Otherwise, it is \code{FALSE}.
#'       \item \code{colocalized_outcomes}: Colocalized outcomes for colocalization confidence set (CoS)
#'       \item \code{cos_id}: Unique identifier for colocalization confidence set (CoS)
#'       \item \code{purity}: Minimum absolute correlation of variables within colocalization confidence set (CoS)
#'       \item \code{top_variable}: The variable with highest variant colocalization probability (VCP)
#'       \item \code{top_variable_vcp}: Variant colocalization probability for the top variable
#'       \item \code{cos_npc}: Normalized probability of colocalization
#'       \item \code{min_npc_outcome}: Minimum normalized probability of colocalized traits
#'       \item \code{n_variables}: Number of variables in colocalization confidence set (CoS)
#'       \item \code{colocalized_index}: Indices of colocalized variables
#'       \item \code{colocalized_variables}: List of colocalized variables
#'       \item \code{colocalized_variables_vcp}: Variant colocalization probabilities for all colocalized variables
#'     }
#'   }
#'   \item When \code{summary_level = 2}:
#'   \itemize{
#'     \item \code{cos_summary}: As described above
#'     \item \code{ucos_summary}: A summary table for trait-specific (uncolocalized) effects
#'   }
#'   \item When \code{summary_level = 3}:
#'   \itemize{
#'     \item \code{cos_summary}: As described above
#'     \item \code{ucos_summary}: A summary table for trait-specific (uncolocalized) effects
#'     \item \code{ambiguous_cos_summary}: A summary table for ambiguous colocalization events from trait-specific effects
#'   }
#' }
#' @details When \code{summary_level = 1}, additional details and examples are introduced in \code{\link{get_cos_summary}}. 
#' When \code{summary_level = 2} or \code{summary_level = 3}, additional details for trait-specific effects and ambiguous 
#' colocalization events are included. See \code{\link{get_ucos_summary}} for details on these tables.
#' 
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
#' true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
#' true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
#' true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' get_colocboost_summary(res)
#'
#' @source See detailed instructions in our tutorial portal:
#'  \url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
#' 
#' @family colocboost_inference
#' @export
#' 
get_colocboost_summary <- function(cb_output,
                                   summary_level = 1,
                                   outcome_names = NULL,
                                   interest_outcome = NULL,
                                   region_name = NULL,
                                   min_abs_corr_between_ucos = 0.5,
                                   median_abs_corr_between_ucos = 0.8) {

    if (!inherits(cb_output, "colocboost")) {
      stop("Input must from colocboost output!")
    }
    cos_summary <- get_cos_summary(cb_output, outcome_names, interest_outcome, region_name) 
    if (summary_level == 1) {
      return(list("cos_summary" = cos_summary))
    }

    if (summary_level == 2){
      ucos_summary <- get_ucos_summary(cb_output, outcome_names, region_name)
      return(list("cos_summary" = cos_summary, 
                  "ucos_summary" = ucos_summary))
    }

    if (summary_level == 3){
      ucos_summary <- get_ucos_summary(
        cb_output, outcome_names, region_name,
         ambiguous_cos = TRUE,
         min_abs_corr_between_ucos = min_abs_corr_between_ucos,
         median_abs_corr_between_ucos = median_abs_corr_between_ucos
      )
      return(list("cos_summary" = cos_summary, 
                  "ucos_summary" = ucos_summary$ucos_summary,
                  "ambiguous_cos_summary" = ucos_summary$ambiguous_cos_summary))
    }
                  
}


#' @rdname get_robust_colocalization
#'
#' @title Recalibrate and summarize robust colocalization events.
#'
#' @description `get_robust_colocalization` filters and recalibrates colocalization results by discarding weaker colocalization events or colocalized outcomes that do not meet evidence thresholds
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param cos_npc_cutoff Minimum threshold of normalized probability of colocalization (NPC) for CoS.
#' @param npc_outcome_cutoff Minimum threshold of normalized probability of colocalized traits in each CoS.
#' @param pvalue_cutoff Maximum threshold of marginal p-values of colocalized variants on colocalized traits in each CoS.
#' @param weight_fudge_factor The strength to integrate weight from different outcomes, default is 1.5
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated colocalization confidence sets (CoS) (default is 0.95).
#'
#' @return A \code{"colocboost"} object with some or all of the following elements:
#'
#' \item{cos_summary}{A summary table for colocalization events.}
#' \item{vcp}{The variable colocalized probability for each variable.}
#' \item{cos_details}{A object with all information for colocalization results.}
#' \item{data_info}{A object with detailed information from input data}
#' \item{model_info}{A object with detailed information for colocboost model}
#' \item{ucos_details}{Updated with trait-specific effects converted from filtered-out outcomes and single-trait CoS that pass the UCoS criterion. If `ucos_details` did not exist, it will be created.}
#'
#' @details 
#' This function filters colocalization confidence sets (CoS) based on evidence thresholds.
#' 
#' For single-trait CoS and filtered-out outcomes from multi-trait CoS, this function uses the stored 
#' `cos_cs_change` (profile likelihood change) and `check_null_max_ucos` (UCoS threshold) from the 
#' original colocboost output to determine if they meet the original UCoS criterion (cs_change >= check_null_max_ucos). 
#' Qualified UCoS candidates are merged into `ucos_details` (or `ucos_details` is created if it doesn't exist).
#' Single-trait CoS are removed from `cos_details` after conversion since they represent trait-specific 
#' rather than colocalized signals.
#'
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
#' true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
#' true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
#' true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' res$cos_details$cos$cos_index
#' filter_res <- get_robust_colocalization(res, cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2)
#' filter_res$cos_details$cos$cos_index
#' 
#' @source See detailed instructions in our tutorial portal: 
#' \url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
#'
#' @family colocboost_inference
#' @export
get_robust_colocalization <- function(cb_output,
                                      cos_npc_cutoff = 0.5,
                                      npc_outcome_cutoff = 0.2,
                                      pvalue_cutoff = NULL,
                                      weight_fudge_factor = 1.5,
                                      coverage = 0.95) {
  if (!inherits(cb_output, "colocboost")) {
    stop("Input must from colocboost object!")
  }


  if (is.null(cb_output$cos_details)) {
    warning("No colocalization results in this region!")
    return(cb_output)
  }

  if (npc_outcome_cutoff == 0 && cos_npc_cutoff == 0 && is.null(pvalue_cutoff)) {
    message("All possible colocalization events are reported regardless of their relative evidence compared to uncolocalized events (cos_npc_cutoff = 0 and npc_outcome_cutoff = 0).")
    return(cb_output)
  } else {
    if (is.null(pvalue_cutoff)) {
      message(paste0(
        "Extracting colocalization results with cos_npc_cutoff = ", cos_npc_cutoff, " and npc_outcome_cutoff = ", npc_outcome_cutoff, ".\n",
        "Keep only CoS with cos_npc >= ", cos_npc_cutoff, ". ",
        "For each CoS, keep the outcomes configurations that the npc_outcome >= ", npc_outcome_cutoff, "."
      ))
    } else {
      if (pvalue_cutoff > 1 | pvalue_cutoff < 0) {
        warning("Please check the pvalue cutoff in [0,1].")
        return(cb_output)
      }
      if (npc_outcome_cutoff == 0 && cos_npc_cutoff == 0) {
        message(paste0(
          "Extracting colocalization results with pvalue_cutoff = ", pvalue_cutoff, ".\n",
          "For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < ", pvalue_cutoff, "."
        ))
      } else {
        message(paste0(
          "Extracting colocalization results with pvalue_cutoff = ", pvalue_cutoff, ", cos_npc_cutoff = ", cos_npc_cutoff, ", and npc_outcome_cutoff = ", npc_outcome_cutoff, ".\n",
          "Keep only CoS with cos_npc >= ", cos_npc_cutoff, ". ",
          "For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < ", pvalue_cutoff, " and npc_outcome >", npc_outcome_cutoff, "."
        ))
      }
    }
  }

  remove_cos <- function(cb_output, remove_idx = NULL) {
    if (length(remove_idx) == 0) {
      return(cb_output)
    }
    ncos <- length(cb_output$cos_details$cos$cos_index)
    if (ncos == 0 | length(remove_idx) == ncos) {
      cb_output$vcp <- NULL
      cb_output$cos_details <- NULL
      cb_output <- c(cb_output, list(vcp = NULL, cos_details = NULL))
      return(cb_output)
    }
    cos_details <- cb_output$cos_details
    cos_details$cos_top_variables <- cos_details$cos_top_variables[-remove_idx, , drop = FALSE]
    cos_details$cos$cos_index <- cos_details$cos$cos_index[-remove_idx]
    cos_details$cos$cos_variables <- cos_details$cos$cos_variables[-remove_idx]
    cos_details$cos_outcomes$outcome_index <- cos_details$cos_outcomes$outcome_index[-remove_idx]
    cos_details$cos_outcomes$outcome_name <- cos_details$cos_outcomes$outcome_name[-remove_idx]
    cos_details$cos_vcp <- cos_details$cos_vcp[-remove_idx]
    cos_details$cos_weights <- cos_details$cos_weights[-remove_idx]
    cos_details$cos_npc <- cos_details$cos_npc[-remove_idx]
    cos_details$cos_min_npc_outcome <- cos_details$cos_min_npc_outcome[-remove_idx]
    cos_details$cos_purity$min_abs_cor <- as.matrix(cos_details$cos_purity$min_abs_cor)[-remove_idx, -remove_idx, drop = FALSE]
    cos_details$cos_purity$median_abs_cor <- as.matrix(cos_details$cos_purity$median_abs_cor)[-remove_idx, -remove_idx, drop = FALSE]
    cos_details$cos_purity$max_abs_cor <- as.matrix(cos_details$cos_purity$max_abs_cor)[-remove_idx, -remove_idx, drop = FALSE]
    # Also update cos_cs_change if it exists
    if (!is.null(cos_details$cos_cs_change)) {
      if (is.data.frame(cos_details$cos_cs_change)) {
        cos_details$cos_cs_change <- cos_details$cos_cs_change[-remove_idx, , drop = FALSE]
      }
    }
    # Ensure names align with remaining CoS structure after removal
    remaining_cos_names <- names(cos_details$cos$cos_index)
    if (length(remaining_cos_names) == length(cos_details$cos_npc)) {
      names(cos_details$cos_npc) <- remaining_cos_names
      names(cos_details$cos_min_npc_outcome) <- remaining_cos_names
    }
    vcp <- as.vector(1 - apply(1 - do.call(cbind, cos_details$cos_vcp), 1, prod))
    names(vcp) <- cb_output$data_info$variables
    cb_output$vcp <- vcp
    cb_output$cos_details <- cos_details
    return(cb_output)
  }
  get_npuc <- function(npc_outcome) {
    max_idx <- which.max(npc_outcome)
    npc_outcome[max_idx] * prod(1 - npc_outcome[-max_idx])
  }
  
  check_ucos_candidates <- function(cb_output, cos_details, ucos_candidates) {
    # Check all UCoS candidates (filtered-out outcomes and single-trait CoS) to see if they qualify as UCoS
    # using the original UCoS criterion: cs_change >= check_null_max_ucos
    # Returns list with: ucos_entries (list of qualified UCoS), n_checked, n_passed
    if (length(ucos_candidates) == 0) {
      return(list(ucos_entries = list(), n_checked = 0, n_passed = 0))
    }
    
    if (is.null(cos_details$cos_cs_change) || is.null(cos_details$check_null_max_ucos)) {
      warning(sprintf("WARNING: cos_cs_change or check_null_max_ucos not available - cannot check UCoS candidates. Only CoS filtering will work; UCoS results may be inaccurate."), call. = FALSE)
      return(list(ucos_entries = list(), n_checked = length(ucos_candidates), n_passed = 0))
    }
    
    cs_change_df <- cos_details$cos_cs_change
    check_null_max_ucos <- cos_details$check_null_max_ucos
    ucos_entries <- list()
    n_passed <- 0
    for (cand in ucos_candidates) {
      original_cos_name <- cand$original_cos_name
      outcome_idx <- cand$outcome_idx
      outcome_name <- cand$outcome_name
      
      # Try to get delta (cs_change) for this outcome from the original CoS
      delta <- NULL
      if (original_cos_name %in% rownames(cs_change_df)) {
        if (outcome_name %in% colnames(cs_change_df)) {
          delta <- cs_change_df[original_cos_name, outcome_name]
        } else if (as.character(outcome_idx) %in% colnames(cs_change_df)) {
          delta <- cs_change_df[original_cos_name, as.character(outcome_idx)]
        } else if (outcome_idx <= ncol(cs_change_df)) {
          delta <- cs_change_df[original_cos_name, outcome_idx]
        }
      }
      
      # Get threshold
      threshold <- NULL
      if (outcome_name %in% names(check_null_max_ucos)) {
        threshold <- check_null_max_ucos[outcome_name]
      } else if (as.character(outcome_idx) %in% names(check_null_max_ucos)) {
        threshold <- check_null_max_ucos[as.character(outcome_idx)]
      } else if (outcome_idx <= length(check_null_max_ucos)) {
        threshold <- check_null_max_ucos[outcome_idx]
      }
      
      # Check if it passes UCoS criterion (cs_change >= check_null_max_ucos)
      if (!is.null(delta) && !is.null(threshold) && !is.na(delta) && !is.na(threshold) && delta >= threshold) {
        n_passed <- n_passed + 1
        
        ucos_entries[[length(ucos_entries) + 1]] <- list(
          ucos_index = cand$cos_index,
          ucos_variables = cand$cos_variables,
          outcome_idx = outcome_idx,
          outcome_name = outcome_name,
          cos_weight = cand$cos_weight,
          purity_min = cand$purity_min,
          purity_median = cand$purity_median,
          purity_max = cand$purity_max,
          top_variable_idx = cand$top_variable_idx,
          top_variable_name = cand$top_variable_name,
          original_cos_name = original_cos_name
        )
      }
    }
    
    
    return(list(ucos_entries = ucos_entries, n_checked = length(ucos_candidates), n_passed = n_passed))
  }

  cos_details <- cb_output$cos_details
  # Preserve original cos_cs_change and check_null_max_ucos before filtering
  # (these contain information about all original CoS, including ones that will be removed)
  original_cos_cs_change <- cos_details$cos_cs_change
  original_check_null_max_ucos <- cos_details$check_null_max_ucos
  
  # Capture initial counts for summary
  initial_cos_count <- length(cos_details$cos$cos_index)
  initial_ucos_count <- if (!is.null(cb_output$ucos_details) && length(cb_output$ucos_details$ucos$ucos_index) > 0) {
    length(cb_output$ucos_details$ucos$ucos_index)
  } else {
    0
  }
  
  coloc_outcome_index <- coloc_outcome <- list()
  colocset_names <- cos_min_npc_outcome <- cos_npc <- c()
  
  # Track filtered-out outcomes from multi-trait CoS for potential UCoS conversion
  filtered_outcome_candidates <- list()
  
  for (i in 1:length(cos_details$cos$cos_index)) {
    cos_npc_config <- cos_details$cos_outcomes_npc[[i]]
    npc_outcome <- cos_npc_config$npc_outcome
    pos_pass <- which(npc_outcome >= npc_outcome_cutoff)
    if (!is.null(pvalue_cutoff)) {
      cos_tmp <- cos_details$cos$cos_index[[i]]
      cos_trait <- cos_details$cos_outcomes$outcome_index[[i]]
      minPV <- sapply(cos_trait, function(tmp) {
        z <- cb_output$data_info$z[[tmp]][cos_tmp]
        pv <- pchisq(z^2, 1, lower.tail = FALSE)
        min(pv)
      })
      pos_pass_pvalue <- which(minPV <= pvalue_cutoff)
      if (length(pos_pass_pvalue) == 0) {
        pos_pass <- NULL
      } else {
        pos_pass_pvalue <- match(cos_trait[pos_pass_pvalue], cos_npc_config$outcomes_index)
        pos_pass <- intersect(pos_pass_pvalue, pos_pass)
      }
    }
    if (length(pos_pass) == 0) {
      # CoS completely removed - add ALL outcomes from this CoS as UCOS candidates
      original_cos_name <- names(cos_details$cos$cos_index)[i]
      original_outcomes <- cos_npc_config$outcomes_index
      outcome_names_all <- rownames(cos_npc_config)
      
      # Get original weight for each outcome
      original_weight <- cos_details$cos_weights[[i]]
      
      # Get purity (same for all outcomes in the CoS)
      cos_purity_min <- if (!is.null(cos_details$cos_purity$min_abs_cor)) {
        as.matrix(cos_details$cos_purity$min_abs_cor)[i, i]
      } else { NA }
      cos_purity_median <- if (!is.null(cos_details$cos_purity$median_abs_cor)) {
        as.matrix(cos_details$cos_purity$median_abs_cor)[i, i]
      } else { NA }
      cos_purity_max <- if (!is.null(cos_details$cos_purity$max_abs_cor)) {
        as.matrix(cos_details$cos_purity$max_abs_cor)[i, i]
      } else { NA }
      
      # Get top variable
      top_var_idx <- if (!is.null(cos_details$cos_top_variables)) {
        cos_details$cos_top_variables$top_index[i]
      } else { 0 }
      top_var_name <- if (!is.null(cos_details$cos_top_variables) && top_var_idx > 0) {
        cos_details$cos_top_variables$top_variables[i]
      } else { NA }
      
      # Add each outcome from the removed CoS as a UCoS candidate
      for (j in seq_along(original_outcomes)) {
        outcome_idx <- original_outcomes[j]
        outcome_name <- outcome_names_all[j]
        
        # Get outcome-specific weight
        outcome_colname <- paste0("outcome", outcome_idx)
        if (outcome_colname %in% colnames(original_weight)) {
          outcome_weight <- original_weight[, outcome_colname, drop = FALSE]
        } else {
          # Fallback: use first column if specific outcome column not found
          outcome_weight <- original_weight[, 1, drop = FALSE]
        }
        
        filtered_outcome_candidates[[length(filtered_outcome_candidates) + 1]] <- list(
          original_cos_idx = i,
          original_cos_name = original_cos_name,
          outcome_idx = outcome_idx,
          outcome_name = outcome_name,
          cos_index = cos_details$cos$cos_index[[i]],
          cos_variables = cos_details$cos$cos_variables[[i]],
          cos_weight = outcome_weight,
          purity_min = cos_purity_min,
          purity_median = cos_purity_median,
          purity_max = cos_purity_max,
          top_variable_idx = top_var_idx,
          top_variable_name = top_var_name,
          is_removed_cos = TRUE  # Mark as from a completely removed CoS
        )
      }
      
      coloc_outcome_index[[i]] <- 0
      coloc_outcome[[i]] <- 0
      cos_min_npc_outcome[i] <- 0
      cos_npc[i] <- 0
      colocset_names[i] <- paste0("remove", i)
    } else {
      # Track outcomes that were filtered out (if original CoS had multiple outcomes)
      original_outcomes <- cos_npc_config$outcomes_index
      outcomes_that_passed <- cos_npc_config$outcomes_index[pos_pass]
      if (length(original_outcomes) > length(pos_pass) && length(pos_pass) > 0) {
        # Some outcomes were filtered out from a multi-trait CoS
        filtered_outcomes <- setdiff(original_outcomes, outcomes_that_passed)
        original_cos_name <- names(cos_details$cos$cos_index)[i]
        for (outcome_filt in filtered_outcomes) {
          # Get outcome name by matching outcome index in the original CoS outcomes
          # This ensures we get the correct outcome name for this specific filtered outcome
          pos_in_original <- which(cos_details$cos_outcomes$outcome_index[[i]] == outcome_filt)
          if (length(pos_in_original) > 0) {
            outcome_name_filt <- cos_details$cos_outcomes$outcome_name[[i]][pos_in_original[1]]
          } else {
            # Fallback: use cos_npc_config (should match, but this is backup)
            outcome_name_filt <- rownames(cos_npc_config)[which(cos_npc_config$outcomes_index == outcome_filt)[1]]
          }
          # Store original weight and purity info for this filtered outcome
          # Get outcome-specific weight from original CoS weights
          original_weight <- cos_details$cos_weights[[i]]
          outcome_colname <- paste0("outcome", outcome_filt)
          if (outcome_colname %in% colnames(original_weight)) {
            outcome_weight <- original_weight[, outcome_colname, drop = FALSE]
          } else {
            # Fallback: use first column if specific outcome column not found
            outcome_weight <- original_weight[, 1, drop = FALSE]
          }
          
          # Get purity (will be same for all outcomes in the CoS)
          cos_purity_min <- if (!is.null(cos_details$cos_purity$min_abs_cor)) {
            as.matrix(cos_details$cos_purity$min_abs_cor)[i, i]
          } else { NA }
          cos_purity_median <- if (!is.null(cos_details$cos_purity$median_abs_cor)) {
            as.matrix(cos_details$cos_purity$median_abs_cor)[i, i]
          } else { NA }
          cos_purity_max <- if (!is.null(cos_details$cos_purity$max_abs_cor)) {
            as.matrix(cos_details$cos_purity$max_abs_cor)[i, i]
          } else { NA }
          
          # Get top variable
          top_var_idx <- if (!is.null(cos_details$cos_top_variables)) {
            cos_details$cos_top_variables$top_index[i]
          } else { 0 }
          top_var_name <- if (!is.null(cos_details$cos_top_variables) && top_var_idx > 0) {
            cos_details$cos_top_variables$top_variables[i]
          } else { NA }
          
          filtered_outcome_candidates[[length(filtered_outcome_candidates) + 1]] <- list(
            original_cos_idx = i,
            original_cos_name = original_cos_name,
            outcome_idx = outcome_filt,
            outcome_name = outcome_name_filt,
            cos_index = cos_details$cos$cos_index[[i]],
            cos_variables = cos_details$cos$cos_variables[[i]],
            cos_weight = outcome_weight,
            purity_min = cos_purity_min,
            purity_median = cos_purity_median,
            purity_max = cos_purity_max,
            top_variable_idx = top_var_idx,
            top_variable_name = top_var_name
          )
        }
      }
      
      cos_min_npc_outcome[i] <- min(npc_outcome[pos_pass])
      if (length(pos_pass) > 1){
        cos_npc[i] <- 1 - get_npuc(npc_outcome[pos_pass])
      } else {
        cos_npc[i] <- 0 # Single-trait CoS get cos_npc=0 (they're not colocalized)
      }
      coloc_outcome_index[[i]] <- sort(cos_npc_config$outcomes_index[pos_pass])
      coloc_outcome[[i]] <- rownames(cos_npc_config)[pos_pass]
      colocset_names[i] <- paste0("cos", i, ":", paste0(paste0("y", coloc_outcome_index[[i]]), collapse = "_"))
      if (grepl("merged", names(cos_details$cos$cos_index)[i])) {
        colocset_names[i] <- paste0(colocset_names[i], ":merged")
      }
    }
  }
  names(coloc_outcome) <- names(coloc_outcome_index) <- names(cos_min_npc_outcome) <- names(cos_npc) <- colocset_names
  cos_details$cos_outcomes <- list("outcome_index" = coloc_outcome_index, "outcome_name" = coloc_outcome)
  cos_details$cos_min_npc_outcome <- cos_min_npc_outcome
  cos_details$cos_npc <- cos_npc

  # - VCP
  cos_weights <- lapply(1:length(cos_details$cos_outcomes$outcome_index), function(idx) {
    w <- cos_details$cos_weights[[idx]]
    config_idx <- cos_details$cos_outcomes$outcome_index[[idx]]
    if (length(config_idx) == 1) {
      if (config_idx == 0) {
        return(matrix(1, nrow = nrow(cos_details$cos_weights[[idx]])))
      }
    }
    w_outcome <- colnames(w)
    config_outcome <- paste0("outcome", config_idx)
    pos <- which(w_outcome %in% config_outcome)
    w[, pos, drop = FALSE]
  })
  cos_details$cos_weights <- cos_weights
  int_weight <- lapply(cos_weights, get_integrated_weight, weight_fudge_factor = weight_fudge_factor)
  # Assign names from colocset_names (entries marked "remove" will be filtered out later)
  names(int_weight) <- names(cos_weights) <- colocset_names
  vcp <- as.vector(1 - apply(1 - do.call(cbind, int_weight), 1, prod))
  names(vcp) <- cb_output$data_info$variables
  cb_output$vcp <- vcp
  cos_details$cos_vcp <- int_weight

  # - resummary results
  cos_re_idx <- lapply(int_weight, function(w) {
    unlist(get_in_cos(w, coverage = coverage))
  })
  cos_re_var <- lapply(cos_re_idx, function(idx) {
    cb_output$data_info$variables[idx]
  })
  # Preserve names from int_weight
  names(cos_re_idx) <- names(cos_re_var) <- names(int_weight)
  coloc_csets <- list("cos_index" = cos_re_idx, "cos_variables" = cos_re_var)
  cos_details$cos <- coloc_csets

  # - hits variables in each csets
  coloc_hits <- coloc_hits_variablenames <- coloc_hits_names <- c()
  for (i in 1:length(int_weight)) {
    inw <- int_weight[[i]]
    if (length(unique(inw)) == 1) {
      coloc_hits <- c(coloc_hits, 0)
      coloc_hits_variablenames <- c(coloc_hits_variablenames, 0)
      coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
    } else {
      pp <- which(inw == max(inw))
      coloc_hits <- c(coloc_hits, pp)
      coloc_hits_variablenames <- c(coloc_hits_variablenames, cb_output$data_info$variables[pp])
      if (length(pp) == 1) {
        coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
      } else {
        coloc_hits_names <- c(coloc_hits_names, paste0(names(int_weight)[i], ".", 1:length(pp)))
      }
    }
  }
  coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_variablenames)
  rownames(coloc_hits) <- coloc_hits_names
  cos_details$cos_top_variables <- coloc_hits
  cb_output$cos_details <- cos_details

  # Remove CoS that have no outcomes passing npc_outcome_cutoff
  remove <- grep("remove", colocset_names)
  cb_output <- remove_cos(cb_output, remove_idx = remove)
  cos_details <- cb_output$cos_details

  # Remove multi-trait CoS that do not pass cos_npc_cutoff
  # Note: Single-trait CoS get cos_npc=0 by design (they're not colocalized).
  # We skip cos_npc_cutoff for single-trait CoS since they already passed npc_outcome_cutoff
  # which is the appropriate filter for single-trait evidence strength.
  if (length(cos_details$cos_npc) > 0) {
    cos_npc_vals <- cos_details$cos_npc
    # Ensure names align with current cos_details structure
    cos_names_actual <- names(cos_details$cos$cos_index)
    if (length(cos_npc_vals) == length(cos_names_actual)) {
      names(cos_npc_vals) <- cos_names_actual
    }
    n_outcome_before <- sapply(cos_details$cos_outcomes$outcome_index, length)
    
    # Only apply cos_npc_cutoff to multi-trait CoS (single-trait have cos_npc=0 and shouldn't be filtered this way)
    is_multi <- n_outcome_before > 1
    if (any(is_multi)) {
      remove_multi <- which(is_multi & cos_details$cos_npc < cos_npc_cutoff)
      remove <- remove_multi
      if (length(remove) > 0) {
        # Collect all outcomes from removed multi-trait CoS as UCoS candidates
        # (before we remove them from cos_details)
        for (rm_idx in remove) {
          cos_name_rm <- names(cos_npc_vals)[rm_idx]
          outcome_indices_rm <- cos_details$cos_outcomes$outcome_index[[rm_idx]]
          outcome_names_rm <- cos_details$cos_outcomes$outcome_name[[rm_idx]]
          
          # Each outcome from this removed CoS becomes a candidate
          for (j in seq_along(outcome_indices_rm)) {
            outcome_idx_rm <- outcome_indices_rm[j]
            outcome_name_rm <- if (length(outcome_names_rm) >= j) outcome_names_rm[j] else NA
            
            # Get outcome-specific weight
            original_weight_rm <- cos_details$cos_weights[[rm_idx]]
            outcome_colname_rm <- paste0("outcome", outcome_idx_rm)
            if (outcome_colname_rm %in% colnames(original_weight_rm)) {
              outcome_weight_rm <- original_weight_rm[, outcome_colname_rm, drop = FALSE]
            } else {
              outcome_weight_rm <- original_weight_rm[, 1, drop = FALSE]
            }
            
            # Get purity
            cos_purity_min_rm <- if (!is.null(cos_details$cos_purity$min_abs_cor)) {
              as.matrix(cos_details$cos_purity$min_abs_cor)[rm_idx, rm_idx]
            } else { NA }
            cos_purity_median_rm <- if (!is.null(cos_details$cos_purity$median_abs_cor)) {
              as.matrix(cos_details$cos_purity$median_abs_cor)[rm_idx, rm_idx]
            } else { NA }
            cos_purity_max_rm <- if (!is.null(cos_details$cos_purity$max_abs_cor)) {
              as.matrix(cos_details$cos_purity$max_abs_cor)[rm_idx, rm_idx]
            } else { NA }
            
            # Get top variable
            top_var_idx_rm <- if (!is.null(cos_details$cos_top_variables)) {
              cos_details$cos_top_variables$top_index[rm_idx]
            } else { 0 }
            top_var_name_rm <- if (!is.null(cos_details$cos_top_variables) && top_var_idx_rm > 0) {
              cos_details$cos_top_variables$top_variables[rm_idx]
            } else { NA }
            
            filtered_outcome_candidates[[length(filtered_outcome_candidates) + 1]] <- list(
              original_cos_idx = rm_idx,
              original_cos_name = cos_name_rm,
              outcome_idx = outcome_idx_rm,
              outcome_name = outcome_name_rm,
              cos_index = cos_details$cos$cos_index[[rm_idx]],
              cos_variables = cos_details$cos$cos_variables[[rm_idx]],
              cos_weight = outcome_weight_rm,
              purity_min = cos_purity_min_rm,
              purity_median = cos_purity_median_rm,
              purity_max = cos_purity_max_rm,
              top_variable_idx = top_var_idx_rm,
              top_variable_name = top_var_name_rm,
              is_removed_cos = TRUE  # Mark as coming from a completely removed CoS
            )
          }
        }
      }
      
    } else {
      # All remaining are single-trait
      remove <- integer(0)
    }
  } else {
    remove <- integer(0)
  }
  cb_output <- remove_cos(cb_output, remove_idx = remove)
  cos_details <- cb_output$cos_details

  # Collect ALL single-trait CoS as UCoS candidates and remove them from cos_details
  # Single-trait signals cannot be colocalization sets by definition
  if (!is.null(cos_details)) {
    n_outcome_check <- sapply(cos_details$cos_outcomes$outcome_index, length)
    single_check <- which(n_outcome_check == 1)
    if (length(single_check) > 0) {
      # Extract single-trait CoS info and add to candidates list BEFORE removal
      cos_names <- names(cos_details$cos$cos_index)
      for (i in seq_along(single_check)) {
        idx <- single_check[i]
        cos_name <- cos_names[idx]
        outcome_idx <- cos_details$cos_outcomes$outcome_index[[idx]][[1]]
        outcome_name <- cos_details$cos_outcomes$outcome_name[[idx]][[1]]
        
        # Get weight (should already be filtered to single outcome)
        outcome_weight <- cos_details$cos_weights[[idx]]
        
        # Get purity
        cos_purity_min <- as.matrix(cos_details$cos_purity$min_abs_cor)[idx, idx]
        cos_purity_median <- as.matrix(cos_details$cos_purity$median_abs_cor)[idx, idx]
        cos_purity_max <- as.matrix(cos_details$cos_purity$max_abs_cor)[idx, idx]
        
        # Get top variable
        top_var_idx <- cos_details$cos_top_variables$top_index[idx]
        top_var_name <- cos_details$cos_top_variables$top_variables[idx]
        
        filtered_outcome_candidates[[length(filtered_outcome_candidates) + 1]] <- list(
          original_cos_idx = idx,
          original_cos_name = cos_name,
          outcome_idx = outcome_idx,
          outcome_name = outcome_name,
          cos_index = cos_details$cos$cos_index[[idx]],
          cos_variables = cos_details$cos$cos_variables[[idx]],
          cos_weight = outcome_weight,
          purity_min = cos_purity_min,
          purity_median = cos_purity_median,
          purity_max = cos_purity_max,
          top_variable_idx = top_var_idx,
          top_variable_name = top_var_name,
          is_single_trait = TRUE  # Mark as single-trait CoS
        )
      }
      
      # Remove ALL single-trait CoS from cos_details (they cannot be colocalization sets)
      cb_output <- remove_cos(cb_output, remove_idx = single_check)
      cos_details <- cb_output$cos_details
    }
  }

  # Check all UCoS candidates (filtered-out outcomes, single-trait CoS, and outcomes from removed CoS)
  # for UCoS eligibility using the original criterion (cs_change >= check_null_max_ucos)
  ucos_result <- if (length(filtered_outcome_candidates) > 0) {
    # Create a temporary cos_details object with original cs_change info for UCoS checking
    # (cos_details might be NULL or filtered, but we need the original cs_change data)
    temp_cos_details <- if (!is.null(cos_details)) cos_details else list()
    if (!is.null(original_cos_cs_change)) {
      temp_cos_details$cos_cs_change <- original_cos_cs_change
    }
    if (!is.null(original_check_null_max_ucos)) {
      temp_cos_details$check_null_max_ucos <- original_check_null_max_ucos
    }
    check_ucos_candidates(cb_output, temp_cos_details, filtered_outcome_candidates)
  } else {
    list(ucos_entries = list(), n_checked = 0, n_passed = 0)
  }
  all_qualified_ucos <- ucos_result$ucos_entries

  # Merge all qualified UCoS candidates into ucos_details
  if (length(all_qualified_ucos) > 0) {
    total_ucos <- length(all_qualified_ucos)
    
    # Get existing ucos_details or create new structure
    if (is.null(cb_output$ucos_details)) {
      # Create new ucos_details structure
      existing_count <- 0
      cb_output$ucos_details <- list(
        "ucos" = list("ucos_index" = list(), "ucos_variables" = list()),
        "ucos_outcomes" = list("outcome_index" = list(), "outcome_name" = list()),
        "ucos_weight" = list(),
        "ucos_purity" = list("min_abs_cor" = matrix(0,0,0), "median_abs_cor" = matrix(0,0,0), "max_abs_cor" = matrix(0,0,0)),
        "ucos_top_variables" = data.frame("top_index" = integer(0), "top_variables" = character(0))
      )
    } else {
      existing_count <- length(cb_output$ucos_details$ucos$ucos_index)
    }
    
    # Append new UCOS entries
    new_ucos_indices <- lapply(all_qualified_ucos, function(e) e$ucos_index)
    new_ucos_variables <- lapply(all_qualified_ucos, function(e) e$ucos_variables)
    new_outcome_indices <- lapply(all_qualified_ucos, function(e) list(e$outcome_idx))
    new_outcome_names <- lapply(all_qualified_ucos, function(e) list(e$outcome_name))
    new_weights <- lapply(all_qualified_ucos, function(e) e$cos_weight)
    new_purity_min <- sapply(all_qualified_ucos, function(e) ifelse(is.na(e$purity_min), 0, e$purity_min))
    new_purity_median <- sapply(all_qualified_ucos, function(e) ifelse(is.na(e$purity_median), 0, e$purity_median))
    new_purity_max <- sapply(all_qualified_ucos, function(e) ifelse(is.na(e$purity_max), 0, e$purity_max))
    new_top_indices <- sapply(all_qualified_ucos, function(e) e$top_variable_idx)
    new_top_variables <- sapply(all_qualified_ucos, function(e) e$top_variable_name)
    
    # Generate names for new UCOS
    new_ucos_names <- paste0("ucos", (existing_count + 1):(existing_count + total_ucos))
    
    # Merge into existing structure
    cb_output$ucos_details$ucos$ucos_index <- c(cb_output$ucos_details$ucos$ucos_index, new_ucos_indices)
    cb_output$ucos_details$ucos$ucos_variables <- c(cb_output$ucos_details$ucos$ucos_variables, new_ucos_variables)
    cb_output$ucos_details$ucos_outcomes$outcome_index <- c(cb_output$ucos_details$ucos_outcomes$outcome_index, new_outcome_indices)
    cb_output$ucos_details$ucos_outcomes$outcome_name <- c(cb_output$ucos_details$ucos_outcomes$outcome_name, new_outcome_names)
    cb_output$ucos_details$ucos_weight <- c(cb_output$ucos_details$ucos_weight, new_weights)
    
    # Update purity matrices (append as diagonal blocks)
    old_n <- nrow(cb_output$ucos_details$ucos_purity$min_abs_cor)
    if (old_n == 0) {
      # Create new matrices
      cb_output$ucos_details$ucos_purity$min_abs_cor <- diag(new_purity_min, nrow = total_ucos, ncol = total_ucos)
      cb_output$ucos_details$ucos_purity$median_abs_cor <- diag(new_purity_median, nrow = total_ucos, ncol = total_ucos)
      cb_output$ucos_details$ucos_purity$max_abs_cor <- diag(new_purity_max, nrow = total_ucos, ncol = total_ucos)
    } else {
      # Append to existing matrices
      total_n <- old_n + total_ucos
      new_min_mat <- diag(new_purity_min, nrow = total_ucos, ncol = total_ucos)
      new_median_mat <- diag(new_purity_median, nrow = total_ucos, ncol = total_ucos)
      new_max_mat <- diag(new_purity_max, nrow = total_ucos, ncol = total_ucos)
      
      # Combine old and new matrices
      min_combined <- matrix(0, nrow = total_n, ncol = total_n)
      min_combined[1:old_n, 1:old_n] <- cb_output$ucos_details$ucos_purity$min_abs_cor
      min_combined[(old_n+1):total_n, (old_n+1):total_n] <- new_min_mat
      
      median_combined <- matrix(0, nrow = total_n, ncol = total_n)
      median_combined[1:old_n, 1:old_n] <- cb_output$ucos_details$ucos_purity$median_abs_cor
      median_combined[(old_n+1):total_n, (old_n+1):total_n] <- new_median_mat
      
      max_combined <- matrix(0, nrow = total_n, ncol = total_n)
      max_combined[1:old_n, 1:old_n] <- cb_output$ucos_details$ucos_purity$max_abs_cor
      max_combined[(old_n+1):total_n, (old_n+1):total_n] <- new_max_mat
      
      cb_output$ucos_details$ucos_purity$min_abs_cor <- min_combined
      cb_output$ucos_details$ucos_purity$median_abs_cor <- median_combined
      cb_output$ucos_details$ucos_purity$max_abs_cor <- max_combined
    }
    
    # Append to top_variables data.frame
    new_top_vars_df <- data.frame(
      "top_index" = new_top_indices,
      "top_variables" = new_top_variables,
      row.names = new_ucos_names
    )
    cb_output$ucos_details$ucos_top_variables <- rbind(cb_output$ucos_details$ucos_top_variables, new_top_vars_df)
    
    # Update names for all lists
    if (existing_count > 0) {
      existing_names <- names(cb_output$ucos_details$ucos_weight)[1:existing_count]
      all_ucos_names <- c(existing_names, new_ucos_names)
    } else {
      all_ucos_names <- new_ucos_names
    }
    names(cb_output$ucos_details$ucos$ucos_index) <- names(cb_output$ucos_details$ucos$ucos_variables) <-
      names(cb_output$ucos_details$ucos_outcomes$outcome_index) <- names(cb_output$ucos_details$ucos_outcomes$outcome_name) <-
      names(cb_output$ucos_details$ucos_weight) <- all_ucos_names
  }
  
  # Validate that only multi-trait CoS remain in cos_details
  if (!is.null(cos_details) && length(cos_details$cos$cos_index) > 0) {
    n_outcome_final <- sapply(cos_details$cos_outcomes$outcome_index, length)
    n_single_remaining <- sum(n_outcome_final == 1)
    if (n_single_remaining > 0) {
      stop("Internal error: Single-trait CoS should have been removed from cos_details")
    }
  }

  # - refine and output
  class(cb_output) <- "colocboost"
  cb_output$cos_summary <- get_cos_summary(cb_output)

  # Calculate final counts
  final_cos_count <- if (!is.null(cb_output$cos_details) && length(cb_output$cos_details$cos$cos_index) > 0) {
    length(cb_output$cos_details$cos$cos_index)
  } else {
    0
  }
  final_ucos_count <- if (!is.null(cb_output$ucos_details) && length(cb_output$ucos_details$ucos$ucos_index) > 0) {
    length(cb_output$ucos_details$ucos$ucos_index)
  } else {
    0
  }
  
  # Summary message
  message(sprintf(
    "Robust filtering: Started with %d CoS, %d UCoS. Finished with %d CoS, %d UCoS.",
    initial_cos_count, initial_ucos_count, final_cos_count, final_ucos_count
  ))

  return(cb_output)
}

#' @rdname get_ambiguous_colocalization
#'
#' @title Get ambiguous colocalization events from trait-specific (uncolocalized) effects.
#'
#' @description `get_ambiguous_colocalization` get the colocalization by discarding the weaker colocalization events or colocalized outcomes
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param min_abs_corr_between_ucos Minimum absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.5.
#' @param median_abs_corr_between_ucos Median absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.8.
#' @param tol A small, non-negative number specifying the convergence tolerance for checking the overlap of the variables in different sets.
#'
#' @return A \code{"colocboost"} object of colocboost output with additional elements:
#' \item{ambiguous_cos}{If exists, a list of ambiguous trait-specific (uncolocalized) effects.}
#'
#' @examples
#' data(Ambiguous_Colocalization)
#' test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
#' res <- get_ambiguous_colocalization(test_colocboost_results)
#' names(res$ambiguous_cos)
#' 
#' @source See detailed instructions in our tutorial portal: 
#' \url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
#'
#' @family colocboost_inference
#' @export
get_ambiguous_colocalization <- function(cb_output, 
                                         min_abs_corr_between_ucos = 0.5, 
                                         median_abs_corr_between_ucos = 0.8,
                                         tol = 1e-9) {

    if (!inherits(cb_output, "colocboost")) {
      stop("Input must from colocboost output!")
    }

    if (!("ucos_details" %in% names(cb_output))) {
      warning(
        "Since you want to extract ambiguous colocalization from trait-specific (uncolocalized) sets,",
        " but there is no output of ucos_details from colocboost.\n",
        " Please run colocboost model with output_level=2!"
      )
      return(cb_output)
    }

    if (is.null(cb_output$ucos_details)){
      message("No trait-specific (uncolocalized) effects in this region!")
      return(cb_output)
    }

    ucos_details <- cb_output$ucos_details
    nucos <- length(ucos_details$ucos$ucos_index)
    if (nucos == 1) {
      message("Only one trait-specific (uncolocalized) effect in this region!")
      return(cb_output)
    }

    # Function to merge ambiguous ucos
    merge_sets <- function(vec) {
      split_lists <- lapply(vec, function(x) as.numeric(unlist(strsplit(x, ";"))))
      result <- list()
      while (length(split_lists) > 0) {
        current <- split_lists[[1]]
        split_lists <- split_lists[-1]
        repeat {
          overlap_index <- NULL
          for (i in seq_along(split_lists)) {
            if (length(intersect(current, split_lists[[i]])) > 0) {
              overlap_index <- i
              break
            }
          }
          if (!is.null(overlap_index)) {
            current <- union(current, split_lists[[overlap_index]])
            split_lists <- split_lists[-overlap_index]
          } else {
            break
          }
        }
        result <- c(result, list(paste(sort(current), collapse = ";")))
      }
      return(result)
    }


    purity <- ucos_details$ucos_purity
    min_abs_cor <- purity$min_abs_cor
    median_abs_cor <- purity$median_abs_cor
    max_abs_cor <- purity$max_abs_cor
    is_ambiguous <- (min_abs_cor > min_abs_corr_between_ucos) * 
      (abs(max_abs_cor - 1) < tol) * 
      (median_abs_cor > median_abs_corr_between_ucos)
    diag(is_ambiguous) <- 0 # no need to check within ucos

    if (sum(is_ambiguous) == 0){
      message("No ambiguous colocalization events!")
      return(cb_output)
    } else {
      message("There exists the ambiguous colocalization events from trait-specific effects. Extracting!")
    }

    temp <- sapply(1:nrow(is_ambiguous), function(x) {
      tt <- c(x, which(is_ambiguous[x, ] != 0))
      return(paste0(sort(tt), collapse = ";"))
    })
    temp <- merge_sets(temp)
    potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
    potential_merged <- potential_merged[which(sapply(potential_merged, length) >= 2)]

   ambiguous_events <- list()
    ambiguous_ucos_names <- c()
    for (i in 1:length(potential_merged)) {
      idx <- potential_merged[[i]]
      test_outcome <- unique(unlist(ucos_details$ucos_outcomes$outcome_index[idx]))
      if (length(test_outcome) == 1) next
      ambiguous_ucos_names[i] <- paste0(names(ucos_details$ucos$ucos_index)[idx], collapse = ";")
      tmp <- list(
        ambiguous_cos = list(
          ucos_index = ucos_details$ucos$ucos_index[idx],
          ucos_variables = ucos_details$ucos$ucos_variables[idx]
        ),
        ambiguous_cos_overlap = list(
          ucos_index = Reduce(intersect, ucos_details$ucos$ucos_index[idx]),
          ucos_variables = Reduce(intersect, ucos_details$ucos$ucos_variables[idx])
        ),
        ambiguous_cos_union = list(
          ucos_index = Reduce(union, ucos_details$ucos$ucos_index[idx]),
          ucos_variables = Reduce(union, ucos_details$ucos$ucos_variables[idx])
        ),
        ambiguous_cos_outcomes = list(
          outcome_idx = unique(unlist(ucos_details$ucos_outcomes$outcome_index[idx])),
          outcome_name = unique(unlist(ucos_details$ucos_outcomes$outcome_name[idx]))
        ),
        ambigous_cos_weight = ucos_details$ucos_weight[idx],     
        ambigous_cos_purity = list(            
          min_abs_cor = min_abs_cor[idx, idx],
          median_abs_cor = median_abs_cor[idx, idx],
          max_abs_cor = max_abs_cor[idx, idx]
        )
      )
      w <- tmp$ambigous_cos_weight
      w <- do.call(cbind, w)
      tmp$recalibrated_cos_vcp <- get_integrated_weight(w)
      tmp$recalibrated_cos <- list(
        "cos_index" = unlist(get_in_cos(tmp$recalibrated_cos_vcp)),
        "cos_variables" = lapply(unlist(get_in_cos(tmp$recalibrated_cos_vcp)), function(idx) cb_output$data_info$variables[idx])
      )
      ambiguous_events[[i]] <- tmp
    }
    names(ambiguous_events) <- ambiguous_ucos_names
    message(paste("There are", length(ambiguous_events), "ambiguous trait-specific effects."))

    cb_output$ambiguous_cos <- ambiguous_events
    return(cb_output)
          
}




#' @rdname get_cos_summary
#'
#' @title Get colocalization summary table from a ColocBoost output.
#'
#' @description `get_cos_summary` get the colocalization summary table with or without the outcomes of interest.
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param outcome_names Optional vector of names of outcomes, which has the same order as Y in the original analysis.
#' @param interest_outcome Optional vector specifying a subset of outcomes from \code{outcome_names} to focus on. When provided, only colocalization events that include at least one of these outcomes will be returned.
#' @param region_name Optional character string. When provided, adds a column with this gene name to the output table for easier filtering in downstream analyses.
#'
#' @return A summary table for colocalization events with the following columns:
#' \item{focal_outcome}{The focal outcome being analyzed if exists. Otherwise, it is \code{FALSE}.}
#' \item{colocalized_outcomes}{Colocalized outcomes for colocalization confidence set (CoS) }
#' \item{cos_id}{Unique identifier for colocalization confidence set (CoS) }
#' \item{purity}{Minimum absolute correlation of variables with in colocalization confidence set (CoS) }
#' \item{top_variable}{The variable with highest variant colocalization probability (VCP) }
#' \item{top_variable_vcp}{Variant colocalization probability for the top variable}
#' \item{cos_npc}{Normalized probability of colocalization}
#' \item{min_npc_outcome}{Minimum normalized probability of colocalized traits}
#' \item{n_variables}{Number of variables in colocalization confidence set (CoS)}
#' \item{colocalized_index}{Indices of colocalized variables}
#' \item{colocalized_variables}{List of colocalized variables}
#' \item{colocalized_variables_vcp}{Variant colocalization probabilities for all colocalized variables}
#'
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
#' true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
#' true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
#' true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' get_cos_summary(res)
#'
#' @source See detailed instructions in our tutorial portal:
#'  \url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
#' 
#' @keywords colocboost_inference
#' @family colocboost_utilities
#' @export
get_cos_summary <- function(cb_output,
                            outcome_names = NULL,
                            interest_outcome = NULL,
                            region_name = NULL) {
  if (!inherits(cb_output, "colocboost")) {
    stop("Input must from colocboost output!")
  }

  coloc_csets <- cb_output$cos_details$cos$cos_index
  if (length(coloc_csets) != 0) {
    analysis_outcome <- cb_output$data_info$outcome_info$outcome_names
    if (!is.null(outcome_names)) {
      analysis_outcome <- outcome_names
    }
    coloc_outcome <- lapply(cb_output$cos_details$cos_outcomes$outcome_index, function(idx) analysis_outcome[idx])
    coloc_sets <- cb_output$cos_details$cos$cos_index
    if (!is.null(cb_output$cos_warnings)) {
      message(cb_output$cos_warnings$warning_message)
    }
    vcp <- as.numeric(cb_output$vcp)

    summary_table <- matrix(NA, nrow = length(coloc_sets), ncol = 12)
    colnames(summary_table) <- c(
      "focal_outcome", "colocalized_outcomes", "cos_id", "purity",
      "top_variable", "top_variable_vcp", "cos_npc", "min_npc_outcome", "n_variables",
      "colocalized_index", "colocalized_variables", "colocalized_variables_vcp"
    )
    summary_table <- as.data.frame(summary_table)
    summary_table[, 1] <- FALSE
    summary_table[, 2] <- unlist(sapply(coloc_outcome, function(tmp) paste0(tmp, collapse = "; ")))
    summary_table[, 3] <- names(coloc_sets)
    summary_table[, 4] <- as.numeric(diag(as.matrix(cb_output$cos_details$cos_purity$min_abs_cor)))
    summary_table[, 5] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) tmp[1]))
    summary_table[, 6] <- sapply(coloc_sets, function(tmp) max(vcp[tmp]))
    summary_table[, 7] <- round(as.numeric(cb_output$cos_details$cos_npc), 4)
    summary_table[, 8] <- round(as.numeric(cb_output$cos_details$cos_min_npc_outcome), 4)
    summary_table[, 9] <- as.numeric(sapply(coloc_sets, length))
    summary_table[, 10] <- unlist(sapply(coloc_sets, function(tmp) paste0(tmp, collapse = "; ")))
    summary_table[, 11] <- unlist(sapply(cb_output$cos_details$cos$cos_variables, function(tmp) paste0(tmp, collapse = "; ")))
    summary_table[, 12] <- unlist(sapply(coloc_sets, function(tmp) paste0(vcp[tmp], collapse = "; ")))
    if (!is.null(region_name)) {
      summary_table$region_name <- region_name
    }
    # - if focal colocalization
    focal_outcome_idx <- which(cb_output$data_info$outcome_info$is_focal)
    if (length(focal_outcome_idx) != 0) {
      focal_outcome <- analysis_outcome[focal_outcome_idx]
      tmp <- sapply(focal_outcome, function(tmp) grep(paste0(tmp, "\\b"), analysis_outcome))
      if.focal <- sapply(coloc_outcome, function(cp) {
        tt <- sapply(focal_outcome, function(tmp) grep(paste0(tmp, "\\b"), cp))
        all(sapply(tt, length) != 0)
      })
      summary_table$focal_outcome <- ifelse(if.focal, focal_outcome, FALSE)
      summary_table <- summary_table[order(summary_table$focal_outcome == "FALSE"), ]
      if (sum(if.focal) == 0) {
        warning("No colocalization with focal outcomes.")
      }
    }
    # - if extract only interest outcome colocalization
    if (!is.null(interest_outcome)) {
      tmp <- sapply(interest_outcome, function(tmp) grep(paste0(tmp, "\\b"), analysis_outcome))
      if (all(sapply(tmp, length) != 0)) {
        if.interest <- sapply(coloc_outcome, function(cp) {
          tt <- sapply(interest_outcome, function(tmp) grep(paste0(tmp, "\\b"), cp))
          all(sapply(tt, length) != 0)
        })
        summary_table$interest_outcome <- paste0(interest_outcome, collapse = "; ")
        summary_table <- summary_table[which(if.interest), ]
        if (sum(if.interest) == 0) {
          warning("No colocalization with interest outcomes.")
        }
      } else {
        warning("Interest outcome is not in the analysis outcomes, please check.")
      }
    }
  } else {
    summary_table <- NULL
  }
  return(summary_table)
}


#' @rdname get_ucos_summary
#'
#' @title Get trait-specific summary table from a ColocBoost output.
#'
#' @description `get_ucos_summary` produces a trait-specific summary table for uncolocalized (single-trait) 
#' associations from ColocBoost results. This is particularly useful for examining trait-specific signals 
#' or for summarizing results from single-trait FineBoost analyses.
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param outcome_names Optional vector of names of outcomes, which has the same order as Y in the original analysis.
#' @param region_name Optional character string. When provided, adds a column with this gene name to the output table for easier filtering in downstream analyses.
#' @param ambiguous_cos Logical indicating whether to include ambiguous colocalization events. The default is FALSE.
#' @param min_abs_corr_between_ucos Minimum absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.5.
#' @param median_abs_corr_between_ucos Median absolute correlation for variants across two trait-specific (uncolocalized) effects to be considered colocalized. The default is 0.8.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{ucos_summary}: A summary table for trait-specific, uncolocalized associations with the following columns:
#'   \itemize{
#'     \item \code{outcomes}: Outcome being analyzed
#'     \item \code{ucos_id}: Unique identifier for trait-specific confidence sets
#'     \item \code{purity}: Minimum absolute correlation of variables within trait-specific confidence sets
#'     \item \code{top_variable}: The variable with highest variant-level probability of association (VPA)
#'     \item \code{top_variable_vpa}: Variant-level probability of association (VPA) for the top variable
#'     \item \code{ucos_npc}: Normalized probability of causal association for the trait-specific confidence set
#'     \item \code{n_variables}: Number of variables in trait-specific confidence set
#'     \item \code{ucos_index}: Indices of variables in the trait-specific confidence set
#'     \item \code{ucos_variables}: List of variables in the trait-specific confidence set
#'     \item \code{ucos_variables_vpa}: Variant-level probability of association (VPA) for all variables in the confidence set
#'     \item \code{region_name}: Region name if provided through the region_name parameter
#'   }
#'   \item \code{ambiguous_cos_summary}: A summary table for ambiguous colocalization events with the following columns:
#'   \itemize{
#'     \item \code{outcomes}: Outcome in the ambiguous colocalization event
#'     \item \code{ucos_id}: Unique identifiers for the ambiguous event 
#'     \item \code{min_between_purity}: Minimum absolute correlation between variables across trait-specific sets in the ambiguous event
#'     \item \code{median_between_purity}: Median absolute correlation between variables across trait-specific sets in the ambiguous event
#'     \item \code{overlap_idx}: Indices of variables that overlap between ambiguous trait-specific sets
#'     \item \code{overlap_variables}: Names of variables that overlap between ambiguous trait-specific sets
#'     \item \code{n_recalibrated_variables}: Number of variables in the recalibrated colocalization set from an ambiguous event
#'     \item \code{recalibrated_index}: Indices of variables in the recalibrated colocalization set from an ambiguous event
#'     \item \code{recalibrated_variables}: Names of variables in the recalibrated colocalization set from an ambiguous event
#'     \item \code{recalibrated_variables_vcp}: Variant colocalization probabilities for recalibrated variables from an ambiguous event
#'     \item \code{region_name}: Region name if provided through the region_name parameter
#'   }
#' }
#'
#' @examples
#' # colocboost example with single trait analysis
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 1  # Only one trait for single-trait analysis
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects the trait
#' true_beta[80, 1] <- 0.2 # SNP11 also affects the trait but with lower effect
#' Y <- X %*% true_beta + rnorm(N, 0, 1)
#' res <- colocboost(X = X, Y = Y, output_level = 2)
#' # Get the trait-specifc effect summary
#' get_ucos_summary(res)
#' 
#' @source See detailed instructions in our tutorial portal: 
#' \url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
#'
#' @keywords colocboost_inference
#' @family colocboost_utilities
#' @export
get_ucos_summary <- function(cb_output, outcome_names = NULL, region_name = NULL, 
                             ambiguous_cos = FALSE,
                             min_abs_corr_between_ucos = 0.5, 
                             median_abs_corr_between_ucos = 0.8) {

  # - check input
  if (!inherits(cb_output, "colocboost")) {
    stop("Input must from colocboost object!")
  }

  if (!("ucos_details" %in% names(cb_output))) {
    warning(
      "Since you want to extract trait-specific (uncolocalized) sets,",
      " but there is no output of ucos_details from colocboost.\n",
      " Please run colocboost model with output_level=2!"
    )
    return(NULL)
  }

  if (is.null(cb_output$ucos_details)){
    message("No trait-specific (uncolocalized) effects in this region!")
    return(NULL)
  }

  specific_cs <- cb_output$ucos_details
  cs_outcome <- cb_output$data_info$outcome_info$outcome_names
  if (!is.null(outcome_names)) {
    cs_outcome <- outcome_names
  }
  vpa <- as.numeric(cb_output$vpa)
  if (length(vpa) == 0){
    w <- do.call(cbind, specific_cs$ucos_weight)
    vpa <- as.vector(1 - apply(1 - w, 1, prod))
  }

  summary_table <- matrix(NA, nrow = length(specific_cs$ucos$ucos_index), ncol = 9)
  colnames(summary_table) <- c(
    "outcomes", "ucos_id", "purity",
    "top_variable", "top_variable_vpa", "n_variables", "ucos_index",
    "ucos_variables", "ucos_variables_vpa"
  )
  summary_table <- as.data.frame(summary_table)
  summary_table[, 1] <- cs_outcome[unlist(specific_cs$ucos_outcomes$outcome_index)]
  summary_table[, 2] <- names(specific_cs$ucos$ucos_index)
  summary_table[, 3] <- as.numeric(diag(as.matrix(specific_cs$ucos_purity$min_abs_cor)))
  summary_table[, 4] <- unlist(sapply(specific_cs$ucos$ucos_variables, function(tmp) tmp[1]))
  summary_table[, 5] <- sapply(specific_cs$ucos$ucos_index, function(tmp) max(vpa[tmp]))
  summary_table[, 6] <- as.numeric(sapply(specific_cs$ucos$ucos_index, length))
  summary_table[, 7] <- unlist(sapply(specific_cs$ucos$ucos_index, function(tmp) paste0(tmp, collapse = "; ")))
  summary_table[, 8] <- unlist(sapply(specific_cs$ucos$ucos_variables, function(tmp) paste0(tmp, collapse = "; ")))
  summary_table[, 9] <- unlist(sapply(specific_cs$ucos$ucos_index, function(tmp) paste0(vpa[tmp], collapse = "; ")))
  if (!is.null(region_name)) {
    summary_table$region_name <- region_name
  }
  summary_table <- as.data.frame(summary_table)
  if (!ambiguous_cos) return(summary_table)


  # advanced summary for ambiguous colocalization at post-processing
  output_summary <- list(
    ucos_summary = summary_table,
    ambiguous_cos_summary = NULL
  )
  test <- get_ambiguous_colocalization(
    cb_output,
    min_abs_corr_between_ucos = min_abs_corr_between_ucos,
    median_abs_corr_between_ucos = median_abs_corr_between_ucos
  )
  if (length(test$ambiguous_cos) == 0) return(output_summary)

  ambiguous_results <- test$ambiguous_cos
  ambiguous_summary <- matrix(NA, nrow = length(ambiguous_results), ncol = 10)
  colnames(ambiguous_summary) <- c(
    "outcomes", "ucos_id", "min_between_purity", "median_between_purity",
    "overlap_idx", "overlap_variables", "n_recalibrated_variables",
    "recalibrated_index", "recalibrated_variables", "recalibrated_variables_vcp"
  )
  ambiguous_summary <- as.data.frame(ambiguous_summary)
  ambiguous_summary[, 1] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$ambiguous_cos_outcomes$outcome_name, collapse = "; ")))
  ambiguous_summary[, 2] <- names(ambiguous_results)
  ambiguous_summary[, 3] <- sapply(ambiguous_results, function(tmp) max(tmp$ambigous_cos_purity$min_abs_cor[lower.tri(tmp$ambigous_cos_purity$min_abs_cor)]) )
  ambiguous_summary[, 4] <- sapply(ambiguous_results, function(tmp) max(tmp$ambigous_cos_purity$median_abs_cor[lower.tri(tmp$ambigous_cos_purity$median_abs_cor)]) )
  ambiguous_summary[, 5] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$ambiguous_cos_overlap$ucos_index, collapse = "; ")))
  ambiguous_summary[, 6] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$ambiguous_cos_overlap$ucos_variables, collapse = "; ")))
  ambiguous_summary[, 7] <- as.numeric(sapply(ambiguous_results, function(tmp) length(tmp$recalibrated_cos$cos_index)))
  ambiguous_summary[, 8] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$recalibrated_cos$cos_index, collapse = "; ")))
  ambiguous_summary[, 9] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$recalibrated_cos$cos_variables, collapse = "; ")))
  ambiguous_summary[, 10] <- unlist(sapply(ambiguous_results, function(tmp) paste0(tmp$recalibrated_cos_vcp[tmp$recalibrated_cos$cos_index], collapse = "; ")))
  if (!is.null(region_name)) {
    ambiguous_summary$region_name <- region_name
  }
  output_summary$ambiguous_cos_summary <- as.data.frame(ambiguous_summary)

  return(output_summary)
}

#' Extract CoS at different coverage
#' 
#' @description `get_cos` extracts colocalization confidence sets (CoS) at different coverage levels 
#' from ColocBoost results. When genotype data (X) or correlation matrix (Xcorr) is provided, it 
#' can also calculate and filter CoS based on purity statistics, ensuring that variants within 
#' each CoS are sufficiently correlated.
#' 
#' @param cb_output Output object from `colocboost` analysis
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated colocalization confidence sets (CoS) (default is 0.95).
#' @param X Genotype matrix of values of the p variables. Used to compute correlations if Xcorr is not provided.
#' @param Xcorr Correlation matrix of correlations between variables. Alternative to X.
#' @param n_purity The maximum number of CoS variables used in calculating the correlation (\dQuote{purity}) statistics. 
#' @param min_abs_corr The minimum absolute correlation value of variants in a CoS to be considered pass (\dQuote{purity}) statistics.
#' @param median_abs_corr The median absolute correlation value of variants in a CoS to be considered pass (\dQuote{purity}) statistics.
#' When the number of variables included in the CoS is greater than this number, the CoS variables are randomly subsampled.
#'
#' @return A list of indices of variables in each CoS.
#'
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
#' true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
#' true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
#' true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' get_cos(res, coverage = 0.99, X = X)
#' get_cos(res, coverage = 0.99, X = X, min_abs_corr = 0.95)
#'
#' @family colocboost_utilities
#' @export
get_cos <- function(cb_output, coverage = 0.95, X = NULL, Xcorr = NULL, n_purity = 100, min_abs_corr = 0.5, median_abs_corr = NULL) {

  if (is.null(cb_output$cos_details$cos)){
    warning("No colocalization results in this region!")
    return(list("cos" = NULL, "cos_purity" = NULL))
  }

  # Refine CoS based on different coverage
  cos_vcp <- cb_output$cos_details$cos_vcp
  cos_diff_coverage <- lapply(cos_vcp, function(w) {
    unlist(get_in_cos(w, coverage = coverage))
  })
  names(cos_diff_coverage) <- paste0(names(cos_diff_coverage), "_coverage_", coverage)

  # Check purity if X or Xcorr exist
  cos_purity <- NULL
  if (!is.null(X) || !is.null(Xcorr)){
    cos_purity <- get_cos_purity(cos_diff_coverage, X = X, Xcorr = Xcorr, n_purity = n_purity)
    # Extract within-CoS purity
    ncos <- nrow(cos_purity[[1]])
    within_purity <- matrix(NA, nrow = ncos, ncol = 3)
    colnames(within_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
    rownames(within_purity) <- names(cos_diff_coverage)
    within_purity[, 1] <- diag(cos_purity$min_abs_cor)
    within_purity[, 2] <- diag(cos_purity$max_abs_cor)
    within_purity[, 3] <- diag(cos_purity$median_abs_cor)
    # Check purity
    if (is.null(median_abs_corr)) {
      is_pure <- which(within_purity[, 1] >= min_abs_corr)
    } else {
      is_pure <- which(within_purity[, 1] >= min_abs_corr | within_purity[, 3] >= median_abs_corr)
    }
    # Filter impured CoS
    if (length(is_pure) == 0) {
      cos_diff_coverage <- NULL
      cos_purity <- NULL
    } else if (length(is_pure) < ncos){
      cos_diff_coverage <- cos_diff_coverage[is_pure]
      cos_purity <- lapply(cos_purity, function(tmp) tmp[is_pure, is_pure, drop = FALSE] )
    }
  }
  cos_refined <- list("cos" = cos_diff_coverage, "cos_purity" = cos_purity)
  return(cos_refined)
}


#' Get integrated weight from different outcomes
#' @keywords cb_get_functions
#' @noRd
get_integrated_weight <- function(weights, weight_fudge_factor = 1.5) {
  av <- apply(weights, 1, function(w) prod(w^(weight_fudge_factor / ncol(weights))))
  return(av / sum(av))
}

#' Extract CoS based on coverage
#' @keywords cb_get_functions
#' @noRd
get_in_cos <- function(weights, coverage = 0.95) {
  # Get indices in decreasing weight order
  temp <- order(weights, decreasing = TRUE)
  weights_ordered <- weights[temp]
  # Find position where cumsum exceeds coverage
  pos <- min(which(cumsum(weights_ordered) >= coverage))
  # Get the actual threshold weight value
  weight_thresh <- weights_ordered[pos]
  # Find all indices with weights >= threshold
  indices <- which(weights >= weight_thresh)
  # Order these indices by their weights in decreasing order
  csets <- indices[order(weights[indices], decreasing = TRUE)]
  return(list(csets))
}


#' Calculate purity within and in-between CoS 
#'
#' @description Calculate purity statistics between all pairs of colocalization confidence sets (CoS)
#'
#' @param cos List of variables in CoS
#' @param X Genotype matrix of values of the p variables. Used to compute correlations if Xcorr is not provided.
#' @param Xcorr Correlation matrix of correlations between variables. Alternative to X.
#' @param n_purity The maximum number of CoS variables used in calculating the correlation (\dQuote{purity}) statistics. 
#' When the number of variables included in the CoS is greater than this number, the CoS variables are randomly subsampled.
#' 
#' @return A list containing three matrices (min_abs_cor, max_abs_cor, median_abs_cor) with
#'   purity statistics for all pairs of CoS. Diagonal elements represent within-CoS purity.
#'
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 
#' true_beta[10, 2] <- 0.4 
#' true_beta[50, 2] <- 0.3 
#' true_beta[80, 3] <- 0.6 
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' cos_res <- get_cos(res, coverage = 0.8)
#' get_cos_purity(cos_res$cos, X = X)
#'
#' @family colocboost_utilities
#' @export
get_cos_purity <- function(cos, X = NULL, Xcorr = NULL, n_purity = 100) {
  # Check inputs
  if (is.null(cos) || length(cos) == 0) return(NULL)
  if (is.null(X) && is.null(Xcorr)) stop("Either X or Xcorr must be provided")
  if (is.numeric(cos)) cos <- list(cos)

  # Get CoS names
  cos_names <- names(cos)
  ncos <- length(cos)
  if (is.null(cos_names)) cos_names <- paste0("cos_", 1:ncos)
  
  # If only one CoS, just return the purity as a matrix
  if (ncos == 1) {
    purity_stats <- get_purity(cos[[1]], X = X, Xcorr = Xcorr, n = n_purity)
    cos_purity <- lapply(1:3, function(ii) {
      mat <- matrix(purity_stats[ii], 1, 1)
      rownames(mat) <- colnames(mat) <- cos_names
      return(mat)
    })
    names(cos_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
    return(cos_purity)
  }
  
  # Initialize empty matrices for purity values
  empty_matrix <- matrix(NA, ncos, ncos)
  colnames(empty_matrix) <- rownames(empty_matrix) <- cos_names
  # Initialize purity matrices with diagonal values (within-CoS purity)
  cos_purity <- lapply(1:3, function(ii) empty_matrix )
  # Compute within-CoS purity
  for (i in 1:ncos){
    purity_stats <- get_purity(cos[[i]], X = X, Xcorr = Xcorr, n = n_purity)
    for (ii in 1:3) {
        cos_purity[[ii]][i, i] <- purity_stats[ii]
      }
  }
  # Compute between-CoS purity for every pair
  for (i in 1:(ncos - 1)) {
    for (j in (i + 1):ncos) {
      cos1 <- cos[[i]]
      cos2 <- cos[[j]]
      purity_stats <- get_between_purity(cos1, cos2, X = X, Xcorr = Xcorr)
      # Update all three purity matrices
      for (ii in 1:3) {
        cos_purity[[ii]][i, j] <- cos_purity[[ii]][j, i] <- purity_stats[ii]
      }
    }
  }
  names(cos_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
  
  return(cos_purity)
}


