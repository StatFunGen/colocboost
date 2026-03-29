#' @title Set of internal functions to refrine colocalization confidence sets
#'
#' @description
#' The `colocboost_refine_cos` functions serves as a summary for the following two refine functions.
#'
#' @details
#' The following functions are included in this set:
#' `merge_cos_ucos` merge a trait-specific confidence set CS into a colocalization set CoS.
#' `merge_ucos` merge two trait-specific confidence sets.
#'
#' These functions are not exported individually and are accessed via `colocboost_refine_cos`.
#'
#' @rdname colocboost_refine_cos
#' @keywords cb_refine_cos
#' @noRd
merge_cos_ucos <- function(cb_obj, out_cos, out_ucos, coverage = 0.95,
                           min_abs_corr = 0.5, tol = 1e-9,
                           median_cos_abs_corr = 0.8) {
  change_obj_each <- out_ucos$change_obj_each
  coloc_sets <- out_cos$cos$cos
  ucos_each <- out_ucos$ucos_each

  # - remove overlap between coloc_sets and single_sets
  is_overlap <- is_highLD <- c()
  for (i in 1:length(coloc_sets)) {
    # cset1 <- coloc_sets[[i]]
    coloc_outcome <- out_cos$cos$coloc_outcomes[[i]]
    ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_cos$cos$avWeight[[i]])))
    pos_name <- order(ttmp)
    out_cos$cos$avWeight[[i]] <- out_cos$cos$avWeight[[i]][, pos_name]
    for (j in 1:length(ucos_each)) {
      cset2 <- ucos_each[[j]]
      fine_outcome <- out_ucos$ucos_outcome[j]

      if (fine_outcome %in% coloc_outcome) {
        # - if fine_y in coloc_y, we check only the overlap
        pp <- which(coloc_outcome == fine_outcome)
        w_tmp <- out_cos$cos$avWeight[[i]][, pp]
        cset1 <- get_in_cos(w_tmp, coverage = coverage)[[1]]
        # - addhoc: if median_cos_abs_corr > 0.8, remove single sets
        X_dict <- cb_obj$cb_data$dict[fine_outcome]
        res <- get_between_purity(cset1, cset2,
          X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
          Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
          miss_idx = cb_obj$cb_data$data[[fine_outcome]]$variable_miss,
          P = cb_obj$cb_model_para$P
        )
        # is.between <- length(intersect(cset1, cset2)) != 0
        is.between <- (abs(res[2] - 1) < tol)
        if (is.between) {
          is_overlap <- c(is_overlap, j)
        } else {
          is.higLD <- res[2] > median_cos_abs_corr
          if (is.higLD) {
            is_highLD <- c(is_highLD, j)
          }
        }
      } else if (!(fine_outcome %in% coloc_outcome)) {
        # - if fine_y not in coloc_y, we check overlap and also min_purity
        change_obj_coloc <- out_cos$cos$cs_change
        cset1 <- coloc_sets[[i]]
        res <- list()
        for (ii in 1:cb_obj$cb_model_para$L) {
          X_dict <- cb_obj$cb_data$dict[ii]
          res[[ii]] <- get_between_purity(cset1, cset2,
            X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
            Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
            miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
            P = cb_obj$cb_model_para$P
          )
        }
        res <- Reduce(pmax, res)
        min_between <- res[1]
        max_between <- res[2]
        ave_between <- res[3]
        is.between <- ((min_between > median_cos_abs_corr) & (abs(max_between - 1) < tol))
        # is.between <- (min_between>min_abs_corr) * (abs(max_between-1)<tol) * (ave_between>median_cos_abs_corr)
        if (is.between) {
          is_overlap <- c(is_overlap, j)
          # -- add weight
          add_avW <- out_ucos$avW_ucos_each[, j]
          out_cos$cos$avWeight[[i]] <- cbind(add_avW, out_cos$cos$avWeight[[i]])
          colnames(out_cos$cos$avWeight[[i]])[1] <- paste0("outcome", fine_outcome)
          ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames(out_cos$cos$avWeight[[i]])))
          pos_name <- order(ttmp)
          out_cos$cos$avWeight[[i]] <- out_cos$cos$avWeight[[i]][, pos_name]
          change_obj_coloc_i <- change_obj_coloc[i, ]
          change_obj_each_j <- change_obj_each[j, ]
          out_cos$cos$cs_change[i, ] <- pmax(change_obj_coloc_i, change_obj_each_j)
          coloc_outcome <- sort(c(coloc_outcome, fine_outcome))
          out_cos$cos$coloc_outcomes[[i]] <- coloc_outcome
        }
      }
    }
  }

  is_overlap <- unique(c(is_overlap, is_highLD))
  if (length(is_overlap) != 0) {
    if (length(is_overlap) == length(ucos_each)) {
      out_ucos$ucos_each <- NULL
      out_ucos$avW_ucos_each <- NULL
      out_ucos$change_obj_each <- NULL
      out_ucos$purity_each <- NULL
    } else {
      out_ucos$ucos_each <- ucos_each[-is_overlap]
      out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[, -is_overlap, drop = FALSE]
      out_ucos$change_obj_each <- change_obj_each[-is_overlap, , drop = FALSE]
      out_ucos$purity_each <- out_ucos$purity_each[-is_overlap, , drop = FALSE]
    }
  }
  if (length(out_ucos$ucos_each) == 0) {
    out_ucos$ucos_each <- NULL
    out_ucos$avW_ucos_each <- NULL
    out_ucos$change_obj_each <- NULL
    out_ucos$purity_each <- NULL
  }
  ll <- list("ucos" = out_ucos, "cos" = out_cos)
  return(ll)
}

#' @importFrom stats na.omit
merge_ucos <- function(cb_obj, past_out,
                       min_abs_corr = 0.5,
                       median_abs_corr = NULL,
                       n_purity = 100,
                       median_cos_abs_corr = 0.8,
                       tol = 1e-9) {
  out_ucos <- past_out$ucos
  out_cos <- past_out$cos
  ucos_each <- out_ucos$ucos_each
  change_obj_each <- out_ucos$change_obj_each

  # calculate between purity
  ncsets <- length(ucos_each)
  min_between <- max_between <- ave_between <- matrix(0, nrow = ncsets, ncol = ncsets)
  for (i.between in 1:(ncsets - 1)) {
    for (j.between in (i.between + 1):ncsets) {
      cset1 <- ucos_each[[i.between]]
      cset2 <- ucos_each[[j.between]]
      y.i <- out_ucos$ucos_outcome[i.between]
      y.j <- out_ucos$ucos_outcome[j.between]
      if (y.i == y.j) {
        next
      }
      yy <- c(y.i, y.j)
      res <- list()
      flag <- 1
      for (ii in yy) {
        X_dict <- cb_obj$cb_data$dict[ii]
        res[[flag]] <- get_between_purity(cset1, cset2,
          X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
          Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
          miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
          P = cb_obj$cb_model_para$P
        )
        flag <- flag + 1
      }
      if (min_abs_corr == 0) {
        res <- Reduce(pmin, res)
      } else {
        res <- Reduce(pmax, res)
      }
      min_between[i.between, j.between] <- min_between[j.between, i.between] <- res[1]
      max_between[i.between, j.between] <- max_between[j.between, i.between] <- res[2]
      ave_between[i.between, j.between] <- ave_between[j.between, i.between] <- res[3]
    }
  }
  # is.between <- (min_between>min_abs_corr) * (abs(max_between-1)<tol) * (ave_between>median_cos_abs_corr)
  is.between <- ((min_between > median_cos_abs_corr) & (abs(max_between - 1) < tol))
  if (sum(is.between) != 0) {
    temp <- sapply(1:nrow(is.between), function(x) {
      tt <- c(x, which(is.between[x, ] != 0))
      return(paste0(sort(tt), collapse = ";"))
    })
    temp <- merge_sets(temp)
    potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
    potential_merged <- potential_merged[which(sapply(potential_merged, length) >= 2)]
    coloc_sets_merged <- avWeight_merged <-
      cs_change_merged <- coloc_outcomes_merged <- list()
    is_merged <- c()
    for (i.m in 1:length(potential_merged)) {
      temp_set <- as.numeric(potential_merged[[i.m]])
      # refine avWeight
      merged <- out_ucos$avW_ucos_each[, temp_set]
      unique_coloc_outcomes <- as.numeric(gsub(".*Y(\\d+).*", "\\1", colnames(merged)))
      if (length(unique(unique_coloc_outcomes)) == 1) next
      # define merged set
      coloc_sets_merged <- c(coloc_sets_merged, list(unique(unlist(ucos_each[temp_set]))))
      colnames(merged) <- paste0("outcome", unique_coloc_outcomes)
      coloc_outcomes_merged <- c(
        coloc_outcomes_merged,
        list(unique(sort(unique_coloc_outcomes)))
      )
      # colnames(temp) <- unique_coloc_outcomes
      avWeight_merged <- c(avWeight_merged, list(merged))
      # refine cs_change
      change_cs_tmp <- change_obj_each[temp_set, , drop = FALSE]
      cs_change_merged <- c(
        cs_change_merged,
        list(apply(change_cs_tmp, 2, max))
      )
      is_merged <- c(is_merged, temp_set)
    }

    if (length(is_merged) != 0) {
      # --- check merged coloc set purity
      purity <- vector(mode = "list", length = length(coloc_sets_merged))
      for (ee in 1:length(coloc_sets_merged)) {
        coloc_t <- coloc_outcomes_merged[[ee]]
        p_tmp <- c()
        for (i3 in coloc_t) {
          pos <- coloc_sets_merged[[ee]]
          X_dict <- cb_obj$cb_data$dict[i3]
          if (!is.null(cb_obj$cb_data$data[[X_dict]]$XtX)) {
            pos <- match(pos, setdiff(1:cb_obj$cb_model_para$P, cb_obj$cb_data$data[[i3]]$variable_miss))
            # - it could happen since merge_cos
            if (sum(is.na(pos)) != 0) {
              pos <- as.numeric(na.omit(pos))
            }
          }
          tmp <- matrix(get_purity(pos,
            X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
            Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
            N = cb_obj$cb_data$data[[i3]]$N, n = n_purity
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
        # add coloc
        merged_colocsets <- coloc_sets_merged[is_pure]
        names(merged_colocsets) <- paste0("merged_cos", 1:length(is_pure))
        out_cos$cos$cos <- c(
          out_cos$cos$cos,
          merged_colocsets
        )
        out_cos$cos$purity <- rbind(
          out_cos$cos$purity,
          purity_all[is_pure, ]
        )
        out_cos$cos$evidence_strength <- c(
          out_cos$cos$evidence_strength,
          rep(0.95, length(is_pure))
        )
        cs_change_tmp <- do.call(rbind, cs_change_merged)[is_pure, , drop = FALSE]
        colnames(cs_change_tmp) <- paste0("change_obj_", 1:cb_obj$cb_model_para$L)
        rownames(cs_change_tmp) <- names(merged_colocsets)
        out_cos$cos$cs_change <- rbind(
          out_cos$cos$cs_change,
          cs_change_tmp
        )
        out_cos$cos$avWeight <- c(
          out_cos$cos$avWeight,
          avWeight_merged[is_pure]
        )
        out_cos$cos$coloc_outcomes <- c(
          out_cos$cos$coloc_outcomes,
          coloc_outcomes_merged[is_pure]
        )
      }

      # - remove single
      if (length(is_merged) == length(ucos_each)) {
        out_ucos$ucos_each <- NULL
        out_ucos$avW_ucos_each <- NULL
        out_ucos$change_obj_each <- NULL
        out_ucos$purity_each <- NULL
      } else {
        out_ucos$ucos_each <- ucos_each[-is_merged]
        out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[, -is_merged, drop = FALSE]
        out_ucos$change_obj_each <- change_obj_each[-is_merged, , drop = FALSE]
        out_ucos$purity_each <- out_ucos$purity_each[-is_merged, , drop = FALSE]
      }
    }
  }
  ll <- list("ucos" = out_ucos, "cos" = out_cos)
  return(ll)
}

#' @title Set of internal utils functions
#'
#' @description
#' The `colocboost_utils` functions serves as a summary for the following two refine functions.
#'
#' @details
#' The following functions are included in this set:
#' `merge_cos_ucos` merge a trait-specific confidence set CS into a colocalization set CoS.
#' `merge_ucos` merge two trait-specific confidence sets.
#'
#' These functions are not exported individually and are accessed via `colocboost_utils`.
#'
#' @rdname colocboost_utils
#' @keywords cb_utils
#' @noRd
get_vcp <- function(past_out, P) {
  if (!is.null(past_out$cos$cos$cos)) {
    avW_coloc_vcp <- sapply(past_out$cos$cos$avWeight, get_integrated_weight)
  } else {
    avW_coloc_vcp <- NULL
  }
  all_weight <- avW_coloc_vcp
  if (length(all_weight) == P) {
    all_weight <- as.matrix(unlist(all_weight))
  }
  if (!is.null(all_weight)) {
    all_weight <- apply(all_weight, 2, as.numeric)
    all_weight <- as.matrix(all_weight)
    vcp <- as.vector(1 - apply(1 - all_weight, 1, prod))
  } else {
    vcp <- rep(0, P)
  }
  return(vcp)
}


get_pip <- function(past_out, R, P) {
  if (length(past_out$cos$cos$cos) != 0) {
    av_coloc <- do.call(cbind, past_out$cos$cos$avWeight)
  } else {
    av_coloc <- NULL
  }
  if (length(past_out$ucos$ucos_each) != 0) {
    av_noncoloc <- past_out$ucos$avW_ucos_each
    tmp <- do.call(rbind, strsplit(colnames(av_noncoloc), ":"))
    colnames(av_noncoloc) <- paste0("outcome", gsub("[^0-9.]+", "", tmp[, 2]))
  } else {
    av_noncoloc <- NULL
  }
  av_all <- cbind(av_coloc, av_noncoloc)
  pip <- vector(mode = "list", length = R)
  if (!is.null(av_all)) {
    av_name <- colnames(av_all)
    for (i in 1:R) {
      pos <- grep(i, av_name)
      if (length(pos) != 0) {
        av_i <- as.matrix(av_all[, pos])
        pip[[i]] <- as.vector(1 - apply(1 - av_i, 1, prod))
      } else {
        pip[[i]] <- rep(0, P)
      }
    }
  }
  return(pip)
}

check_two_overlap_sets <- function(total, i, j) {
  t1 <- total[[i]]
  t2 <- total[[j]]
  if (t1 != "ONE" & t2 != "ONE") {
    return(ifelse(t1 > t2, i, j))
  } else if (t1 == "ONE" & t2 != "ONE") {
    return(i)
  } else if (t1 != "ONE" & t2 == "ONE") {
    return(j)
  } else if (t1 == "ONE" & t2 == "ONE") {
    return(sample(c(i, j), 1))
  }
}

# Function to merge overlapping sets
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


get_avWeigth <- function(cb_model, coloc_outcomes, update, pos.coloc, name_weight = FALSE) {
  avWeight <- lapply(coloc_outcomes, function(i) {
    pos <- which(update[i, ] != 0)
    weight <- cb_model[[i]]$weights_path[match(pos.coloc, pos), ]
  })
  avWeight <- do.call(cbind, avWeight)
  if (name_weight) {
    colnames(avWeight) <- paste0("outcome", coloc_outcomes)
  }
  return(avWeight)
}


get_max_profile <- function(cb_obj, check_null_max = 0.025, 
                            check_null_max_ucos = 0.015,
                            check_null_method = "profile") {
  for (i in 1:cb_obj$cb_model_para$L) {
    cb <- cb_obj$cb_model[[i]]
    scaling_factor <- if (!is.null(cb_obj$cb_data$data[[i]]$N)) cb_obj$cb_data$data[[i]]$N - 1 else 1
    if (check_null_method == "profile") {
      cb$check_null_max <- 1000 * check_null_max / scaling_factor
      cb$check_null_max_ucos <- 1000 * check_null_max_ucos / scaling_factor
    } else {
      cb$check_null_max <- check_null_max
      cb$check_null_max_ucos <- check_null_max_ucos
    }
    cb_obj$cb_model[[i]] <- cb
  }
  return(cb_obj)
}


# ### Function for check cs for each weight
# w_cs <- function(weights, coverage = 0.95) {
#   indices <- unlist(get_in_cos(weights, coverage = coverage))
#   result <- rep(0, length(weights))
#   result[indices] <- 1
#   return(result)
# }

### Function for check cs for each weight
w_cs <- function(weights, coverage = 0.95) {
  indices <- unlist(get_in_cos(weights, coverage = coverage))
  return(indices)
}

#' Pure R implementation (fallback)
#' @noRd
get_merge_ordered_with_indices <- function(vector_list) {
  
  # Fallback function for priority-based merging (embedded)
  get_merge_ordered_with_priority_fallback <- function(vector_list) {
    # Convert all vectors to character
    vector_list <- lapply(vector_list, as.character)
    
    # Start with empty result
    result <- character()
    processed_elements <- character()
    
    # Process vectors in priority order (1st, 2nd, 3rd, etc.)
    for (i in seq_along(vector_list)) {
      current_vec <- vector_list[[i]]
      
      # Add elements from current vector that haven't been processed yet
      new_elements <- setdiff(current_vec, processed_elements)
      
      if (length(new_elements) > 0) {
        # Find positions of new elements in current vector to maintain relative order
        new_positions <- which(current_vec %in% new_elements)
        ordered_new_elements <- current_vec[new_positions]
        
        # Add to result
        result <- c(result, ordered_new_elements)
        processed_elements <- c(processed_elements, ordered_new_elements)
      }
    }
    
    return(result)
  }
  
  # Convert all vectors to character
  vector_list <- lapply(vector_list, as.character)
  n_vectors <- length(vector_list)
  # Step 1: Get all unique elements
  all_elements <- unique(unlist(vector_list))
  n_elements <- length(all_elements)
  # Step 2: Build a graph of ordering constraints
  # Use an adjacency list: for each element, store elements that must come after it
  graph <- new.env(hash = TRUE, parent = emptyenv(), size = n_elements)
  for (elem in all_elements) {
    graph[[elem]] <- character()
  }
  # Add edges based on consecutive pairs in each vector
  for (vec in vector_list) {
    for (i in seq_len(length(vec) - 1)) {
      from_elem <- vec[i]
      to_elem <- vec[i + 1]
      if (from_elem != to_elem) { # Avoid self-loops
        # Add to_elem to the list of elements that must come after from_elem
        graph[[from_elem]] <- unique(c(graph[[from_elem]], to_elem))
      }
    }
  }
  # Step 3: Compute in-degrees (number of incoming edges for each node)
  in_degree <- new.env(hash = TRUE, parent = emptyenv(), size = n_elements)
  for (elem in all_elements) {
    in_degree[[elem]] <- 0
  }
  for (from_elem in all_elements) {
    for (to_elem in graph[[from_elem]]) {
      in_degree[[to_elem]] <- in_degree[[to_elem]] + 1
    }
  }
  # Step 4: Topological sort using Kahn's algorithm
  # Start with nodes that have no incoming edges
  queue <- all_elements[sapply(all_elements, function(elem) in_degree[[elem]] == 0)]
  result <- list()
  result_idx <- 1
  while (length(queue) > 0) {
    # Take the first element from the queue
    current <- queue[1]
    queue <- queue[-1]
    result[[result_idx]] <- current  # List assignment is fast
    result_idx <- result_idx + 1
    # Process all neighbors (elements that must come after current)
    neighbors <- graph[[current]]
    for (next_elem in neighbors) {
      in_degree[[next_elem]] <- in_degree[[next_elem]] - 1
      if (in_degree[[next_elem]] == 0) {
        queue <- c(queue, next_elem)
      }
    }
  }
  result <- unlist(result)
  # Step 5: Check for cycles and use fallback if needed
  if (length(result) != n_elements) {
    # Different variable orders detected - use priority-based fallback
    warning("Variable names in different datasets have different orders (e.g., dataset1 has A then B, but dataset2 has B then A). ",
            "Using priority-based fallback to prioritize the variable order in the earlier datasets.")
    result <- get_merge_ordered_with_priority_fallback(vector_list)
  }
  
  result
}



#' @title Set of internal functions to re-organize ColocBoost output format
#'
#' @description
#' The `colocboost_output_reorganization` functions access basic properties inferences from a fitted ColocBoost model. This documentation serves as a summary for all related post-inference functions.
#'
#' @details
#' The following functions are included in this set:
#' `get_data_info` get formatted \code{data_info} in ColocBoost output
#' `get_cos_details` get formatted \code{cos_details} in ColocBoost output
#' `get_model_info` get formatted \code{model_info} in ColocBoost output
#' `get_full_output` get formatted additional information in ColocBoost output when \code{output_level!=1}
#'
#' These functions are not exported individually and are accessed via `colocboost_output_reorganization`.
#'
#' @rdname colocboost_output_reorganization
#' @keywords cb_reorganization
#' @noRd
colocboost_output_reorganization <- function() {
  message("This function re-formats colocboost output as internal used. See details for more information.")
}

#' Get formatted \code{data_info} in ColocBoost output
#' @noRd
#' @keywords cb_reorganization
get_data_info <- function(cb_obj) {
  ## - analysis data information
  n_outcome <- cb_obj$cb_model_para$L
  n_variables <- cb_obj$cb_model_para$P
  analysis_outcome <- cb_obj$cb_model_para$outcome_names
  variables <- cb_obj$cb_data$variable.names
  focal_outcome <- NULL
  is_focal <- rep(FALSE, n_outcome)
  if (!is.null(cb_obj$cb_model_para$focal_outcome_idx)) {
    focal_outcome <- analysis_outcome[cb_obj$cb_model_para$focal_outcome_idx]
    is_focal[cb_obj$cb_model_para$focal_outcome_idx] <- TRUE
  }
  is_sumstat <- grepl("sumstat_outcome", names(cb_obj$cb_data$data))
  N <- cb_obj$cb_model_para$N
  check_no_N <- sapply(cb_obj$cb_model_para$N, is.null)
  if (sum(check_no_N)!=0){
    N[which(check_no_N)] <- "NA"
    N <- unlist(N)
  }
  outcome_info <- data.frame(
    "outcome_names" = analysis_outcome, "sample_size" = N,
    "is_sumstats" = is_sumstat, "is_focal" = is_focal
  )
  rownames(outcome_info) <- paste0("y", 1:n_outcome)

  ## - marginal associations
  z_scores <- lapply(cb_obj$cb_model, function(cb) {
    as.numeric(cb$z_univariate)
  })
  betas <- lapply(cb_obj$cb_model, function(cb) {
    as.numeric(cb$beta)
  })
  names(z_scores) <- names(betas) <- analysis_outcome
  ## - output data info
  data.info <- list(
    "n_outcomes" = n_outcome,
    "n_variables" = n_variables,
    "outcome_info" = outcome_info,
    "variables" = variables,
    "coef" = betas,
    "z" = z_scores
  )
  return(data.info)
}

#' Get formatted \code{cos_details} in ColocBoost output
#' @noRd
#' @keywords cb_reorganization
get_cos_details <- function(cb_obj, coloc_out, data_info = NULL) {
  if (is.null(data_info)) {
    data_info <- get_data_info(cb_obj)
  }


  ### ----- Define the colocalization results
  coloc_sets <- coloc_out$cos
  if (length(coloc_sets) != 0) {
    # - colocalization outcome configurations
    tmp <- get_cos_evidence(cb_obj, coloc_out, data_info)
    normalization_evidence <- tmp$normalization_evidence
    npc <- tmp$npc
    cos_min_npc_outcome <- sapply(normalization_evidence, function(cp) min(cp$npc_outcome))

    # - colocalized outcomes
    analysis_outcome <- cb_obj$cb_model_para$outcome_names
    coloc_outcome_index <- coloc_outcome <- list()
    colocset_names <- c()
    for (i in 1:length(coloc_out$cos)) {
      coloc_outcome_index[[i]] <- coloc_out$coloc_outcomes[[i]]
      coloc_outcome[[i]] <- analysis_outcome[coloc_outcome_index[[i]]]
      colocset_names[i] <- paste0("cos", i, ":", paste0(paste0("y", coloc_outcome_index[[i]]), collapse = "_"))
      if (grepl("merged", names(coloc_sets)[i])) {
        colocset_names[i] <- paste0(colocset_names[i], ":merged")
      }
    }
    names(coloc_outcome) <- names(coloc_outcome_index) <- colocset_names
    names(npc) <- names(normalization_evidence) <- names(cos_min_npc_outcome) <- colocset_names
    coloc_outcomes <- list("outcome_index" = coloc_outcome_index, "outcome_name" = coloc_outcome)

    # - colocalized sets for variables
    coloc_csets_variableidx <- coloc_out$cos
    coloc_csets_variablenames <- lapply(coloc_csets_variableidx, function(coloc_tmp) {
      cb_obj$cb_model_para$variables[coloc_tmp]
    })
    coloc_csets_variableidx <- lapply(coloc_csets_variablenames, function(variable) match(variable, data_info$variables))
    names(coloc_csets_variableidx) <- names(coloc_csets_variablenames) <- colocset_names
    coloc_csets_original <- list("cos_index" = coloc_csets_variableidx, "cos_variables" = coloc_csets_variablenames)

    # - colocalized set cs_change
    cs_change <- coloc_out$cs_change
    rownames(cs_change) <- colocset_names
    colnames(cs_change) <- analysis_outcome

    # - VCP
    cos_weights <- lapply(coloc_out$avWeight, function(w) {
      pos <- match(data_info$variables, cb_obj$cb_model_para$variables)
      return(w[pos, , drop = FALSE])
    })
    int_weight <- lapply(cos_weights, get_integrated_weight, weight_fudge_factor = cb_obj$cb_model_para$weight_fudge_factor)
    names(int_weight) <- names(cos_weights) <- colocset_names

    # - resummary results
    cos_re_idx <- lapply(int_weight, function(w) {
      unlist(get_in_cos(w, coverage = cb_obj$cb_model_para$coverage))
    })
    cos_re_var <- lapply(cos_re_idx, function(idx) {
      data_info$variables[idx]
    })
    coloc_csets <- list("cos_index" = cos_re_idx, "cos_variables" = cos_re_var)
    
    # - recalculate purity
    purity <- vector(mode = "list", length = length(cos_re_idx))
    for (ee in 1:length(cos_re_idx)) {
      coloc_t <- coloc_outcome_index[[ee]]
      p_tmp <- c()
      for (i3 in coloc_t) {
        pos <- cos_re_idx[[ee]]
        X_dict <- cb_obj$cb_data$dict[i3]
        if (!is.null(cb_obj$cb_data$data[[X_dict]]$XtX)) {
          pos <- match(pos, setdiff(1:cb_obj$cb_model_para$P, cb_obj$cb_data$data[[i3]]$variable_miss))
        }
        tmp <- matrix(get_purity(pos,
                                 X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
                                 Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                 N = cb_obj$cb_data$data[[i3]]$N, n = cb_obj$cb_model_para$n_purity
        ), 1, 3)
        p_tmp <- rbind(p_tmp, tmp)
      }
      purity[[ee]] <- matrix(apply(p_tmp, 2, max), 1, 3)
    }
    purity_all <- do.call(rbind, purity)
    purity_all <- as.data.frame(purity_all)
    colnames(purity_all) <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
    coloc_out$purity <- purity_all
    if (is.null(cb_obj$cb_model_para$median_abs_corr)) {
      is_pure <- which(purity_all[, 1] >= cb_obj$cb_model_para$min_abs_corr)
    } else {
      is_pure <- which(purity_all[, 1] >= cb_obj$cb_model_para$min_abs_corr | purity_all[, 3] >= cb_obj$cb_model_para$median_abs_corr)
    }
    if (length(is_pure)==0){
      coloc_results <- NULL
      vcp <- NULL
      return(list("cos_results" = coloc_results, "vcp" = vcp))
    } else if (length(is_pure)!=length(cos_re_idx)){
      int_weight = int_weight[is_pure]
      coloc_csets <- lapply(coloc_csets, function(cs) cs[is_pure])
      coloc_outcomes <- lapply(coloc_outcomes, function(cs) cs[is_pure])
      normalization_evidence <- normalization_evidence[is_pure]
      npc <- npc[is_pure]
      cos_min_npc_outcome <- cos_min_npc_outcome[is_pure]
      cos_weights <- cos_weights[is_pure]
      coloc_out$purity <- purity_all[is_pure,,drop = FALSE]
      colocset_names <- colocset_names[is_pure]
    }
    vcp <- as.vector(1 - apply(1 - do.call(cbind, int_weight), 1, prod))
    names(vcp) <- data_info$variables
    
    
    # - hits variables in each csets
    coloc_hits <- coloc_hits_variablenames <- coloc_hits_names <- c()
    for (i in 1:length(int_weight)) {
      inw <- int_weight[[i]]
      pp <- which(inw == max(inw))
      coloc_hits <- c(coloc_hits, pp)
      coloc_hits_variablenames <- c(coloc_hits_variablenames, data_info$variables[pp])
      if (length(pp) == 1) {
        coloc_hits_names <- c(coloc_hits_names, names(int_weight)[i])
      } else {
        coloc_hits_names <- c(coloc_hits_names, paste0(names(int_weight)[i], ".", 1:length(pp)))
      }
    }
    coloc_hits <- data.frame("top_index" = coloc_hits, "top_variables" = coloc_hits_variablenames)
    rownames(coloc_hits) <- coloc_hits_names

    # - purity
    ncos <- length(coloc_csets$cos_index)
    if (ncos >= 2) {
      empty_matrix <- matrix(NA, ncos, ncos)
      colnames(empty_matrix) <- rownames(empty_matrix) <- colocset_names
      csets_purity <- lapply(1:3, function(ii) {
        diag(empty_matrix) <- coloc_out$purity[, ii]
        return(empty_matrix)
      })
      for (i in 1:(ncos - 1)) {
        for (j in (i + 1):ncos) {
          cset1 <- coloc_csets$cos_index[[i]]
          cset2 <- coloc_csets$cos_index[[j]]
          y.i <- coloc_outcomes$outcome_index[[i]]
          y.j <- coloc_outcomes$outcome_index[[j]]
          yy <- unique(c(y.i, y.j))
          res <- list()
          flag <- 1
          for (ii in yy) {
            X_dict <- cb_obj$cb_data$dict[ii]
            res[[flag]] <- get_between_purity(cset1, cset2,
              X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
              Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
              miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
              P = cb_obj$cb_model_para$P
            )
            flag <- flag + 1
          }
          res <- Reduce(pmax, res)
          csets_purity <- lapply(1:3, function(ii) {
            csets_purity[[ii]][i, j] <- csets_purity[[ii]][j, i] <- res[ii]
            return(csets_purity[[ii]])
          })
        }
      }
      names(csets_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
    } else {
      csets_purity <- lapply(1:length(coloc_out$purity), function(i) {
        tmp <- as.matrix(coloc_out$purity[i])
        rownames(tmp) <- colnames(tmp) <- colocset_names
        tmp
      })
      names(csets_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
    }

    # - save coloc_results
    coloc_results <- list(
      "cos" = coloc_csets,
      "cos_outcomes" = coloc_outcomes,
      "cos_vcp" = int_weight,
      "cos_outcomes_npc" = normalization_evidence,
      "cos_npc" = npc,
      "cos_min_npc_outcome" = cos_min_npc_outcome,
      "cos_purity" = csets_purity,
      "cos_top_variables" = coloc_hits,
      "cos_weights" = cos_weights
    )


    # - missing variable and warning message
    missing_variables_idx <- Reduce(union, lapply(cb_obj$cb_data$data, function(cb) cb$variable_miss))
    missing_variables <- cb_obj$cb_model_para$variables[missing_variables_idx]
    cos_missing_variables_idx <- lapply(coloc_csets_original$cos_variables, function(variable) {
      missing <- intersect(variable, missing_variables)
      if (length(missing) != 0) {
        match(missing, data_info$variables)
      } else {
        NULL
      }
    })
    cos_missing_variables <- lapply(cos_missing_variables_idx, function(variable) {
      if (!is.null(variable)) {
        data_info$variables[variable]
      } else {
        NULL
      }
    })
    warning_needed <- any(!sapply(cos_missing_variables, is.null))
    if (warning_needed) {
      is_missing <- which(!sapply(cos_missing_variables, is.null))
      cos_missing_variables_idx <- cos_missing_variables_idx[is_missing]
      cos_missing_variables <- cos_missing_variables[is_missing]
      cos_missing_vcp <- lapply(cos_missing_variables_idx, function(idx) {
        vcp[idx]
      })
      warning_message <- paste(
        "CoS", paste(names(cos_missing_variables_idx), collapse = ","),
        "contains missing variables in at least one outcome.",
        "The missing variables will cause the ~0 VCP scores."
      )
      cos_warnings <- list(
        "cos_missing_info" = list(
          "index" = cos_missing_variables_idx,
          "variables" = cos_missing_variables,
          "vcp" = cos_missing_vcp
        ),
        "warning_message" = warning_message
      )
      coloc_results$cos_warnings <- cos_warnings
    }
  } else {
    coloc_results <- NULL
    vcp <- NULL
  }
  return(list("cos_results" = coloc_results, "vcp" = vcp))
}

#' Get formatted \code{model_info} in ColocBoost output
#' @noRd
#' @keywords cb_reorganization
get_model_info <- function(cb_obj, outcome_names = NULL) {
  if (is.null(outcome_names)) {
    data_info <- get_data_info(cb_obj)
    outcome_names <- data_info$outcome_info$outcome_names
  }

  profile_loglik <- cb_obj$cb_model_para$profile_loglike
  n_updates <- cb_obj$cb_model_para$num_updates
  n_updates_outcome <- cb_obj$cb_model_para$num_updates_outcome
  model_coveraged <- cb_obj$cb_model_para$coveraged
  model_coveraged_outcome <- cb_obj$cb_model_para$coveraged_outcome
  jk_update <- cb_obj$cb_model_para$real_update_jk
  if (!is.null(jk_update)){
    rownames(jk_update) <- paste0("jk_star_", 1:nrow(jk_update))
    colnames(jk_update) <- outcome_names
  }
  outcome_proximity_obj <- lapply(cb_obj$cb_model, function(cb) cb$obj_path)
  outcome_coupled_best_update_obj <- lapply(cb_obj$cb_model, function(cb) cb$obj_single)
  outcome_profile_loglik <- lapply(cb_obj$cb_model, function(cb) cb$profile_loglike_each)
  names(outcome_proximity_obj) <- names(outcome_coupled_best_update_obj) <-
    names(outcome_profile_loglik) <- names(n_updates_outcome) <- 
    names(model_coveraged_outcome) <- outcome_names
  ll <- list(
    "model_coveraged" = model_coveraged,
    "n_updates" = n_updates,
    "profile_loglik" = profile_loglik,
    "outcome_profile_loglik" = outcome_profile_loglik,
    "outcome_proximity_obj" = outcome_proximity_obj,
    "outcome_coupled_best_update_obj" = outcome_coupled_best_update_obj,
    "outcome_model_coveraged" = model_coveraged_outcome,
    "outcome_n_updates" = n_updates_outcome,
    "jk_star" = jk_update
  )
  return(ll)
}

#' Get formatted additional information in ColocBoost output when \code{output_level!=1}
#' @noRd
#' @keywords cb_reorganization
get_full_output <- function(cb_obj, past_out = NULL, variables = NULL, cb_output = NULL, weaker_ucos = FALSE) {
  cb_model <- cb_obj$cb_model
  cb_model_para <- cb_obj$cb_model_para

  ## - obtain the order of variables based on the variables names if it has position information
  if (is.null(variables)) {
    ordered <- 1:length(cb_obj$cb_model_para$variables)
  } else {
    ordered <- match(variables, cb_obj$cb_model_para$variables)
  }

  ## - reorder all output
  # - cb_model
  tmp <- lapply(cb_model, function(cb) {
    cb$beta <- cb$beta[ordered]
    cb$weights_path <- cb$weights_path[, ordered]
    cb$change_loglike <- cb$change_loglike[ordered]
    cb$correlation <- as.numeric(cb$correlation[ordered])
    cb$z <- as.numeric(cb$z[ordered])
    cb$ld_jk <- cb$ld_jk[, ordered]
    cb$z_univariate <- as.numeric(cb$z_univariate[ordered])
    cb$beta_hat <- as.numeric(cb$beta_hat[ordered])
    cb$multi_correction <- as.numeric(cb$multi_correction[ordered])
    cb$multi_correction_univariate <- as.numeric(cb$multi_correction_univariate[ordered])
    return(cb)
  })
  cb_model <- tmp

  # - sets
  if (!is.null(past_out)) {
    out_ucos <- past_out$ucos
    # - single sets
    if (!is.null(out_ucos$ucos_each)) {
      out_ucos$ucos_each <- lapply(out_ucos$ucos_each, function(cs) {
        match(cb_model_para$variables[cs], variables)
      })
      out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[ordered, , drop = FALSE]

      # - re-orginize specific results
      analysis_outcome <- cb_obj$cb_model_para$outcome_names
      specific_outcome_index <- specific_outcome <- list()
      specific_cs_names <- c()
      for (i in 1:length(out_ucos$ucos_each)) {
        cc <- out_ucos$avW_ucos_each[, i, drop = FALSE]
        tmp_names <- colnames(cc)
        specific_outcome_index[[i]] <- as.numeric(gsub(".*Y([0-9]+).*", "\\1", tmp_names))
        specific_outcome[[i]] <- analysis_outcome[specific_outcome_index[[i]]]
        specific_cs_names[i] <- paste0("ucos", i, ":y", specific_outcome_index[[i]])
      }
      names(specific_outcome) <- names(specific_outcome_index) <- specific_cs_names
      specific_outcomes <- list("outcome_index" = specific_outcome_index, "outcome_name" = specific_outcome)

      # - specific sets for variables
      specific_cs_variableidx <- out_ucos$ucos_each
      specific_cs_variablenames <- lapply(specific_cs_variableidx, function(specific_tmp) {
        cb_obj$cb_model_para$variables[specific_tmp]
      })
      specific_cs_variableidx <- lapply(specific_cs_variablenames, function(variable) match(variable, variables))
      names(specific_cs_variableidx) <- names(specific_cs_variablenames) <- specific_cs_names
      specific_css <- list("ucos_index" = specific_cs_variableidx, "ucos_variables" = specific_cs_variablenames)

      # - specific set cs_change 
      cs_change <- out_ucos$change_obj_each
      rownames(cs_change) <- specific_cs_names
      colnames(cs_change) <- analysis_outcome
      index_change <- as.data.frame(which(cs_change != 0, arr.ind = TRUE))
      change_outcomes <- analysis_outcome[index_change$col]
      change_values <- diag(as.matrix(cs_change[index_change$row, index_change$col]))
      cs_change <- data.frame("ucos_outcome" = change_outcomes, "ucos_delta" = change_values)

      # - filter weak ucos
      check_null_max_ucos <- sapply(cb_model, function(cb) cb$check_null_max_ucos)
      remove_weak <- sapply(1:nrow(cs_change), function(ic) {
        outcome_tmp <- cs_change$ucos_outcome[ic]
        delta_tmp <- cs_change$ucos_delta[ic]
        pp <- which(cb_obj$cb_model_para$outcome_names == outcome_tmp)
        check_tmp <- check_null_max_ucos[pp]
        delta_tmp >= check_tmp
      })
      keep_ucos <- which(remove_weak)
      if (length(keep_ucos) == 0) {
        specific_results <- NULL
      } else {
        specific_outcomes$outcome_index <- specific_outcomes$outcome_index[keep_ucos]
        specific_outcomes$outcome_name <- specific_outcomes$outcome_name[keep_ucos]
        specific_css$ucos_index <- specific_css$ucos_index[keep_ucos]
        specific_css$ucos_variables <- specific_css$ucos_variables[keep_ucos]
        cs_change <- cs_change[keep_ucos, , drop = FALSE]
        out_ucos$avW_ucos_each <- out_ucos$avW_ucos_each[, keep_ucos, drop = FALSE]
        specific_cs_names <- specific_cs_names[keep_ucos]
        out_ucos$purity_each <- out_ucos$purity_each[keep_ucos, , drop = FALSE]

        # - ucos_weight
        specific_w <- lapply(1:ncol(out_ucos$avW_ucos_each), function(ii) out_ucos$avW_ucos_each[, ii, drop = FALSE])
        names(specific_w) <- specific_cs_names

        # - hits variables in each csets
        cs_hits <- sapply(1:length(specific_w), function(jj) {
          inw <- specific_w[[jj]]
          which(inw == max(inw))[1]
        })
        cs_hits_variablenames <- sapply(cs_hits, function(ch) variables[ch])
        specific_cs_hits <- data.frame("top_index" = cs_hits, "top_variables" = cs_hits_variablenames) # save
        rownames(specific_cs_hits) <- specific_cs_names

        # - add npc_outcome for ucos (Dec 08, 2025)
        ucoloc_info <- list(
            ucos = specific_css$ucos_index,
            outcome = specific_outcomes$outcome_index,
            outcome_name = specific_outcomes$outcome_name
        )
        npc_outcome <- get_ucos_evidence(cb_obj, ucoloc_info)

        # - purity
        nucos <- length(specific_css$ucos_index)
        if (nucos >= 2) {
          empty_matrix <- matrix(NA, nucos, nucos)
          colnames(empty_matrix) <- rownames(empty_matrix) <- specific_cs_names
          specific_cs_purity <- lapply(1:3, function(ii) {
            diag(empty_matrix) <- out_ucos$purity_each[, ii]
            return(empty_matrix)
          })
          for (i in 1:(nucos - 1)) {
            for (j in (i + 1):nucos) {
              cset1 <- specific_css$ucos_index[[i]]
              cset2 <- specific_css$ucos_index[[j]]
              y.i <- specific_outcomes$outcome_index[[i]]
              y.j <- specific_outcomes$outcome_index[[j]]
              yy <- unique(c(y.i, y.j))
              res <- list()
              flag <- 1
              for (ii in yy) {
                X_dict <- cb_obj$cb_data$dict[ii]
                res[[flag]] <- get_between_purity(cset1, cset2,
                  X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
                  Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                  miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                  P = cb_obj$cb_model_para$P
                )
                flag <- flag + 1
              }
              res <- Reduce(pmax, res)
              specific_cs_purity <- lapply(1:3, function(ii) {
                specific_cs_purity[[ii]][i, j] <- specific_cs_purity[[ii]][j, i] <- res[ii]
                return(specific_cs_purity[[ii]])
              })
            }
          }
          names(specific_cs_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
        } else {
          specific_cs_purity <- lapply(1:3, function(ii) {
            mm <- as.matrix(out_ucos$purity_each[,ii])
            rownames(mm) <- colnames(mm) <- specific_cs_names
            return(mm)
          })
          names(specific_cs_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
        }

        # - cos&ucos purity
        cos <- cb_output$cos_details$cos$cos_index
        ncos <- length(cos)
        if (ncos != 0) {
          empty_matrix <- matrix(NA, ncos, nucos)
          colnames(empty_matrix) <- specific_cs_names
          rownames(empty_matrix) <- names(cos)
          cos_ucos_purity <- lapply(1:3, function(ii) empty_matrix)
          for (i in 1:ncos) {
            for (j in 1:nucos) {
              cset1 <- cos[[i]]
              cset2 <- specific_css$ucos_index[[j]]
              y.i <- cb_output$cos_details$cos_outcomes$outcome_index[[i]]
              y.j <- specific_outcomes$outcome_index[[j]]
              yy <- unique(c(y.i, y.j))
              res <- list()
              flag <- 1
              for (ii in yy) {
                X_dict <- cb_obj$cb_data$dict[ii]
                res[[flag]] <- get_between_purity(cset1, cset2,
                  X = get_genotype_matrix(cb_obj$cb_data$data[[X_dict]]),
                  Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                  miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                  P = cb_obj$cb_model_para$P
                )
                flag <- flag + 1
              }
              res <- Reduce(pmax, res)
              cos_ucos_purity <- lapply(1:3, function(ii) {
                cos_ucos_purity[[ii]][i, j] <- res[ii]
                return(cos_ucos_purity[[ii]])
              })
            }
          }
          names(cos_ucos_purity) <- c("min_abs_cor", "max_abs_cor", "median_abs_cor")
        } else {
          cos_ucos_purity <- NULL
        }


        # - save coloc_results
        specific_results <- list(
          "ucos" = specific_css,
          "ucos_outcomes" = specific_outcomes,
          "ucos_weight" = specific_w,
          "ucos_top_variables" = specific_cs_hits,
          "ucos_purity" = specific_cs_purity,
          "cos_ucos_purity" = cos_ucos_purity,
          "ucos_outcomes_npc" = npc_outcome
        )
      }
    } else {
      specific_results <- NULL
    }

    # - cb_model_para
    cb_model_para$N <- as.numeric(unlist(cb_model_para$N))
    cb_model_para$variables <- variables

    ll <- list(
      "ucos_details" = specific_results,
      "cb_model" = cb_model,
      "cb_model_para" = cb_model_para
    )
  } else {
    # - cb_model_para
    cb_model_para$N <- as.numeric(unlist(cb_model_para$N))
    cb_model_para$variables <- variables
    ll <- list(
      "ucos_detials" = NULL,
      "cb_model" = cb_model,
      "cb_model_para" = cb_model_para
    )
  }

  return(ll)
}

