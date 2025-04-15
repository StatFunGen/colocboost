merge_cos_ucos <- function(cb_obj, out_cos, out_ucos, coverage = 0.95,
                               min_abs_corr = 0.5, tol = 1e-9,
                               median_cos_abs_corr = 0.8){


    change_obj_each <-  out_ucos$change_obj_each
    coloc_sets <-  out_cos$cos$cos
    ucos_each <-  out_ucos$ucos_each

    # - remove overlap between coloc_sets and single_sets
    is_overlap <- is_highLD <- c()
    for (i in 1:length(coloc_sets)){
        # cset1 <- coloc_sets[[i]]
        coloc_outcome <-  out_cos$cos$coloc_outcomes[[i]]
        ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames( out_cos$cos$avWeight[[i]])))
        pos_name <- order(ttmp)
        out_cos$cos$avWeight[[i]] <-  out_cos$cos$avWeight[[i]][,pos_name]
        for (j in 1:length(ucos_each)){
            cset2 <- ucos_each[[j]]
            fine_outcome <-  out_ucos$ucos_outcome[j]

            if (fine_outcome %in% coloc_outcome){

                # - if fine_y in coloc_y, we check only the overlap
                pp <- which(coloc_outcome == fine_outcome)
                w_tmp <-  out_cos$cos$avWeight[[i]][,pp]
                cset1 <- get_in_cos(w_tmp, coverage = coverage)[[1]]
                # - addhoc: if median_cos_abs_corr > 0.8, remove single sets
                X_dict <- cb_obj$cb_data$dict[fine_outcome]
                res <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                          Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                          N = cb_obj$cb_data$data[[fine_outcome]]$N,
                                          miss_idx = cb_obj$cb_data$data[[fine_outcome]]$variable_miss,
                                          P = cb_obj$cb_model_para$P)
                # is.between <- length(intersect(cset1, cset2)) != 0
                is.between <- (abs(res[2]-1) < tol)
                if (is.between){
                    is_overlap = c(is_overlap, j)
                } else {
                    is.higLD <- res[2] > median_cos_abs_corr
                    if (is.higLD){ is_highLD = c(is_highLD, j)}
                }

            } else if (!(fine_outcome %in% coloc_outcome)){

                # - if fine_y not in coloc_y, we check overlap and also min_purity
                change_obj_coloc <-  out_cos$cos$cs_change
                cset1 <- coloc_sets[[i]]
                res <- list()
                for (ii in 1:cb_obj$cb_model_para$L){
                    X_dict <- cb_obj$cb_data$dict[ii]
                    res[[ii]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                    Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                    N = cb_obj$cb_data$data[[ii]]$N,
                                                    miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                                                    P = cb_obj$cb_model_para$P)
                }
                res <- Reduce(pmax, res)
                min_between <- res[1]
                max_between <- res[2]
                ave_between <- res[3]
                is.between <- ((min_between>median_cos_abs_corr) & (abs(max_between-1) < tol))
                # is.between <- (min_between>min_abs_corr) * (abs(max_between-1)<tol) * (ave_between>median_cos_abs_corr)
                if (is.between){
                    is_overlap = c(is_overlap, j)
                    # -- add weight
                    add_avW <-  out_ucos$avW_ucos_each[,j]
                    out_cos$cos$avWeight[[i]] <- cbind(add_avW,  out_cos$cos$avWeight[[i]])
                    colnames( out_cos$cos$avWeight[[i]])[1] <- paste0("outcome", fine_outcome)
                    ttmp <- as.numeric(gsub("[^0-9.]+", "", colnames( out_cos$cos$avWeight[[i]])))
                    pos_name <- order(ttmp)
                    out_cos$cos$avWeight[[i]] <-  out_cos$cos$avWeight[[i]][,pos_name]
                    change_obj_coloc_i <- change_obj_coloc[i,]
                    change_obj_each_j <- change_obj_each[j,]
                    out_cos$cos$cs_change[i,] <- pmax(change_obj_coloc_i, change_obj_each_j)
                    coloc_outcome <- sort(c(coloc_outcome, fine_outcome))
                    out_cos$cos$coloc_outcomes[[i]] <- coloc_outcome
                }

            }
        }
    }

    is_overlap <- unique(c(is_overlap,is_highLD))
    if (length(is_overlap) != 0){
        if (length(is_overlap) == length(ucos_each)){
             out_ucos$ucos_each = NULL
             out_ucos$avW_ucos_each <- NULL
             out_ucos$change_obj_each <- NULL
             out_ucos$purity_each <- NULL
        } else {
             out_ucos$ucos_each = ucos_each[-is_overlap]
             out_ucos$avW_ucos_each <-  out_ucos$avW_ucos_each[,-is_overlap,drop = FALSE]
             out_ucos$change_obj_each <- change_obj_each[-is_overlap,,drop = FALSE]
             out_ucos$purity_each <-  out_ucos$purity_each[-is_overlap,,drop = FALSE]
        }
    }
    if (length( out_ucos$ucos_each) == 0){
         out_ucos$ucos_each = NULL
         out_ucos$avW_ucos_each <- NULL
         out_ucos$change_obj_each <- NULL
         out_ucos$purity_each <- NULL
    }
    ll <- list("ucos" =  out_ucos, "cos" =  out_cos)
    return(ll)

}

#' @importFrom stats na.omit
merge_ucos <- function(cb_obj, past_out,
                       min_abs_corr = 0.5,
                       median_abs_corr = NULL,
                       n_purity = 100,
                       median_cos_abs_corr = 0.8,
                       tol = 1e-9){


    out_ucos <- past_out$ucos
    out_cos <- past_out$cos
    ucos_each <-  out_ucos$ucos_each
    change_obj_each <-  out_ucos$change_obj_each

    # calculate between purity
    ncsets <- length(ucos_each)
    min_between <- max_between <- ave_between <- matrix(0, nrow =  ncsets, ncol =  ncsets)
    for (i.between in 1:(ncsets-1)){
        for (j.between in (i.between+1):ncsets){
            cset1 <- ucos_each[[i.between]]
            cset2 <- ucos_each[[j.between]]
            y.i <-  out_ucos$ucos_outcome[i.between]
            y.j <-  out_ucos$ucos_outcome[j.between]
            if (y.i == y.j){ next }
            yy <- c(y.i, y.j)
            res <- list()
            flag <- 1
            for (ii in yy){
                X_dict <- cb_obj$cb_data$dict[ii]
                res[[flag]] <- get_between_purity(cset1, cset2, X = cb_obj$cb_data$data[[X_dict]]$X,
                                                  Xcorr = cb_obj$cb_data$data[[X_dict]]$XtX,
                                                  N = cb_obj$cb_data$data[[ii]]$N,
                                                  miss_idx = cb_obj$cb_data$data[[ii]]$variable_miss,
                                                  P = cb_obj$cb_model_para$P)
                flag <- flag + 1
            }
            if (min_abs_corr == 0){
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
    is.between <- ((min_between>median_cos_abs_corr) & (abs(max_between-1) < tol))
    if (sum(is.between) != 0){
        temp <- sapply(1:nrow(is.between), function(x){
            tt <- c(x, which(is.between[x,] != 0))
            return(paste0(sort(tt), collapse = ";"))
        })
        temp <- merge_sets(temp)
        potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
        potential_merged <- potential_merged[which(sapply(potential_merged, length) >= 2)]
        coloc_sets_merged <- avWeight_merged <-
            cs_change_merged <- coloc_outcomes_merged <- list()
        is_merged <- c()
        for (i.m in 1:length(potential_merged)){
            temp_set <- as.numeric(potential_merged[[i.m]])
            # refine avWeight
            merged <-  out_ucos$avW_ucos_each[, temp_set]
            unique_coloc_outcomes <- as.numeric(gsub(".*Y(\\d+).*", "\\1", colnames(merged)))
            if (length(unique(unique_coloc_outcomes))==1) next
            # define merged set
            coloc_sets_merged <- c(coloc_sets_merged, list( unique(unlist(ucos_each[temp_set])) ))
            colnames(merged) <- paste0("outcome", unique_coloc_outcomes)
            coloc_outcomes_merged <- c(coloc_outcomes_merged,
                                       list(unique(sort(unique_coloc_outcomes))))
            # colnames(temp) <- unique_coloc_outcomes
            avWeight_merged <- c(avWeight_merged, list(merged))
            # refine cs_change
            change_cs_tmp <- change_obj_each[temp_set, , drop = FALSE]
            cs_change_merged <- c(cs_change_merged,
                                  list(apply(change_cs_tmp, 2, max)))
            is_merged <- c(is_merged, temp_set)
        }

        if (length(is_merged) != 0){

            # --- check merged coloc set purity
            purity = vector(mode='list', length=length(coloc_sets_merged))
            for(ee in 1:length(coloc_sets_merged)){
                coloc_t <- coloc_outcomes_merged[[ee]]
                p_tmp <- c()
                for (i3 in coloc_t){
                    pos <- coloc_sets_merged[[ee]]
                    X_dict <- cb_obj$cb_data$dict[i3]
                    if (!is.null(cb_obj$cb_data$data[[X_dict]]$XtX)){
                        pos <- match(pos, setdiff(1:cb_obj$cb_model_para$P, cb_obj$cb_data$data[[i3]]$variable_miss))
                        # - it could happen since merge_cos
                        if(sum(is.na(pos)) != 0){pos <- as.numeric(na.omit(pos))}
                    }
                    tmp <- matrix(get_purity(pos,X=cb_obj$cb_data$data[[X_dict]]$X,
                                             Xcorr=cb_obj$cb_data$data[[X_dict]]$XtX,
                                             N=cb_obj$cb_data$data[[i3]]$N,n=n_purity),1,3)
                    p_tmp <- rbind(p_tmp, tmp)
                }
                purity[[ee]] <- matrix(apply(p_tmp, 2, max),1,3)

            }
            purity_all <- do.call(rbind, purity)
            purity_all = as.data.frame(purity_all)
            colnames(purity_all) = c("min_abs_corr","mean_abs_corr","median_abs_corr")
            if (is.null(median_abs_corr)) {
                is_pure = which(purity_all[,1] >= min_abs_corr)
            } else {
                is_pure = which(purity_all[,1] >= min_abs_corr | purity_all[,3] >= median_abs_corr)
            }
            if (length(is_pure) > 0){
                # add coloc
                merged_colocsets <- coloc_sets_merged[is_pure]
                names(merged_colocsets) <- paste0("merged_cos", 1:length(is_pure))
                 out_cos$cos$cos <- c( out_cos$cos$cos,
                                           merged_colocsets)
                 out_cos$cos$purity <- rbind( out_cos$cos$purity,
                                                purity_all[is_pure, ])
                 out_cos$cos$evidence_strength <- c( out_cos$cos$evidence_strength,
                                                       rep(0.95,length(is_pure)))
                cs_change_tmp <- do.call(rbind, cs_change_merged)[is_pure, , drop=FALSE]
                colnames(cs_change_tmp) <- paste0("change_obj_", 1:cb_obj$cb_model_para$L)
                rownames(cs_change_tmp) <- names(merged_colocsets)
                 out_cos$cos$cs_change <- rbind( out_cos$cos$cs_change,
                                                   cs_change_tmp)
                 out_cos$cos$avWeight <- c( out_cos$cos$avWeight,
                                              avWeight_merged[is_pure])
                 out_cos$cos$coloc_outcomes <- c( out_cos$cos$coloc_outcomes,
                                                      coloc_outcomes_merged[is_pure])
            }

            # - remove single
            if (length(is_merged) == length(ucos_each)){
                 out_ucos$ucos_each <- NULL
                 out_ucos$avW_ucos_each <- NULL
                 out_ucos$change_obj_each <- NULL
                 out_ucos$purity_each <- NULL
            } else {
                 out_ucos$ucos_each <- ucos_each[-is_merged]
                 out_ucos$avW_ucos_each <-  out_ucos$avW_ucos_each[, -is_merged, drop = FALSE]
                 out_ucos$change_obj_each <- change_obj_each[-is_merged, , drop = FALSE]
                 out_ucos$purity_each <-  out_ucos$purity_each[-is_merged, , drop = FALSE]
            }
        }
    }
    ll <- list("ucos" =  out_ucos, "cos" =  out_cos)
    return(ll)
}


get_vcp <- function(past_out, P){
    if (!is.null(past_out$cos$cos$cos)){
        avW_coloc_vcp <- sapply(past_out$cos$cos$avWeight, get_integrated_weight)
    } else {avW_coloc_vcp <- NULL}
    all_weight <- avW_coloc_vcp
    if (length(all_weight) == P){
        all_weight <- as.matrix(unlist(all_weight))
    }
    if (!is.null(all_weight)){
        all_weight <- apply(all_weight, 2, as.numeric)
        all_weight <- as.matrix(all_weight)
        vcp <- as.vector(1 - apply(1 - all_weight,1,prod))
    } else {
        vcp <- rep(0, P)
    }
    return(vcp)
}


get_pip <- function(past_out, R, P){
    if (length(past_out$cos$cos$cos)!=0){
        av_coloc <- do.call(cbind, past_out$cos$cos$avWeight)
    } else {av_coloc = NULL}
    if (length(past_out$ucos$ucos_each)!=0){
        av_noncoloc <- past_out$ucos$avW_ucos_each
        tmp <- do.call(rbind, strsplit(colnames(av_noncoloc), ":"))
        colnames(av_noncoloc) <- paste0("outcome", gsub("[^0-9.]+", "", tmp[,2]))
    } else {av_noncoloc = NULL}
    av_all <- cbind(av_coloc, av_noncoloc)
    pip <- vector(mode='list', length=R)
    if (!is.null(av_all)){
        av_name <- colnames(av_all)
        for (i in 1:R){
            pos <- grep(i, av_name)
            if (length(pos) != 0){
                av_i <- as.matrix(av_all[, pos])
                pip[[i]] <- as.vector(1 - apply(1-av_i,1,prod))
            } else {
                pip[[i]] <- rep(0,P)
            }
        }
    }
    return(pip)
}

check_two_overlap_sets <- function(total, i, j){
  t1 <- total[[i]]
  t2 <- total[[j]]
  if (t1 != "ONE" & t2 != "ONE"){
    return(ifelse(t1>t2, i, j))
  } else if (t1 == "ONE" & t2 != "ONE"){
    return(i)
  } else if (t1 != "ONE" & t2 == "ONE"){
    return(j)
  } else if (t1 == "ONE" & t2 == "ONE"){
    return(sample(c(i,j), 1))
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



get_avWeigth <- function(cb_model, coloc_outcomes, update, pos.coloc, name_weight=FALSE){
  
  avWeight <- lapply(coloc_outcomes, function(i){
    pos <- which(update[i,] != 0)
    weight <- cb_model[[i]]$weights_path[match(pos.coloc, pos), ]
  })
  avWeight <- do.call(cbind, avWeight)
  if (name_weight) {
    colnames(avWeight) <- paste0("outcome", coloc_outcomes)
  }
  return(avWeight)
  
}



get_max_profile <- function(cb_obj, check_null_max=0.02, check_null_method = "profile"){
  for (i in 1:cb_obj$cb_model_para$L){
    cb <- cb_obj$cb_model[[i]]
    scaling_factor <- if(!is.null(cb_obj$cb_data$data[[i]]$N)) cb_obj$cb_data$data[[i]]$N-1 else 1
    if (check_null_method == "profile"){
      cb$check_null_max <- 1000*check_null_max / scaling_factor
    } else {
      cb$check_null_max <- check_null_max
    }
    cb_obj$cb_model[[i]] <- cb
  }
  return(cb_obj)
}


### Function for check cs for each weight
w_cs <- function(w, coverage = 0.95){
  n <- sum(cumsum(sort(w,decreasing = TRUE)) < coverage) + 1
  o <- order(w,decreasing = TRUE)
  result = rep(0,length(w))
  result[o[1:n]] = 1
  return(result)
}

get_integrated_weight <- function(avWeight, weight_fudge_factor = 1.5){
  av <- apply(avWeight, 1, function(w) prod(w^(weight_fudge_factor/ncol(avWeight))))
  return(av / sum(av))
}

get_in_cos <- function(weights, coverage = 0.95){
  
  temp <- order(weights, decreasing=T)
  csets <- temp[1:min(which(cumsum(weights[temp]) >= coverage))] # 95%
  return(list(csets))
  
}


#' Pure R implementation (fallback)
#' @noRd
get_merge_ordered_with_indices <- function(vector_list) {
  # Quick validation
  if (!is.list(vector_list) || length(vector_list) == 0) {
    stop("Input must be a non-empty list of vectors")
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
      if (from_elem != to_elem) {  # Avoid self-loops
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
  result <- character()
  while (length(queue) > 0) {
    # Take the first element from the queue
    current <- queue[1]
    queue <- queue[-1]
    result <- c(result, current)

    # Process all neighbors (elements that must come after current)
    neighbors <- graph[[current]]
    for (next_elem in neighbors) {
      in_degree[[next_elem]] <- in_degree[[next_elem]] - 1
      if (in_degree[[next_elem]] == 0) {
        queue <- c(queue, next_elem)
      }
    }
  }

  # Step 5: Check for cycles (if result doesn't include all elements, there’s a cycle)
  if (length(result) != n_elements) {
    stop("Cycle detected in ordering constraints; cannot produce a valid merged order")
  }

  result
}