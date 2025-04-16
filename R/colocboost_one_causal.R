#' @title Set of functions for ColocBoost based on per-trait-one-causal assumption, including LD-free mode and one iteration mode.
#'
#' @description
#' The `colocboost_one_causal` functions access the set of functions for ColocBoost based on per-trait-one-causal assumption.
#'
#'
#' @details
#' The following functions are included in this set:
#' `get_abc` get the colocalization summary table with or without the specific outcomes.
#'
#' These functions are not exported individually and are accessed via `colocboost_one_causal`.
#'
#' @rdname colocboost_one_causal
#' @keywords cb_one_causal
#' @noRd
colocboost_one_causal <- function(cb_model, cb_model_para, cb_data,
                                  jk_equiv_corr = 0.8,
                                  jk_equiv_loglik = 1,
                                  tau = 0.01,
                                  learning_rate_decay = 1,
                                  func_simplex = "z2z",
                                  lambda = 0.5,
                                  lambda_target_outcome = 1,
                                  LD_free = FALSE){
  
  if (jk_equiv_corr != 0){
    cb_obj <- colocboost_one_iteration(cb_model, cb_model_para, cb_data,
                                       jk_equiv_corr = jk_equiv_corr,
                                       jk_equiv_loglik = jk_equiv_loglik,
                                       tau = tau,
                                       learning_rate_decay = learning_rate_decay,
                                       func_simplex = func_simplex,
                                       lambda = lambda,
                                       lambda_target_outcome = lambda_target_outcome,
                                       LD_free = LD_free)
  } else {
    cb_obj <- colocboost_diagLD(cb_model, cb_model_para, cb_data,
                                jk_equiv_corr = jk_equiv_corr,
                                jk_equiv_loglik = jk_equiv_loglik,
                                tau = tau,
                                learning_rate_decay = learning_rate_decay,
                                func_simplex = func_simplex,
                                lambda = lambda,
                                lambda_target_outcome = lambda_target_outcome,
                                LD_free = LD_free)
  }
  return(cb_obj)
  
}



#' Under one causal per trait assumption with one iteration mode
#'
#' @description
#' Executes one iteration of the ColocBoost algorithm under the one-causal-variant-per-trait
#' assumption. This function provides a simplified approach to examine colocalization patterns
#' for traits with strong marginal associations.
#'
#' @details
#' This version implements an LD-dependent mode that can handle scenarios where LD structures
#' may not perfectly match between datasets but still contain informative signals. The function
#' identifies the most strongly associated variants and evaluates their colocalization evidence.
#'
#' @keywords cb_one_causal
#' @noRd
colocboost_one_iteration <- function(cb_model, cb_model_para, cb_data,
                                     jk_equiv_corr = 0.8,
                                     jk_equiv_loglik = 1,
                                     tau = 0.01,
                                     learning_rate_decay = 1,
                                     func_simplex = "z2z",
                                     lambda = 0.5,
                                     lambda_target_outcome = 1,
                                     LD_free = FALSE){
    
    
    if (sum(cb_model_para$update_y == 1) != 0){
      
      ######## - some traits updated
      # - step 1: check update clusters
      real_update <- boost_check_update_jk_one_causal(cb_model, cb_model_para, cb_data,
                                                      jk_equiv_corr = jk_equiv_corr,
                                                      jk_equiv_loglik = jk_equiv_loglik)
      
      # - step 2: boost update
      for (i_update in 1:length(real_update)){
          # - update cb_model_parameter and report results
          cb_model_para$update_temp <- real_update[[i_update]]
          cb_model_para$jk <- rbind(cb_model_para$jk, real_update[[i_update]]$update_jk)
          cb_model_para$update_status <- cbind(cb_model_para$update_status, as.matrix(real_update[[i_update]]$update_status))
          cb_model_para$real_update_jk <- rbind(cb_model_para$real_update_jk, real_update[[i_update]]$real_update_jk)
          # - update cb_model
          cb_model <- colocboost_update(cb_model, cb_model_para, cb_data,
                                        tau = tau,
                                        learning_rate_decay = learning_rate_decay,
                                        func_simplex = func_simplex,
                                        lambda = lambda,
                                        lambda_target_outcome = lambda_target_outcome,
                                        LD_free = LD_free)
      }
    }
    # -- remove redundant parameters
    cb_model_para$num_updates <- 0
    rm_elements <- c("update_temp", "update_y")
    if (!is.null(cb_model_para$need_more)) {rm_elements <- c(rm_elements, "need_more")}
    if (!is.null(cb_model_para$true_stop)) {rm_elements <- c(rm_elements, "true_stop")}
    cb_model_para[rm_elements] <- NULL
    for (i in 1:length(cb_model)){ cb_model[[i]]$obj_path <- as.numeric(unlist(cb_model[[i]]$obj_path[-1]))}
    for (i in 1:length(cb_model)){ cb_model[[i]]$obj_single <- as.numeric(unlist(cb_model[[i]]$obj_single[-1]))}
    
    cb_obj <- list("cb_data" = cb_data, "cb_model" = cb_model, "cb_model_para" = cb_model_para)
    class(cb_obj) <- "colocboost"
    return(cb_obj)
    
}


#' Identify and best update trait group for one causal variant - one iteration mode
#'
#' @details
#' The function operates by partitioning traits into different equivalence groups based on
#' correlation patterns and log-likelihood differences. 
#'
#' @keywords cb_one_causal
#' @noRd
boost_check_update_jk_one_causal <- function(cb_model, cb_model_para, cb_data,
                                             prioritize_jkstar = TRUE, 
                                             jk_equiv_corr = 0.8,
                                             jk_equiv_loglik = 1){
  
    pos.update <- which(cb_model_para$update_y == 1)
    update_jk <- rep(NA, cb_model_para$L+1)
    # - check pairwise equivalent of jk
    if (length(pos.update) == 1){
      
      update_status <- rep(0, cb_model_para$L)
      real_update_jk <- rep(NA, cb_model_para$L)
      update_status[pos.update] <- -1
      correlation <- abs(cb_model[[pos.update]]$correlation)
      jk <- which(correlation == max(correlation))
      jk <- ifelse(length(jk) == 1, jk, sample(jk,1))
      update_jk[c(1,pos.update+1)] <- c(jk, jk)
      real_update_jk[pos.update] <- jk
      real_update <- list("update_status" = update_status, 
                          "real_update_jk" = real_update_jk,
                          "update_jk" = update_jk)
      real_update <- list(real_update)
      
    } else {
      # -- extract data and model based on pos.update
      model_update = cb_model[pos.update]
      X_dict <- cb_data$dict[pos.update]
      adj_dep <- sapply(pos.update, function(ii) cb_data$data[[ii]]$dependency)
      
      # -- define jk and jk_each
      cor_vals_each <- do.call(cbind, lapply(model_update, function(cc) cc$correlation))
      abs_cor_vals_each <- sweep(abs(cor_vals_each), 2, adj_dep, `*`)
      jk_each <- apply(abs_cor_vals_each, 2, function(temp){
        jk_temp <- which(temp == max(temp))
        return(ifelse(length(jk_temp) == 1, jk_temp, sample(jk_temp,1)))
      })
      change_each_pair <- check_pair_jkeach(jk_each, pos.update,
                                            model_update = model_update,
                                            cb_data = cb_data,
                                            X_dict = X_dict,
                                            jk_equiv_corr = jk_equiv_corr,
                                            jk_equiv_loglik = jk_equiv_loglik)
      # define category with same jk
      temp <- sapply(1:nrow(change_each_pair), function(x){
        tt <- c(x, which(change_each_pair[x,] != 0))
        return(paste0(sort(tt), collapse = ";"))
      })
      temp <- merge_sets(temp)
      cate <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
      
      real_update <- lapply(cate, function(cate_each){
          update_status <- rep(0, cb_model_para$L)
          real_update_jk <- rep(NA, cb_model_para$L)
          true_update <- as.numeric(cate_each)
          if (length(true_update)>1){
            update_status[pos.update[true_update]] <- 1
          } else {
            update_status[pos.update[true_update]] <- -1
          }
          abs_cor_vals_redefine <- rowSums(abs_cor_vals_each[,true_update,drop=FALSE])
          jk_redefine <- which(abs_cor_vals_redefine == max(abs_cor_vals_redefine))
          jk_redefine <- ifelse(length(jk_redefine) == 1, jk_redefine, sample(jk_redefine,1))
          real_update_jk[pos.update[true_update]] <- jk_redefine
          update_jk[c(1,pos.update+1)] <- c(jk_redefine, jk_each)
          return(list("update_status" = update_status, 
                      "real_update_jk" = real_update_jk,
                      "update_jk" = update_jk))
      })
    }
    return(real_update)
}



#' ColocBoost under one causal per trait assumption - LD-free mode
#' @description
#' Executes one iteration of the ColocBoost algorithm under the one-causal-variant-per-trait
#' assumption. This function provides a simplified approach to examine colocalization patterns
#' for traits without using any LD information.
#' 
#' @keywords cb_one_causal
#' @noRd
colocboost_diagLD <- function(cb_model, cb_model_para, cb_data,
                              jk_equiv_corr = 0,
                              jk_equiv_loglik = 0.1,
                              tau = 0.01,
                              learning_rate_decay = 1,
                              func_simplex = "z2z",
                              lambda = 0.5,
                              lambda_target_outcome = 1,
                              LD_free = FALSE){
  
  if (sum(cb_model_para$update_y == 1) == 1){
    
    pos.update <- which(cb_model_para$update_y == 1)
    update_jk <- rep(NA, cb_model_para$L+1)
    update_status <- rep(0, cb_model_para$L)
    real_update_jk <- rep(NA, cb_model_para$L)
    update_status[pos.update] <- -1
    correlation <- abs(cb_model[[pos.update]]$correlation)
    jk <- which(correlation == max(correlation))
    jk <- ifelse(length(jk) == 1, jk, sample(jk,1))
    update_jk[c(1,pos.update+1)] <- c(jk, jk)
    real_update_jk[pos.update] <- jk
    # - update cb_model_parameter and report results
    cb_model_para$update_temp <- list("update_status" = update_status, 
                                      "real_update_jk" = real_update_jk)
    cb_model_para$jk <- rbind(cb_model_para$jk, update_jk)
    cb_model_para$update_status <- cbind(cb_model_para$update_status, as.matrix(update_status))
    cb_model_para$real_update_jk <- rbind(cb_model_para$real_update_jk, real_update_jk)
    cb_model <- colocboost_update(cb_model, cb_model_para, cb_data,
                                  tau = tau, learning_rate_decay = learning_rate_decay, func_simplex = func_simplex,
                                  lambda = lambda, lambda_target_outcome = lambda_target_outcome, LD_free = LD_free)
    
  }
  
  if (sum(cb_model_para$update_y == 1) > 1){
    
    cb_model_tmp <- cb_model
    ######## - initial update to get weight
    pos.update <- which(cb_model_para$update_y == 1)
    model_update = cb_model[pos.update]
    X_dict <- cb_data$dict[pos.update]
    adj_dep <- sapply(pos.update, function(ii) cb_data$data[[ii]]$dependency)
    cor_vals_each <- do.call(cbind, lapply(model_update, function(cc) cc$correlation))
    abs_cor_vals_each <- sweep(abs(cor_vals_each), 2, adj_dep, `*`)
    jk_each <- apply(abs_cor_vals_each, 2, function(temp){
      jk_temp <- which(temp == max(temp))
      return(ifelse(length(jk_temp) == 1, jk_temp, sample(jk_temp,1)))
    })
    weights <- c()
    for (iy in pos.update){
        update_status <- rep(0, cb_model_para$L)
        real_update_jk <- rep(NA, cb_model_para$L)
        update_status[iy] <- -1
        jk <- jk_each[which(pos.update==iy)]
        real_update_jk[iy] <- jk
        cb_model_para$update_temp <- list("update_status" = update_status, 
                                          "real_update_jk" = real_update_jk)
        # - update cb_model
        cb_model_tmp <- colocboost_update(cb_model_tmp, cb_model_para, cb_data,
                                          tau = tau, learning_rate_decay = learning_rate_decay, func_simplex = func_simplex,
                                          lambda = lambda, lambda_target_outcome = lambda_target_outcome, LD_free = LD_free)
        weights <- rbind(weights, cb_model_tmp[[iy]]$weights_path)
    }
    ###### overlap weights
    overlap_pair <- check_pair_overlap(weights, coverage = 0.95)
    ###### change change loglikelihood
    change_each_pair <- check_pair_jkeach(jk_each, pos.update,
                                          model_update = model_update,
                                          cb_data = cb_data,
                                          X_dict = X_dict,
                                          jk_equiv_corr = jk_equiv_corr,
                                          jk_equiv_loglik = jk_equiv_loglik)
    change_each_pair <- change_each_pair * overlap_pair
    # define category with same jk
    temp <- sapply(1:nrow(change_each_pair), function(x){
      tt <- c(x, which(change_each_pair[x,] != 0))
      return(paste0(sort(tt), collapse = ";"))
    })
    temp <- merge_sets(temp)
    cate <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
    update_jk <- rep(NA, cb_model_para$L+1)
    real_update <- lapply(cate, function(cate_each){
      update_status <- rep(0, cb_model_para$L)
      real_update_jk <- rep(NA, cb_model_para$L)
      true_update <- as.numeric(cate_each)
      if (length(true_update)>1){
        update_status[pos.update[true_update]] <- 1
      } else {
        update_status[pos.update[true_update]] <- -1
      }
      abs_cor_vals_redefine <- rowSums(abs_cor_vals_each[,true_update,drop=FALSE])
      jk_redefine <- which(abs_cor_vals_redefine == max(abs_cor_vals_redefine))
      jk_redefine <- ifelse(length(jk_redefine) == 1, jk_redefine, sample(jk_redefine,1))
      real_update_jk[pos.update[true_update]] <- jk_redefine
      update_jk[c(1,pos.update+1)] <- c(jk_redefine, jk_each)
      return(list("update_status" = update_status, 
                  "real_update_jk" = real_update_jk,
                  "update_jk" = update_jk))
    })
    
    # - step 2: boost update
    for (i_update in 1:length(real_update)){
      # - update cb_model_parameter and report results
      cb_model_para$update_temp <- real_update[[i_update]]
      cb_model_para$jk <- rbind(cb_model_para$jk, real_update[[i_update]]$update_jk)
      cb_model_para$update_status <- cbind(cb_model_para$update_status, as.matrix(real_update[[i_update]]$update_status))
      cb_model_para$real_update_jk <- rbind(cb_model_para$real_update_jk, real_update[[i_update]]$real_update_jk)
      # - update cb_model
      cb_model <- colocboost_update(cb_model, cb_model_para, cb_data,
                                    tau = tau,
                                    learning_rate_decay = learning_rate_decay,
                                    func_simplex = func_simplex,
                                    lambda = lambda,
                                    lambda_target_outcome = lambda_target_outcome,
                                    LD_free = LD_free)
    }
  }
  # -- remove redundant parameters
  cb_model_para$num_updates <- 0
  rm_elements <- c("update_temp", "update_y")
  if (!is.null(cb_model_para$need_more)) {rm_elements <- c(rm_elements, "need_more")}
  if (!is.null(cb_model_para$true_stop)) {rm_elements <- c(rm_elements, "true_stop")}
  cb_model_para[rm_elements] <- NULL
  for (i in 1:length(cb_model)){ cb_model[[i]]$obj_path <- as.numeric(unlist(cb_model[[i]]$obj_path[-1]))}
  for (i in 1:length(cb_model)){ cb_model[[i]]$obj_single <- as.numeric(unlist(cb_model[[i]]$obj_single[-1]))}
  
  cb_obj <- list("cb_data" = cb_data, "cb_model" = cb_model, "cb_model_para" = cb_model_para)
  class(cb_obj) <- "colocboost"
  return(cb_obj)
  
}


#' Check for overlapping causal effects between trait pairs
#'
#' @description
#' This function evaluates pairs of traits to determine if they share overlapping
#' causal effects. It will used to create the trait group when LD information is not available.
#'
#' @keywords cb_one_causal
#' @noRd
check_pair_overlap <- function(weights, coverage = 0.95){
  
  overlap_pair <- matrix(NA, nrow = nrow(weights), ncol =  nrow(weights))
  for (i in 1:nrow(weights)){
    w1 <- weights[i,]
    cos1 <- get_in_cos(w1, coverage = coverage)[[1]]
    for (j in 1:nrow(weights)){
      if (j != i){
        w2 <- weights[j,]
        cos2 <- get_in_cos(w2, coverage = coverage)[[1]]
        
        overlap <- intersect(cos1, cos2)
        if ( length(overlap)!=0 ){
          cumsum1 <- sum(w1[overlap])
          cumsum2 <- sum(w2[overlap])
          overlap_pair[i,j] <- (cumsum1>0.5) & (cumsum2>0.5)
        } else {
          overlap_pair[i,j] <- FALSE
        }
        
      } else {
        overlap_pair[i,j] <- FALSE
      }
    }
  }
  overlap_pair <- overlap_pair + t(overlap_pair)
  return(overlap_pair)
  
}
