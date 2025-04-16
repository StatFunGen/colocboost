#' @title Set of functions for post inferences from ColocBoost
#'
#' @description
#' The `colocboost_post_inference` functions access basic properties inferences from a fitted ColocBoost model. This documentation serves as a summary for all related post-inference functions.
#'
#'
#' @details
#' The following functions are included in this set:
#' `get_cormat` a fast function to calulate correlation matrix from individual level data.
#' `check_null_post` a function to remove the spurious signals.
#' `get_purity` a function to calculate within-CoS purity.
#' `get_modularity` a function to calculate modularity for a given number of clusters.
#' `get_n_cluster` a function to get the number of clusters based on modularity-based hierarchical clustering method.
#' `get_between_purity` a function to calculate purity between two CoS.
#' `get_cos_evidence` a function to get the evidence of colocalization.
#' `w_purity` a function to calculate within-CoS purity for each single weight from one SEC.
#'
#' These functions are not exported individually and are accessed via `colocboost_post_inference`.
#'
#' @rdname colocboost_post_inference
#' @keywords cb_post_inference
#' @noRd
colocboost_post_inference <- function() {
  message("This function post inferences of colocboost output. See details for more information.")
}


#' @title A fast function to calculate correlation matrix (LD matrix) from individual level data
#' @description
#' This function calculates the correlation matrix (LD matrix) from individual level data.
#' 
#' @param X A matrix of individual level data.
#' @param intercepte A logical value indicating whether to include an intercept in the model. Default is FALSE.
#' 
#' @return A correlation matrix (LD matrix).
#' 
#' @rdname colocboost_post_inference
#' @keywords cb_post_inference
#' @export
get_cormat <- function(X, intercepte = FALSE){
  X = t(X)
  # Center each variable
  if (!intercepte){
    X = X - rowMeans(X)
  }
  # Standardize each variable
  X = X / sqrt(rowSums(X^2))
  # Calculate correlations
  cr = tcrossprod(X)
  return(cr)
}


#' Function to remove the spurious signals
#' @importFrom utils head tail
#' @keywords cb_post_inference
#' @noRd
check_null_post <- function(cb_obj, 
                            coloc_sets_temp,
                            coloc_outcomes,
                            check_null = 0.1,
                            check_null_method = "profile",
                            weaker_effect = TRUE){
  
  extract_last <- function(lst) { tail(lst, n = 1)}
  extract_first <- function(lst) { head(lst, n = 1)}
  
  get_profile <- function(cs_beta, X = NULL, Y = NULL, N = NULL,
                          XtX = NULL, YtY = NULL, XtY = NULL, miss_idx, adj_dep = 1){
    
    if (!is.null(X)){
      mean((Y-X%*%as.matrix(cs_beta))^2) *N/(N-1)
    } else if (!is.null(XtY)){
      scaling_factor <- if (!is.null(N)) (N - 1) else 1
      yty <- YtY / scaling_factor
      xtx <- XtX
      if (length(miss_idx)!=0){
        xty <- XtY[-miss_idx] / scaling_factor
        cs_beta <- cs_beta[-miss_idx]
      } else { xty <- XtY / scaling_factor  }
      
      ( yty - 2*sum(cs_beta*xty) + sum( (xtx %*% as.matrix(cs_beta)) * cs_beta ) ) * adj_dep
    }
    
  }
  
  get_cs_obj <- function(cs_beta, res, tau, func_simplex, lambda, adj_dep, LD_free,
                         X = NULL, Y = NULL, N = NULL,
                         XtX = NULL, YtY = NULL, XtY = NULL, miss_idx = NULL){
    
    correlation <- get_correlation(X = X, res = res, XtY = XtY, N = N, YtY = YtY, 
                                   XtX = XtX, beta_k = cs_beta, miss_idx = miss_idx)
    z <- get_z(correlation, n=N, res)
    abs_cor <- abs(correlation)
    jk <- which(abs_cor == max(abs_cor))
    jk <- ifelse(length(jk) == 1, jk, sample(jk,1))
    P <- length(z)
    ld_jk <- get_LD_jk(jk, X = X, XtX = XtX, N = N,
                       remain_idx = setdiff(1:P, miss_idx), P = P)
    ld_feature <- sqrt(abs(ld_jk))
    # - calculate delta
    delta <- boost_KL_delta(z = z, ld_feature = ld_feature, adj_dep=adj_dep,
                            func_simplex = func_simplex, lambda = lambda)
    scaling_factor <- if (!is.null(N)) (N - 1) else 1
    cov_Xtr <- if(!is.null(X)){
      t(X)%*%res / scaling_factor
    } else {
      res / scaling_factor
    }
    obj_ld <- if (LD_free) ld_feature else rep(1, length(ld_feature))
    obj_ld[miss_idx] <- 0
    exp_term <- adj_dep * obj_ld * (abs(cov_Xtr))
    return(tau * matrixStats::logSumExp(exp_term / tau + log(delta)))
    
  }
  
  update_res <- function(X = NULL, Y = NULL, XtX = NULL, XtY = NULL, N = NULL, cs_beta, miss_idx){
    if (!is.null(X)){
      return(Y - X%*%cs_beta)
    } else if (!is.null(XtX)){
      scaling.factor <- if (!is.null(N)) N-1 else 1
      xtx <- XtX / scaling.factor
      if (length(miss_idx)!=0){
        xty <- XtY[-miss_idx] / scaling.factor
        res.tmp <- rep(0, length(XtY))
        res.tmp[-miss_idx] <- xty - xtx%*%cs_beta[-miss_idx]
      } else {
        xty <- XtY / scaling.factor
        res.tmp <- xty - xtx%*%cs_beta
      }
      return(res.tmp)
    }
  }
  # - add hoc
  cut <- if(length(cb_obj$cb_data)==1) 0.2 else 1
  
  # ----- null filtering
  cb_data <- cb_obj$cb_data
  cs_change <- check_cs_change <- matrix(0, nrow = length(coloc_sets_temp), ncol = cb_obj$cb_model_para$L)
  colnames(cs_change) <- colnames(check_cs_change) <- paste0("change_obj_", 1:cb_obj$cb_model_para$L)
  max_change <- matrix(0, nrow = length(coloc_sets_temp), ncol = cb_obj$cb_model_para$L)
  for (i in 1:length(coloc_sets_temp)){
    cs_variants <- as.numeric(unlist(coloc_sets_temp[[i]]))
    for (j in coloc_outcomes){
      cs_beta <- cb_obj$cb_model[[j]]$beta
      cs_beta[cs_variants] <- 0
      X_dict <- cb_data$dict[j]
      adj_dep <- cb_data$data[[j]]$dependency
      if (check_null_method == "profile"){
        cs_profile <- get_profile(cs_beta, X = cb_data$data[[X_dict]]$X, Y = cb_data$data[[j]]$Y, 
                                  XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[j]]$XtY, 
                                  YtY = cb_data$data[[j]]$YtY, N = cb_data$data[[j]]$N, 
                                  miss_idx = cb_data$data[[j]]$variable_miss, adj_dep = adj_dep)
        last_profile <- extract_last(cb_obj$cb_model[[j]]$profile_loglike_each)
        change <- abs(cs_profile - last_profile) 
        # - add hoc
        if (min(cb_obj$cb_model[[j]]$multi_correction_univariate[cs_variants]) >= cut){ change <- 0 }
        # ------
        # - check_null
        check_cs_change[i,j] <- change  / diff(range(cb_obj$cb_model[[j]]$profile_loglike_each))
        cs_change[i,j] <- change # / diff(range(cb_obj$cb_model[[j]]$profile_loglike_each))
        first_profile <- extract_first(cb_obj$cb_model[[j]]$profile_loglike_each)
        max_change[i,j] <- (first_profile - last_profile) >= cb_obj$cb_model[[j]]$check_null_max
      } else if (check_null_method == "obj"){
        res <- update_res(X = cb_data$data[[X_dict]]$X, Y = cb_data$data[[j]]$Y, 
                          XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[j]]$XtY, 
                          N = cb_data$data[[j]]$N, cs_beta,
                          miss_idx = cb_data$data[[j]]$variable_miss)
        cs_obj <- get_cs_obj(cs_beta, res, cb_obj$cb_model_para$tau, cb_obj$cb_model_para$func_simplex, 
                             cb_obj$cb_model_para$lambda, adj_dep = adj_dep,
                             LD_free = cb_obj$cb_model_para$LD_free,
                             X = cb_data$data[[X_dict]]$X, N = cb_data$data[[j]]$N, 
                             XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[X_dict]]$XtY, 
                             YtY = cb_data$data[[X_dict]]$YtY,
                             miss_idx = cb_data$data[[j]]$variable_miss)
        last_obj <- min(cb_obj$cb_model[[j]]$obj_path)
        change <- abs(cs_obj - last_obj) 
        if (length(cb_obj$cb_model[[j]]$obj_path)==1){
          total_obj <- 1
        } else {
          total_obj <- diff(range(cb_obj$cb_model[[j]]$obj_path))
        }
        check_cs_change[i,j] <- change / total_obj
        cs_change[i,j] <- change
        max_change[i,j] <- total_obj >= cb_obj$cb_model[[j]]$check_null_max
      }
    }
  }
  if (!weaker_effect){ 
    check_cs_change <- cs_change 
    check_null_tmp <- sapply(1:cb_obj$cb_model_para$L, function(j) cb_obj$cb_model[[j]]$check_null_max)
  } else { check_null_tmp <- rep(check_null, cb_obj$cb_model_para$L) }
  cs_change <- as.data.frame(cs_change * max_change) # changed
  for (ii in 1:nrow(check_cs_change)){
    cs_tmp <- check_cs_change[ii,]
    check_cs_change[ii,] <- (cs_tmp>check_null_tmp)
  }
  is_non_null <- which(rowSums( check_cs_change * max_change ) != 0)
  
  ll = list("cs_change" = cs_change,
            "is_non_null" = is_non_null)
  return(ll)
}

#' Function to calculate within-CoS purity
#' @keywords cb_post_inference
#' @noRd
#' @importFrom stats na.omit
get_purity <- function(pos, X=NULL, Xcorr=NULL, N = NULL, n = 100) {
    get_upper_tri = Rfast::upper_tri
    get_median    = Rfast::med
    
    if (sum(is.na(pos))!=0){
        pos <- as.numeric(na.omit(pos))
    }
    
    if (length(pos) == 1)
        return(c(1,1,1))
    else {
        
        # Subsample the columns if necessary.
        if (length(pos) > n)
            pos = sample(pos,n)
        
        if (is.null(Xcorr)) {
            X_sub = X[,pos]
            X_sub = as.matrix(X_sub)
            corr <- suppressWarnings({get_cormat(X_sub)})
            corr[which(is.na(corr))] = 0
            value = abs(get_upper_tri(corr))
        } else {
            Xcorr <- Xcorr # if (!is.null(N)) Xcorr/(N-1) else Xcorr
            value = abs(get_upper_tri(Xcorr[pos,pos]))
        }
        return(c(min(value),
                 sum(value)/length(value),
                 get_median(value)))
    }
}

# ------ Calculate modularity ----------
#' Function to calculate modularity for a given number of clusters
#' @param Weight A matrix of weights
#' @param B A matrix of binary values indicating the clusters
#' @return The modularity value
#' @keywords cb_post_inference
#' @noRd
get_modularity <- function(Weight, B){
    if (dim(Weight)[1] == 1){
        Q <- 0
    } else {
        W_pos <- Weight * (Weight > 0)
        W_neg <- Weight * (Weight < 0)
        N <- dim(Weight)[1]
        K_pos <- colSums(W_pos)
        K_neg <- colSums(W_neg)
        m_pos <- sum(K_pos)
        m_neg <- sum(K_neg)
        m <- m_pos + m_neg
        cate <- B %*% t(B)
        if (m_pos == 0 & m_neg == 0){
            Q <- 0
        } else {
            if (m_pos == 0){
                Q_positive <- 0
                Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
            } else if (m_neg == 0){
                Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
                Q_negative <- 0
            } else {
                Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
                Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
            }
        }
        Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
    }
}

#' Function to get the number of clusters based on modularity-based hierarchical clustering method
#' @param hc A hierarchical clustering object
#' @param Sigma A matrix of weights
#' @param m The number of clusters
#' @param min_cluster_corr The minimum correlation threshold for clustering within the same cluster
#' @return A list containing the number of clusters and modularity values
#' @keywords cb_post_inference
#' @noRd
#' @importFrom stats cutree
get_n_cluster <- function(hc, Sigma, m=ncol(Sigma), min_cluster_corr = 0.8){
    if (min(Sigma) > min_cluster_corr){
        IND = 1
        Q = 1
    } else {
        Q <- c()
        if (ncol(Sigma) < 10){m = ncol(Sigma)}
        for(i in 1:m)
        {
            index=cutree(hc,i)
            B=sapply(1:i, function(t) as.numeric(index==t))
            Q[i] <- get_modularity(Sigma, B)
        }
        
        IND=which(Q==max(Q))
        L=length(IND)
        if (L>1) IND=IND[1]
    }
    return(list("n_cluster" = IND,
                "Qmodularity" = Q))
}

#' Function to calculate within-CoS purity for each single weight from one SEC
#' @keywords cb_post_inference
#' @noRd
w_purity <- function(weights, X=NULL, Xcorr=NULL, N = NULL, n = 100, coverage = 0.95, 
                     min_abs_corr = 0.5, median_abs_corr = NULL, miss_idx = NULL){
    
    tmp <- apply(weights, 2, w_cs, coverage = coverage)
    tmp_purity <- apply(tmp, 2, function(tt){
        pos <- which(tt == 1)
        # deal with missing snp here
        if (!is.null(Xcorr)){
            pos <- match(pos, setdiff(1:length(tmp), miss_idx))
        }
        get_purity(pos, X=X, Xcorr=Xcorr, N=N, n = n)
    })
    if (is.null(median_abs_corr)) {
        is_pure = which(tmp_purity[1,] >= min_abs_corr)
    } else {
        is_pure = which(tmp_purity[1,] >= min_abs_corr | tmp_purity[3,] >= median_abs_corr)
    }
    return(is_pure)
}

#' Calculate purity between two CoS
#' @keywords cb_post_inference
#' @noRd
#' @importFrom stats na.omit
get_between_purity <- function(pos1, pos2, X=NULL, Xcorr=NULL, N = NULL, miss_idx = NULL, P = NULL){
    
    get_matrix_mult <- function(X_sub1, X_sub2){
        
        X_sub1 = t(X_sub1)
        X_sub2 = t(X_sub2)
        # Standardize each variable
        X_sub1 = X_sub1 / sqrt(rowSums(X_sub1^2))
        X_sub1[is.nan(X_sub1)] = 0
        X_sub2 = X_sub2 / sqrt(rowSums(X_sub2^2))
        X_sub2[is.nan(X_sub2)] = 0
        # Calculate correlations
        cr = tcrossprod(X_sub1, X_sub2)
        return(cr)
        
    }
    get_median = Rfast::med
    
    if (is.null(Xcorr)){
        X_sub1 = X[,pos1,drop=FALSE]
        X_sub2 = X[,pos2,drop=FALSE]
        value <- abs(get_matrix_mult(X_sub1, X_sub2))
        
    } else {
        pos1 <- na.omit(match(pos1, setdiff(1:P, miss_idx)))
        pos2 <- na.omit(match(pos2, setdiff(1:P, miss_idx)))
        # scaling_factor <- if (!is.null(N)) N-1 else 1
        if (length(pos1)!=0 & length(pos2)!=0){
            value = abs(Xcorr[pos1,pos2]) # / scaling_factor
        } else {
            value = 0
        }
        
    }
    return(c(min(value), max(value), get_median(value)))
}

#' Function to get the evidence of colocalization
#' @keywords cb_post_inference
#' @noRd
#' @importFrom stats var
#' @importFrom utils tail
get_cos_evidence <- function(cb_obj, coloc_out, data_info){

  get_cos_config <- function(w, config_idx, weight_fudge_factor = 1.5, coverage = 0.95){
    
    w_outcome <- colnames(w)
    config_outcome <- paste0("outcome", config_idx)
    pos <- which(w_outcome %in% config_outcome)
    w_config <- w[, pos, drop=FALSE]
    int_w <- get_integrated_weight(w_config, weight_fudge_factor = weight_fudge_factor)
    unlist(get_in_cos(int_w, coverage = coverage))
    
  }
  
  get_cos_profile <- function(cs_beta, outcome_idx, X = NULL, Y = NULL, N = NULL,
                          XtX = NULL, YtY = NULL, XtY = NULL, miss_idx = NULL, adj_dep = 1){
    
    if (!is.null(X)){
      cos_profile <- mean((Y-X%*%as.matrix(cs_beta))^2) *N/(N-1)
      yty <- var(Y)
    } else if (!is.null(XtY)){
      scaling_factor <- if (!is.null(N)) (N - 1) else 1
      yty <- YtY / scaling_factor
      xtx <- XtX
      if (length(miss_idx)!=0){
        xty <- XtY[-miss_idx] / scaling_factor
        cs_beta <- cs_beta[-miss_idx]
      } else { xty <- XtY / scaling_factor  }
      
      cos_profile <- ( yty - 2*sum(cs_beta*xty) + sum( (xtx %*% as.matrix(cs_beta)) * cs_beta ) ) * adj_dep
    }
    delta <- yty - cos_profile
    if (delta <= 0){
      warning(paste("Warning message: potential sumstat & LD mismatch may happens for outcome", outcome_idx,
                    ". Using logLR  = CoS(profile) - max(profile). Please check our website."))
    }
    cos_profile
  }
  
  get_outcome_profile_change <- function(outcome_idx, cos, cb_obj, check_null_max){
    extract_last <- function(lst) { tail(lst, n = 1)}
    cb_data <- cb_obj$cb_data
    cs_beta <- rep(0, cb_obj$cb_model_para$P)
    cs_beta[cos] <- cb_obj$cb_model[[outcome_idx]]$beta[cos]
    X_dict <- cb_data$dict[outcome_idx]
    cos_profile <- get_cos_profile(cs_beta, outcome_idx, 
                                   X = cb_data$data[[X_dict]]$X, Y = cb_data$data[[outcome_idx]]$Y,
                                   XtX = cb_data$data[[X_dict]]$XtX, XtY = cb_data$data[[outcome_idx]]$XtY,
                                   YtY = cb_data$data[[outcome_idx]]$YtY, N = cb_data$data[[outcome_idx]]$N,
                                   miss_idx = cb_data$data[[outcome_idx]]$variable_miss,
                                   adj_dep = cb_data$data[[outcome_idx]]$dependency)
    max_profile <- max(cb_obj$cb_model[[outcome_idx]]$profile_loglike_each)
    ifelse(max_profile < cos_profile, 0, max_profile - cos_profile)
  }
  
  # - Calculate best configuration likelihood explained by minimal configuration
  get_normalization_evidence <- function(profile_change, null_max, outcomes, outcome_names) {
    # Define the exponential likelihood ratio normalization (ELRN)
    logLR_normalization <- function(ratio) { 1 - exp( - 2*ratio ) }
    
    ratio <- profile_change / null_max
    prob <- logLR_normalization(ratio)
    df <- data.frame(outcomes_index = outcomes, relative_logLR = ratio, npc_outcome = prob)
    rownames(df) <- outcome_names[outcomes]
    sorted_df <- df[order(-df$relative_logLR), ]
    return(sorted_df)
  }
  
  get_npuc <- function(npc_outcome){
    max_idx <- which.max(npc_outcome)
    npc_outcome[max_idx] * prod(1 - npc_outcome[-max_idx])
  }
  
  avWeight <- coloc_out$avWeight
  cs_change <- coloc_out$cs_change
  check_null_max <- sapply(cb_obj$cb_model, function(cb) cb$check_null_max)
  outcome_names <- data_info$outcome_info$outcome_names
  n_cos <- length(avWeight)
  npc <- c() 
  normalization_evidence <- list()
  for (i in 1:n_cos){
    w <- avWeight[[i]]
    outcomes <- coloc_out$coloc_outcomes[[i]]
    # most likely cos
    cos <- get_cos_config(w, outcomes,weight_fudge_factor = cb_obj$cb_model_para$weight_fudge_factor,coverage = cb_obj$cb_model_para$coverage)
    profile_change_outcome <- sapply(outcomes, function(outcome_idx){
      get_outcome_profile_change(outcome_idx, cos, cb_obj, check_null_max)
    })
    normalization_evidence[[i]] <- get_normalization_evidence(profile_change = profile_change_outcome, 
                                                              null_max = check_null_max[outcomes], 
                                                              outcomes, outcome_names)
    
    # - calcualte normalization probability of colocalization (NPC)
    npuc <- get_npuc(normalization_evidence[[i]]$npc_outcome)
    npc[i] <- 1 - npuc
  }
  names(normalization_evidence) <- names(npc) <- names(coloc_out$cos)
  return(list(normalization_evidence=normalization_evidence, npc=npc))
}


