#' @rdname colocboost
#'
#' @title ColocBoost: A gradient boosting informed multi-omics xQTL colocalization method
#'
#' @description `colocboost` implements a proximity adaptive smoothing gradient boosting approach for multi-trait colocalization at gene loci,
#'              accommodating multiple causal variants. This method, introduced by Cao et al. (2025), is particularly suited for scaling
#'              to large datasets involving numerous molecular quantitative traits and disease traits.
#'              In brief, this function fits a multiple linear regression model \eqn{Y = XB + E} in matrix form. 
#'              ColocBoost can be generally used in multi-task variable selection regression problem.
#'
#' @details The function \code{colocboost} implements the proximity smoothed gradient boosting method from Cao et al (2025).
#' There is an additional step to help merge the confidence sets with small \code{between_putiry}
#' (default is 0.8) but within the same locus. This step addresses potential instabilities in linkage disequilibrium (LD) estimation
#' that may arise from small sample sizes or discrepancies in minor allele frequencies (MAF) across different confidence sets.
#'
#' @param X A list of genotype matrices for different outcomes, or a single matrix if all outcomes share the same genotypes.
#'          Each matrix should have column names, if sample sizes and variables possibly differing across matrices.
#' @param Y A list of vectors of outcomes or an N by L matrix if it is considered for the same X and multiple outcomes.
#' @param sumstat A list of data.frames of summary statistics.
#'                  The coloumns of data.frame should include either \code{z} or \code{beta}/\code{sebeta}.
#'                  \code{n} is the sample size for the summary statistics, it is highly recommendation to provide.
#'                  \code{variant} is required if sumstat for different outcomes do not have the same number of variables.
#'                  \code{var_y} is the variance of phenotype (default is 1 meaning that the Y is in the \dQuote{standarized} scale).
#' @param LD A list of correlation matrix indicating the LD matrix for each genotype. It also could be a single matrix if all sumstats were
#'           obtained from the same gentoypes.
#' @param dict_YX A L by 2 matrix of dictionary for \code{X} and \code{Y} if there exist subsets of outcomes corresponding to the same X matrix.
#'                  The first column should be 1:L for L outcomes. The second column should be the index of \code{X} corresponding to the outcome.
#'                  The innovation: do not provide the same matrix in \code{X} to reduce the computational burden.
#' @param dict_sumstatLD A L by 2 matrix of dictionary for \code{sumstat} and \code{LD} if there exist subsets of outcomes corresponding to the same sumstat.
#'                  The first column should be 1:L for L sumstat The second column should be the index of \code{LD} corresponding to the sumstat.
#'                  The innovation: do not provide the same matrix in \code{LD} to reduce the computational burden.
#' @param outcome_names The names of outcomes, which has the same order for Y.
#' @param target_outcome_idx The index of the target outcome if perform GWAS-xQTL ColocBoost
#' @param target_outcome_variables If \code{target_outcome_variables = TRUE}, only consider the variables exist in the target outcome.
#' @param overlap_variables If \code{overlap_variables = TRUE}, only perform colocalization in the overlapped region.
#' @param intercept If \code{intercept = TRUE}, the intercept is fitted. Setting \code{intercept = FALSE} is generally not recommended.
#' @param standardize If \code{standardize = TRUE}, standardize the columns of genotype and outcomes to unit variance.
#' @param effect_est Matrix of variable regression coefficients (i.e. regression beta values) in the genomic region
#' @param effect_se Matrix of standard errors associated with the beta values
#' @param effect_n A scalar or a vector of sample sizes for estimating regression coefficients. Highly recommended!
#' 
#' @param M The maximum number of gradient boosting rounds. If the number of outcomes are large, it will be automatically increased to a larger number.
#' @param stop_thresh The stop criterion for overall profile loglikelihood function.
#' @param tau The smooth parameter for proximity adaptive smoothing weights for the best update jk-star.
#' @param learning_rate_init The minimum learning rate for updating in each iteration.
#' @param learning_rate_decay The decayrate for learning rate. If the objective function is large at the early iterations,
#'                  we need to have the higher learning rate to improve the computational efficiency.
#' @param dynamic_learning_rate If \code{dynamic_learning_rate = TRUE}, the dynamic learning rate based on \code{learning_rate_init} and \code{learning_rate_decay} will be used in SEC.
#' @param prioritize_jkstar When \code{prioritize_jkstar = TRUE}, the selected outcomes will prioritize best update j_k^star in SEC.
#' @param func_compare The criterion when we update jk-star in SEC (default is "min_max").
#' @param jk_equiv_corr The LD cutoff between overall best update jk-star and marginal best update jk-l for lth outcome
#' @param jk_equiv_loglik The change of loglikelihood cutoff between overall best update jk-star and marginal best update jk-l for lth outcome
#' @param coloc_thresh The cutoff of checking if the best update jk-star is the potential causal variable for outcome l if jk-l is not similar to jk-star (used in Delayed SEC).
#' @param lambda The ratio \[0,1\] for z^2 and z in fun_prior simplex, defult is 0.5
#' @param lambda_target_outcome The ratio for z^2 and z in fun_prior simplex for the target outcome, default is 1
#' @param func_simplex The data-driven local association simplex \eqn{\delta} for smoothing the weights. Default is "LD_z2z" is the elastic net for z-score and also weighted by LD.
#' @param func_multi_test The alternative method to check the stop criteria. When \code{func_multi_test = "lfdr"}, boosting iterations will be stopped
#'                      if the local FDR for all variables are greater than \code{lfsr_max}.
#' @param stop_null The cutoff of nominal p-value when \code{func_multi_test = "Z"}.
#' @param multi_test_max The cutoff of the smallest FDR for pre-filtering the outcomes when \code{func_multi_test = "lfdr"} or \code{func_multi_test = "lfsr"}.
#' @param multi_test_thresh The cutoff of the smallest FDR for stop criteria when \code{func_multi_test = "lfdr"} or \code{func_multi_test = "lfsr"}.
#' @param ash_prior The prior distribution for calculating lfsr when \code{func_multi_test = "lfsr"}.
#' @param p.adjust.methods The adjusted pvalue method in stats:p.adj  when \code{func_multi_test = "fdr"}
#' @param residual_correlation The residual correlation based on the sample overlap, it is diagonal if it is NULL.
#' 
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated colocalization confidence sets (CoS) (default is 0.95).
#' @param min_cluster_corr The small correlation for the weights distributions across different iterations to be decided having only one cluster.
#' @param dedup If \code{dedup = TRUE}, the duplicate confidence sets will be removed in the post-processing.
#' @param overlap If \code{overlap = TRUE}, the overlapped confidence sets will be removed in the post-processing.
#' @param n_purity The maximum number of confidence set (CS) variables used in calculating the correlation (\dQuote{purity}) statistics.
#'                  When the number of variables included in the CS is greater than this number, the CS variables are randomly subsampled.
#' @param min_abs_corr Minimum absolute correlation allowed in a confidence set. The default is 0.5 corresponding to a squared correlation of 0.25,
#'                      which is a commonly used threshold for genotype data in genetic studies.
#' @param median_abs_corr An alternative "purity" threshold for the CS. Median correlation between pairs of variables in a CS less than this
#'                          threshold will be filtered out and not reported. When both min_abs_corr and median_abs_corr are set,
#'                          a CS will only be removed if it fails both filters. Default set to NULL but it is recommended to set it to 0.8 in practice.
#' @param median_cos_abs_corr Median absolute correlation between variants allowed to merge multiple colocalized sets. The default is 0.8 corresponding to a stringent threshold
#'                          to merge colocalized sets, which may resulting in a huge set.
#' @param tol A small, non-negative number specifying the convergence tolerance for checking the overlap of the variables in different sets.
#' @param merge_cos When \code{merge_cos = TRUE}, the sets for only one outcome will be merged if passed the \code{median_cos_abs_corr}.
#' @param sec_coverage_thresh A number between 0 and 1 specifying the weight in each SEC (default is 0.8).
#' @param weight_fudge_factor The strenght to integrate weight from differnt outcomes, default is 1.5
#' @param check_null The cut off value for change conditional objective function. Default is 0.1.
#' @param check_null_method The metric to check the null sets. Default is "profile"
#' @param check_null_max The smallest value of change of profile loglikelihood for each outcome.
#' @param weaker_effect If \code{weaker_effect = TRUE}, consider the weaker single effect due to coupling effects
#' @param LD_free When \code{LD_free = FALSE}, objective function doesn't include LD information.
#' @param output_level When \code{output_level = 2}, return the ucos details for the single specific effects. 
#'                      When \code{output_level = 3}, return the entire Colocboost model to diagnostic results (more space).
#'
#' @return A \code{"colocboost"} object with some or all of the following elements:
#'
#' \item{cos_summary}{A summary table for colocalization events.}
#' \item{vcp}{The variable colocalized probability for each variable.}
#' \item{cos_details}{A object with all information for colocalization results.}
#' \item{data_info}{A object with detailed information from input data}
#' \item{model_info}{A object with detailed information for colocboost model}
#'
#' @importFrom stats na.omit
#' @export

colocboost <- function(X = NULL, Y = NULL, # individual data
                       sumstat = NULL, LD = NULL, # summary statistics: either Z, bhat, sebhat, N, var_Y,
                       ###### - index dict for X match multiple Y / LD match multiple sumstat
                       dict_YX = NULL, # Y index for 1st column, X index for 2nd column
                       dict_sumstatLD = NULL, # sumstat index for 1st column, LD index for 2nd column
                       outcome_names = NULL, # the names of outcomes
                       ###### - GWAS-xQTL verison
                       target_outcome_idx = NULL, 
                       target_outcome_variables = TRUE,
                       overlap_variables = FALSE,
                       intercept=TRUE, # centered genotype and phenotype
                       standardize=TRUE, # standardized genotype and phenotype
                       ###### - HyPrColoc input
                       effect_est = NULL, # same as HyPrColoc, beta hat matrix: with rowname of variable names
                       effect_se = NULL, # same as HyPrColoc, sebeta hat matrix with rowname of variable names
                       effect_n = NULL,

                       ###### - Model Parameters
                       M = NULL, # maximum iteration time
                       stop_thresh = 1e-06, # stop criterion for profile_log and objective functions
                       tau = 0.01, # kernal_tau parameter for proximity smoothing weight
                       learning_rate_init = 0.01, # minimum step size for updating in each boosting round
                       learning_rate_decay = 1, # decay of learning rate
                       dynamic_learning_rate = TRUE, # if dynamic learning rate
                       prioritize_jkstar = TRUE, # if prioritize j_k^star for updating
                       func_compare = "min_max", # the criterion when we update j_k^star in Logic 3
                       jk_equiv_corr = 0.8, # check if jk_star ~ jk_r
                       jk_equiv_loglik = 1, # check if jk_star ~ jk_r.
                       coloc_thresh = 0.1,
                       lambda = 0.5, # the ratio for z^2 and z in weight penalty
                       lambda_target_outcome = 1, # the ratio for z^2 and z in weight penalty for target outcome
                       func_simplex = "LD_z2z", # data-driven association simplex
                       func_multi_test = "lfdr", 
                       stop_null = 1,
                       multi_test_max = 1,
                       multi_test_thresh = 1,
                       ash_prior = "normal", # only applicable if func_multi_test = lfsr
                       p.adjust.methods = "fdr",
                       residual_correlation = NULL, # phenotypic correlation, it is diagonal if it is NULL
                       
                       ###### - Post Inference Parameters
                       coverage = 0.95, # coverage of CoS
                       min_cluster_corr = 0.8, # only one cluster if all weights correlations bigger than this cut off
                       dedup = TRUE, # if remove the duplicate CoS in the post-processing
                       overlap = TRUE, # if remove the overlapped CoS
                       n_purity = 100, # the number of variables in purity
                       min_abs_corr = 0.5, # the cut off value of purity in each CoS
                       median_abs_corr = NULL,
                       median_cos_abs_corr = 0.8, # minimum LD between two CoS
                       tol = 1e-9, # tol for LD
                       merge_cos = TRUE, # if merge two sets for one outcome
                       sec_coverage_thresh = 0.8,
                       weight_fudge_factor = 1.5, # integrative weight
                       check_null = 0.1, # the cut off value for change conditional objective function
                       check_null_method = "profile", # the metric to check the null sets.
                       check_null_max = 0.02, # the smallest value of change of profile loglikelihood for each outcome.
                       weaker_effect = TRUE,
                       LD_free = FALSE,
                       output_level = 1){


    ###################### ---- one module for data object
    message("Starting checking the input data.")
    # - check if all missing
    check_individual <- ( is.null(X) & is.null(Y) )
    check_sumstat <- ( is.null(sumstat) & (is.null(effect_est) & is.null(effect_se)) )
    if (check_individual & check_sumstat){
        warning("Error: No individual data (X, Y) or summary statistics (sumstat) or (effect_est, effect_se) are provided! Please check!")
        return(NULL)
    }
    
    # - check individual level data
    if (!is.null(X) & !is.null(Y)){

        # --- check input
        if (is.data.frame(X)){ X <- as.matrix(X) }
        if (is.data.frame(Y)){ Y <- as.matrix(Y) }
        if (is.matrix(X)){ X <- list(X) }
        if (is.atomic(Y) && !is.list(Y)){
            Y <- as.matrix(Y)
            if (ncol(Y) == 1){
                Y <- list(Y)
            } else {
                n <- nrow(X[[1]])
                if (nrow(Y) == n){
                    Y <- lapply(1:ncol(Y), function(i) Y[,i,drop=FALSE])
                } else if (ncol(Y) == n){
                    Y <- lapply(1:nrow(Y), function(i) t(Y[i,,drop=FALSE]))
                } else {
                    stop("X and Y do not have the same sample size!")
                }
            }
        } else {
            Y <- lapply(1:length(Y), function(ii){
                if (is.null(dict_YX)){
                    idx <- ii
                } else {
                    idx <- dict_YX[ii,2]
                }
                n <- nrow(X[[idx]])
                y <- Y[[ii]]
                y <- as.matrix(y)
                if (nrow(y) == n) {
                    return(y)
                } else if (ncol(y) == n){
                    return(t(y))
                } else {
                    stop("X and Y do not have the same sample size!")
                }
            })
        }
        # --- check if variables in individual data
        p.ind <- unique(sapply(X, ncol))
        if (length(p.ind) != 1){
            variable.tmp <- sapply(X, function(xx){if (!is.null(colnames(xx))) TRUE else FALSE})
            if (!all(variable.tmp)){
                warning("Error: X matrices do not have the same number of variables. Provide variable names to the colnames of X matrix.")
                return(NULL)
            }
        } else {
            # --- check if there is only one X, default variable names as X_1, ..., X_p
            X <- lapply(X, function(xx){
              if (is.null(colnames(xx))){ colnames(xx) <- paste0("X_", 1:ncol(xx)) }
              return(xx)
            })
        }
        
        keep_variable_individual <- lapply(X, colnames)
        if (!is.list(X) & !is.list(Y)){
            warning("Error: Input X and Y must be the list containing genotype matrics and all phenotype vectors!")
            return(NULL)
        } else {
            if (length(X) == 1){
                yx_dict <- rep(1, length(Y))
            } else if (length(X) == length(Y)){
                yx_dict <- 1:length(Y)
            } else {
                if (is.null(dict_YX)){
                    warning("Error: Please provide the dict_YX since you have multiple Y but only few X!")
                    return(NULL)
                } else {
                    # - dict for Y to X mapping
                    yx_dict <- rep(NA, length(Y))
                    for (i in 1:length(Y)){
                        tmp <- unique(dict_YX[dict_YX[,1]==i,2])
                        if (length(tmp) == 0){
                            warning(paste("Error: You don't provide matched X for outcome", i))
                            return(NULL)
                        } else if (length(tmp) != 1){
                            warning(paste("Error: You provide different matched X for outcome", i))
                            return(NULL)
                        } else {
                            yx_dict[i] <- tmp
                        }
                    }
                    # - check if Y and X all mapped
                    if (max(yx_dict) > length(X)){
                        warning("Error: You don't provide enough X matrices!")
                        return(NULL)
                    }
                }
            }
        }
        # keep_variable_individual <- lapply(yx_dict, function(i) keep_variable_individual[[i]] )
        if (any(sapply(X, anyNA))) {
            warning("Error: Input X must not contain missing values.")
            return(NULL)
        }
        if (any(sapply(Y, anyNA))) {
            pos <- which(sapply(Y, anyNA))
            if (length(pos) == 1){
                Y_miss <- Y[[pos]]
                samples_kept = which(!is.na(Y_miss))
                Y[[pos]] = as.matrix(Y_miss[samples_kept])
                # - check if only y for x
                x_dict <- yx_dict[pos]
                X_miss <- X[[x_dict]]
                X_miss <- X_miss[samples_kept,]
                if (sum(yx_dict == x_dict) != 1){
                    x_dict <- max(yx_dict)+1
                    yx_dict[x_dict] = x_dict
                }
                X[[x_dict]] = X_miss
            } else {
                for (i in pos){
                    Y_miss <- Y[[i]]
                    samples_kept = which(!is.na(Y_miss))
                    Y[[i]] = as.matrix(Y_miss[samples_kept])
                    # - check if only y for x
                    x_dict <- yx_dict[i]
                    X_miss <- X[[x_dict]]
                    X_miss <- X_miss[samples_kept,]
                    if (sum(yx_dict == x_dict) != 1){
                        x_dict <- max(yx_dict)+1
                        yx_dict[x_dict] = x_dict
                    }
                    X[[x_dict]] = X_miss
                }
            }
        }
    } else {yx_dict <- keep_variable_individual <- NULL}

    # - check summary-level data
    if ( (!is.null(sumstat)) | (!is.null(effect_est) & !is.null(effect_se)) ){
        
        # --- check input of (effect_est, effect_se)
        if ( (!is.null(effect_est) & !is.null(effect_se)) ){
            if ( !all(dim(effect_est)==dim(effect_se)) ){
                warning("Error: effect_est and effect_se should be the same dimension! Please check!")
                return(NULL)
            }
            variables <- rownames(effect_est)
            effect_est <- as.matrix(effect_est)
            effect_se <- as.matrix(effect_se)
            if (is.null(variables)){
              variables <- paste0("variant_", 1:nrow(effect_est))
            }
            sumstat <- list()
            for (sum_iy in 1:ncol(effect_est)){
                sumstat_y <- data.frame(
                    "beta" = effect_est[, sum_iy],
                    "sebeta" = effect_se[, sum_iy],
                    "variant" = variables
                )
                if (!is.null(effect_n)){
                    if (length(effect_n)==1){
                        sumstat_y$n <- effect_n
                    } else if (length(effect_n) == ncol(effect_est)){
                        sumstat_y$n <- effect_n[sum_iy]
                    } else {
                        warning("effect_n does not match effect_est, ignoring sample size!")
                    }
                }
                sumstat[[sum_iy]] <- sumstat_y
            }
        }
        
        if ( !is.null(sumstat) ){
            if (is.data.frame(sumstat)){ sumstat <- list(sumstat) }
            if (!is.list(sumstat)){
              warning("Error: Input sumstat must be the list containing summary level data for all outcomes!")
              return(NULL)
            }
            # --- check if variables names in summary data
            variable.tmp <- sapply(sumstat, function(xx){ if (("variant" %in% colnames(xx))){ TRUE } else { FALSE } })
            if (!all(variable.tmp)){
                warning("sumstat do not contain variables names. Assigning default names as 'variant_1, ..., variant_P'. ",
                        "Recommend to provide variables names to the column named 'variant'.")
                sumstat <- lapply(sumstat, function(xx){
                    if (!("variant" %in% colnames(xx))){ xx$variant <- paste0("variant_", 1:nrow(xx)) }
                    return(xx)
                })
            }
        }
        keep_variable_sumstat <- lapply(sumstat, function(xx){ xx$variant})
        
        # --- check input of LD
        if (is.null(LD)){
          
          # if no LD input, set diagonal matrix to LD
          warning("Providing the LD for summary statistics data is highly recommended. ",
                  "Without LD, only a single iteration will be performed under the assumption of one causal variable per outcome. ",
                  "Additionally, the purity of CoS cannot be evaluated!")
          
          p.sumstat <- sapply(keep_variable_sumstat, length)
          p.unique <- unique(p.sumstat)
          if (length(p.unique)==1){
              ld <- diag(1, nrow = p.unique)
              colnames(ld) <- rownames(ld) <- keep_variable_sumstat[[1]]
              LD <- list(ld)
              sumstatLD_dict <- rep(1, length(sumstat))
          } else {
              LD <- lapply(keep_variable_sumstat, function(sn){
                  ld <- diag(1, nrow = length(sn))
                  colnames(ld) <- rownames(ld) <- sn
                  return(ld)
              })
              sumstatLD_dict <- 1:length(sumstat)
          }
          
          # change some algorithm parameters
          M <- 1 # one iteration
          min_abs_corr = 0 # remove purity checking
          jk_equiv_corr = 0
          jk_equiv_loglik = 0.1
          func_simplex = "only_z2z"
          
        } else {
          
          if (is.data.frame(LD)){ LD <- as.matrix(LD) }
          if (is.matrix(LD)){  LD <- list(LD) }
          if (length(LD) == 1){
            sumstatLD_dict <- rep(1, length(sumstat))
          } else if (length(LD) == length(sumstat)){
            sumstatLD_dict <- 1:length(sumstat)
          } else {
            if (is.null(dict_sumstatLD)){
              warning("Error: Please provide the dict_sumstatLD since you have multiple sumstat but only few LD!")
              return(NULL)
            } else {
              # - dict for sumstat to LD mapping
              sumstatLD_dict <- rep(NA, length(sumstat))
              for (i in 1:length(sumstat)){
                tmp <- unique(dict_sumstatLD[dict_sumstatLD[,1]==i,2])
                if (length(tmp) == 0){
                  warning(paste("Error: You don't provide matched LD for sumstat", i))
                  return(NULL)
                } else if (length(tmp) != 1){
                  warning(paste("Error: You provide multiple matched LD for sumstat", i))
                  return(NULL)
                } else {
                  sumstatLD_dict[i] <- tmp
                }
              }
              if (max(sumstatLD_dict) > length(LD)){
                warning("Error: You don't provide enough LD matrices!")
                return(NULL)
              }
            }
          }
          
        }
        
        # - checking sample size existency
        n_exist <- sapply(sumstat, function(ss){"n" %in% colnames(ss)})
        if (!all(n_exist)){
            # The sample size (n) is not provided, so use unadjusted z-scores.
            # The choice of n=2, yty=1 is mostly arbitrary except in that it
            # ensures var(y) = yty/(n-1) = 1
            p_no <- which(!n_exist)
            warning("Providing the sample size (n), or even a rough estimate of n, ",
                    "is highly recommended. Without n, the implicit assumption is ",
                    "n is large (Inf) and the effect sizes are small (close to zero).",
                    "outcome ", paste(p_no, collapse = ","), " in sumstat don't contain 'n'!")
        }
          
        Z <- N_sumstat <- Var_y <- SeBhat <- vector(mode='list', length=length(sumstat))
        for (i.summstat in 1:length(sumstat)){
          
          summstat_tmp <- sumstat[[i.summstat]]
          colVar <- colnames(summstat_tmp)
          # - check p-value
          if ( !("z" %in% colVar) ){
            if (!("beta" %in% colVar) || !("sebeta" %in% colVar)){
              warning("Error: Please provide either z or (bhat, sebhat)!")
              return(NULL)}
            bhat <- summstat_tmp[,"beta"]
            sebhat <- summstat_tmp[,"sebeta"]
            if (anyNA(bhat) || anyNA(sebhat)){
              warning("Error: bhat, sebhat cannot have missing values")
              return(NULL)}
            if (any(sebhat<=0)){
              warning("Error: sebhat cannot have zero or negative elements")
              return(NULL)}
            z <- bhat/sebhat
            SeBhat[[i.summstat]] <- sebhat
          } else {
            z <- summstat_tmp[,"z"]
          }
          if (anyNA(z)){
            warning(paste("summary statistic dataset", i.sumstat, "contains NA values that are replaced with 0"))
            z[is.na(z)] = 0
          }
          
          # - check N
          if ( !("n" %in% colVar) ){
            z <- z
          } else {
            
            n = summstat_tmp[,"n"]
            n <- as.numeric(na.omit(n)) # add hoc
            if (!all(n>1)){
              warning("Error: Sample size N must be greater than 1!") 
              return(NULL)}
            n = median(n)
            N_sumstat[[i.summstat]] <- n
            
            # When n is provided, compute the adjusted z-scores.
            z <- z * sqrt( (n-1)/(z^2+n-2) )
            if ("var_y" %in% colVar){
              var_y = unique(summstat_tmp[,"var_y"])
              # --- need to check
              if (length(var_y) != 1){ var_y = 1}
              if (var_y <= 0){
                warning("Error: var_y cannot have zero or negative elements")
                return(NULL)}
              Var_y[[i.summstat]] <- var_y
            }
            
          }
          Z[[i.summstat]] <- z
        }
    } else {
        Z <- N_sumstat <- Var_y <- SeBhat <- sumstatLD_dict <- keep_variable_sumstat <- NULL
    }
    # - initial colocboost object
    keep_variables <- c(keep_variable_individual, keep_variable_sumstat)
    overlapped_variables <- Reduce("intersect", keep_variables)
    mean_variables <- mean(sapply(keep_variables, length))
    min_variables <- min(sapply(keep_variables, length))
    if (min_variables < 100){
      warning("Warning message about the number of variables.\n",
              "The smallest number of variables across outcomes is ", min_variables, "<100.",
              " If this is what you expected, this is not a problem.",
              " If this is not what you expected, please check input data.")
    }
    if (length(overlapped_variables)<=1){
      warning("Error: No or only 1 overlapping variables were found across all outcomes, colocalization cannot be performed. ",
           "Please verify the variable names across different outcomes.")
      return(NULL)
    } else if ( (length(overlapped_variables)/mean_variables)<0.1 ){
      warning("Warning message about the overlapped variables.\n",
              "The average number of variables across outcomes is ", mean_variables, 
              ". But only ", length(overlapped_variables), " number of variables overlapped (<10%).\n",
              "If this is what you expected, this is not a problem.\n",
              "If this is not you expected, please check if the variable name matched across outcomes.")
    }
    cb_data <- colocboost_init_data(X = X, Y = Y, dict_YX = yx_dict,
                                    Z = Z, LD = LD, N_sumstat = N_sumstat, dict_sumstatLD = sumstatLD_dict,
                                    Var_y = Var_y, SeBhat = SeBhat,
                                    keep_variables = keep_variables,
                                    target_outcome_idx = target_outcome_idx,
                                    target_outcome_variables = target_outcome_variables,
                                    overlap_variables = overlap_variables,
                                    intercept = intercept,
                                    standardize = standardize,
                                    residual_correlation = residual_correlation)

    ##################  colocboost updates   ###################################
    message("Starting gradient boosting algorithm.")
    cb_obj <- colocboost_workhorse(cb_data,
                                   M = M,
                                   prioritize_jkstar = prioritize_jkstar,
                                   tau = tau,
                                   learning_rate_init = learning_rate_init,
                                   learning_rate_decay = learning_rate_decay,
                                   func_simplex = func_simplex,
                                   lambda = lambda,
                                   lambda_target_outcome = lambda_target_outcome,
                                   stop_thresh = stop_thresh,
                                   func_multi_test = func_multi_test,
                                   stop_null = stop_null,
                                   multi_test_max = multi_test_max,
                                   multi_test_thresh = multi_test_thresh,
                                   ash_prior = ash_prior,
                                   p.adjust.methods = p.adjust.methods,
                                   jk_equiv_corr = jk_equiv_corr,
                                   jk_equiv_loglik = jk_equiv_loglik,
                                   func_compare = func_compare,
                                   coloc_thresh = coloc_thresh,
                                   LD_free = LD_free,
                                   dynamic_learning_rate = dynamic_learning_rate,
                                   target_outcome_idx = target_outcome_idx,
                                   outcome_names = outcome_names)

    # --- post-processing of the colocboost updates
    message("Starting assemble analyses and results summary.")
    cb_output <- colocboost_assemble(cb_obj, 
                                     coverage = coverage,
                                     weight_fudge_factor = weight_fudge_factor,
                                     check_null = check_null,
                                     check_null_method = check_null_method,
                                     check_null_max = check_null_max,
                                     dedup = dedup,
                                     overlap = overlap,
                                     n_purity = n_purity,
                                     min_abs_corr = min_abs_corr,
                                     sec_coverage_thresh = sec_coverage_thresh,
                                     median_abs_corr = median_abs_corr,
                                     min_cluster_corr = min_cluster_corr,
                                     median_cos_abs_corr = median_cos_abs_corr,
                                     weaker_effect = weaker_effect,
                                     merge_cos =  merge_cos,
                                     tol = tol,
                                     output_level = output_level)
    class(cb_output) <- "colocboost"
    return(cb_output)
}

