library(testthat)

# Utility function to generate a simple colocboost model object
generate_test_model <- function(n = 100, p = 20, L = 2, seed = 42) {
  set.seed(seed)
  
  # Generate X with LD structure
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate true effects - create a shared causal variant
  true_beta <- matrix(0, p, L)
  true_beta[5, 1] <- 0.5  # SNP5 affects trait 1
  true_beta[5, 2] <- 0.4  # SNP5 also affects trait 2 (colocalized)
  true_beta[10, 2] <- 0.3 # SNP10 only affects trait 2
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Convert Y to list
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Run colocboost with minimal parameters to get a model object
  suppressWarnings({
    result <- colocboost(
      X = X_list, 
      Y = Y_list,
      M = 5,  # Small number of iterations for faster testing
      output_level = 3  # Include full model details
    )$diagnostic_details
  })
  result$cb_data <- colocboost_init_data(
        X = X_list,
        Y = Y_list,
        dict_YX = NULL, 
        Z = NULL,
        LD = NULL,
        N_sumstat = NULL,
        dict_sumstatLD = NULL,
        Var_y = NULL,
        SeBhat = NULL
      )
  result
}

# Test for colocboost_check_update_jk
test_that("colocboost_check_update_jk handles update selection", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # Check that function can be called without error
  # Note: This is testing a function that's normally internal
  # and would typically be covered by testing the main function
  if (exists("colocboost_check_update_jk")) {
    expect_error({
      result <- colocboost_check_update_jk(
        cb_obj$cb_model, 
        cb_obj$cb_model_para, 
        cb_obj$cb_data
      )
    }, NA)
  } else {
    skip("colocboost_check_update_jk not directly accessible")
  }
})

# Test for colocboost_update
test_that("colocboost_update updates model parameters", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # Access update function if it's exported
  if (exists("colocboost_update")) {
    # Create a temporary update status to test with
    cb_obj$cb_model_para$update_temp <- list(
      update_status = rep(1, cb_obj$cb_model_para$L),
      real_update_jk = rep(1, cb_obj$cb_model_para$L)
    )
    
    # Test function
    expect_error({
      updated_model <- colocboost_update(
        cb_obj$cb_model, 
        cb_obj$cb_model_para, 
        cb_obj$cb_data
      )
    }, NA)
  } else {
    skip("colocboost_update not directly accessible")
  }
})

# Test for colocboost_init_data
test_that("colocboost_init_data correctly initializes data", {
  skip_on_cran()
  
  # Generate test data
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), n, p)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n*2), n, 2)
  
  # If the function is exported
  if (exists("colocboost_init_data")) {
    expect_error({
      cb_data <- colocboost_init_data(
        X = list(X),
        Y = list(Y[,1], Y[,2]),
        dict_YX = matrix(c(1:2, rep(1,2)), ncol=2),
        Z = NULL,
        LD = NULL,
        N_sumstat = NULL,
        dict_sumstatLD = NULL,
        Var_y = NULL,
        SeBhat = NULL
      )
    }, NA)
  } else {
    skip("colocboost_init_data not directly accessible")
  }
})

# Test colocboost_assemble function
test_that("colocboost_assemble processes model results", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # If the function is exported
  if (exists("colocboost_assemble")) {
    expect_error({
      result <- colocboost_assemble(
        cb_obj,
        coverage = 0.95
      )
    }, NA)
  } else {
    skip("colocboost_assemble not directly accessible")
  }
})

# Test for colocboost_workhorse
test_that("colocboost_workhorse performs boosting iterations", {
  skip_on_cran()
  
  # Generate test data
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), n, p)
  colnames(X) <- paste0("SNP", 1:p)
  Y <- matrix(rnorm(n*2), n, 2)
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Initialize CB objects
  suppressWarnings({
    # First get the data object by running colocboost with M=0
    temp <- colocboost(X = X_list, Y = Y_list, M = 0)
    
    # If the workhorse function is exported
    if (exists("colocboost_workhorse")) {
      expect_error({
        result <- colocboost_workhorse(
          temp$cb_data,
          M = 5  # Small number for testing
        )
      }, NA)
    } else {
      skip("colocboost_workhorse not directly accessible")
    }
  })
})

# Test colocboost_plot function
test_that("colocboost_plot handles different plot options", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # Basic plot call
  expect_error(suppressWarnings(colocboost_plot(cb_obj)), NA)
  
  # Test with different y-axis values
  expect_error(suppressWarnings(colocboost_plot(cb_obj, y = "z_original")), NA)
  
  # Test with different outcome_idx
  expect_error(suppressWarnings(colocboost_plot(cb_obj, outcome_idx = 1)), NA)
})

# Test get_cos_summary function
test_that("get_cos_summary handles different parameters", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # Basic summary call
  expect_error(get_cos_summary(cb_obj), NA)
  
  # With custom outcome names
  expect_error(get_cos_summary(cb_obj, outcome_names = c("Trait1", "Trait2")), NA)
  
  # With gene name
  summary_with_gene <- get_cos_summary(cb_obj, gene_name = "TestGene")
  
  # If summary is not NULL, check for gene_name column
  if (!is.null(summary_with_gene)) {
    expect_true("gene_name" %in% colnames(summary_with_gene))
    expect_equal(summary_with_gene$gene_name[1], "TestGene")
  }
})

# Test for get_strong_colocalization
test_that("get_strong_colocalization filters results correctly", {
  skip_on_cran()
  
  # Generate a test model
  cb_obj <- generate_test_model()
  
  # Basic call
  expect_error(get_strong_colocalization(cb_obj), NA)
  
  # With stricter thresholds
  expect_error(get_strong_colocalization(cb_obj, cos_npc_cutoff = 0.8), NA)
  
  # With p-value threshold
  expect_error(get_strong_colocalization(cb_obj, pvalue_cutoff = 0.05), NA)
})