# Validate and Process All Input Data for ColocBoost

Internal function to validate and process both individual-level and
summary-level input data

## Usage

``` r
colocboost_validate_input_data(
  X = NULL,
  Y = NULL,
  sumstat = NULL,
  LD = NULL,
  dict_YX = NULL,
  dict_sumstatLD = NULL,
  effect_est = NULL,
  effect_se = NULL,
  effect_n = NULL,
  overlap_variables = FALSE,
  M = 500,
  min_abs_corr = 0.5,
  jk_equiv_corr = 0.8,
  jk_equiv_loglik = 1,
  func_simplex = "LD_z2z",
  cos_npc_cutoff = 0.2,
  npc_outcome_cutoff = 0.2
)
```

## Arguments

- X:

  A list of genotype matrices for different outcomes, or a single matrix
  if all outcomes share the same genotypes.

- Y:

  A list of vectors of outcomes or an N by L matrix if it is considered
  for the same X and multiple outcomes.

- sumstat:

  A list of data.frames of summary statistics.

- LD:

  A list of correlation matrices indicating the LD matrix for each
  genotype.

- dict_YX:

  A L by 2 matrix of dictionary for X and Y if there exist subsets of
  outcomes corresponding to the same X matrix.

- dict_sumstatLD:

  A L by 2 matrix of dictionary for sumstat and LD if there exist
  subsets of outcomes corresponding to the same sumstat.

- effect_est:

  Matrix of variable regression coefficients (i.e. regression beta
  values) in the genomic region

- effect_se:

  Matrix of standard errors associated with the beta values

- effect_n:

  A scalar or a vector of sample sizes for estimating regression
  coefficients.

- overlap_variables:

  If overlap_variables = TRUE, only perform colocalization in the
  overlapped region.

- M:

  The maximum number of gradient boosting rounds for each outcome
  (default is 500).

- min_abs_corr:

  Minimum absolute correlation allowed in a confidence set.

## Value

A list containing:

- X:

  Processed list of genotype matrices

- Y:

  Processed list of phenotype vectors

- yx_dict:

  Dictionary mapping Y to X

- keep_variable_individual:

  List of variable names for each X matrix

- sumstat:

  Processed list of summary statistics data.frames

- LD:

  Processed list of LD matrices

- sumstatLD_dict:

  Dictionary mapping sumstat to LD

- keep_variable_sumstat:

  List of variant names for each sumstat

- Z:

  List of z-scores for each outcome

- N_sumstat:

  List of sample sizes for each outcome

- Var_y:

  List of phenotype variances for each outcome

- SeBhat:

  List of standard errors for each outcome

- M_updated:

  Updated M value (may be changed if LD not provided)

- min_abs_corr_updated:

  Updated min_abs_corr value (may be changed if LD not provided)

- jk_equiv_corr_updated:

  Updated jk_equiv_corr value

- jk_equiv_loglik_updated:

  Updated jk_equiv_loglik value

- func_simplex_updated:

  Updated func_simplex value
