# Extract CoS at different coverage

`get_cos` extracts colocalization confidence sets (CoS) at different
coverage levels from ColocBoost results. When genotype data (X) or
correlation matrix (Xcorr) is provided, it can also calculate and filter
CoS based on purity statistics, ensuring that variants within each CoS
are sufficiently correlated.

## Usage

``` r
get_cos(
  cb_output,
  coverage = 0.95,
  X = NULL,
  Xcorr = NULL,
  n_purity = 100,
  min_abs_corr = 0.5,
  median_abs_corr = NULL
)
```

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- coverage:

  A number between 0 and 1 specifying the “coverage” of the estimated
  colocalization confidence sets (CoS) (default is 0.95).

- X:

  Genotype matrix of values of the p variables. Used to compute
  correlations if Xcorr is not provided.

- Xcorr:

  Correlation matrix of correlations between variables. Alternative to
  X.

- n_purity:

  The maximum number of CoS variables used in calculating the
  correlation (“purity”) statistics.

- min_abs_corr:

  The minimum absolute correlation value of variants in a CoS to be
  considered pass (“purity”) statistics.

- median_abs_corr:

  The median absolute correlation value of variants in a CoS to be
  considered pass (“purity”) statistics. When the number of variables
  included in the CoS is greater than this number, the CoS variables are
  randomly subsampled.

## Value

A list of indices of variables in each CoS.

## See also

Other colocboost_utilities:
[`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md),
[`get_cos_purity()`](https://statfungen.github.io/colocboost/reference/get_cos_purity.md),
[`get_cos_summary()`](https://statfungen.github.io/colocboost/reference/get_cos_summary.md),
[`get_hierarchical_clusters()`](https://statfungen.github.io/colocboost/reference/get_hierarchical_clusters.md),
[`get_ucos_summary()`](https://statfungen.github.io/colocboost/reference/get_ucos_summary.md)

## Examples

``` r
# colocboost example
set.seed(1)
N <- 1000
P <- 100
# Generate X with LD structure
sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
X <- MASS::mvrnorm(N, rep(0, P), sigma)
colnames(X) <- paste0("SNP", 1:P)
L <- 3
true_beta <- matrix(0, P, L)
true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
Y <- matrix(0, N, L)
for (l in 1:L) {
  Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
}
res <- colocboost(X = X, Y = Y)
#> Validating input data.
#> Starting gradient boosting algorithm.
#> Gradient boosting for outcome 1 converged after 98 iterations!
#> Gradient boosting for outcome 3 converged after 106 iterations!
#> Gradient boosting for outcome 2 converged after 107 iterations!
#> Performing inference on colocalization events.
#> Extracting colocalization results with pvalue_cutoff = 0.001, cos_npc_cutoff = 0.2, and npc_outcome_cutoff = 0.2.
#> Keep only CoS with cos_npc >= 0.2. For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < 0.001 and npc_outcome >0.2.
get_cos(res, coverage = 0.99, X = X)
#> $cos
#> $cos$`cos1:y1_y2_coverage_0.99`
#> [1] 10  9 11
#> 
#> 
#> $cos_purity
#> $cos_purity$min_abs_cor
#>                          cos1:y1_y2_coverage_0.99
#> cos1:y1_y2_coverage_0.99                0.8333458
#> 
#> $cos_purity$max_abs_cor
#>                          cos1:y1_y2_coverage_0.99
#> cos1:y1_y2_coverage_0.99                0.8816448
#> 
#> $cos_purity$median_abs_cor
#>                          cos1:y1_y2_coverage_0.99
#> cos1:y1_y2_coverage_0.99                0.9047605
#> 
#> 
get_cos(res, coverage = 0.99, X = X, min_abs_corr = 0.95)
#> $cos
#> NULL
#> 
#> $cos_purity
#> NULL
#> 
```
