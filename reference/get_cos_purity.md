# Calculate purity within and in-between CoS

Calculate purity statistics between all pairs of colocalization
confidence sets (CoS)

## Usage

``` r
get_cos_purity(cos, X = NULL, Xcorr = NULL, n_purity = 100)
```

## Arguments

- cos:

  List of variables in CoS

- X:

  Genotype matrix of values of the p variables. Used to compute
  correlations if Xcorr is not provided.

- Xcorr:

  Correlation matrix of correlations between variables. Alternative to
  X.

- n_purity:

  The maximum number of CoS variables used in calculating the
  correlation (“purity”) statistics. When the number of variables
  included in the CoS is greater than this number, the CoS variables are
  randomly subsampled.

## Value

A list containing three matrices (min_abs_cor, max_abs_cor,
median_abs_cor) with purity statistics for all pairs of CoS. Diagonal
elements represent within-CoS purity.

## See also

Other colocboost_utilities:
[`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md),
[`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md),
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
true_beta[10, 1] <- 0.5 
true_beta[10, 2] <- 0.4 
true_beta[50, 2] <- 0.3 
true_beta[80, 3] <- 0.6 
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
cos_res <- get_cos(res, coverage = 0.8)
get_cos_purity(cos_res$cos, X = X)
#> $min_abs_cor
#>                         cos1:y1_y2_coverage_0.8
#> cos1:y1_y2_coverage_0.8                       1
#> 
#> $max_abs_cor
#>                         cos1:y1_y2_coverage_0.8
#> cos1:y1_y2_coverage_0.8                       1
#> 
#> $median_abs_cor
#>                         cos1:y1_y2_coverage_0.8
#> cos1:y1_y2_coverage_0.8                       1
#> 
```
