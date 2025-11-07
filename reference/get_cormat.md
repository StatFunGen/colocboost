# A fast function to calculate correlation matrix (LD matrix) from individual level data

This function calculates the correlation matrix (LD matrix) from
individual level data.

## Usage

``` r
get_cormat(X, intercepte = TRUE)
```

## Arguments

- X:

  A matrix of individual level data.

- intercepte:

  A logical value indicating whether to include an intercept in the
  model. Default is FALSE.

## Value

A correlation matrix (LD matrix).

## See also

Other colocboost_utilities:
[`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md),
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
cormat <- get_cormat(X)
```
