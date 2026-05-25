# Perform modularity-based hierarchical clustering for a correlation matrix

This function performs a modularity-based hierarchical clustering
approach to identify clusters from a correlation matrix.

## Usage

``` r
get_hierarchical_clusters(cormat, min_cluster_corr = 0.8)
```

## Arguments

- cormat:

  A correlation matrix.

- min_cluster_corr:

  The small correlation for the weights distributions across different
  iterations to be decided having only one cluster. Default is 0.8.

## Value

A list containing:

- cluster:

  A binary matrix indicating the cluster membership of each variable.

- Q_modularity:

  The modularity values for the identified clusters.

## See also

Other colocboost_utilities:
[`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md),
[`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md),
[`get_cos_purity()`](https://statfungen.github.io/colocboost/reference/get_cos_purity.md),
[`get_cos_summary()`](https://statfungen.github.io/colocboost/reference/get_cos_summary.md),
[`get_ucos_summary()`](https://statfungen.github.io/colocboost/reference/get_ucos_summary.md)

## Examples

``` r
# Example usage
set.seed(1)
N <- 100
P <- 4
sigma <- matrix(0.2, nrow = P, ncol = P)
diag(sigma) <- 1
sigma[1:2, 1:2] <- 0.9
sigma[3:4, 3:4] <- 0.9
X <- MASS::mvrnorm(N, rep(0, P), sigma)
cormat <- get_cormat(X)
clusters <- get_hierarchical_clusters(cormat)
clusters$cluster
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> [3,]    0    1
#> [4,]    0    1
clusters$Q_modularity
#> [1] 0.0000000 0.3618688 0.2714016 0.1809344
```
