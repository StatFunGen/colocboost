# Get colocalization summary table from a ColocBoost output.

`get_cos_summary` get the colocalization summary table with or without
the outcomes of interest.

## Usage

``` r
get_cos_summary(
  cb_output,
  outcome_names = NULL,
  interest_outcome = NULL,
  region_name = NULL
)
```

## Source

See detailed instructions in our tutorial portal:
<https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html>

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- outcome_names:

  Optional vector of names of outcomes, which has the same order as Y in
  the original analysis.

- interest_outcome:

  Optional vector specifying a subset of outcomes from `outcome_names`
  to focus on. When provided, only colocalization events that include at
  least one of these outcomes will be returned.

- region_name:

  Optional character string. When provided, adds a column with this gene
  name to the output table for easier filtering in downstream analyses.

## Value

A summary table for colocalization events with the following columns:

- focal_outcome:

  The focal outcome being analyzed if exists. Otherwise, it is `FALSE`.

- colocalized_outcomes:

  Colocalized outcomes for colocalization confidence set (CoS)

- cos_id:

  Unique identifier for colocalization confidence set (CoS)

- purity:

  Minimum absolute correlation of variables with in colocalization
  confidence set (CoS)

- top_variable:

  The variable with highest variant colocalization probability (VCP)

- top_variable_vcp:

  Variant colocalization probability for the top variable

- cos_npc:

  Normalized probability of colocalization

- min_npc_outcome:

  Minimum normalized probability of colocalized traits

- n_variables:

  Number of variables in colocalization confidence set (CoS)

- colocalized_index:

  Indices of colocalized variables

- colocalized_variables:

  List of colocalized variables

- colocalized_variables_vcp:

  Variant colocalization probabilities for all colocalized variables

## See also

Other colocboost_utilities:
[`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md),
[`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md),
[`get_cos_purity()`](https://statfungen.github.io/colocboost/reference/get_cos_purity.md),
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
#> Gradient boosting for outcome 1 converged after 106 iterations!
#> Gradient boosting for outcome 3 converged after 113 iterations!
#> Gradient boosting for outcome 2 converged after 117 iterations!
#> Performing inference on colocalization events.
#> Extracting colocalization results with pvalue_cutoff = 0.001, cos_npc_cutoff = 0.2, and npc_outcome_cutoff = 0.2.
#> Keep only CoS with cos_npc >= 0.2. For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < 0.001 and npc_outcome >0.2.
get_cos_summary(res)
#>   focal_outcome colocalized_outcomes     cos_id    purity top_variable
#> 1         FALSE               Y1; Y2 cos1:y1_y2 0.9067592        SNP10
#>   top_variable_vcp cos_npc min_npc_outcome n_variables colocalized_index
#> 1        0.8884754       1               1           2             10; 9
#>   colocalized_variables            colocalized_variables_vcp
#> 1           SNP10; SNP9 0.888475419309202; 0.101770958665919
```
