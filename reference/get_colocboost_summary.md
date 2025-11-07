# Get summary tables from a ColocBoost output.

`get_colocboost_summary` get colocalization and trait-specific summary
table with or without the outcomes of interest.

## Usage

``` r
get_colocboost_summary(
  cb_output,
  summary_level = 1,
  outcome_names = NULL,
  interest_outcome = NULL,
  region_name = NULL,
  min_abs_corr_between_ucos = 0.5,
  median_abs_corr_between_ucos = 0.8
)
```

## Source

See detailed instructions in our tutorial portal:
<https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html>

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- summary_level:

  When `summary_level = 1`, return basic summary table for
  colocalization results. See details in `get_ucos_summary` function
  when `summary_level = 2`.

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

- min_abs_corr_between_ucos:

  Minimum absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.5.

- median_abs_corr_between_ucos:

  Median absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.8.

## Value

A list containing results from the ColocBoost analysis:

- When `summary_level = 1` (default):

  - `cos_summary`: A summary table for colocalization events with the
    following columns:

    - `focal_outcome`: The focal outcome being analyzed if exists.
      Otherwise, it is `FALSE`.

    - `colocalized_outcomes`: Colocalized outcomes for colocalization
      confidence set (CoS)

    - `cos_id`: Unique identifier for colocalization confidence set
      (CoS)

    - `purity`: Minimum absolute correlation of variables within
      colocalization confidence set (CoS)

    - `top_variable`: The variable with highest variant colocalization
      probability (VCP)

    - `top_variable_vcp`: Variant colocalization probability for the top
      variable

    - `cos_npc`: Normalized probability of colocalization

    - `min_npc_outcome`: Minimum normalized probability of colocalized
      traits

    - `n_variables`: Number of variables in colocalization confidence
      set (CoS)

    - `colocalized_index`: Indices of colocalized variables

    - `colocalized_variables`: List of colocalized variables

    - `colocalized_variables_vcp`: Variant colocalization probabilities
      for all colocalized variables

- When `summary_level = 2`:

  - `cos_summary`: As described above

  - `ucos_summary`: A summary table for trait-specific (uncolocalized)
    effects

- When `summary_level = 3`:

  - `cos_summary`: As described above

  - `ucos_summary`: A summary table for trait-specific (uncolocalized)
    effects

  - `ambiguous_cos_summary`: A summary table for ambiguous
    colocalization events from trait-specific effects

## Details

When `summary_level = 1`, additional details and examples are introduced
in
[`get_cos_summary`](https://statfungen.github.io/colocboost/reference/get_cos_summary.md).
When `summary_level = 2` or `summary_level = 3`, additional details for
trait-specific effects and ambiguous colocalization events are included.
See
[`get_ucos_summary`](https://statfungen.github.io/colocboost/reference/get_ucos_summary.md)
for details on these tables.

## See also

Other colocboost_inference:
[`get_ambiguous_colocalization()`](https://statfungen.github.io/colocboost/reference/get_ambiguous_colocalization.md),
[`get_robust_colocalization()`](https://statfungen.github.io/colocboost/reference/get_robust_colocalization.md)

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
get_colocboost_summary(res)
#> $cos_summary
#>   focal_outcome colocalized_outcomes     cos_id    purity top_variable
#> 1         FALSE               Y1; Y2 cos1:y1_y2 0.9047605        SNP10
#>   top_variable_vcp cos_npc min_npc_outcome n_variables colocalized_index
#> 1        0.8519757  0.9999          0.9999           2             10; 9
#>   colocalized_variables            colocalized_variables_vcp
#> 1           SNP10; SNP9 0.851975728205876; 0.118928231146593
#> 
```
