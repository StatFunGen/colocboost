# Get trait-specific summary table from a ColocBoost output.

`get_ucos_summary` produces a trait-specific summary table for
uncolocalized (single-trait) associations from ColocBoost results. This
is particularly useful for examining trait-specific signals or for
summarizing results from single-trait FineBoost analyses.

## Usage

``` r
get_ucos_summary(
  cb_output,
  outcome_names = NULL,
  region_name = NULL,
  ambiguous_cos = FALSE,
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

- outcome_names:

  Optional vector of names of outcomes, which has the same order as Y in
  the original analysis.

- region_name:

  Optional character string. When provided, adds a column with this gene
  name to the output table for easier filtering in downstream analyses.

- ambiguous_cos:

  Logical indicating whether to include ambiguous colocalization events.
  The default is FALSE.

- min_abs_corr_between_ucos:

  Minimum absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.5.

- median_abs_corr_between_ucos:

  Median absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.8.

## Value

A list containing:

- `ucos_summary`: A summary table for trait-specific, uncolocalized
  associations with the following columns:

  - `outcomes`: Outcome being analyzed

  - `ucos_id`: Unique identifier for trait-specific confidence sets

  - `purity`: Minimum absolute correlation of variables within
    trait-specific confidence sets

  - `top_variable`: The variable with highest variant-level probability
    of association (VPA)

  - `top_variable_vpa`: Variant-level probability of association (VPA)
    for the top variable

  - `ucos_npc`: Normalized probability of causal association for the
    trait-specific confidence set

  - `n_variables`: Number of variables in trait-specific confidence set

  - `ucos_index`: Indices of variables in the trait-specific confidence
    set

  - `ucos_variables`: List of variables in the trait-specific confidence
    set

  - `ucos_variables_vpa`: Variant-level probability of association (VPA)
    for all variables in the confidence set

  - `region_name`: Region name if provided through the region_name
    parameter

- `ambiguous_cos_summary`: A summary table for ambiguous colocalization
  events with the following columns:

  - `outcomes`: Outcome in the ambiguous colocalization event

  - `ucos_id`: Unique identifiers for the ambiguous event

  - `min_between_purity`: Minimum absolute correlation between variables
    across trait-specific sets in the ambiguous event

  - `median_between_purity`: Median absolute correlation between
    variables across trait-specific sets in the ambiguous event

  - `overlap_idx`: Indices of variables that overlap between ambiguous
    trait-specific sets

  - `overlap_variables`: Names of variables that overlap between
    ambiguous trait-specific sets

  - `n_recalibrated_variables`: Number of variables in the recalibrated
    colocalization set from an ambiguous event

  - `recalibrated_index`: Indices of variables in the recalibrated
    colocalization set from an ambiguous event

  - `recalibrated_variables`: Names of variables in the recalibrated
    colocalization set from an ambiguous event

  - `recalibrated_variables_vcp`: Variant colocalization probabilities
    for recalibrated variables from an ambiguous event

  - `region_name`: Region name if provided through the region_name
    parameter

## See also

Other colocboost_utilities:
[`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md),
[`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md),
[`get_cos_purity()`](https://statfungen.github.io/colocboost/reference/get_cos_purity.md),
[`get_cos_summary()`](https://statfungen.github.io/colocboost/reference/get_cos_summary.md),
[`get_hierarchical_clusters()`](https://statfungen.github.io/colocboost/reference/get_hierarchical_clusters.md)

## Examples

``` r
# colocboost example with single trait analysis
set.seed(1)
N <- 1000
P <- 100
# Generate X with LD structure
sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
X <- MASS::mvrnorm(N, rep(0, P), sigma)
colnames(X) <- paste0("SNP", 1:P)
L <- 1  # Only one trait for single-trait analysis
true_beta <- matrix(0, P, L)
true_beta[10, 1] <- 0.5 # SNP10 affects the trait
true_beta[80, 1] <- 0.2 # SNP11 also affects the trait but with lower effect
Y <- X %*% true_beta + rnorm(N, 0, 1)
res <- colocboost(X = X, Y = Y, output_level = 2)
#> Validating input data.
#> Starting gradient boosting algorithm.
#> Gradient boosting for outcome 1 converged after 58 iterations!
#> Performing inference on colocalization events.
#> No colocalization results in this region!
#> Extracting outcome-specific (uncolocalized) results with pvalue_cutoff = 1e-05, and npc_outcome_cutoff = 0.2.
#> For each uCoS, keep the outcome-specific (uncolocalized) events that pvalue of variants for the outcome < 1e-05 and npc_outcome >0.2.
# Get the trait-specifc effect summary
get_ucos_summary(res)
#>   outcomes  ucos_id    purity top_variable top_variable_vpa n_variables
#> 1       Y1 ucos1:y1 0.9067592        SNP10        0.7702073           2
#> 2       Y1 ucos2:y1 0.6623469        SNP80        0.5205935           4
#>       ucos_index             ucos_variables
#> 1          10; 9                SNP10; SNP9
#> 2 80; 81; 77; 79 SNP80; SNP81; SNP77; SNP79
#>                                                             ucos_variables_vpa
#> 1                                         0.770207273736694; 0.227530535623682
#> 2 0.520593464150777; 0.344977939462854; 0.0617682672655711; 0.0355614888280541
```
