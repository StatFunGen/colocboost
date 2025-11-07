# Recalibrate and summarize robust colocalization events.

`get_robust_colocalization` get the colocalization by discarding the
weaker colocalization events or colocalized outcomes

## Usage

``` r
get_robust_colocalization(
  cb_output,
  cos_npc_cutoff = 0.5,
  npc_outcome_cutoff = 0.2,
  pvalue_cutoff = NULL,
  weight_fudge_factor = 1.5,
  coverage = 0.95
)
```

## Source

See detailed instructions in our tutorial portal:
<https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html>

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- cos_npc_cutoff:

  Minimum threshold of normalized probability of colocalization (NPC)
  for CoS.

- npc_outcome_cutoff:

  Minimum threshold of normalized probability of colocalized traits in
  each CoS.

- pvalue_cutoff:

  Maximum threshold of marginal p-values of colocalized variants on
  colocalized traits in each CoS.

- weight_fudge_factor:

  The strength to integrate weight from different outcomes, default is
  1.5

- coverage:

  A number between 0 and 1 specifying the “coverage” of the estimated
  colocalization confidence sets (CoS) (default is 0.95).

## Value

A `"colocboost"` object with some or all of the following elements:

- cos_summary:

  A summary table for colocalization events.

- vcp:

  The variable colocalized probability for each variable.

- cos_details:

  A object with all information for colocalization results.

- data_info:

  A object with detailed information from input data

- model_info:

  A object with detailed information for colocboost model

- ucos_from_cos:

  A object with information for trait-specific effects if exists after
  removing weaker signals.

## See also

Other colocboost_inference:
[`get_ambiguous_colocalization()`](https://statfungen.github.io/colocboost/reference/get_ambiguous_colocalization.md),
[`get_colocboost_summary()`](https://statfungen.github.io/colocboost/reference/get_colocboost_summary.md)

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
res$cos_details$cos$cos_index
#> $`cos1:y1_y2`
#> [1] 10  9
#> 
filter_res <- get_robust_colocalization(res, cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2)
#> Extracting colocalization results with cos_npc_cutoff = 0.5 and npc_outcome_cutoff = 0.2.
#> Keep only CoS with cos_npc >= 0.5. For each CoS, keep the outcomes configurations that the npc_outcome >= 0.2.
filter_res$cos_details$cos$cos_index
#> $`cos1:y1_y2`
#> [1] 10  9
#> 
```
