# Recalibrate and summarize robust uncolocalized events.

`get_robust_ucos` get the uncolocalized events by discarding the weaker
signals if "ucos_details" with "output_level = 2" exist.

## Usage

``` r
get_robust_ucos(cb_output, npc_outcome_cutoff = 0.2, pvalue_cutoff = NULL)
```

## Source

See detailed instructions in our tutorial portal:
<https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html>

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- npc_outcome_cutoff:

  Minimum threshold of normalized probability of colocalized traits in
  each CoS.

- pvalue_cutoff:

  Maximum threshold of marginal p-values of colocalized variants on
  colocalized traits in each CoS.

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

- ucos_details:

  A object with information for trait-specific effects if exists after
  removing weaker signals.

## See also

Other colocboost_inference:
[`get_ambiguous_colocalization()`](https://statfungen.github.io/colocboost/reference/get_ambiguous_colocalization.md),
[`get_colocboost_summary()`](https://statfungen.github.io/colocboost/reference/get_colocboost_summary.md),
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
res <- colocboost(X = X, Y = Y, output_level = 2)
#> Validating input data.
#> Starting gradient boosting algorithm.
#> Gradient boosting for outcome 1 converged after 98 iterations!
#> Gradient boosting for outcome 3 converged after 106 iterations!
#> Gradient boosting for outcome 2 converged after 107 iterations!
#> Performing inference on colocalization events.
#> Extracting colocalization results with pvalue_cutoff = 0.001, cos_npc_cutoff = 0.2, and npc_outcome_cutoff = 0.2.
#> Keep only CoS with cos_npc >= 0.2. For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < 0.001 and npc_outcome >0.2.
#> Extracting outcome-specific (uncolocalized) results with pvalue_cutoff = 1e-05, and npc_outcome_cutoff = 0.2.
#> For each uCoS, keep the outcome-specific (uncolocalized) events that pvalue of variants for the outcome < 1e-05 and npc_outcome >0.2.
res$ucos_details$ucos$ucos_index
#> $`ucos1:y2`
#> [1] 50
#> 
#> $`ucos2:y3`
#> [1] 80
#> 
filter_res <- get_robust_ucos(res, npc_outcome_cutoff = 0.2, pvalue_cutoff = 1e-6)
#> Extracting outcome-specific (uncolocalized) results with pvalue_cutoff = 1e-06, and npc_outcome_cutoff = 0.2.
#> For each uCoS, keep the outcome-specific (uncolocalized) events that pvalue of variants for the outcome < 1e-06 and npc_outcome >0.2.
filter_res$ucos_details$ucos$ucos_index
#> $`ucos1:y2`
#> [1] 50
#> 
#> $`ucos2:y3`
#> [1] 80
#> 
```
