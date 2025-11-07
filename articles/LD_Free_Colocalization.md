# LD mismatch and LD-free Colocalization

This vignette demonstrates LD mismatch diagnosis in the `colocboost`
package and how to perform LD-mismatch and LD-free colocalization
analysis, when some traits completely lack LD information or share only
partial variant coverage with other traits.

``` r

library(colocboost)
```

## 1. LD mismatch diagnosis

The `colocboost` assumes that the LD matrix accurately estimates the
correlations among variants from the original GWAS genotype data.
Typically, the LD matrix comes from some public databases of genotypes
in a suitable reference population. An inaccurate LD matrix may lead to
unreliable colocalization results, especially if the LD matrix is
significantly different from the one estimated from the original
genotype data.

#### Why LD Mismatch Matters

An inaccurate LD matrix can cause inconsistencies between the summary
statistics and the reference LD matrix, leading to:

- Biased estimates of causal variants.
- Increased computational time due to slower algorithm convergence.
- Potentially misleading colocalization results.

ColocBoost provides diagnostic warnings for assessing the consistency of
the summary statistics with the reference LD matrix.

- Estimated residual variance of the model is negative or greater than
  phenotypic variance (`rtr < 0` or `rtr > var_y`; see details in
  Supplementary Note S3.5.2).
- The trait-specific gradient boosting model fails to converge.
- Change in trait-specific profile log-likelihood according to a CoS is
  negative (see details in Supplementary Note S3.5.3).

#### Example of including LD mismatch

In this example, we create a simulated dataset with LD mismatch by
changing the sign of Z-scores for 1% of variants for each trait.

``` r

# Create a simulated dataset with LD mismatch
data("Sumstat_5traits")
data("Ind_5traits")
LD <- get_cormat(Ind_5traits$X[[1]])

# Change sign of Z-score for 1% of variants for each trait by including mismatched LD
set.seed(123)
miss_prop <- 0.005 
sumstat <- lapply(Sumstat_5traits$sumstat, function(ss){
  p <- nrow(ss)
  pos_miss <- sample(1:p, ceiling(miss_prop * p))
  ss$z[pos_miss] <- -ss$z[pos_miss]
  return(ss)
})
```

#### Running ColocBoost with LD Mismatch

When running `colocboost` with an LD mismatch, you may encounter
diagnostic warnings. These warnings are not errors, and the analysis
will still proceed. However, the results may be less reliable due to the
mismatch, and the computational time may increase as the algorithm takes
longer to converge.

``` r

res <- colocboost(sumstat = sumstat, LD = LD)
#> Validating input data.
#> Starting gradient boosting algorithm.
#> Warning in colocboost_workhorse(cb_data, M = M, prioritize_jkstar =
#> prioritize_jkstar, : ColocBoost gradient boosting for outcome 5 did not
#> coverage in 500 iterations! Please check consistency between summary statistics
#> and LD matrix. See details in tutorial website
#> https://statfungen.github.io/colocboost/articles/.
#> Gradient boosting for outcome 4 converged after 523 iterations!
#> Gradient boosting at 1000 iterations, still updating.
#> Warning in colocboost_workhorse(cb_data, M = M, prioritize_jkstar =
#> prioritize_jkstar, : ColocBoost gradient boosting for outcome 1 did not
#> coverage in 500 iterations! Please check consistency between summary statistics
#> and LD matrix. See details in tutorial website
#> https://statfungen.github.io/colocboost/articles/.
#> Gradient boosting for outcome 2  stop since rtr < 0 or max(correlation) > 1 after 1213 iterations! Results for this locus are not stable, please check if mismatch between sumstat and LD! See details in tutorial website https://statfungen.github.io/colocboost/articles/.
#> Gradient boosting for outcome 3  stop since rtr < 0 or max(correlation) > 1 after 1475 iterations! Results for this locus are not stable, please check if mismatch between sumstat and LD! See details in tutorial website https://statfungen.github.io/colocboost/articles/.
#> Performing inference on colocalization events.
#> Warning in get_cos_profile(cs_beta, outcome_idx, X = cb_data$data[[X_dict]]$X,
#> : Warning message: potential sumstat & LD mismatch may happens for outcome 2 .
#> Using logLR = CoS(profile) - max(profile). Please check our website
#> https://statfungen.github.io/colocboost/articles/.
res$cos_details$cos$cos_index
#> $`cos1:y1_y2_y3_y4`
#> [1] 229 186 194 205 168
```

These warnings serve as diagnostic tools to alert users about potential
inconsistencies in the input data.

``` r

res$cos_details$cos_outcomes_npc
#> $`cos1:y1_y2_y3_y4`
#>    outcomes_index relative_logLR npc_outcome
#> Y3              3      2.3111393   0.9901696
#> Y1              1      1.9710178   0.9805913
#> Y4              4      0.9895934   0.8618184
#> Y2              2      0.0000000   0.0000000
```

**Note**: In the above example, the normalized probability of trait 2 is
0, indicating that colocalization with trait 2 may be less reliable due
to the LD mismatch. This is a warning, not an error, and the
colocalization analysis will still proceed. Therefore, in this case, we
suggest treating the colocalization of trait 2 with caution.

Potential solutions include:

- Using a more accurate LD matrix for trait 2.
- Excluding trait 2 from the analysis or colocalization results.
- Using the LD-mismatch or LD-free approach for colocalization analysis
  (below).
- Check the visualization of colocalization results including trait 2
  using `colocboost_plot(res)`. Remove the potential spurious signals
  when LD mismatch is detected using
  - `get_robust_colocalization(res, cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2)`
    to exclude the trait 2 in the above example if the signals are not
    reasonable.
  - `get_robust_colocalization(res, pvalue_cutoff = 1e-5, cos_npc_cutoff = 0, npc_outcome_cutoff = 0)`
    to include all colocalized traits with the larger marginal evidence,
    but the mismatch is detected.

## 2. LD-free and LD-mismatch colocalization analysis

When there is substantial discordance between the LD matrix and summary
statistics, the reliability of colocalization analysis may be
compromised. Such discordance can arise when the LD matrix and summary
statistics are derived from different populations or when the LD matrix
is estimated from a smaller or less representative reference sample.
This can lead to unexpected results, such as biased causal variant
identification or reduced accuracy in the analysis.

To address these challenges, ColocBoost provides two alternative
approaches for colocalization analysis with the assumption of one causal
variant per trait per region:

- **One iteration approach** (recommended): performing only 1 iteration
  of gradient boosting with the LD matrix ensures that:

  - The LD matrix is only used to check the equivalence among
    trait-specific best update variants.
  - The accuracy of the results is improved compared to completely
    ignoring the LD matrix.

This method is particularly useful when the LD matrix is mismatched but
still provides valuable insights into variant correlations.

``` r

# Perform only 1 iteration of gradient boosting with LD matrix
res_mismatch <- colocboost(sumstat = sumstat, LD = LD, M = 1)
#> Validating input data.
#> Starting gradient boosting algorithm.
#> Running ColocBoost with assumption of one causal per outcome per region!
#> Performing inference on colocalization events.
```

- **LD-free**: when the mismatch between the LD matrix and summary
  statistics is too large to be useful or when LD information is
  completely unavailable, ColocBoost provides an LD-free approach.

``` r

res_free <- colocboost(sumstat = sumstat)
#> Validating input data.
#> Warning in colocboost(sumstat = sumstat): Providing the LD for summary
#> statistics data is highly recommended. Without LD, only a single iteration will
#> be performed under the assumption of one causal variable per outcome.
#> Additionally, the purity of CoS cannot be evaluated!
#> Starting gradient boosting algorithm.
#> Running ColocBoost with assumption of one causal per outcome per region!
#> Performing inference on colocalization events.
```

While this method is computationally efficient, it has limitations due
to the strong assumption of a single causal variant per trait per
region. Users should interpret the results with caution, especially in
regions with complex LD structures or multiple causal variants.

ColocBoost also provides a flexibility to use HyPrColoc compatible
format for summary statistics without LD matrix.

``` r

# Loading the Dataset
data(Ind_5traits)
X <- Ind_5traits$X
Y <- Ind_5traits$Y

# Coverting to HyPrColoc compatible format
effect_est <- effect_se <- effect_n <- c()
for (i in 1:length(X)){
  x <- X[[i]]
  y <- Y[[i]]
  effect_n[i] <- length(y)
  output <- susieR::univariate_regression(X = x, y = y)
  effect_est <- cbind(effect_est, output$beta)
  effect_se <- cbind(effect_se, output$sebeta)
}
colnames(effect_est) <- colnames(effect_se) <- c("Y1", "Y2", "Y3", "Y4", "Y5")
rownames(effect_est) <- rownames(effect_se) <- colnames(X[[1]])


# Run colocboost
res <- colocboost(effect_est = effect_est, effect_se = effect_se, effect_n = effect_n)
#> Validating input data.
#> Warning in colocboost(effect_est = effect_est, effect_se = effect_se, effect_n
#> = effect_n): Providing the LD for summary statistics data is highly
#> recommended. Without LD, only a single iteration will be performed under the
#> assumption of one causal variable per outcome. Additionally, the purity of CoS
#> cannot be evaluated!
#> Starting gradient boosting algorithm.
#> Running ColocBoost with assumption of one causal per outcome per region!
#> Performing inference on colocalization events.

# Identified CoS
res$cos_details$cos$cos_index
#> $`cos1:y1_y3_y4`
#> [1] 186 205 194 168
```
