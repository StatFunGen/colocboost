# Get ambiguous colocalization events from trait-specific (uncolocalized) effects.

`get_ambiguous_colocalization` get the colocalization by discarding the
weaker colocalization events or colocalized outcomes

## Usage

``` r
get_ambiguous_colocalization(
  cb_output,
  min_abs_corr_between_ucos = 0.5,
  median_abs_corr_between_ucos = 0.8,
  tol = 1e-09
)
```

## Source

See detailed instructions in our tutorial portal:
<https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html>

## Arguments

- cb_output:

  Output object from `colocboost` analysis

- min_abs_corr_between_ucos:

  Minimum absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.5.

- median_abs_corr_between_ucos:

  Median absolute correlation for variants across two trait-specific
  (uncolocalized) effects to be considered colocalized. The default is
  0.8.

- tol:

  A small, non-negative number specifying the convergence tolerance for
  checking the overlap of the variables in different sets.

## Value

A `"colocboost"` object of colocboost output with additional elements:

- ambiguous_cos:

  If exists, a list of ambiguous trait-specific (uncolocalized) effects.

## See also

Other colocboost_inference:
[`get_colocboost_summary()`](https://statfungen.github.io/colocboost/reference/get_colocboost_summary.md),
[`get_robust_colocalization()`](https://statfungen.github.io/colocboost/reference/get_robust_colocalization.md)

## Examples

``` r
data(Ambiguous_Colocalization)
test_colocboost_results <- Ambiguous_Colocalization$ColocBoost_Results
res <- get_ambiguous_colocalization(test_colocboost_results)
#> There exists the ambiguous colocalization events from trait-specific effects. Extracting!
#> There are 1 ambiguous trait-specific effects.
names(res$ambiguous_cos)
#> [1] "ucos1:y1;ucos2:y2"
```
