# Package index

## Example Data

Example datasets for demonstration and testing.

- [`Ambiguous_Colocalization`](https://statfungen.github.io/colocboost/reference/Ambiguous_Colocalization.md)
  : A real data example includes an ambiguous colocalization between
  eQTL and GWAS
- [`Heterogeneous_Effect`](https://statfungen.github.io/colocboost/reference/Heterogeneous_Effect.md)
  : Individual level data for 2 traits and 2 causal variants with
  heterogeneous effects
- [`Ind_5traits`](https://statfungen.github.io/colocboost/reference/Ind_5traits.md)
  : Individual level data for 5 traits
- [`Non_Causal_Strongest_Marginal`](https://statfungen.github.io/colocboost/reference/Non_Causal_Strongest_Marginal.md)
  : Individual level data for 2 traits and 2 causal variants, but the
  strongest marginal association is not causal
- [`Sumstat_5traits`](https://statfungen.github.io/colocboost/reference/Sumstat_5traits.md)
  : Summary level data for 5 traits
- [`Weaker_GWAS_Effect`](https://statfungen.github.io/colocboost/reference/Weaker_GWAS_Effect.md)
  : Individual level data for 2 traits and 2 causal variants with weaker
  effects for focal trait

## Model fitting

Main interface function for fitting multi-trait colocalization model.

- [`colocboost()`](https://statfungen.github.io/colocboost/reference/colocboost.md)
  : ColocBoost: A gradient boosting informed multi-omics xQTL
  colocalization method

## Inference and summary

Functions for inference and summary from fitted model.

- [`get_robust_colocalization()`](https://statfungen.github.io/colocboost/reference/get_robust_colocalization.md)
  : Recalibrate and summarize robust colocalization events.
- [`get_colocboost_summary()`](https://statfungen.github.io/colocboost/reference/get_colocboost_summary.md)
  : Get summary tables from a ColocBoost output.
- [`get_ambiguous_colocalization()`](https://statfungen.github.io/colocboost/reference/get_ambiguous_colocalization.md)
  : Get ambiguous colocalization events from trait-specific
  (uncolocalized) effects.
- [`get_robust_ucos()`](https://statfungen.github.io/colocboost/reference/get_robust_ucos.md)
  : Recalibrate and summarize robust uncolocalized events.

## Visualization

Functions for visualizing ColocBoost result.

- [`colocboost_plot()`](https://statfungen.github.io/colocboost/reference/colocboost_plot.md)
  : Plot visualization plot from a ColocBoost output.

## Utilities

Helper functions and utilities

- [`get_cormat()`](https://statfungen.github.io/colocboost/reference/get_cormat.md)
  : A fast function to calculate correlation matrix (LD matrix) from
  individual level data
- [`get_cos()`](https://statfungen.github.io/colocboost/reference/get_cos.md)
  : Extract CoS at different coverage
- [`get_cos_purity()`](https://statfungen.github.io/colocboost/reference/get_cos_purity.md)
  : Calculate purity within and in-between CoS
- [`get_cos_summary()`](https://statfungen.github.io/colocboost/reference/get_cos_summary.md)
  : Get colocalization summary table from a ColocBoost output.
- [`get_hierarchical_clusters()`](https://statfungen.github.io/colocboost/reference/get_hierarchical_clusters.md)
  : Perform modularity-based hierarchical clustering for a correlation
  matrix
- [`get_ucos_summary()`](https://statfungen.github.io/colocboost/reference/get_ucos_summary.md)
  : Get trait-specific summary table from a ColocBoost output.
