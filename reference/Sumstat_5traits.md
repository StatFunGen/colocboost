# Summary level data for 5 traits

An example dataset with simulated statistics for 5 traits

## Usage

``` r
Sumstat_5traits
```

## Format

### `Sumstat_5traits`

A list with 2 elements

- sumstat:

  Summary statistics for 5 traits

- true_effect_variants:

  List of causal variants

## Source

The Sumstat_5traits dataset contains 5 simulated summary statistics,
where it is directly derived from the Ind_5traits dataset using marginal
association. The dataset is specifically designed for evaluating and
demonstrating the capabilities of ColocBoost in multi-trait
colocalization analysis with summary association data. See Cao etc. 2025
for details. Due to the file size limitation of CRAN release, this is a
subset of simulated data. See full dataset in colocboost paper repo
<https://github.com/StatFunGen/colocboost-paper>.

## See also

Other colocboost_data:
[`Ambiguous_Colocalization`](https://statfungen.github.io/colocboost/reference/Ambiguous_Colocalization.md),
[`Heterogeneous_Effect`](https://statfungen.github.io/colocboost/reference/Heterogeneous_Effect.md),
[`Ind_5traits`](https://statfungen.github.io/colocboost/reference/Ind_5traits.md),
[`Non_Causal_Strongest_Marginal`](https://statfungen.github.io/colocboost/reference/Non_Causal_Strongest_Marginal.md),
[`Weaker_GWAS_Effect`](https://statfungen.github.io/colocboost/reference/Weaker_GWAS_Effect.md)
