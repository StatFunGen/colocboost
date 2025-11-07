# Individual level data for 5 traits

An example dataset with simulated genotypes and traits for 5 traits

## Usage

``` r
Ind_5traits
```

## Format

### `Ind_5traits`

A list with 3 elements

- X:

  List of genotype matrices

- Y:

  List of traits

- true_effect_variants:

  List of causal variants

## Source

The Ind_5traits dataset contains 5 simulated phenotypes alongside
corresponding genotype matrices. The dataset is specifically designed
for evaluating and demonstrating the capabilities of ColocBoost in
multi-trait colocalization analysis with individual-level data. See Cao
etc. 2025 for details. Due to the file size limitation of CRAN release,
this is a subset of simulated data. See full dataset in colocboost paper
repo <https://github.com/StatFunGen/colocboost-paper>.

## See also

Other colocboost_data:
[`Ambiguous_Colocalization`](https://statfungen.github.io/colocboost/reference/Ambiguous_Colocalization.md),
[`Heterogeneous_Effect`](https://statfungen.github.io/colocboost/reference/Heterogeneous_Effect.md),
[`Non_Causal_Strongest_Marginal`](https://statfungen.github.io/colocboost/reference/Non_Causal_Strongest_Marginal.md),
[`Sumstat_5traits`](https://statfungen.github.io/colocboost/reference/Sumstat_5traits.md),
[`Weaker_GWAS_Effect`](https://statfungen.github.io/colocboost/reference/Weaker_GWAS_Effect.md)
