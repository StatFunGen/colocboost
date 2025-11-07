# Individual level data for 2 traits and 2 causal variants, but the strongest marginal association is not causal

An example dataset with simulated genotypes and traits for 2 traits and
2 common causal variants, but the strongest marginal association is not
causal variant.

## Usage

``` r
Non_Causal_Strongest_Marginal
```

## Format

### `Non_Causal_Strongest_Marginal`

A list with 3 elements

- X:

  List of genotype matrices

- Y:

  List of traits

- variant:

  indices of two causal variants

## Source

The Non_Causal_Strongest_Marginal dataset contains 2 simulated
phenotypes alongside corresponding genotype matrices. There are two
causal variants, but the strongest marginal association is not a causal
variant. Due to the file size limitation of CRAN release, this is a
subset of simulated data to generate Figure 2b in Cao etc. 2025. See
full dataset in colocboost paper repo
<https://github.com/StatFunGen/colocboost-paper>.

## See also

Other colocboost_data:
[`Ambiguous_Colocalization`](https://statfungen.github.io/colocboost/reference/Ambiguous_Colocalization.md),
[`Heterogeneous_Effect`](https://statfungen.github.io/colocboost/reference/Heterogeneous_Effect.md),
[`Ind_5traits`](https://statfungen.github.io/colocboost/reference/Ind_5traits.md),
[`Sumstat_5traits`](https://statfungen.github.io/colocboost/reference/Sumstat_5traits.md),
[`Weaker_GWAS_Effect`](https://statfungen.github.io/colocboost/reference/Weaker_GWAS_Effect.md)
