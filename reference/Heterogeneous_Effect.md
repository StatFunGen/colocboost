# Individual level data for 2 traits and 2 causal variants with heterogeneous effects

An example dataset with simulated genotypes and traits for 2 traits and
2 common causal variants with heterogeneous effects

## Usage

``` r
Heterogeneous_Effect
```

## Format

### `Heterogeneous_Effect`

A list with 3 elements

- X:

  List of genotype matrices

- Y:

  List of traits

- variant:

  indices of two causal variants

## Source

The Heterogeneous_Effect dataset contains 2 simulated phenotypes
alongside corresponding genotype matrices. There are two causal
variants, both of which have heterogeneous effects on two traits. Due to
the file size limitation of CRAN release, this is a subset of simulated
data to generate Figure 2b in Cao etc. 2025. See full dataset in
colocboost paper repo <https://github.com/StatFunGen/colocboost-paper>.

## See also

Other colocboost_data:
[`Ambiguous_Colocalization`](https://statfungen.github.io/colocboost/reference/Ambiguous_Colocalization.md),
[`Ind_5traits`](https://statfungen.github.io/colocboost/reference/Ind_5traits.md),
[`Non_Causal_Strongest_Marginal`](https://statfungen.github.io/colocboost/reference/Non_Causal_Strongest_Marginal.md),
[`Sumstat_5traits`](https://statfungen.github.io/colocboost/reference/Sumstat_5traits.md),
[`Weaker_GWAS_Effect`](https://statfungen.github.io/colocboost/reference/Weaker_GWAS_Effect.md)
