# Input Data Format

``` r

library(colocboost)
```

This vignette documents the standard input data formats of `colocboost`.

## 1. Individual Level Data

For analyses using individual-level data, the basic format for single
trait is as follows:

- `X` is an $`N \times P`$ matrix with $`N`$ individuals and $`P`$
  variants. Including variant names as column names is highly
  recommended, especially when working with multiple $`X`$ matrices and
  $`Y`$ vectors.
- `Y` is a length $`N`$ vector containing phenotype values for the same
  $`N`$ individuals as $`X`$.

The input format for multiple traits is similar, but `X` should be a
list of genotype matrices, each corresponding to a different trait. `Y`
should also be a list of phenotype vectors. For example:

- `X = list(X1, X2, X3, X4, X5)` where each `Xi` is a matrix for trait
  `i` - with the dimension of $`N_i \times P_i`$, where $`N_i`$ and
  $`P_i`$ do not need to be the same for different traits.
- `Y = list(Y1, Y2, Y3, Y4, Y5)` where each `Yi` is a vector for trait
  `i` - with $`N_i`$ individuals.

`colocboost` also offers flexible input options (see detailed usage with
different input formats, refer to [Individual Level Data
Colocalization](https://statfungen.github.io/colocboost/articles/Individual_Level_Colocalization.html)):

- Single $`X`$ matrix with $`N \times P`$, and $`Y`$ matrix with
  $`N \times L`$ for $`L`$ traits.
- Multiple $`X`$ matrices and unmatched $`Y`$ vectors with a mapping
  dictionary (example shown in section 3 below).

## 2. Summary Statistics

For analyses using summary statistics, the basic format for single trait
is as follows:

- `sumstat` is a data frame with required columns `z` or (`beta`,
  `sebeta`), and optional columns but highly recommended `n` and
  `variant`.

``` r

data(Sumstat_5traits)
head(Sumstat_5traits$sumstat[[1]])
#>              z    n variant
#> 451 -1.0945531 1153    rs_1
#> 452 -0.4113347 1153    rs_2
#> 453 -0.4113347 1153    rs_3
#> 454 -0.7467923 1153    rs_4
#> 455 -0.3018575 1153    rs_5
#> 456 -0.5256479 1153    rs_6
```

    - `z` or (`beta`, `sebeta`) - required: either z-score or (effect size and standard error)
    - `n` - highly recommended: sample size for the summary statistics, it is highly recommendation to provide.
    - `variant` - highly recommended: required if sumstat for different outcomes do not have the same number of variables (multiple sumstat and multiple LD).

- `LD` is a matrix of LD. This matrix does not need to contain the exact
  same variants as in `sumstat`, but the `colnames` and `rownames` of
  `LD` should include the `variant` names for proper alignment.

The input format for multiple traits is similar, but `sumstat` should be
a list of data frames `sumstat = list(sumstat1, sumstat2, sumstat3)`.
The flexibility of input format for multiple traits is as follows (see
detailed usage with different input formats, refer to [Summary
Statistics
Colocalization](https://statfungen.github.io/colocboost/articles/Summary_Level_Colocalization.html)):

- One LD matrix with a superset of variants in `sumstat` for all traits
  is allowed.
- Multiple LD matrices, each corresponding to a different trait, are
  also allowed for the trait-specific LD structure.
- Multiple LD matrices and unmatched `sumstat` data frames with a
  mapping dictionary are also allowed (example shown in section 3
  below).

## 3. Optional: mapping between arbitrary input $`X`$ and $`Y`$

For analysis when including multiple genotype matrices `X` with
unmatched arbitrary phenotype vectors `Y`, a mapping dictionary
`dict_YX` is required to indicate the relationship between `X` and `Y`.
Similarly, when multiple LD matrices with unmatched arbitrary multiple
summary statistics `sumstat` are used, a mapping dictionary
`dict_sumstatLD` is required to indicate the relationship between
`sumstat` and `LD`.

For example, considering three genotype matrices `X = list(X1, X2, X3)`
and 6 phenotype vectors `Y = list(Y1, Y2, Y3, Y4, Y5, Y6)`, where

- `X1` is for trait 1, trait 2, trait 3
- `X2` is for trait 4, trait 5
- `X3` is for trait 6

Then, you need to define a 6 by 2 matrix mapping dictionary `dict_YX` as
follows:

- The first column should be `c(1,2,3,4,5,6)` for 6 traits.
- The second column should be `c(1,1,1,2,2,3)` for 3 genotype matrices.

Here, each row indicates the trait index and the corresponding genotype
matrix index.

``` r

dict_YX <- cbind(c(1,2,3,4,5,6), c(1,1,1,2,2,3))
dict_YX
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    2    1
#> [3,]    3    1
#> [4,]    4    2
#> [5,]    5    2
#> [6,]    6    3
```

## 4. HyPrColoc compatible format: effect size and standard error matrices

ColocBoost also provides a flexibility to use HyPrColoc compatible
format for summary statistics with and without LD matrix. For example,
when analyze $`L`$ traits for the same $`P`$ variants with the specified
effect size and standard error matrices:

- `effect_est` (required) is $`P \times L`$ matrix of variable
  regression coefficients (i.e. regression beta values) in the genomic
  region.
- `effect_se` (required) is $`P \times L`$ matrix of standard errors for
  the regression coefficients.
- `effect_n` (highly recommended) is either a scalar or a vector of
  sample sizes for estimating regression coefficients.
- `LD` (optional) is LD matrix for the $`P`$ variants. If it is not
  provided, it will apply LD-free ColocBoost.

See more details about HyPrColoc compatible format in [Summary
Statistics
Colocalization](https://statfungen.github.io/colocboost/articles/Summary_Level_Colocalization.html)).

See more details about data format to implement LD-free ColocBoost and
LD-mismatch diagnosis in [LD mismatch and LD-free
Colocalization](https://statfungen.github.io/colocboost/articles/LD_Free_Colocalization.html)).
