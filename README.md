# ColocBoost for multi-trait colocalization in molecular QTL and GWAS studies
[![Codecov test coverage](https://codecov.io/gh/StatFunGen/colocboost/branch/main/graph/badge.svg)](https://codecov.io/gh/StatFunGen/colocboost?branch=main)
[![CRAN Version](https://www.r-pkg.org/badges/version/colocboost)](https://cran.r-project.org/package=colocboost)

This R package implements ColocBoost --- motivated and designed for colocalization analysis ([first formulated here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)) of multiple genetic association studies --- as a multi-task learning approach to variable selection regression with highly correlated predictors and sparse effects, based on frequentist statistical inference. It provides statistical evidence to identify which subsets of predictors have non-zero effects on which subsets of response variables.

## Installation

### Conda
Install major releases from conda (recommended)

```bash
conda install -c dnachun r-colocboost
```

### CRAN
Install released versions from cran

```r
install.packages("colocboost")
```

### GitHub
Install the development version from GitHub

```r
devtools::install_github("StatFunGen/colocboost")
```

## Usage

### Single-trait Fine-mapping (FineBoost)
Run FineBoost for single-trait fine-mapping (similar interface to SuSiE)
```r
result <- colocboost(X=X, Y=y)
```

### Multi-trait Colocalization
```r
# Basic multi-trait analysis
result <- colocboost(X=list(X), Y=list(y1, y2, y3))

# Using summary statistics
result <- colocboost(sumstat=list(sumstat1, sumstat2), LD=LD_matrix)

# View colocalization summary
summary <- get_cos_summary(result)

# Visualize results
colocboost_plot(result)

# Filter for stronger colocalization evidence
filtered <- get_strong_colocalization(result, cos_npc_cutoff = 0.5)
```

For more complex analyses involving multiple datasets mixing individual level and summary statistics data, we recommend using [this pipeline wrapper](https://github.com/StatFunGen/pecotmr/blob/main/R/colocboost_pipeline.R) from the `pecotmr` package. The `pecotmr` package can be installed either from source or from our conda package at https://anaconda.org/dnachun/r-pecotmr.

## Citation

If you use ColocBoost in your research, please cite:

Cao X, Sun H, Feng R, Mazumder R, Najar CFB, Li YI, de Jager PL, Bennett D, The Alzheimer's Disease Functional Genomics Consortium, Dey KK, Wang G. (2025+). Integrative multi-omics QTL colocalization maps regulatory architecture in aging human brain. bioRxiv. [https://doi.org/](https://doi.org/)

## Documentation

For detailed documentation, use the R help system:

```r
?colocboost
?colocboost_plot
?get_cos_summary
?get_strong_colocalization
```

## License

This package is released under the MIT License.
