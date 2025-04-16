# ColocBoost for multi-trait colocalization in molecular QTL and GWAS studies
[![Codecov test coverage](https://codecov.io/gh/StatFunGen/colocboost/branch/main/graph/badge.svg)](https://codecov.io/gh/StatFunGen/colocboost?branch=main)
[![CRAN Version](https://www.r-pkg.org/badges/version/colocboost)](https://cran.r-project.org/package=colocboost)

This R package implements ColocBoost --- motivated and designed for colocalization analysis ([first formulated here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)) of multiple genetic association studies --- as a multi-task learning approach to variable selection regression with highly correlated predictors and sparse effects, based on frequentist statistical inference. It provides statistical evidence to identify which subsets of predictors have non-zero effects on which subsets of response variables.

## Installation

### CRAN
Install released versions from CRAN - pre-built packages are available on macOS and Windows

```r
install.packages("colocboost")
```

### GitHub
Install the development version from GitHub

```r
devtools::install_github("StatFunGen/colocboost")
```

### Conda
Install major releases using pre-built conda package with a conda-compatible package manager (recommended)

Global pixi installation is the easiest way to use the conda package
```bash
pixi global install r-base # Install r-base as a global package if not already installed
pixi global install --environment r-base r-colocboost # Inject r-colocboost into r-base global environment
```
The package can also be added to a local pixi environment
```bash
pixi workspace channel add dnachun # Add the dnachun channel to the workspace
pixi add r-colocboost # Add r-colocboost as a dependency to the environment
```
Micromamba is recommended instead of conda or mamba for traditional conda environments
```bash
micromamba install -c dnachun r-colocboost
mamba install -c dnachun r-colocboost
conda install -c dnachun r-colocboost
```
## Usage

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

### Single-trait Fine-mapping (FineBoost) - Special Case
Run FineBoost for single-trait fine-mapping (similar interface to SuSiE)
```r
result <- colocboost(X=X, Y=y)
```

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
