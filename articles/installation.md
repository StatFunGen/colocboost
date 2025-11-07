# Installation

## Installation

### CRAN (Stable Release)

Install major and stable releases from CRAN (Linux, macOS and Windows)

``` r

install.packages("colocboost")
```

### GitHub

Install the development version from GitHub

``` r

devtools::install_github("StatFunGen/colocboost")
```

### Conda

Install minor releases using pre-built conda package with a
conda-compatible package manager (recommended)

Global pixi installation is the easiest way to use the conda package

``` bash
pixi global install -c dnachun -c conda-forge r-base # Install r-base as a global package if not already installed
pixi global install -c dnachun -c conda-forge --environment r-base r-colocboost # Inject r-colocboost into r-base global environment
```

The package can also be added to a local pixi environment

``` bash
pixi workspace channel add dnachun # Add the dnachun channel to the workspace
pixi add r-colocboost # Add r-colocboost as a dependency to the environment
```

Micromamba is recommended instead of conda or mamba for traditional
conda environments

``` bash
micromamba install -c dnachun -c conda-forge r-colocboost
mamba install -c dnachun -c conda-forge r-colocboost
conda install -c dnachun -c conda-forge r-colocboost
```
