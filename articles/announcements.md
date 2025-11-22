# News

## ColocBoost on the media

- *April 20, 2025*: Manuscript describing ColocBoost methods with
  applications is posted on
  [medRxiv](https://www.medrxiv.org/content/10.1101/2025.04.17.25326042v1).
- *May 2, 2025*: `colocboost` R package is available on
  [CRAN](https://CRAN.R-project.org/package=colocboost).

## Software updates

- `v1.0.7` (**Important update**) Improvements to ColocBoost (check out
  the full details in
  [PR](https://github.com/StatFunGen/colocboost/pull/116) and
  [PR](https://github.com/StatFunGen/colocboost/pull/121)).
  - Enhanced `colocboost` main function with post-filtering and only
    keep the robust colocalization events.
  - Enhanced `colocboost_plot` function with flexible highlighting
    options and new visualization styles.
  - Optimized performance and computational efficiency
  - Improved documentation and examples for the wrapper pipeline
  - Minor bug fixes for increased stability
- `v1.0.6` Memory optimization and visualization improvements with bug
  fixes [CRAN](https://CRAN.R-project.org/package=colocboost).
  - Optimized LD-free version to reduce memory usage by eliminating
    large identity LD matrix generation
  - Enhanced `colocboost_plot` function with improved horizontal and
    vertical spacing labels
  - Resolved minor bugs for improved stability For detailed changelog,
    see [PR](https://github.com/StatFunGen/colocboost/pull/102).
- `v1.0.4` Initial public stable release on
  [CRAN](https://CRAN.R-project.org/package=colocboost). Due to the file
  size limitation of CRAN release, all example data are the subset of
  the simulated data. See full dataset in [colocboost paper
  repo](https://github.com/StatFunGen/colocboost-paper).
- `v1.0.0` Initial public release. Pre-release versions `0.1.*` for
  reproducing our earlier work, including the first version of
  manuscript on medRxiv, are still available on
  [GitHub](https://github.com/StatFunGen/colocboost/tags).
