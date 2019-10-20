
# mtree <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/caravagn/mtree.svg?branch=master)](https://travis-ci.org/caravagn/mtree)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

Thee `mtree` package provides a simple implementation of mutation trees
that can be build of binary data reporting the presence or absence of
specific mutations in multiple tumour regions. This type of models can
be used to study the evolutionary trajectories of a tumour from bulk
sequencing data generated via targeted sequencing.

`mtree` is part of the `evoverse`, a package that gathers multiple R
packages to implement Cancer Evolution analyses; see more [about
evoverse](https://caravagn.github.io/evoverse).

The package provides methods to create the trees as S3 objects, as well
as a sampler to generate them with a Monte Carlo strategy that have been
first described in the [revolver](https://caravagn.github.io/revolver)
algorithm. The package provides also functions to plot and analyze the
trees.

#### Help and support

`mtree` has its own webpage at [GitHub
pages](https://caravagn.github.io/mtree/).

-----

### Installation

You can install the released version of `mtree` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagn/mtree")
```

-----

#### Copyright and contacts

Giulio Caravagna, PhD. *Institute of Cancer Research, London, UK*.

  - Personal webpage:
    [https://bit.ly/2kc9E6Y](https://sites.google.com/site/giuliocaravagna/),
  - Email address: <giulio.caravagna@icr.ac.uk> and
    <gcaravagn@gmail.com>
  - Twitter feed: [@gcaravagna](https://twitter.com/gcaravagna)
  - GitHub space: [caravagn](https://github.com/caravagn)
