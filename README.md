
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.888409.svg)](https://doi.org/10.5281/zenodo.888409)

<!-- README.md is generated from README.Rmd. Please edit that file -->
DistRegLMERR - Logisitic Mean Embeddings Ridge Regression (LMERR) for the Distribution Regression problem.
============

The goal of DistRegLMERR is to ...

### SAA 2018 Abstract
A Bayesian model of Distribution Regression using a Mean Embedding Ridge Regression (MERR) algorithm is developed to address two primary shortcomings of current Archaeological Predictive Modeling (APM) practice; 1) neglecting the richness of archaeological landforms by collapsing a site to a single point or observation; and 2) disregarding the implicit and explicit uncertainty of archaeological data, predictions, and model parameters. This research addresses the first hurdle by developing a Logistic MERR approach to Distribution Regression. This method first samples a distribution of variable measurements from the spatial area of each site, then uses a kernel to project the distributions into a non-geographical feature space to calculate mean embeddings, finally Kernel Ridge Regression estimates similarity coefficients for inference and prediction. The primary benefits of the MERR approach to APM are the consideration of archaeological landform richness and variation, explicitly modeling similarity between sites and the environment, and allowing for similarity metrics specific to archaeological research questions. The second hurdle is addressed by applying the MERR method within a Bayesian framework for probabilistic modeling. As such, the uncertainty of data and parameters can be explicitly modeled with priors resulting in a posterior predictive distribution useful for quantifying and visualizing risk.

### Citation

Please cite this compendium as:

> Harris, Matthew D., (2017). *Logisitic Mean Embeddings Ridge Regression (LMERR) for the Distribution Regression problem.*. Accessed 10 Sep 2017. Online at <https://doi.org/10.5281/zenodo.888409>

### Installation

You can install DistRegLMERR from github with:

``` r
# install.packages("devtools")
devtools::install_github("mrecos/DistRegLMERR")
```

### Licenses

**Text and figures:** [CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code:** See the [DESCRIPTION](DESCRIPTION) file

**Data:** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/) attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please see our [contributor guidelines](CONTRIBUTING.md). Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
