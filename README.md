
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.888409.svg)](https://doi.org/10.5281/zenodo.888409) [![Build Status](https://travis-ci.org/mrecos/DistRegLMERR.svg?branch=master)](https://travis-ci.org/mrecos/DistRegLMERR)

PRE-RELEASE
===========

DistRegLMERR - Distribution Regression with Kernel Logistic Regression on Mean Embeddings.
------------------------------------------------------------------------------------------

### SAA 2018 Abstract

A model of Distribution Regression using Kernel Logistic Regression (KLR) on Mean Feature Space Embeddings is developed to address two primary shortcomings of current Archaeological Predictive Modeling (APM) practice; 1) neglecting the richness of archaeological landforms by collapsing a site to a single point or observation; and 2) disregarding the implicit and explicit uncertainty of archaeological data, predictions, and model parameters. This research addresses the first hurdle by developing a KLR on embeddings approach to Distribution Regression. This method first samples a distribution of variable measurements from the spatial area of each site, then uses a kernel to project the distributions into a non-geographical feature space to calculate mean embeddings, finally Kernel Logistic Regression estimates similarity coefficients for inference and prediction. The primary benefits of this approach to APM are the consideration of archaeological landform richness and variation, explicitly modeling similarity between sites and the environment, and allowing for similarity metrics specific to archaeological research questions. The second hurdle is addressed by applying the MERR method within a Bayesian framework for probabilistic modeling. As such, the uncertainty of data and parameters can be explicitly modeled with priors resulting in a posterior predictive distribution useful for quantifying and visualizing risk.

-   Note: Changes from original abstract; 1) change of terms to Kernel Logistic Regression (KLR); and 2) deemphasized Bayesian approach as it is a work in progress. The repo contains working Stan code for probabilistic models, but they are not ready for prime-time.

![](https://github.com/mrecos/DistRegLMERR/blob/master/analysis/images/KLR_map.jpg?raw=true)

### Package

This packaage contains the functions neccesary to compute Kernel Linear Regression (KLR) on mean kernel embeddings, functions for preparing site and background data, and a function for simulate archaeological site location data.

-   build\_K - Function takes list of traning data, scalar value for `sigma` hyperparameter, and a distance method to compute a mean embedding similarity kernel. This kernel is a pair-wise (N x N) matrix of the mean similarity between the attirbutes describing each site location and background group.
-   KLR - Function takes the similarity kenrel matrix `K`, a vector of presence/absence coded as `1` ro `0`, and a scalar value for the `lambda` regularizing hyperparameter; optionally values for maximum iterations and threshold. This function performs Kernel Logistic Regression (KLR) via Iterative Reweighted Least Squares (IRLS). The objective is to approximate a set of parameters that minimize the negative likelihood of the parameters given the data and response. This function returns a list of `pred`, the estimated response (probability of site-presence) for the training data, and `alphas`, the approximated parameters from the IRLS algortithm.
-   KLR\_predict - Function takes a list of the training data, a list of the testing data, a vector of the approximated `alphas` parameters, a scalar value for the `sigma` kernel hyperparameter, and a distance method. This function predicts the probability of site presence for new oberservations based on the training data and `alphas` parameters. This is accomplished by building the `k*k` kernel matrix as the similarity between the training test data then computing the inverse logit of `k*k %*% alphas`. The output is the predicted probability of site presence for each traning data example.

### Citation

Please cite this compendium as:

> Harris, Matthew D., (2017). *Distribution Regression with Kernel Logistic Regression on Mean Embeddings*. Accessed 10 Sep 2017. Online at <https://doi.org/10.5281/zenodo.888409>

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
