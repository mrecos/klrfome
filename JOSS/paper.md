---
title: "KLRfome - Kernel Logistic Regression on Focal Mean Embeddings"
authors:
- affiliation: 1
  name: Matthew D. Harris
  orcid: 0000-0002-4627-7692
date: "14 April 2018"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- Archaeology
- Archaeological Science
- Kernel Methods
- Machine Learning
affiliations:
- index: 1
  name: AECOM Technologies
---

# Summary

- KLRfome is a package for predicting a geographic area's realtive sensitivity for the presence of archaeological sites. This is achieved by fitting, predicting, and visualizing a Kernel Logistic Regression model on mean feature embeddings. This regression algorithm and package are created to improve upon the current methods in archaeological predictive modeling. These improvements include 1) modeling rich descriptions of archaeological landforms by mitigating undesiarable spatial corraltion between samples, 2) explicitly modeling the similarity between archaeological sites as characterized by rich features, 3) the ability to define research specific similairty measures, and 4) focal window prediction that can be modified based on theory or managment goals. This work is based on Szabó's work on feature mean embeddings for distribution regression [@Szabo] and the Kernel Logistic Regression algortihm in [@Zhu]. The KLRfome package is available at [https://github.com/mrecos/klrfome](https://github.com/mrecos/klrfome) and DOI [https://doi.org/10.5281/zenodo.1218403](https://doi.org/10.5281/zenodo.1218403)


![Data ingest](https://github.com/mrecos/klrfome/blob/master/README_images/KLRfome_dataflow.png?raw=true)
![Focal prediction](https://github.com/mrecos/klrfome/blob/master/README_images/KLRfome_prediction.png?raw=true)

# References