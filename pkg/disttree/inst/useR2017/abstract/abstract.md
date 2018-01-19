---
title: "Distributional Trees and Forests"
author: |
  | Lisa Schlosser^1^, Torsten Hothorn^2^, Achim Zeileis^1^
  |
  | 1. University of Innsbruck
  | 2. University of Zurich
institute:
- $^1$Department of Statistics
- $^2$Epidemiology, Biostatistics and Prevention Institute
output:
  html_document: default
  pdf_document: default
  word_document: default
nocite: null
bibliography: ref.bib
---

**Keywords**: Distributional regression, recursive partitioning, decision trees, random forests

**Webpages**: https://R-Forge.R-project.org/projects/partykit/

In regression analysis one is interested in the relationship between a dependent variable and one or more explanatory variables. Various methods to fit statistical models to the data set have been developed, starting from ordinary linear models considering only the mean of the response variable and ranging to probabilistic models where all parameters of a distribution are fit to the given data set. 

If there is a strong variation within the data it might be advantageous to split the data first into more homogeneous subgroups based on given covariates and then fit a local model in each subgroup rather than fitting one global model to the whole data set. This can be done by applying regression trees and forests.

Both of these two concepts, parametric modeling and algorithmic trees, have been investigated and developed further, however, mostly separated from each other. Therefore, our goal is to embed the progress made in the field of probabilistic modeling in the idea of algorithmic tree and forest models. In particular, more flexible models such as GAMLSS [@stasinopoulos2005] should be fitted in the nodes of a tree in order to capture location, scale, shape as well as censoring, tail behavior etc. while non-additive effects of the explanatory variables can be detected by the splitting algorithm used to build the tree.

The corresponding implementation is provided in an *R* package **disttree** which is available on R-Forge and includes the two main functions `disttree` and `distforest`. Next to the data set and a formula the user only has to specify a distribution family and receives a tree/forest model with a set of distribution parameters for each final node. One possible way to specify a distribution family is to hand over a **gamlss.dist** family object [@stasinopoulos2007generalized].
In `disttree` and `distforest` the fitting function `distfit` is applied within a tree building algorithm chosen by the user. Either the MOB algorithm, an algorithm for model-based recursive partitioning [@zeileis2008model], or the ctree algorithm [@hothorn2006unbiased] can be used as a framework. These algorithms are both implemented in the **partykit** package [@hothorn2015package].

# References
