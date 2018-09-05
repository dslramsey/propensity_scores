Using propensity scores for causal inference in ecology: options, considerations and a case study
================

Notes
-----

This repository contains data and code from:

Ramsey, D.S.L., Forsyth, D.M., Wright, E., McKay, M., and Westbrooke, I. (2018). "Using propensity scores for causal inference in ecology: options, considerations and a case study" *Methods in Ecology and Evolution*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1403922.svg)](https://doi.org/10.5281/zenodo.1403922)

Getting started
---------------

**File descriptions**:

*Propensity\_simulations.R* – R code used to conduct Monte Carlo simulations of propensity score methods. Uses functions defined in *propensity\_simulation\_functions.R*

*Propensity\_simulation\_functions.R* – R functions used to perform propensity score simulations and helper functions to calculate IPTW weights calc.pswts() and weights from full matching get.match.weights() that target either the ATE or ATT.

*Propensity\_TreeCover\_Analysis.R* – R code used to estimate effects of possum control on canopy tree condition using propensity scores. Uses data provided in *Tree\_Cover\_data.csv*.

Prerequisites
-------------

The simulation script require packages *optmatch*, *tidyverse*, *survey*, *forcats* and *RColorBrewer* while the tree cover analysis requires in addition, *cobalt* , *readr* and *treatsens*.
