# R package - Indirect methods for estimating age-specific death rates
-----------------------------------
[![Build Status](https://travis-ci.org/mpascariu/MortalityEstimate.svg?branch=master)](https://travis-ci.org/mpascariu/MortalityEstimate)
[![issues](https://img.shields.io/github/issues-raw/mpascariu/MortalityEstimate.svg)]()
[![license](https://img.shields.io/github/license/mpascariu/MortalityEstimate.svg)]()

This repository includes R code for estimating mortality curves and full life tables
based on given one or two pieces of information: life expectancy, child mortality, or child and adult mortality. The implemented models are better suited to the practical needs of mortality estimation, since the input parameters are continuous yet the second one is optional.


Installation
============

1. Make sure you have the most recent version of R
2. Run the following code in your R console 

```r
# install.packages("devtools")

library(devtools)
install_github("mpascariu/MortalityEstimate")
```

Help
===============
All functions are documented in the standard way, which means that 
once you load the package using ```library(MortalityEstimate)```
you can just type ```?LinearLink``` to see the help file. 

Check the examples provided in the `Kannisto` 
`lifetable` and `LinearLink` functions.

