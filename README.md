# The Linear-Link model
-----------------------------------

This repository includes R code for estimating mortality curves and life tables
based on a value of life expectancy at birth. An initial version the model was 
presented at The Twelfth International Longevity Risk and 
Capital Markets Solutions Conference in Chicago (September 30, 2016).
`The Linear Link: Deriving Age-Specific Death Rates from Life Expectancy`
by [Marius Pascariu](http://findresearcher.sdu.dk:8080/portal/da/person/mpascariu) 
and [Vladimir Canudas-Romo](http://www.sdu.dk/ansat/vcanudas). Since then 
[Jose Manuel Aburto](https://twitter.com/jm_aburto) and 
[Ugofilippo Basellini](https://twitter.com/ugobas) joined the project and
made important contributions.

Installation
============

1. Make sure you have the most recent version of R
2. Run the following code in your R console 

```r
# install.packages("devtools")

library(devtools)
install_github("mpascariu/LinearLink")
```

Help
===============
All functions are documented in the standard way, which means that 
once you load the package using ```library(LinearLink)```
you can just type ```?LinearLink``` to see the help file. 

Check the examples provided in the `Kannisto` 
`lifetable` and `LinearLink` functions.

For now the model was tested to work for life expectancy at birth. For the other
ages further changes need to be implemented.

Abstract
========

Predicting the human longevity level in the future by directly 
forecasting life expectancy offers numerous advantages compared 
with methods based on extrapolation of age-specific death rates. 
However, the reconstruction of accurate life tables starting from 
a given level of life expectancy at birth or any other age is not 
straightforward. Model life tables were extensively used in the past 
for estimating age patterns of mortality in data-poor countries.
We propose a new model inspired by indirect estimation techniques used 
in demography that can be used to estimate full life tables 
given a predicted life expectancy.