# EXPRESSO
EXpression PREdiction with Summary Statistics Only 

## Table of contents
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Quick tutorial](#Quick_tutorial)
* [Input files](#Input_files)
* [Contact](#Contact)

## Introduction
EXPRESSO (EXpression PREdiction with Summary Statistics Only) could bulid gene expression model with eQTL summary statistics and reference panel only. It also integrates 3D genomic data to define cis-regulatory regions properly and uses epigenetic annotation to prioritize causal variants. It is developed and maintained by Lida Wang at [Dajiang Liu's Group](https://dajiangliu.blog).

## Installation
The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the MASS, data.table, BEDMatrix and caret packages installed.

```
install.packages("devtools")
library(devtools)
```
And also, you need the latest version of [fast.lasso](https://github.com/zhanxw/fast.lasso) and [rareGWAMA](https://github.com/dajiangliu/rareGWAMA) to be installed.

```
devtools::install_github("zhanxw/fast.lasso")
devtools::install_github("dajiangliu/rareGWAMA")
```
Then you could install EXPRESSO from the repository here.

```
devtools::install_github("LidaWangPSU/EXPRESSO/EXPRESSO")
library(EXPRESSO)
```
Here we go.

## Quick tutorial
### Bulid gene expression prediction model
```
res.tmp <- EXPRESSO(sumstatFile,annoFile,windowFile,refFile,out_path,minMaf,maxIter,gene.vec,append=F)
```
Input includes
* sumstatFile: summary statistics file
* annoFile: epigenetic annotation file 
* windowFile: 3D genomic windows file
* refFile: reference panel file 
* out_path: pre-specified output path
* minMaf: filter by minimum allele frequency
* maxIter: maximum iteration of gradient descent algorithm
* gene.vec: a list of input gene id

### Output results
We perform EXPRESSO by three different tunning parameter methods, including pseudo variable selection, BIC and MSE.
The weight output includes:
* gene: gene id
* snp: snp id
* weight: corresponding non-zero weight


## Usage
We provided example input data [here](https://github.com/LidaWangPSU/EXPRESSO/tree/main/example_data).

Data were subsetted from GTEx whole blood tissue as an example to run the script.

Example of R script used to run EXPRESSO can be found [here](https://github.com/LidaWangPSU/EXPRESSO/blob/main/example_data/example.code.R).

## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)
