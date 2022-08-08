# WGCNA-shinyApp

[![R version](https://img.shields.io/badge/R-v4.1.1-salmon)](https://www.r-project.org) [![TBtools version](https://img.shields.io/badge/TBtools-%3Ev1.09-greenyellow)](https://www.yuque.com/cjchen/hirv8i/fzc4g9) ![lifecycle](https://img.shields.io/badge/lifecycle-Experimental-lightcyan) [![license](https://img.shields.io/badge/license-MIT-red)](https://opensource.org/licenses/MIT) [![Myblog](https://img.shields.io/badge/Blog-ShanwLearnBioinfo-purple)](http://www.shawnlearnbioinfo.top/)

A shiny app for WGCNA...

# Getting started

R version: `>4.1.1`

OS: `MacOS > 10.10`, `Win 7-11`, `linux must have a graphic interface`

``` bash
# clone this repo to your machine
git clone git@github.com:ShawnWx2019/WGCNA-shinyApp.git WGCNAshiny

cd WGCNAshiny

## Method 1.

Rscript WGCNAbyClick.v1.R

## Method 2. open WGCNAbyClick.v1.R by Rstudio or other IDE you perfer and run this script.
```

![](images/paste-E878B059.png)

![](images/paste-9FA8BCAD.png)

![](images/paste-F0179D86.png)

![](images/paste-52ABF860.png)

# Input data prepare

## datExpr

**Data source:**

-   transcriptomics

    -   readcount.

    -   normalized readcount (FPKM, RPKM, TPM, CPM)

    -   microarray data

-   metabolomics

    -   peak area.

-   proteomics,

    -protein abundance.

-   ...

**Format:**

-   Gene/metabolite/protein ID in row and sample ID in column.

-   The sample ID should not contain spaces, special symbols, and should not start with numbers.

-   DO NOT use pure numbers as gene/metabolite/protein ID.

-   Only accepted tab-delimited file, such as `.txt` or `.tsv`.
