---
title: "How to install the necessary packages"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This file will explain how to install the deconvolution tools used in this project, as well as other relevant tools used.

## R Based Tools

The following deconvolution tools are availale as R packages:

| Tool      | How to install                                  | Latest Version Tested to Work |           |             |
| --------  | ----------------------------------------------- | ----------------------------- | --------- | ----------- |
| MuSiC     | `devtools::install_github('xuranw/MuSiC')`      | 0.1.1                         | [code](https://xuranw.github.io/MuSiC/index.html)  | [article](https://doi.org/10.1038/s41467-018-08023-x) |
| BSeq-SC   | `devtools::install_github('shenorrlab/bseqsc')` | 1.0                           | [code](https://shenorrlab.github.io/bseqsc/index.html)  | [article](https://doi.org/10.1016/j.cels.2016.08.011) |
| SCDC      | `devtools::install_github("meichendong/SCDC")`  | 0.0.0.9000                    | [code](https://meichendong.github.io/SCDC/index.html)  | [article](https://doi.org/10.1093/bib/bbz166) |
| BisqueRNA | `install.packages('BisqueRNA')`                 | 1.0.4                         | [code](https://github.com/cozygene/bisque)  | [article](https://doi.org/10.1038/s41467-020-15816-6) |
| MOMF      | `devtools::install_github('sqsun/MOMF')`        | 0.2.0                         | [code](https://github.com/sqsun/MOMF)  | [article](https://doi.org/10.3390/cells8101161) |
| CPM       | `install.packages('scBio')`                     | 0.1.6                         | [code](https://github.com/amitfrish/scBio)  | [article](https://doi.org/10.1038/s41592-019-0355-5) |

Many of the above tools need the following packages to be installed: *devtools* (`install.packages('devtools')`) and *BiocManager* (`install.packages('BiocManager')`).


## Python Based Tools

The following deconvolution tools are available as python modules (**python3!**):

| Tool      | How to install                     | Latest Version Tested to Work |           |             |
| --------  | ---------------------------------- | ----------------------------- | --------- | ----------- |
| AutoGeneS | `python3 -m pip install autogenes` | 1.0.4                         | [code](https://github.com/theislab/AutoGeneS)  | [article](https://doi.org/10.1101/2020.02.21.940650) |
| Scaden    | `python3 -m pip install scaden`    | 0.9.4                         | [code](https://github.com/KevinMenden/scaden)  | [article](https://doi.org/10.1126/sciadv.aba2619) |

An 'interface' to these tools via R was implemented, so that as many tools as possible would be ran in the same place. For that the R package *reticulate* must be installed: `install.packages('reticulate')`.

## Web Based Tools

[CIBERSORTx](https://cibersortx.stanford.edu) was the only web based deconvolution tool used in this project. 

## Modified Scripts

[DWLS](https://github.com/dtsoucas/DWLS) was downloaded from its GitHub page and minor changes were done to the script in order to run it along with the other tools. The modified version of DWLS is present in the *R/* folder.






