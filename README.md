# Kasumi <img src="https://www.dropbox.com/scl/fi/gpnvjqsumo5nxeo0kzfmz/kasumi_badge.png?rlkey=6gda70ivmy5biml7qnaujx1ai&raw=1" align="right" height = "139">

<!-- badges: start -->
<!-- badges: end -->

## Overview

Kasumi is a method for the identification of spatially localized  neighborhoods of intra- and intercellular relationships, persistent across  samples and conditions. Kasumi learns compressed explainable representations of spatial omics samples while preserving relevant biological signals that are readily deployable for data exploration and hypothesis generation, facilitating translational tasks.

## System Requirements

**kasumi** requires a standard configuration and enough RAM to store the analyzed dataset and to support in-memory operations.

The package requires R version 4.1 or higher. This package is developed on macOS Sonoma. The package is compatible with Linux and maxOS operating systems and  should be compatible with Windows.

## Installation

You can install the latest version from GitHub with `remotes`:

```r
# install.packages("remotes")
remotes::install_github("jtanevski/kasumi")
```

**kasumi** currently depends on **mistyR** version 1.99.10 or higher. You can install the latest version from GitHub:

```r
remotes::install_github("jtanevski/mistyR")
```

## Citation
If you use **kasumi** for your research please cite the [following publication](https://doi.org/10.1101/2024.03.06.583691): 

> Jovan Tanevski, Loan Vulliard, Felix Hartmann, Julio Saez-Rodriguez. Learning tissue representation by identification of persistent local patterns in spatial omics data. bioRxiv 2024.03.06.58369 (2024). https://doi.org/10.1101/2024.03.06.583691
