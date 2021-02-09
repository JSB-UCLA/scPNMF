# scPNMF

`scPNMF` is a method to facilitate gene selection for targeted gene profiling by learning a sparse gene encoding of single cells. Compared with existing gene selection methods, `scPNMF` has two advantages. First, its selected informative genes can better distinguish cell types, with a small number, e.g., < 200 gene. Second, it enables the alignment of new targeted gene profiling data with reference data in a low-dimensional space to help the prediction of cell types in the new data.

## Installation

`scPNMF` can be installed from Github with the following code in `R`:

``` r
install.packages("devtools")
library(devtools)

install_github("JSB-UCLA/scPNMF", build_vignettes = TRUE)
```

## Usage

For detailed info on `scPNMF` method and applications, please check out the package vignettes with the following code in `R`: 

``` r
browseVignettes("scPNMF")
```

## Contact

Any questions or suggestions on `scPNMF` are welcomed! Please report it on [issues](https://github.com/JSB-UCLA/scPNMF/issues), or contact Dongyuan Song (<dongyuansong@ucla.edu>) or Kexin Li (<aileenlikexin@outlook.com>).

