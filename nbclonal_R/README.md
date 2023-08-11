# nbclonal

R package for calculating clonal mutations and estimating transmission bottleneck.

Some of the functions in this package uses `parallrl:mclapply()` to reduce the running time. Unfortunately, `mclapply()` can only run on mac but not on windows. If you are using an operating system other than mac, you may need to adjust the functions using `mclapply()` or simply change them into `lapply()` and modify the arguments.

# How to install

You can install this package using `devtools`:

```r
devtools::install_github("MXTeresa/nbclonal")
```

# Overview

You can inspect all functions included in the package by:

```r
library("nbclonal")
ls("package:nbclonal")
```

# Sample usage and plotting

Please see the file [NbClonal.Rmd](https://github.com/MXTeresa/nbclonal/blob/50ef89de021f88cebcc56960e791efbc2c83b04b/NbClonal.Rmd) and the generated pdf file [NbClonal.pdf](https://github.com/MXTeresa/nbclonal/blob/50ef89de021f88cebcc56960e791efbc2c83b04b/NbClonal.pdf) for a detailed description of each function, how they are developed, and ways to plot the results. Since some of the functions has been updated after generating this pdf, the functions may vary a little bit as the ones included in the package. Please use the final version of the functions listed in the package instead of copy-pasting functions from the pdf and rmd file.
