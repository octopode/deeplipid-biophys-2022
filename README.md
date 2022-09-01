# Phase behavior of plasmalogens defines pressure specialization in deep-sea invertebrates
## J. R. Winnikoff, D. Milshteyn, M. A. Pedraza Joya, A. M. Armando, O. Quehenberger, A. Sodt, E. A. Dennis, E. Lyman, S. H. D. Haddock, and I. Budin

### What's here

This repo contains the raw data and analysis scripts used to generate the figures and statistical claims in the paper.

### How to use it

Start by opening the `.Rproj` file in the top-level directory. This will set the working directory in R so that `here()` calls in the scripts open the correct data files and helper scripts. From within RStudio, navigate to `02-scripts` and open the `.R` files for the analyses/figures you want to run or modify

### Dependencies

The SAXS profiles (`.dat`) files in this repo were generated using [BioXTAS RAW](https://sourceforge.net/projects/bioxtasraw/).
Original HDF5 image files are very large and can be furnished upon request.

High-pressure fluorescence spectroscopy (HPFS) data used in the paper were collected using the [spectackler](https://github.com/octopode/spectackler) suite of Python drivers and scripts, which I also maintain.

All other analyses in the paper were carried out using the R core in [RStudio](https://www.rstudio.com/products/rstudio/download/) with the [tidyverse](https://www.tidyverse.org/) family of libraries plus [phytools](https://cran.r-project.org/web/packages/phytools/index.html).

### Support

If you are unable to reproduce an analysis or otherwise have questions, please raise an issue or email jwinnikoff@gmail.com.