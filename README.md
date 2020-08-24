# PhantomPurgeR
R package for the optimal purging of phantom molecules by the robust estimation of the sample index hopping rate in multiplexed droplet-based single-cell RNA-seq data


Paper
-------

[Farouni, R., Djambazian, H., Ferri, L. E., Ragoussis, J., and Najafabadi, H. S. (2020). Model-based analysis of sample index hopping reveals its widespread artifacts in multiplexed single-cell RNA-sequencing. *Nature Communications*.](https://www.nature.com/articles/s41467-020-16522-z)

The bioRxiv preprint can be accessed [here](https://www.biorxiv.org/content/10.1101/617225v5).



Website
---------

The paper's [website](https://csglab.github.io/PhantomPurgeR/) hosts code and R Markdown notebooks implementing the statistical modeling approach

## Installation

In R run:

```{r}
library("devtools")
devtools::install_github("csglab/PhantomPurgeR")
```

Installation Requirements
----------
The code was tested using R 3.6.2 running on Ubuntu 18.04.2 LTS and R 3.5.2 running on CentOS Linux 7.

To run the code, the following R packages are required.

- [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)
- [DropletUtils](https://www.bioconductor.org/packages/release/bioc/html/DropletUtils.html)
- tidyverse
- matrixStats
- broom
- furrr
- tictoc
- data.table
- cowplot
- scales

Contact
---------

We use GitHub [issues](https://github.com/csglab/PhantomPurgeR/issues) for tracking requests and bugs. Please submit an new issue if you any comment or you would like to report a software bug.
