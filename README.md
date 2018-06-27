### News

15-6-2018: Vignette on simulated single cell epigenomes added.

# cisTopic: Probabilistic modelling of cis-regulatory topics from single cell epigenomics data

cisTopic is an R-package to simultaneously identify cell states and cis-regulatory topics from single cell epigenomics data.

## Installation

For installling and load cisTopic, run:

```r
devtools::install_github("aertslab/cisTopic")
library(cisTopic)
```

## Dependencies

cisTopic requires RcisTarget (version >= 1.0.2). You can install it with:

```r
devtools::install_github("aertslab/RcisTarget")
library(RcisTarget)
```

RcisTarget feather databases are available [here](https://resources.aertslab.org/cistarget/).

For installing other Bioconductor dependencies:

```r
source("https://bioconductor.org/biocLite.R")
suppressMessages(biocLite(c('Rsubread',
                            'GenomicRanges',
                            'destiny',
                            'rGREAT',
                            'ChIPseeker',
			    'AUCell'
                            )))
```
