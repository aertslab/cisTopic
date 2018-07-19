### News

<ul>
<li>19-7-2018: Vignette on cisTopic in a big data set (Lake et al., 2018).</li>
<li>19-7-2018: Update for liftOver.</li>
<li>19-7-2018: Update for big data management.</li>
<li>12-7-2018: Updates for mm9, dm3 and dm6 motif enrichment.</li>
</ul>

# cisTopic: Probabilistic modelling of cis-regulatory topics from single cell epigenomics data

cisTopic is an R-package to simultaneously identify cell states and cis-regulatory topics from single cell epigenomics data.

## Installation

For installling and load cisTopic, run:

```r
devtools::install_github("aertslab/cisTopic")
library(cisTopic)
```

## Dependencies

cisTopic requires AUCell (version >= 1.2.4) and RcisTarget (version >= 1.0.2). You can install them with:

```r
devtools::install_github("aertslab/AUCell")
library(AUCell)
```

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
                            'ChIPseeker'
                            )))
```

Finally, cisTopic requires NMF (version >= 0.23.6). It can be installed with:

```r
devtools::install_github('renozao/NMF@devel')
library(NMF)
```
