### News

<ul>
<li>23-11-2018:  cisTopic v0.2.0: Umap dimensionality reduction, new plotting functions, region clustering, predictive distribution and enrichment of region sets (e.g. ChIP-seq signatures, topic-specific cistromes) in cells. 
<li>19-7-2018: Vignette on cisTopic in a big data set (Lake et al., 2018).</li>
<li>19-7-2018: Update for liftOver.</li>
<li>19-7-2018: Update for big data management.</li>
<li>12-7-2018: Updates for mm9, dm3 and dm6 motif enrichment.</li>
</ul>

# cisTopic: Probabilistic modelling of cis-regulatory topics from single cell epigenomics data

cisTopic is an R-package to simultaneously identify cell states and cis-regulatory topics from single cell epigenomics data.

## Dependencies

The following packages have to be installed manually before installing cisTopic:

```r
install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.4.0.tar.gz", repos=NULL)
install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/RcisTarget_1.2.0.tar.gz", repos=NULL)
```

## Installation

For installling and load cisTopic, run:

```r
devtools::install_github("aertslab/cisTopic")
library(cisTopic)
```

## Databases

RcisTarget feather databases are available [here](https://resources.aertslab.org/cistarget/).
