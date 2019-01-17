### News
  - ***2018-11-23***:
    - cisTopic v0.2.0: Umap dimensionality reduction, new plotting functions, region clustering,
      predictive distribution and enrichment of region sets (e.g. ChIP-seq signatures, topic-specific cistromes) in cells. 
  - ***2018-07-19***:
    - Vignette on cisTopic in a big data set (Lake et al., 2018).
    - Update for liftOver.
    - Update for big data management.
  - ***2018-07-12***:
    - Updates for mm9, dm3 and dm6 motif enrichment.
# cisTopic: Probabilistic modelling of cis-regulatory topics from single cell epigenomics data
cisTopic is an R-package to simultaneously identify cell states and cis-regulatory topics from single cell epigenomics data.
## Dependencies (for R < 3.5)
The following packages have to be installed manually before installing cisTopic:
```r
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/AUCell")
```
## Installation
For installing and loading cisTopic, run:
```r
devtools::install_github("aertslab/cisTopic")
library(cisTopic)
```
## Databases
RcisTarget feather databases are available at [https://resources.aertslab.org/cistarget/](https://resources.aertslab.org/cistarget/).

## Tutorial
You can find the rendered vignette [here](https://rawcdn.githack.com/aertslab/cisTopic/f6091ff94c79dfc927ae8eddf201d10a3f092923/vignettes/CompleteAnalysis.html).
