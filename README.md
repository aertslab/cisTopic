### News
  - ***2019-4-12***:
    - New functions for CellRanger ATAC input and tutorial on the 5k PBMC data set from 10X . 
      
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

## Tutorials
  - [Basic tutorial on simulated single cell epigenomes from melanoma cell lines](https://rawcdn.githack.com/aertslab/cisTopic/f628c6f60918511ba0fa4a85366ebf52db5940f7/vignettes/CompleteAnalysis.html).
  - [10X tutorial on the 5k PBMC data set from 10X](https://rawcdn.githack.com/aertslab/cisTopic/86e54ce1aa8ebd1836ca8566eb361d1141a0462d/vignettes/10X_workflow.html).