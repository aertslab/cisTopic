## News
### ***2019-12-20***:
- cisTopic v3: Faster topic modelling based on WarpLDA (see vignettes for details).
- The function runModels() is deprecated. Use runCGSModels() for modelling based on Collapsed Gibbs Sampling (equivalent to runModels()), or runWarpLDAModels() (for modelling based on WarpLDA).
- Version 2 objects (with or without models) can be used and analyzed with version 3.
      
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
### Version 2
  - [Basic tutorial on simulated single cell epigenomes from melanoma cell lines](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/CompleteAnalysis.html). Data available [here](https://drive.google.com/drive/folders/18ETGIKgXkILo3Xfv9KuysOMqchmSfFX2?usp=sharing).
  - [10X tutorial on the 5k PBMC data set from 10X](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/10X_workflow.html). Data available [here](https://drive.google.com/drive/folders/1QORpLPsXejva3oFhECLrnAh5a7FVmJF1?usp=sharing).
  - [Running GREAT and motif enrichment with the mm10 and hg38 genome assemblies](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/Runningwithmm10andhg38.html).
  
### Version 3
 - [Basic tutorial on simulated single cell epigenomes from melanoma cell lines](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_CompleteAnalysis.html). Data available [here](https://drive.google.com/drive/folders/18ETGIKgXkILo3Xfv9KuysOMqchmSfFX2?usp=sharing).
 - [10X tutorial on the 5k PBMC data set from 10X](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_10X_workflow.html). Data available [here](https://drive.google.com/drive/folders/1QORpLPsXejva3oFhECLrnAh5a7FVmJF1?usp=sharing).
 

