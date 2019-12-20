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
  - [Basic tutorial on simulated single cell epigenomes from melanoma cell lines](https://rawcdn.githack.com/aertslab/cisTopic/f628c6f60918511ba0fa4a85366ebf52db5940f7/vignettes/CompleteAnalysis.html). Data available [here](https://drive.google.com/drive/folders/18ETGIKgXkILo3Xfv9KuysOMqchmSfFX2?usp=sharing).
  - [10X tutorial on the 5k PBMC data set from 10X](https://rawcdn.githack.com/aertslab/cisTopic/8d15fa2813312aa0b20c1042604079558829e947/vignettes/10X_workflow.html). Data available [here](https://drive.google.com/drive/folders/1QORpLPsXejva3oFhECLrnAh5a7FVmJF1?usp=sharing).
  - [Running GREAT and motif enrichment with the mm10 and hg38 genome assemblies](https://rawcdn.githack.com/aertslab/cisTopic/a83e1cd7ab45283890e944ee9754fca401a5c424/vignettes/Runningwithmm10andhg38.html).
  
### Version 3
 - [Basic tutorial on simulated single cell epigenomes from melanoma cell lines](https://rawcdn.githack.com/aertslab/cisTopic/e7d1e4f835bae95d859d92d509f2ac92d5bee732/vignettes/WarpLDA_CompleteAnalysis.html). Data available [here](https://drive.google.com/drive/folders/18ETGIKgXkILo3Xfv9KuysOMqchmSfFX2?usp=sharing).
