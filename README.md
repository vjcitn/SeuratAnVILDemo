# SeuratAnVILDemo
script that runs in Bioc 3.14 environment, except version 4.0.2 (a downgrade relative to what is on CRAN 12/13/2021) of SeuratObject must be used

In a Bioconductor 3.14 cloud environment

- BiocManager::install(c("AnVIL", "LoomExperiment", "cellgeni/sceasy"), update=FALSE, ask=FALSE)
- source the script SeuratRevisedOK.r
