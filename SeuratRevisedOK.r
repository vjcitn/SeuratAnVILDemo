## ----attach--------------------------------------------------------------------------------
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.1-3.14")
# use appropriate methods to install AnVIL, LoomExperiment, reticulate,
# cellgeni/sceasy
#
# this script should be run AFTER those packages are installed

ii = installed.packages()
stopifnot(all(c("AnVIL", "LoomExperiment", "reticulate", "sceasy") %in% rownames(ii)))

system("wget https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.0.2.tar.gz")
install.packages("SeuratObject_4.0.2.tar.gz", repos=NULL, type="source")

suppressPackageStartupMessages({
library(AnVIL)
library(LoomExperiment)
library(reticulate)
library(sceasy)
})


## ----dol-----------------------------------------------------------------------------------
loompy <- reticulate::import('loompy')

# Assign environment variables and view them
project <- avworkspace_namespace()
workspace <- avworkspace()
bucket <- avbucket()
project
workspace
bucket


## ----dosta---------------------------------------------------------------------------------
drs_url <- 'drs://jade-terra.datarepo-prod.broadinstitute.org/v1_4c9d9b2d-a469-40bc-b79b-c3465c154675_57130684-c89c-4b27-9030-298f50cbc95c'
stats <- drs_stat(drs_url)
stats


## ----docop---------------------------------------------------------------------------------
try(drs_cp(drs_url, ".")) # will quietly fail if loom file present


## ----look----------------------------------------------------------------------------------
list.files()


## ----doconv--------------------------------------------------------------------------------
sceasy::convertFormat('sc-landscape-human-liver-10XV2.loom', from="loom", to="anndata",
                       outFile='trial_liver.h5ad')

list.files()


## ----doconv2-------------------------------------------------------------------------------
sceasy::convertFormat('trial_liver.h5ad', from="anndata", to="seurat",
                       outFile='liver_seurat.rds')


## ----dorest--------------------------------------------------------------------------------
# Load libraries needed for using Seurat
library(dplyr)
library(Seurat)
library(patchwork)

list.files()

liver <- readRDS(file = "liver_seurat.rds")

# Look at the Seurat object
liver

liver = subset(liver, subset = input_id == "02f207ae-217b-42ce-9733-c03e61541fcc")
liver

summary(liver$input_id)

# Visualize QC metrics as a violin plot
#VlnPlot(liver, features = "pct_mitochondrial_molecules", ncol = 1)
VlnPlot(liver, features = "n_molecules", ncol = 1)

# Filtering cells
liver <- subset(liver, subset = nFeaturess_RNA > 200 & nFeaturess_RNA < 2500 & pct_mitochondrial_molecules < 5)

liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)


liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(liver), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(liver)

plot1

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot2

# This step takes several minutes
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

liver <- RunPCA(liver, features = VariableFeatures(object = liver))

ElbowPlot(liver)

# FindNeighbors uses the first 10 principle components
liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
liver <- RunUMAP(liver, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(liver, reduction = "umap")

saveRDS(liver, file = "liver_seurat_clusters.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(liver, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# CHANGED
library(magrittr)
library(dplyr)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Example of setting up differential expression testing with ROC method
cluster2.markers <- FindMarkers(liver, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# Save Seurat object
saveRDS(liver, file = "liver_final.rds")

# Copy the rds file generated in the notebook into the workspace bucket using AnVIL package gsutil commands
gsutil_cp("./*.rds", bucket)
# Run list command to see if file is in the bucket
gsutil_ls(bucket)

