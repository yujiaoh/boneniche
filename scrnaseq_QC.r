# Pre-process of scRNA-seq data

#Load package
library(dplyr)
library(Seurat)
library(patchwork)
library(scDblFinder)

#Load data
data.dir <- "path/to/my/file"
obj = Read10X_h5(glue::glue("{data.dir}/filtered_feature_bc_matrix.h5"))

#Create project
Basic.Seurat.Preprocessing = function(obj.data){
    obj = CreateSeuratObject(
        counts = obj.data,
        min.cells = 3,
        min.features = 200)
    obj[['percent.mt']] = PercentageFeatureSet(obj, pattern="^mt-")
    obj}

obj <- Basic.Seurat.Preprocessing (obj)

#QC plot
VlnPlot(obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

#filtering low quality data
obj <- subset(obj, 
              subset = nFeature_RNA > 200 & 
              nFeature_RNA < 8000 & 
              percent.mt < 10)

#Find doublets
Basic.Doublet.Finder = function(obj){
    sce = as.SingleCellExperiment(obj)
    sce = scDblFinder(
        sce,
        samples='orig.ident',
        clusters = FALSE,
        BPPARAM= BiocParallel::MulticoreParam(10))
    sce}
obj.dbl = Basic.Doublet.Finder(obj)
obj@meta.data$scDblFinder.score = obj.dbl$scDblFinder.score
obj@meta.data$scDblFinder.class = obj.dbl$scDblFinder.class

#Remove doublets
obj = subset(obj,subset = (scDblFinder.class == 'singlet'))

#Save .rds for integration
saveRDS(obj, file = "path/to/my/file/obj.rds")