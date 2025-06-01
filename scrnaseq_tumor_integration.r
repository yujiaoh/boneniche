
# Integration of mouse tumor scRNA-seq data using Seurat
# Please note that slight variations may occur between runs due to random seed initialization.

#Load package
library(dplyr)
library(Seurat)
library(patchwork)

#Load processed data
pt<- readRDS("path/to/my/file.rds")
bm<- readRDS("path/to/my/file.rds")

pt$orig.ident = 'PrimaryTumor'
bm$orig.ident = 'BoneMetastases'

#Perform integration
combined<-merge(pt, y =bm, add.cell.ids =c("pt", "bm"), project ="ALL")

combine.list <- SplitObject(combined, split.by = "orig.ident")
combine.list <- lapply(X = combine.list, FUN = function(x) {
         x <- NormalizeData(x)
         x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
        })
features <- SelectIntegrationFeatures(object.list = combine.list, 
                                      nfeatures = 5000)
anchors <- FindIntegrationAnchors(object.list = combine.list, 
                                  anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined)<-"integrated"
combined<-ScaleData(combined, 
                    verbose =FALSE)
combined<-RunPCA(combined, 
                 npcs =50, 
                 verbose =FALSE)
combined<-RunUMAP(combined, 
                  reduction ="pca", 
                  dims =1:30)
combined<-FindNeighbors(combined, 
                        reduction ="pca", 
                        dims =1:30)
combined<-FindClusters(combined, 
                       resolution =c(0.2,0.4,0.6))

# Visualization
DimPlot(combined, 
        reduction ="umap", 
        group.by ="orig.ident")
DimPlot(combined, 
        reduction ="umap", 
        group.by =c("integrated_snn_res.0.2")
DimPlot(combined, 
        reduction ="umap", 
        group.by =c("integrated_snn_res.0.4")
DimPlot(combined, 
        reduction ="umap", 
        group.by =c("integrated_snn_res.0.6")

# Save integrated data for further analysis
saveRDS(combined, "path/to/my/file.rds")

#Identify cell type markers
allmarkers <-FindAllMarkers(combined, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
write.csv(allmarkers, "path/to/markers.csv")