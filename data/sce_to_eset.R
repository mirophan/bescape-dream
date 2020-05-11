library(Seurat)
library(BisqueRNA)

sce_to_seurat <- function(sce_path, sc_anno_path, filename='eset.RDS'){
    sce <- readRDS(sce_path)
    sc_anno <- read.csv(sc_anno_path, header=FALSE)

    seurat <- as.Seurat(sce, counts=NULL, data='X')

# TODO add celltype_dream as param
    seurat@meta.data$celltype_dream <- sc_anno$V1
    Idents(seurat)<-"celltype_dream"
    seurat@assays$RNA@counts<-seurat@assays$RNA@data
    saveRDS(seurat, file=filename)
    return(seurat)

}

seurat_to_eset <- function(seurat, delim='-', idx=2, filename='seurat.RDS'){
    # Idents(liver_h5ad) %>% head() # to check sample names
    out.eset <- BisqueRNA::SeuratToExpressionSet(seurat, delim, idx, version = 'v3')
    # fix GeneID
    out.eset@featureData@data$GeneID <- featureNames(out.eset)
    #assign cluster to just characters
    out.eset@phenoData@data$cluster<-as.character(out.eset@phenoData@data$cellType)
    #update vardata
    varMetadata(out.eset)<-data.frame(labelDescription=colnames(out.eset@phenoData@data),stringsAsFactors=FALSE)

    out.eset@featureData@varMetadata<-data.frame(labelDescription="GeneID",stringsAsFactors=FALSE)

    rownames(out.eset@featureData@varMetadata)<-"GeneID"

    saveRDS(out.eset, file=filename)
    return(out.eset)
}

#martin.seurat <- sce_to_seurat(sce_path='martin_raw_sce.RDS', sc_anno_path='martin_geosketch_scanno.csv', filename='martin_raw_seurat_geo.RDS')
#martin.eset <- seurat_to_eset(martin.seurat, delim='_', idx=1, filename='martin_raw_eset_geo2.RDS')


#brca.seurat <- sce_to_seurat(sce_path='brca_raw_sce.RDS', sc_anno_path='brca_geosketch_scanno.csv', filename='brca_raw_seurat_geo.RDS')
brca.eset <- seurat_to_eset(brca.seurat, delim='-', idx=2, filename='brca_raw_eset_geo.RDS')


