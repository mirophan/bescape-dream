library(tidyverse)
library(Biobase)
#library(dplyr)
#library(readr)
#library(purrr)
#library(magrittr)
#library(tibble)
#library(tidyr)

library(SCDC)

## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets

#Declare parameters needed for scdc(TODO)
#scfilelist <- c('../docs/datasets/scdc/gep/baron_sc.rds', '../docs/datasets/scdc/gep/segerstolpe_eset.rds')
#scfilelist <- c('./gep/martin_raw_eset_geo.RDS', './gep/brca_raw_eset_geo.RDS')
scfilelist <- c('../unfiltered_data/brca_tumor_eset.RDS',
                '../unfiltered_data/crc_tumor_eset_update.RDS',
                '../unfiltered_data/martin_raw_eset.RDS') #TKT: please don't change the order
crcindex<-c(2,3)
brcaindex<-1

celltypevar = 'cluster' #variable name containing the cell type annot in @phenoData of the eset
samplevar = 'SubjectName' # variable name in @phenoData@data$... identifying sample name
#celltypesel = c('alpha','beta','delta','gamma','acinar','ductal')

# intersecting set: BRCA x Martin
#celltypesel = c('NK.cells','endothelial.cells','fibroblasts','macrophages','memory.B.cells','memory.CD4.T.cells','memory.CD8.T.cells','myeloid.dendritic.cells','naive.B.cells','naive.CD4.T.cells','naive.CD8.T.cells','regulatory.T.cells')

#Union of cells in DREAM (the last celltypes are missing in one of the datasets)
celltypesel = c('NK.cells','endothelial.cells','fibroblasts','macrophages','memory.B.cells','memory.CD4.T.cells','memory.CD8.T.cells','myeloid.dendritic.cells','naive.B.cells','naive.CD4.T.cells','naive.CD8.T.cells','regulatory.T.cells',"monocytes","neutrophils")

# intersecting set: BRCA x Smille
#celltypesel = c('NK.cells','endothelial.cells','fibroblasts','macrophages','memory.CD4.T.cells','memory.CD8.T.cells','monocytes','myeloid.dendritic.cells','others','regulatory.T.cells')

input_df <- readr::read_csv("./input/input.csv")

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

indications <- input_df$cancer.type

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("./input/", expression_files)

scdata_all = list(0)
nscdata<-length(scfilelist)
for(index in 1:nscdata){
scdata_all[[index]]<-readRDS(scfilelist[index])
}




do_scdc <- function(expression_path, dataset_name, indication){
    
    # This reads in the input file and converts to a matrix which will be
    # input to SCDC
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("Gene") %>% 
        as.matrix()

    #creating phenoData with just sample ID (TKT)
    pData<-as.data.frame(colnames(expression_matrix))
    colnames(pData)<-"sample"
    metadata<-data.frame(labelDescription=c("sample"),row.names=c("sample"))
    phenoData<-new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
    row.names(phenoData)<-colnames(expression_matrix)

    bulk<-Biobase::ExpressionSet(assayData=expression_matrix,phenoData=phenoData)#no phenotype data available
    
    # scdc deconvolution - doesnt need to run, just need to make sure the bulk.eset is in the correct format
    #ens <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata, ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)

    if(indication=="CRC"){
        ens <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata_all[c(2,3)], ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T) #TKT: order is defined, please don't change the order, i.e. list 2 and 3 are CRC and Martin respectively
    } else if(indication=="BRCA"){
        ens <- SCDC_prop(bulk.eset = bulk, sc.eset = scdata_all[[1]], ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)
    } else {
        ens <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata_all, ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)
    }
    

    #hack on sc datasets with missing cell types. must be appended at the end and follow the order of the celltype list
    if(indication=="CRC"){
    for(index in 1:2){
        t<-as.data.frame(ens$prop.only[[index]],rownames=TRUE)
        if(!("monocytes" %in% colnames(t))){
            t$monocytes<-0}   
        if(!("neutrophils" %in% colnames(t))){
            t$neutrophils<-0}
        ens$prop.only[[index]]<-as.matrix(t)
        rm(t)}
    }

    if(indication=="CRC"){
        result_matrix <-  as.data.frame(wt_prop(ens$w_table[1, 1:2], ens$prop.only)) #only take 2 scRNAseq
    } else if(indication == "BRCA") {
        result_matrix <-  as.data.frame(ens$prop.est) #No ensemble
    } else {
        result_matrix <-  as.data.frame(wt_prop(ens$w_table[1, 1:nscdata], ens$prop.only))
    } #ensemble with 3 RNASeq
    
    
    
    #For missing celltypes: Run based on intercept of cell types across scRNASeq datasets and correct
    recalib<-function(resultmat,cell_sum,cells_add){
        sumind<-which(colnames(resultmat) %in% cells_add)
        resultmat<-resultmat %>% dplyr::mutate(!!cell_sum:=rowSums(resultmat[,sumind]))
        resultmat<-resultmat %>% dplyr::select(-cells_add)
        return(resultmat/rowSums(resultmat))#recalibrate based on new Rowsum to add up to one
    }
    #Collapsing by celltypes
    result_matrix<-recalib(result_matrix,"B.cells",c("memory.B.cells","naive.B.cells"))
    result_matrix<-recalib(result_matrix,"CD4.T.cells",c("memory.CD4.T.cells","naive.CD4.T.cells","regulatory.T.cells"))
    result_matrix<-recalib(result_matrix,"CD8.T.cells",c("memory.CD8.T.cells","naive.CD8.T.cells"))
    result_matrix<-recalib(result_matrix,"monocytic.lineage",c("macrophages","myeloid.dendritic.cells","monocytes"))
    #result_matrix$neutrophils<-0 #rm hard coded version
    #result_matrix$monocytes<-result_matrix$macrophages+result_matrix$myeloid.dendritic.cells#rm hard coded version
    
    
    #For SCDC based on current output/predictions.csv(TKT)
    #1) put the sample into column, use as key
    #2) convert wide to long format using all the remaining variables(cell types)   
    result_df_wide<- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id")
    result_df_wide$sample.id<-as.factor(result_df_wide$sample.id)
    result_df <- tidyr::gather(result_df_wide,cell.type,prediction,-sample.id,factor_key= TRUE) 
    
    gc()
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
   
}

## Run MCP-Counter on each of the expression files
result_dfs <- purrr::pmap(list(expression_paths, dataset_names, indications), do_scdc) #map2 only takes 2 arguments 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)


## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions_unionCells_ds1ds2_coarsed_wMartin_tissueType.csv")

    
