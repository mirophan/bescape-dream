library(tidyverse)
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
scfilelist <- c('./gep/martin_eset.RDS', './gep/brca_eset.RDS')

celltypevar = 'cellType' #variable name containing the cell type annot in @phenoData of the eset
samplevar = 'SubjectName' # variable name in @phenoData@data$... identifying sample name
celltypesel = c('memory.B.cells','naive.B.cells','memory.CD4.T.cells','naive.CD4.T.cells','regulatory.T.cells','memory.CD8.T.cells','naive.CD8.T.cells','NK.cells','neutrophils','monocytes','myeloid.dendritic.cells','macrophages','fibroblasts','endothelial.cells')

# intersecting set: BRCA x Martin
#celltypesel = c('NK.cells','endothelial.cells','fibroblasts','macrophages','memory.B.cells','memory.CD4.T.cells','memory.CD8.T.cells','myeloid.dendritic.cells','naive.B.cells','naive.CD4.T.cells','naive.CD8.T.cells','regulatory.T.cells')

# intersecting set: BRCA x Smille
#celltypesel = c('NK.cells','endothelial.cells','fibroblasts','macrophages','memory.CD4.T.cells','memory.CD8.T.cells','monocytes','myeloid.dendritic.cells','others','regulatory.T.cells')

input_df <- readr::read_csv("./input/input.csv")

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("./input/", expression_files)

## Load sc datasets: martin and brca
scdata = list(0)
nscdata<-length(scfilelist)
for(index in 1:nscdata){
scdata[[index]]<-readRDS(scfilelist[index])
}


do_scdc <- function(expression_path, dataset_name){
    
    # This reads in the input file and converts to a matrix which will be
    # input to SCDC
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("Gene") %>% 
        as.matrix()

    #creating phenoData with just sample ID (TKT)
    pData<-as.data.frame(colnames(expression_matrix))
    colnames(pData)<-"sample.id"
    metadata<-data.frame(labelDescription=c("sample.id"),row.names=c("sample.id"))
    phenoData<-new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
    row.names(phenoData)<-colnames(expression_matrix)

    bulk<-Biobase::ExpressionSet(assayData=expression_matrix,phenoData=phenoData)#no phenotype data available
    
    # scdc deconvolution - doesnt need to run, just need to make sure the bulk.eset is in the correct format
    ens <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata, ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)

    # actually select best sc-dataset for each sample
    result_matrix <-  as.data.frame(wt_prop(ens$w_table[1, 1:2], ens$prop.only))
    
    #For missing celltypes: Run based on intercept of cell types across scRNASeq datasets and correct
    recalib<-function(resultmat,cells_zero,cell_sum,cells_add){
        sumind<-which(colnames(resultmat) %in% cells_add)
        resultmat<-resultmat %>% dplyr::mutate(!!cells_zero:=0) %>%
            dplyr::mutate(!!cell_sum:=rowSums(resultmat[,sumind]))
        return(resultmat/rowSums(resultmat))#recalibrate based on new Rowsum to add up to one
    }
    #Example
    #result_matrix<-recalib(result_matrix,"neutrophils","monocytes",c("macrophages","myeloid.dentritic.cells"))
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
    
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

## Run MCP-Counter on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_scdc) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)


## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

    
