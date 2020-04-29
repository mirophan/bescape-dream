library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(tibble)
library(tidyr)

library(SCDC)

## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets

#Declare parameters needed for scdc(TODO)
#scdata
#celltypevar
#samplevar
#celltypesel

input_df <- readr::read_csv("./input/input.csv")

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("./input/", expression_files)

do_scdc <- function(expression_path, dataset_name){
    
    # This reads in the input file and converts to a matrix which will be
    # input to SCDC
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("Gene") %>% 
        as.matrix()

    #creating phenoData with just sample ID
    pData<-as.data.frame(colnames(expression_matrix))
    colnames(pData)<-"sample.id"
    metadata<-data.frame(labelDescription=c("sample.id"),row.names=c("sample.id"))
    phenoData<-new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
    row.names(phenoData)<-colnames(expression_matrix)

    bulk<-Biobase::ExpressionSet(assayData=expression_matrix,phenoData=phenoData)#no phenotype data available
    
    # scdc deconvolution - doesnt need to run, just need to make sure the bulk.eset is in the correct format
    result_matrix <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata, ct.varname = params$celltypevar, sample = params$samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)
    
    # Convert the result matrix back to a dataframe
    result_df <- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("mcpcounter.cell.type") %>% 
        dplyr::as_tibble()
    result_df <- SCDC
    
    # Stack the predictions into one column
    result_df <- tidyr::gather(
        result_df,
        key = "sample.id", 
        value = "prediction", 
        -mcpcounter.cell.type) 

    #For SCDC based on current output/predictions.csv
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

    
