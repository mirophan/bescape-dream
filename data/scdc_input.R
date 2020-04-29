library(dplyr)
library(readr)
library(purrr)
library(magrittr)
library(tibble)
library(tidyr)

library(SCDC)

## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets

input_df <- readr::read_csv("dream/input/input.csv")

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)

do_scdc <- function(expression_path, dataset_name){
    
    # This reads in the input file and converts to a matrix which will be
    # input to SCDC
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("Gene") %>% 
        as.matrix() 
    
    # scdc deconvolution
    result_matrix <- MCPcounter::MCPcounter.estimate(
        expression_matrix,
        probesets = probesets,
        genes = genes,
        featuresType = 'HUGO_symbols')
    
    # Convert the result matrix back to a dataframe
    result_df <- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("mcpcounter.cell.type") %>% 
        dplyr::as_tibble()
    
    # Stack the predictions into one column
    result_df <- tidyr::gather(
        result_df,
        key = "sample.id", 
        value = "prediction", 
        -mcpcounter.cell.type) 
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

## Run MCP-Counter on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_scdc) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)


## Write result into output directory
readr::write_csv(combined_result_df, "output/predictions.csv")

    
