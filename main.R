library(tercen)
library(tidyverse)
source("R/gsea.R")
library(clusterProfiler)
library(gage)
data(egSymb)

# get the relevant input data
getDataFrame = function(ctx){
  # Check if there is only 1 row factor
  id = ctx %>% 
    rselect()
  
  if (dim(id)[2] != 1){
    stop("Require exactly one row variable containing Kinase Names")
  }
  
  # Check if there is only 1 column factor
  contrast = ctx %>% 
    cselect()
  
  if (dim(contrast)[2] != 1){
    stop("Require exactly one column variable containing contrasts")
  }
  
  # Check if there is a color value
  if (length(ctx$colors) != 1){
    stop("Need 1 color containing Kinase Family")
  }
  
  Comparison = ctx$cselect() %>%
    dplyr::mutate(.ci = 0:(n()-1))
  
  KinaseName = ctx$rselect() %>%
    dplyr::mutate(.ri = 0:(n()-1))
  
  family = ctx$select(ctx$colors) %>% 
    cbind(ctx %>% select(.ri, .ci)) %>% 
    distinct(.ri, `Kinase Family`)
  
  df = ctx %>% 
    select(.ri, .ci, .y) %>%
    left_join(Comparison, by = ".ci") %>%
    left_join(KinaseName, by = ".ri") %>%
    left_join(family, by = ".ri") %>%
    dplyr::rename("FC" = ".y") %>%
    rename_with(~ c("Comparison", 'KinaseName', "Family"), .cols = c(4, 5, 6))
  
  return(df)
  
}

# Parameters
fdr_thr = ctx$op.value("FDRThreshold", as.double, 0.3) 

# Data
idmap_kin <- read_csv("data/idmap_id_entrez_uniname_kinase.csv")

# input_df <- read_csv("UKA_test.csv") %>% clean_tercen_columns() %>% filter(contrast == "T1-15P vs C-N")


# Get Tercen data
ctx = tercenCtx()
input_df = getDataFrame(ctx)


# Do GSEA and make downloadable output
gsea_full_result = input_df %>% 
  group_by(.ci, Comparison) %>% # Do GSEA separately for each comparison (contrast)
  do(gsea_per_supergroup(., idmap_kin))

## Make a downloadable export csv in tercen
exported_table <- gsea_full_result %>% select(-.ci)

## Get filename from workflow, group, datastep
nms <- get_names(ctx)

ts <- format(Sys.time(), "%Y-%m-%d-%H%M%S")

if(!is.null(nms$GRP)) {
  filename <- paste(nms$WF, nms$GRP, nms$DS, ts, sep = "_")
} else {
  filename <- paste(nms$WF, nms$DS, ts, sep = "_")
}

## Bring the downloadable csv to Tercen
file_to_tercen(file_path = exported_table, filename = paste0("GSEA_", filename)) %>%
  ctx$addNamespace() %>%
  as_relation(relation_name = "CSV Export") %>%
  as_join_operator(list(), list()) %>%
  save_relation(ctx)


# As output, make a simple kinase-to-pathway table

result <- make_kin_pw_matrix(gsea_full_result, input_df)

result %>% 
  ctx$addNamespace() %>%
  ctx$save()


