
gsea_per_supergroup <- function(input_df, idmap_kin){
  
  comp <- unique(input_df$Comparison)
  
  
  ctx$log(paste0('Working on comparison: ', comp))
  
  # Parse input to create vector of fold changes with entrezids
  uka_fcv <- parse_uka(input_df, idmap_kin = idmap_kin)
  
  # Do GSEA
  gseaKEGG <- gseKEGG(geneList = uka_fcv,
                      organism = "hsa",
                      pvalueCutoff = fdr_thr,
                      pAdjustMethod = "fdr",
                      verbose = F)
  
  
  gsea_result <- tryCatch(
    {
      gseaKEGG@result
    },
    error=function(e) {
      print(e)
      
      ctx$log(paste0("For comparison: ", comp, " no enrichment was found. Try increasing the FDR cutoff. 
                   Skipping to next comparison."))
      
      return(NULL)
    }
  )
  
  # Make result table
  parse_gsea_clusterprofiler_res(gsea_result) %>% mutate(Comparison = comp)
  
}

parse_uka <- function(input_df, idmap_kin){
  input_df %>%
    left_join(idmap_kin, by = c("KinaseName" = "uniprotname")) %>%
    dplyr::select(entrezid, FC) %>%
    arrange(desc(FC)) %>%
    tibble::deframe()  
}


make_kin_pw_matrix <- function(pw_results, input_df){
  
  df_fc <- pw_results %>%
    # Make a df with pathway to single kinase
    mutate(core_genes = str_split(core_genes, ", ")) %>%
    unnest(core_genes) %>%
    rename("KinaseName" = core_genes)  %>% # Rename for joining
    left_join(input_df, by = c("KinaseName", "Comparison")) %>% # Join input to get FC
    arrange(Comparison, p.val, family) # Would this work in Tercen?
  # %>%
  #   pivot_wider(names_from = "KinaseName", values_from = 'FC', values_fill = 0) 
  
  return(df_fc)
  
}



parse_gsea_clusterprofiler_res <- function(df){
  res <- df %>% dplyr::rename("pathwayname" = "Description",
                              'core_genes_entrezid'= "core_enrichment",
                              "pathwayid" = "ID",
                              "FDR" = "p.adjust",
                              "p.val" = 'pvalue',
                              'set.size' = "setSize") %>%
    dplyr::select(-leading_edge, -qvalue, -rank) %>%
    mutate(core_genes = NA) %>%
    mutate(direction = ifelse(sign(NES)>=0, "UP", "DOWN"))
  res$core_genes_e <- gsub("/", ", ", res$core_genes_e)
  
  for (i in seq_along(res$pathwayid)){
    coregenes_entrez <- strsplit(res[i, "core_genes_e"], ", ")[[1]]
    pwgenes <- eg2sym(coregenes_entrez)
    pwgenes <- paste(pwgenes, collapse = ", ")
    res$core_genes[res$pathwayid == res$pathwayid[i]] <- pwgenes
  }
  rownames(res) <- NULL
  return(res)
}

