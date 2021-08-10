## DESEQ2 Functions
# functions to map entrez ids for pathway analysis
# also to get results for treatment conditions, including up and down regulated genes
# also perform GO and KEGG pathway analysis via clusterProfiler

library(clusterProfiler)
library(AnnotationDbi)

# map IDS starting with ensembld ids

convertIDs <-
  function(ids,
           from,
           to,
           db,
           ifMultiple = c("putNA", "useFirst")) {
    stopifnot(inherits(db, "AnnotationDb"))
    ifMultiple <- match.arg(ifMultiple)
    suppressWarnings(selRes <- AnnotationDbi::select(
      db,
      keys = ids,
      keytype = from,
      columns = c(from, to)
    ))
    if (ifMultiple == "putNA") {
      duplicatedIds <- selRes[duplicated(selRes[, 1]), 1]
      selRes <- selRes[!selRes[, 1] %in% duplicatedIds,]
    }
    return(selRes[match(ids, selRes[, 1]), 2])
  }

# splitter to strip version number from ensembl ids

splitter <- function(x) {
  split_var <- strsplit(x, ".", fixed = TRUE) [[1]][[1]]
  return(split_var)
}


# function to look at upregulated genes & GO & KEGG Pathways using ClusterProfiler functions

genes_and_paths_up <- function(dds, condition) {
  # results
  res <-
    results(dds, name = condition)
  res_old <- res
  res <-
    res[order(res$padj), ]
  res <-
    res[!is.na(res$padj), ]
  res <-
    res[res$padj < 0.05,]
  
  # pathway analysis .. just upregulated
  res_genes <-
    res[res$log2FoldChange > 0.5, ]
  res_genes <-
    res_genes[order(res_genes$log2FoldChange, decreasing = TRUE), ]
  
  ens_genes <-
    lapply(rownames(res_genes), splitter)
  ens_genes <- as.character(ens_genes)
  
  back <-  rownames(res_old)
  back <- as.character(back)
  
  geneList <- res_genes$log2FoldChange
  
  # Entrez ids for clusterprofiler functions
  
  ent <-
    convertIDs(ens_genes,
               "ENSEMBL",
               "ENTREZID",
               org.Hs.eg.db,
               ifMultiple = c("putNA"))
  
  # GO enrichment
  
  rad_enrich <-  enrichGO(
    gene          = ent,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE,
    minGSSize = 2
    
  )
  
  # KEGG enrichment
  
  rad_enrich_KEGG <-
    enrichKEGG(
      ent,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      # ent_back,
      minGSSize = 2,
      maxGSSize = 500,
      qvalueCutoff = 0.05,
      use_internal_data = FALSE
    )
  
  # annotate res genes / the significant output
  
  syms <-
    convertIDs(ens_genes,
               "ENSEMBL",
               "SYMBOL",
               org.Hs.eg.db,
               ifMultiple = c("putNA"))
  
  res_genes$symbols <- syms
  
  # return a list of (1) significantly upregulated genes (2) enriched GO pathways (3) enriched KEGG pathways (4) the entire results including insignificant results
  out <- list(res_genes, rad_enrich, rad_enrich_KEGG, res_old)
  
}


# down regulated genes and pathways


genes_and_paths_down <- function(dds, condition) {
  # results
  # results
  res <-
    results(dds, name = condition)
  res_old <- res
  res <-
    res[order(res$padj), ]
  res <-
    res[!is.na(res$padj), ]
  res <-
    res[res$padj < 0.05,]
  
  # pathway analysis .. just upregulated
  res_genes <-
    res[res$log2FoldChange < -0.5, ]
  
  
  res_genes <-
    res_genes[order(res_genes$log2FoldChange, decreasing = FALSE), ]
  
  ens_genes <-
    lapply(rownames(res_genes), splitter)
  ens_genes <- as.character(ens_genes)
  
  back <-  rownames(res_old)
  back <- as.character(back)
  
  geneList <- res_genes$log2FoldChange
  
  ent <-
    convertIDs(ens_genes,
               "ENSEMBL",
               "ENTREZID",
               org.Hs.eg.db,
               ifMultiple = c("putNA"))
  
  
  
  
  # GO pathways
  
  rad_enrich <-  enrichGO(
    gene          = ent,
    # changed from ent
    minGSSize = 2,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # KEGG pathways
  
  rad_enrich_KEGG <-
    enrichKEGG(
      ent,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      # ent_back,
      minGSSize = 2,
      maxGSSize = 500,
      qvalueCutoff = 0.05,
      use_internal_data = FALSE
    )
  
  # annotate res genes / the significant output
  syms <-
    convertIDs(ens_genes,
               "ENSEMBL",
               "SYMBOL",
               org.Hs.eg.db,
               ifMultiple = c("putNA"))
  
  res_genes$symbols <- syms
  
  
  out <- list(res_genes, rad_enrich, rad_enrich_KEGG, res_old)
  
}
