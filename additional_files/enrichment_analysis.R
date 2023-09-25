### Stavros Giannoukakos ###
### ELLBA - Gene Enrichment analysis v0.0.2

args <- commandArgs(TRUE)
if (length(args) == 3) {
  # Input of the biological rel. matrix
  gene_list <- unlist(strsplit(args[1], ","))
  # Output directory where all analysis will be stored
  main_outdir <- args[2]
  # Feature matrix to use
  feature_mat <- args[3]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


suppressPackageStartupMessages({
  library("gprofiler2")
})


gene_list <- gsub("\\..*","",gene_list)

# Performing functional profiling of the selected gene list
gostlink <- gost(query = gene_list, organism = "hsapiens", ordered_query = F, 
                significant = TRUE, measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", as_short_link = T) # , sources=c('GO','KEGG', 'REAC')


if (feature_mat=="GeneExpression") {
  cat("Use the following links to proceed with the online g:Profiler analysis.\nEach link corresponds to the aforementioned feature matrix:\n\n", 
    file=paste0(main_outdir,'/gProfiler_click_link.txt'), append = TRUE)
}

cat(paste0(feature_mat,": ", gostlink, "\n"), file=paste0(main_outdir,'/gProfiler_click_link.txt'), append = TRUE)
