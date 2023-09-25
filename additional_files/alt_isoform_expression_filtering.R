### Stavros Giannoukakos ###
### ELLBA - Alternative Isoform Expression filtering

args <- commandArgs(TRUE)
if (length(args) == 2) {
  # Input Alt. Isoform Expression Matrix
  matfile <- args[1]
  # Output directory where all analysis will be stored
  main_outdir <- args[2]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


suppressPackageStartupMessages({
  library("dplyr")
  library("reshape")
  library("data.table")
})


alt_isoform_expression <- data.frame(fread(matfile))
row.names(alt_isoform_expression) <- alt_isoform_expression$GeneID; alt_isoform_expression$GeneID <- NULL

# Remove rows with NAs
altIso <- alt_isoform_expression[!(rowSums(is.na(alt_isoform_expression))), ]

# Output the filtered matrix
write.table(data.frame("GeneID"=rownames(altIso), altIso), file=paste(main_outdir,"/alternative_isoform_expression_matrix.filtered.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)
