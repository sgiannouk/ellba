### Stavros Giannoukakos ###
### ELLBA - Gene Fusion filtering v0.0.2

args <- commandArgs(TRUE)
if (length(args) == 2) {
  # Input Fusion Matrix
  matfile <- args[1]
  # Output directory where all analysis will be stored
  main_outdir <- args[2]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


suppressPackageStartupMessages({
  library("data.table")
  library("tidyr")
  library("dplyr")
})


fusions <- fread(matfile, data.table=FALSE)

# Separating gene fusion1 from gene fusion2
fusions <- fusions %>%  extract(fusion, into = c("fusion1", "fusion2"), "(.*)_([^_]+)$")
fusions$fusion2 <- NULL  # Eliminating fusion1 column
# Collapsing by fusion1
fusions <- data.frame(dplyr::group_by(fusions, fusion1) %>% dplyr::summarise_all(sum))
row.names(fusions) <- fusions$fusion1; fusions$fusion1 <- NULL
# Replace all values above 1 with 1
fusions[fusions > 1] <- 1
# Mild filtering where a fusion must be expressed 
# in at least the 5% of the population
population = length(fusions)
fusions <- fusions[rowSums(fusions) >= round( (population * 5)/100, 0), ]

if (nrow(fusions) == 0) {
  print("No fusion genes were detected..")
} else {
  write.table(data.frame("fusion"=rownames(fusions), fusions), file=paste(main_outdir,"/gene_fusion_expression_matrix.filtered.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)
}
