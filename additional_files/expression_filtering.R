### Stavros Giannoukakos ###
### ELLBA - Gene Expression filtering

args <- commandArgs(TRUE)
if (length(args) == 6) {
  # Input Gene Expression Matrix
  matfile <- args[1]
  # Input of the Clinical data
  clin_data <- args[2]
  # Output directory where all analysis will be stored
  main_outdir <- args[3]
  # Number of cores to use
  num_threads <- as.numeric(args[4])
  # Type of analysis
  typeof <- as.character(args[5])
  # Groups
  groups <- unlist(strsplit(args[6], ","))
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


suppressPackageStartupMessages({
  library("edgeR")
  library("DESeq2")
  library("RUVSeq")
  library("reshape")
  library("ggplot2")
  library("data.table")
  library("RColorBrewer")
  library("BiocParallel")
  library("dplyr")
})


plotPCAfunc <- function(normalised_data, norm_type){
  
              dt <- newSeqExpressionSet(normalised_data, phenoData = clinical_data)
              
              jpeg(file.path(outdir, paste0("Fig3.A.", norm_type, ".rle.plot.jpg")), width = 20, height = 6, units = 'in', res = 600)
              par(mar=c(8,3,3,3))
              plotRLE(dt, outline=FALSE, ylim=c(-3, 3), col=colors[factor(dt$Classification)],
                      cex.axis =.8, las=2.5, main = paste0(norm_type," normalisation"))
              dev.off()
              jpeg(file.path(outdir, paste0("Fig3.B.", norm_type, ".pca.plot.jpg")), width = 20, height = 6, units = 'in', res = 600)
              par(mar=c(5,3,3,3))
              plotPCA(dt, col=colors[factor(dt$Classification)], k=2, labels=F, pch=19, cex=1.5,
                      main = paste0(norm_type," normalisation"))
              dev.off()

 }



paletteColors <- c("#F8B323","#656A59","#46B2B5","#8CAA7E","#D36F68","#826276","#FF9DA7","#6F5438","#8A8B79","#800000",
                   "#57826F","#785D37","#62BBA5","#FFB84D","#AAA488","#B2432F","#3A6589","#9B5672","#908150","#168E7F",
                   "#C29365","#233B43","#65ADC2","#6E8E84","#2D6D66","#59A77F","#5D8CA8","#D5695D","#D3BA68","#65A479",
                   "#C17529","#AD84C6","#6F8183","#8784C7","#84ACB6","#5D739A","#A5B592","#F3A447","#E7BC29","##E28394",
                   "#D092A7","#9C85C0","#809EC2","#94B6D2","#DD8047","#A5AB81","#D8B25C","#7BA79D","#968C8C","#7F8FA9",
                   "#F07F09","#9F2936","#1B587C","#4E8542","#604878","#C19859","#6997AF","#645135","#439EB7","#79A8A4")


register(MulticoreParam(num_threads))
outdir <- file.path(main_outdir, "plots")
dir.create(outdir, showWarnings = FALSE)

clinical_data <- read.table(clin_data, sep="\t", header=T)
clinical_data <- clinical_data[!duplicated(clinical_data$SampleName), ]
rownames(clinical_data) <- clinical_data$SampleName

# Loading Gene Expression
expression_mat <- data.frame(fread(matfile, sep='\t', header=T, nThread=num_threads))
row.names(expression_mat) <- expression_mat[ ,1]; expression_mat[ ,1] <- NULL
expression_mat <- round(expression_mat, digits = 0)

# Keep only data that are in the expression file
clinical_data <- clinical_data[clinical_data$SampleName %in% colnames(expression_mat), ]
# Output the filtered clinical data
if (typeof == "gene") {write.table(clinical_data, file=paste(main_outdir,"/clinical_data.filtered.tsv", sep=""), sep="\t", row.names = F, quote=FALSE)}

# Drop unnecessary columns
clinical_data <- clinical_data %>% dplyr::select(SampleName, Classification, Abbreviation, Institution)

# Order clinical data based on Classification and SampleName
clinical_data <- clinical_data[order(clinical_data$Classification, clinical_data$SampleName), ]
clinical_data[2:length(clinical_data)] <- lapply(clinical_data[2:length(clinical_data)], factor)  # Converting to factors

# Ordering the gene expression data based on the clinical data  
expression_mat <- expression_mat[ ,order(match(colnames(expression_mat),  rownames(clinical_data) )), drop = FALSE      ]
rm(clin_data, matfile, num_threads)




################### Exploratory Analysis ################### 
# Colour pallet to use on the next plots
colors <- c('#8FBDD3','#BB6464')


# PCA plot with library size and heatmap
dds <- DESeqDataSetFromMatrix(countData = expression_mat, colData = clinical_data , design = ~ Classification)
rld <- vst(dds, blind=TRUE, fitType='local')
data_pca <- plotPCA(rld, intgroup = "Classification", returnData=TRUE)
percentVar <- round(100 * attr(data_pca, "percentVar"))
data_pca <- merge(data_pca, clinical_data[ ,c(1,4)], by.x=0, by.y=0)
row.names(data_pca) <- data_pca$Row.names; data_pca$Row.names <- NULL;  data_pca$SampleName <- NULL

# Plotting the PCA
p1 <- ggplot(data_pca, aes(PC1, PC2, color=Classification, label=name)) +
      geom_point(size=2) +
      geom_text(size=2, hjust=.5, vjust=-1.1) +
      theme_bw() +
      scale_color_manual(values=c(colors[1],colors[2])) +
      theme(text = element_text(size = 12),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "right",
            legend.justification = "center") +
      labs(title = "PCA plot highlighting the different groups",
           x = paste0("PC1: ",percentVar[1],"% variance"),
           y = paste0("PC2: ",percentVar[2],"% variance"))
ggsave(p1, file=file.path(outdir, paste0("Fig1.A.", typeof, ".pca.unfiltered.condition.jpeg")), device="jpeg", width = 10, height = 8, units = "in", dpi = 350)


p2 <- ggplot(data_pca, aes(PC1, PC2, color=Institution, label=name)) +
      geom_point(size=3, alpha=.7) +
      scale_color_manual(values = paletteColors[1:length(unique(clinical_data$Institution))]) +
      theme_bw() +
      theme(text = element_text(size = 12),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "right",
            legend.justification = "center") +
      labs(title = "PCA plot highlighting the different batch",
           x = paste0("PC1: ",percentVar[1],"% variance"),
           y = paste0("PC2: ",percentVar[2],"% variance"))
ggsave(p2, file=file.path(outdir, paste0("Fig1.B.", typeof, ".pca.unfiltered.batch.jpeg")), device="jpeg", width = 10, height = 8, units = "in", dpi = 350)
rm(p1, p2)


is_outlier <- function(x) {return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))}

libsize <- data.frame(Size = colSums(expression_mat))
libsize <- merge(libsize, clinical_data[,c(1,2,4)], by=0)
row.names(libsize) <- libsize$Row.names; libsize$Row.names <- NULL
libsize <- libsize[,c(2,1,3,4)]

dat <- libsize %>%
       tibble::rownames_to_column(var="outlier") %>%
       group_by(Classification) %>%
       mutate(is_outlier=ifelse(is_outlier(Size), Size, as.numeric(NA)))

dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
potential_outliers <- as.vector(na.omit(dat$outlier))

op <- 0
if (length(potential_outliers)!=0) {
  print("Potential outlies samples:")
  print(potential_outliers)
  
  op <- ggplot(dat, aes(y=Size, x=Classification, fill=Classification)) +
        geom_boxplot(outlier.shape = NA,  width=0.15) +
        geom_jitter(aes(color = factor(is_outlier)), height = 0.5, width = .05, size=1.8) +
        ggrepel::geom_text_repel(aes(label=outlier), na.rm=TRUE, nudge_y=0.05, nudge_x=.05, size=2.5) +
        theme_bw() +
        scale_fill_manual(values=colors) +
        theme(text = element_text(size = 12),
              legend.position = "None") +
        labs(title = "Boxplot plot based on each sample's library size",
             caption=paste0("No sample was found to be potential outliers using the\ninterquartile range criterion (1.5 IQR rule)"),
             y = "Total Read Counts",
             x = "")
  ggsave(op, file=file.path(outdir, paste0("Fig2.", typeof, ".boxplot.outliers.jpeg")), device="jpeg", width = 10, height = 8, units = "in", dpi = 350)
  write.csv(potential_outliers, file=paste(outdir,"/potential_outlier_samples.",typeof,".tsv", sep=""), row.names=FALSE, quote=FALSE)
  
} else {
  
  print(paste0("The 1.5 IQR rule could not detect any potential outlier!"))
  
  op <- ggplot(dat, aes(y=Size, x=Classification, fill=Classification)) +
        geom_boxplot(outlier.shape = NA,  width=0.15) +
        geom_jitter(color = "#7f7f7f", height = 0.5, width = .05, size=1.8) +
        theme_bw() +
        scale_fill_manual(values=colors) +
        theme(text = element_text(size = 12),
              legend.position = "None") +
        labs(title = "Boxplot plot based on each sample's library size",
             caption="The indicated samples might be potential outliers using the\ninterquartile range criterion (1.5 IQR rule)",
             y = "Total Read Counts",
             x = "")
  ggsave(op, file=file.path(outdir, paste0("Fig2.", typeof, ".boxplot.outliers.jpeg")), device="jpeg", width = 10, height = 8, units = "in", dpi = 350)
}
rm(is_outlier, dat, op, rld, percentVar, data_pca, dds, libsize, potential_outliers)





################### Filtering and CPM ###################
# Creating edgeR object with our gene expression data
edgeR_table <- DGEList(counts=expression_mat, group=clinical_data$Classification)
# Filtering out lowly expressed genes
keep <- filterByExpr(edgeR_table, group=clinical_data$Classification)
edgeR_table <- edgeR_table[keep, ,keep.lib.sizes=FALSE]
filt_data <- edgeR_table$counts
rm(keep, edgeR_table)

cpm_matrix <- cpm(filt_data)
write.table(data.frame("featureID"=rownames(cpm_matrix), cpm_matrix), file=file.path(main_outdir, paste0(typeof, "_expression_matrix.filtered.scaled.tsv")), sep="\t", row.names = F, quote=FALSE)
# plotPCAfunc(cpm_matrix, paste0(typeof,".cpm"))
rm(cpm_matrix)
