# load libraries
# R Version 3.4.1 (2017-06-30)
library(DESeq2) # DESeq2 1.16.1
library(pheatmap) # pheatmap 1.0.10

# set options
options(stringsAsFactors = F)

######################################################

# source code
analyze.tissue("2019-04-25_C2C12_RNAseq_kallisto_mapping.txt", 0.05)

######################
###### FUNCTION ######
######################
# This function takes in the C2C12 kallisto-mapped gene counts file and processes it through DESeq2 modeling to find 
# differentially expressed genes between control and MOTSc treated samples. Using DESeq2 results, the function will
# also create a heatmap plot.
# INPUT: counts.file = 2019-04-25_C2C12_RNAseq_kallisto_mapping.txt
#        FDR = false discovery rate (0.05)

analyze.tissue <- function(counts.file, FDR) {

  # read kallisto mappings
  my.data <- read.csv(counts.file, sep = "\t", header = T)

  # sum read over genes (to not have results over transcripts for DEseq2)
  my.data.per.gene <- aggregate(my.data[,5:16],by=list(my.data$tGeneSymbol),FUN=sum)

  # round counts (DESeq needs integers)
  my.data.per.gene[,2:13] <- round(my.data.per.gene[,2:13])
  rownames(my.data.per.gene) <- my.data.per.gene$Group.1

  # get the genes with no reads out
  my.null <- which(apply(my.data.per.gene[,2:13], 1, sum) <= 5) # see deseq2 vignetter
  my.filtered.matrix <- my.data.per.gene[-my.null,2:13]
  
  # get number of mapped reads
  print(apply(my.filtered.matrix,2,sum))
  
  # remove genes that start with "ERCC" because they are housekeeping genes
  remove <- grep("Ercc", row.names(my.filtered.matrix), ignore.case = F, fixed = T)
  my.filtered.matrix <- my.filtered.matrix[-remove,]
  
  # prepare matrix for DESeq2
  my.status <- rep("ALL",dim(my.filtered.matrix)[2])
  my.status[grep("CTL",colnames(my.filtered.matrix))] <- "CTL"
  my.status[grep("MOTSc",colnames(my.filtered.matrix))] <- "MOTSc"
  
  dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                            status = my.status )
    
  # get matrix using treatment as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                colData = dataDesign,
                                design = ~ status)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds) # no outliers reported
  
  res <- results(dds.deseq, contrast=c("status","MOTSc", "CTL")) # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
  
  # parse sample names
  my.sample.names <- unlist(strsplit(colnames( my.filtered.matrix ), c("_abundance.tsv")))
  
  # assign colors
  my.colors <- c(rep("orange",6), rep("dodgerblue",6))
  
  # determine normalized expression value
  tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
  
  colnames(tissue.cts) <- c(paste("CTL",c(1:6),sep="_"),paste("MOTSc",c(1:6),sep="_"))
  
  ### Heatmap Plot
  # exclude NA
  res <- res[!is.na(res$padj),]
  
  my.genes <- rownames(res)[res$padj < FDR]
  my.num <- length(my.genes)
  
  # heatmap drawing - only if there is at least one gene
  my.heatmap.out <- paste(Sys.Date(), "FDR", (FDR*100),"C2C12_Heatmap_significant_genes.pdf", sep = "_")
  
  pdf(my.heatmap.out, width = 8, height = 5, onefile = F)
  my.heatmap.title <- paste("Significant ", "(FDR <", (FDR*100), "%) ", my.num, " genes", sep = "")
  pheatmap(tissue.cts[my.genes,],
           cluster_cols = F,
           cluster_rows = T,
           colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
           show_rownames = F, scale="row",
           main = my.heatmap.title, cellwidth = 25)
  dev.off()
  
  # output result tables to files
  my.outprefix <- paste(Sys.Date(), "C2C12_RNAseq", sep = "_")
  my.out.ct.mat <- paste(my.outprefix,"log2_counts_matrix.txt", sep = "_")
  
  my.out.stats <- paste(my.outprefix,"all_genes_statistics.txt", sep = "_")
  my.out.fdr <- paste(my.outprefix,"FDR", (FDR*100), "genes_statistics.txt", sep = "_")
  my.out.rdata <- paste(my.outprefix,"statistics.RData", sep = "_")
  
  write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
  write.table(res[my.genes,], file = my.out.fdr, sep = "\t" , row.names = T, quote=F)
  save(res,file=my.out.rdata)

}
