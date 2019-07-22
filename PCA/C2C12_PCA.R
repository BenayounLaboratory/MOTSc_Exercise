# This script generates a PCA plot of PC1 and PC2 of gene expression (transcriptomics) for C2C12 normalized counts

# R Version 3.4.1 (2017-06-30)

######################################################

# read DESeq2 normalized counts
my.data <- read.csv("2019-04-25_C2C12_RNAseq_log2_counts_matrix.txt", stringsAsFactors = F, sep = "\t", header = T) 

str(my.data)

######################################################

# PCA using prcomp function
data.pca <- prcomp(t(my.data), scale = T)

# get the precentage for each component
summary(data.pca)

# Importance of components%s:
#                           PC1     PC2     PC3      PC4      PC5      PC6      PC7     PC8      PC9     PC10     PC11
# Standard deviation     51.6561 47.9309 41.8080 39.20053 38.51421 37.37337 36.04521 34.9335 33.48459 32.71347 30.35743
# Proportion of Variance  0.1592  0.1371  0.1043  0.09167  0.08849  0.08332  0.07751  0.0728  0.06689  0.06384  0.05498
# Cumulative Proportion   0.1592  0.2962  0.4005  0.49217  0.58066  0.66399  0.74150  0.8143  0.88118  0.94502  1.00000
#                           PC12
# Standard deviation     2.514e-13
# Proportion of Variance 0.000e+00
# Cumulative Proportion  1.000e+00

# find the color code for orange, dodgerblue
col2rgb("orange")       # CTL 255, 165, 0
col2rgb("dodgerblue")   # MOTSc 30, 144, 255

# make PCA plot based PC1 and PC2
pdf("PCA_PC1-2_C2C12.pdf", width =11, height =8.5)
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
  plot(data.pca$x[,c(1,2)], 
       pch = 17,  
       xlab = 'PC1 (16%)', 
       ylab = 'PC2 (14%)', 
       cex.lab = 1.5,
       col = "white",cex=0.6) 
  points(data.pca$x[grep("CTL", colnames(my.data)) ,c(1,2)],
         pch=17, 
         cex=2.5, 
         col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255))
  points(data.pca$x[grep("MOTSc", colnames(my.data)) ,c(1,2)],
         pch=17,
         cex=2.5, 
         col=rgb(red=30,green=144,blue=255, alpha = 204, maxColorValue = 255))
  legend("topright",c("CTL","MOTSc"), inset=c(0,-0.19),horiz=T,
  col=c("orange","dodgerblue"), text.col=c("orange","dodgerblue"), pch=c(17,17), cex = 1 , pt.cex=3)
dev.off()