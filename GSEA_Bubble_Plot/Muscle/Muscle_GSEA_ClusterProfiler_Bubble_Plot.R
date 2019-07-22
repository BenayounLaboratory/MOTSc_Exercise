# load libraries
# R Version 3.4.1 (2017-06-30)
library('pheatmap') # pheatmap 1.0.10
library(ggplot2) # ggplot2 3.1.0
library(scales) # scales 1.0.0

# set options
options(stringsAsFactors = F)
options(digits = 5)

################################################################

# source code
get_enrich_balloons("GO_BP_FDR15_pathways", 15, "GO-BP")
get_enrich_balloons("KEGG_FDR15_pathways", 15, "KEGG")

####################
##### FUNCTION #####
####################
# This function makes bubble plots from Muscle GSEA statistics (FDR 15)
# INPUT: my.data.name = output file name describing the gene set name and FDR
#        FDR = false discovery rate (0.15)
#        gene.set = name of gene set collection (GO-BP or KEGG)
#        my.threshold = if threshold is 0, the bubble plot will display pathways
#                 that are present in at least 0 or more columns
#        my.colnames = names of columns for the bubble plot

get_enrich_balloons <- function(my.data.name, FDR, gene.set, my.threshold = 0,
                                  my.colnames = c("Muscle")) {
  
  if (FDR == 5) {
    my.enrich.sets <- list.files(".", pattern = "*FDR5.txt")
  }
  
  if (FDR == 10) {
    my.enrich.sets <- list.files(".", pattern = "*FDR10.txt")
  }
  
  if (FDR == 15) {
    my.enrich.sets <- list.files(".", pattern = "*FDR15.txt")
  }
  
  # reorder files based on colnames
  my.columns <<- c()

  for (i in 1:length(my.colnames)) {
     my.columns <<- c(my.columns,
                     grep(paste("_", my.colnames[i], "_", sep = ""), my.enrich.sets))
  }
  
  # get data from significant FDR
  my.tissues <- vector(length=length(my.colnames), mode="list")
  names(my.tissues) <- my.colnames
  my.pathways <- c()
  
  for ( i in 1:length(my.columns)) {
    my.file <- my.enrich.sets[my.columns[i]]
    print(my.file);
    my.tissues[[i]] <- read.csv(paste("./", my.file, sep=""), sep="\t", header=T) 
    row.names(my.tissues[[i]]) <- my.tissues[[i]]$Description
    my.pathways <- unique(c(my.pathways, rownames(my.tissues[[i]])))
    
  }
  print(my.pathways)
  
  ## prepapre output data
  # p-val matrix
  my.matrix <- matrix(0,length(my.pathways),length(my.colnames)) # default: -log10(1) pval == 0 no enrichment
  
  # Enrichment matrix
  my.matrix2 <- matrix(0,length(my.pathways),length(my.colnames)) # initialize with Enrichment = 0 if no enrich
  
  # matrix with record of significance
  my.matrix3 <- matrix(0,length(my.pathways),length(my.colnames)) # to get sigificant pathways
  
  colnames(my.matrix)  <- my.colnames
  colnames(my.matrix2) <- my.colnames
  rownames(my.matrix)  <- my.pathways
  rownames(my.matrix2) <- my.pathways
  colnames(my.matrix3) <- my.colnames
  rownames(my.matrix3) <- my.pathways
  
  # collect data from files
  for (i in 1:length(my.pathways)) {
    #print(my.pathways[i])
    
    for (j in 1:length(my.colnames)) { # tissues 
      
      my.id <- which(rownames(my.tissues[[j]]) %in% my.pathways[i])
      #print(paste(j, length(my.id)))
      if(length(my.id) == 1) { # if was significant in this tissue (and not on both tail ends, which would be 2)
        #print('IN')
        
        my.matrix[i,j] <- -log10(my.tissues[[j]]$p.adjust[my.id]+1e-10) # log(0) is undefined
        
        my.matrix2[i,j] <- my.tissues[[j]]$NES[my.id]
        
        my.matrix3[i,j] <- 1
        
      }
    }
  }
  
  # find pathways significant in threshold (minimum) or more tissues (1 potential from our experiments)
  my.sigs <- apply(my.matrix2,1,sum) > my.threshold
  my.res.enrich <- data.frame(my.matrix2[my.sigs,])
  my.pval.enrich <- data.frame(my.matrix[my.sigs,])
  
  # sort by average change
  my.average <- apply(my.res.enrich,1,mean)
  my.sorted <- sort(my.average,index.return=T,decreasing=T)
  my.res.enrich2 <- my.res.enrich[my.sorted$ix,]
  my.res.enrich2 <- as.data.frame(my.res.enrich2)
  
  my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])

  # format for ggplot
  my.res.enrich2$Pathnames <- names(my.sorted$x)
  my.res.enrich3 <- cbind(my.res.enrich2[,c('Pathnames')],my.res.enrich2[,1],my.pval.enrich2[,1])
  colnames(my.res.enrich3) <- c('Pathnames','signed_enrichment','minusLog10Pval')
  my.res.enrich3 <- as.data.frame(my.res.enrich3)
  my.res.enrich3$signed_enrichment <- as.double(my.res.enrich3$signed_enrichment)
  my.res.enrich3$minusLog10Pval <- as.double(my.res.enrich3$minusLog10Pval)
  
  my.txtname <- paste(Sys.Date(),"Muscle_Enrichment_table",my.data.name,"pathways_significant_in",my.threshold,"or_more.txt", sep="_")
  write.table(my.res.enrich3,file=my.txtname,sep="\t",quote=F)
  
  my.max <<- max(my.res.enrich3$signed_enrichment)
  my.min <<- min(my.res.enrich3$signed_enrichment)

  my.values <- c(0, my.min + (my.max-my.min) /8 * 0, my.min + (my.max-my.min) /8 * 1, my.min + (my.max-my.min) /8 * 2, my.min + (my.max-my.min) /8 * 3,
                                 my.min + (my.max-my.min) /8 * 4, my.min + (my.max-my.min) /8 * 5, my.min + (my.max-my.min) /8 * 6, my.min + (my.max-my.min) /8 * 7,
                                 my.min + (my.max-my.min) /8 * 8)

  my.scaled <<- rescale(my.values, to = c(0,1))
  my.color.vector <- colorRampPalette(c("blanchedalmond", "firebrick1", "firebrick2", "firebrick3", "firebrick4"))(9)
  
  # to preserve the wanted order
  my.res.enrich3$Pathnames <- factor(my.res.enrich3$Pathnames, levels = rev(unique(my.res.enrich3$Pathnames)))
  
  my.pdfname <- paste(Sys.Date(),"Muscle_Enrichment_BALLOON_plot",my.data.name,"pathways_significant_in",my.threshold,"or_more.pdf", sep="_")
  
  if (gene.set == "GO-BP") {
    pdf(my.pdfname, onefile=T, height = 100, width=15)
  }

  if (gene.set == "KEGG") {
    pdf(my.pdfname, onefile=T, height = 25, width=15)
  }
  
  my.plot <- ggplot(my.res.enrich3,aes(x=my.colnames,y=Pathnames,colour=signed_enrichment,size=minusLog10Pval))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- ggplot(my.res.enrich3,aes(x=my.colnames,y=Pathnames,colour=signed_enrichment,size=minusLog10Pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("GSEA") + labs(x = "Tissue/condition", y = "Pathways")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
  print(my.plot)
  dev.off()  
  
}
