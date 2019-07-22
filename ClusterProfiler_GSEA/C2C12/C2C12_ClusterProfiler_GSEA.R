# This script perform GSEA at FDR 5 and 15 using the clusterProfiler package for C2C12 DESeq2 all gene statistics

# load libraries
# R Version 3.5.0 (2018-04-23)
library(clusterProfiler) # clusterProfiler 3.10.1
library(org.Mm.eg.db) # org.Mm.eg.db 3.7.0

# set options
options(stringsAsFactors = F)

######################################################

# load data from DESeq2
res.C2C12 <- read.csv("2019-04-25_C2C12_RNAseq_all_genes_statistics.txt", sep = "\t", header = T)

entrezID.background  <- bitr(rownames(res.C2C12), 
                            fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
# 8536 mapping correcly to annotation for C2C12

######################################################

# prepare GeneLists using DEseq2 t-statistic to rank genes
res.C2C12 <- as.data.frame(res.C2C12)
res.C2C12$symbol <- rownames(res.C2C12)
res.C2C12 <- merge(res.C2C12, entrezID.background, by.x = "symbol", by.y = "SYMBOL")
C2C12.geneList = res.C2C12$stat
names(C2C12.geneList) = res.C2C12$ENTREZID
C2C12.geneList = sort(C2C12.geneList, decreasing = TRUE)

## GO Gene Set Enrichment Analysis
go.bp.gsea <- gseGO(geneList     = C2C12.geneList,
                    OrgDb        = org.Mm.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    nPerm        = 1000,
                    minGSSize    = 100,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.15, # change to 0.05 for FDR 5
                    verbose      = FALSE)

# write results to file
write.table(go.bp.gsea@result, file = paste(Sys.Date(),"C2C12_GOBP_GSEA_Analysis_FDR15.txt", sep = "_"), quote = F, sep = "\t")
head(go.bp.gsea)

# GSEA plot of interesting pathways
pdf(paste(Sys.Date(),"C2C12_GO_BP_FDR15_Selected_GSEA_plots.pdf", sep = "_"))
gseaplot(go.bp.gsea, geneSetID = "GO:0045333", title = "GO:0045333 cellular respiration")
gseaplot(go.bp.gsea, geneSetID = "GO:0030198", title = "GO:0030198 extracellular matrix organization")
gseaplot(go.bp.gsea, geneSetID = "GO:0070848", title = "GO:0070848 response to growth factor")
gseaplot(go.bp.gsea, geneSetID = "GO:0040017", title = "GO:0040017 positive regulation of locomotion")
gseaplot(go.bp.gsea, geneSetID = "GO:0007005", title = "GO:0007005 mitochondrion organization")
gseaplot(go.bp.gsea, geneSetID = "GO:0045087", title = "GO:0045087 innate immune response")
gseaplot(go.bp.gsea, geneSetID = "GO:0008610", title = "GO:0008610 lipid biosynthetic process")
gseaplot(go.bp.gsea, geneSetID = "GO:0006954", title = "GO:0006954 inflammatory response")
gseaplot(go.bp.gsea, geneSetID = "GO:0016236", title = "GO:0016236 macroautophagy")
gseaplot(go.bp.gsea, geneSetID = "GO:0006888", title = "GO:0006888 ER to Golgi vesicle-mediated transport")
gseaplot(go.bp.gsea, geneSetID = "GO:0071804", title = "GO:0071804 cellular potassium ion transport")
gseaplot(go.bp.gsea, geneSetID = "GO:0019827", title = "GO:0019827 stem cell population maintenance")
gseaplot(go.bp.gsea, geneSetID = "GO:0048284", title = "GO:0048284 organelle fusion")
dev.off()

## KEGG Gene Set Enrichment Analysis 
kegg.gsea <- gseKEGG(geneList     = C2C12.geneList,
                     organism     = 'mmu',
                     keyType      = 'ncbi-geneid',
                     nPerm        = 1000,
                     pvalueCutoff = 0.15, # change to 0.05 for FDR 5
                     verbose      = FALSE)

# write results to file
write.table(kegg.gsea@result, file = paste(Sys.Date(),"C2C12_KEGG_GSEA_Analysis_FDR15.txt", sep = "_"), quote = F, sep = "\t")
head(kegg.gsea)

# GSEA plot of interesting pathways
pdf(paste(Sys.Date(),"C2C12_KEGG_FDR15_Selected_GSEA_plots.pdf", sep = "_"))
gseaplot(kegg.gsea, geneSetID = "mmu04218", title = "mmu04218	Cellular senescence")
gseaplot(kegg.gsea, geneSetID = "mmu04152", title = "mmu04152	AMPK signaling pathway")
gseaplot(kegg.gsea, geneSetID = "mmu01210", title = "mmu01210	2-Oxocarboxylic acid metabolism")
gseaplot(kegg.gsea, geneSetID = "mmu00020", title = "mmu00020	Citrate cycle (TCA cycle)")
gseaplot(kegg.gsea, geneSetID = "mmu04211", title = "mmu04211	Longevity regulating pathway")
gseaplot(kegg.gsea, geneSetID = "mmu00190", title = "mmu00190	Oxidative phosphorylation")
dev.off()