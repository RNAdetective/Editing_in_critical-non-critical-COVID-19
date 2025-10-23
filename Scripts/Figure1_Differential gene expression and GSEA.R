# This script performs differential gene and transcript expression, plots DEGs and ADAR isoforms and performs GSEA

### Deseq2 - gene 

#load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(EnhancedVolcano)
library(ggplotify)
library(plotly)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)  # Human
library(ReactomePA)
library(enrichplot)
library(ggpubr)

# Set the working directory
directory <- "G:\\Manuscript_2_rerunning_new_SG\\DESeq2"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Aiswarya_DESeq2__"

# Input count data
data  <- as.matrix(read.csv("gene_count_matrix.csv", header=T, row.names=1, sep=","))
head(data)

# Create experiment labels (two conditions - Non_critical and Critical)
colData <- read.csv("colData.csv", header=T, row.names=1, sep=",")
head(colData)

# Making sure that the row names in col data match the column names in count data 
all(colnames(data) %in% rownames(colData))

# Making sure that the row names in col data are in the same order as the count data 
all(colnames(data) == rownames(colData))

# Create DESeq input matrix  
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData = colData, 
                              design = ~condition)
dds

# Keeping only rows with atleast 10 reads in 90% of samples
keep <- rowSums(counts(dds) >= 10) >= 62
dds_filtered <- dds[keep,]
vsd <- vst(dds_filtered, blind = TRUE)
vsd_matrix <- assay(vsd)

# Plot PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) # round the percentages

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("Critical" = "tomato", "Non_critical" = "darkblue"))
ggsave("PCA_deseq2_gene.png", width = 10, height = 8, dpi = 300)

# Run DESeq 
dds <- DESeq(dds_filtered)
res <- results(dds)
res
res <- results(dds, contrast = c("condition", "Critical","Non_critical"))  # Critical vs non critical 

# Output without filters 
write.csv(res, file = paste0(outputPrefix, "-deseq2gene_results-beforeFilters.csv"))

# Subset the results based on Padj abd log2FC
ressubset= subset(res,abs(res$log2FoldChange) >= 0.58 & padj<0.05)
ressubset <- ressubset[order(ressubset$padj),] 

# Save data results and normalized reads to csv
resdata <- merge(as.data.frame(ressubset), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-deseq2gene_results-with-padj0.05-FC1.5.csv"))

# Upregulated in Critical (Log2Fc > 0.58)
upregulated <- subset(res, log2FoldChange >= 0.58 & padj < 0.05)
dim(upregulated)
# Downregulated in critical (Log2fc < -0.58)
downregulated <- subset(res, log2FoldChange <= -0.58 & padj < 0.05)
dim(downregulated)

write.csv(upregulated, paste0(outputPrefix, "-upregulated.csv"))
write.csv(downregulated, paste0(outputPrefix, "-downregulated.csv"))


#=================== Annotate =====================================================
# Using ensembl as the dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("entrezgene_id", "external_gene_name", "chromosome_name", 
                "start_position", "end_position", "strand", "gene_biotype", 
                "description")

annotations <- getBM(attributes = attributes, mart = ensembl)

res_df <- read.csv("G:\\Manuscript_2_rerunning_new_SG\\DESeq2\\Aiswarya_DESeq2__-deseq2gene_results-with-padj0.05-FC1.5.csv")
res_df$Gene_ID <- sub(".*\\|", "", res_df$gene)

# Merge DESeq2 results with annotations and write
annotated_res <- merge(res_df, annotations, by.x = "Gene_ID", by.y = "external_gene_name", all.x = TRUE)
write.csv(annotated_res, "Deseq2_filtered_annotated.csv")

#===================== Gene set enrichment analysis ============================

# Ranking genes by log2FC
ranked_genes <- res$log2FoldChange
head(ranked_genes)
names(ranked_genes) <- rownames(res)

# Spliting row names ie gene names into ENSEMBL and SYMBOL for annotation later
gene_info <- do.call(rbind, strsplit(rownames(res), "\\|"))
colnames(gene_info) <- c("ENSEMBL", "SYMBOL")

# Assigning SYMBOLs as names to log2FC values
ranked_genes <- res$log2FoldChange
names(ranked_genes) <- gene_info[, "SYMBOL"]

# Removing NAs and duplicates
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# Sorting in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Maping SYMBOL to ENTREZ - using bitr function
gene_df <- bitr(names(ranked_genes), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db) # note: 27.34% of gene IDs are fail to map...

# Matching fold changes to ENTREZ IDs
ranked_genes_entrez <- ranked_genes[gene_df$SYMBOL]
names(ranked_genes_entrez) <- gene_df$ENTREZID

# Cleanup
ranked_genes_entrez <- ranked_genes_entrez[!is.na(names(ranked_genes_entrez))]
ranked_genes_entrez <- sort(ranked_genes_entrez, decreasing = TRUE)

# Reactome - GSEA
reactome_gsea <- gsePathway(geneList = ranked_genes_entrez,
                            organism = "human",
                            exponent = 1, # this is default, this shows the contribution of most highly ranked genes in our gene set to the overall enrichment score
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            verbose = FALSE)
# Dot plot                                                        
?dotplot
# dot plot ordered (pathways displayed from top to bottom by NES)
dotplot(
  reactome_gsea,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 20,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "Reactome Pathways",
  #orderBy = "NES",
  label_format = 30
) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16))
ggsave("Reactome_dotplot_deseq2OP.png", width = 10, height = 11, dpi = 300)



dotplot(
  reactome_gsea,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 20,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "Reactome Pathways",
  label_format = 40
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, angle = 90, hjust = 1),  
    axis.title.x = element_text(size = 16)
  )

ggsave("Reactome_dotplot_deseq2OP_horizontal_tilted_labels.png", width = 18, height = 8, dpi = 300)


# Save as dataframe and reorder according to NES values
reactome_gsea_df <- as.data.frame(reactome_gsea)
ordered_reactome_gsea <- reactome_gsea_df %>%
  arrange(desc(NES))
write.csv(ordered_reactome_gsea, "DESEq2_reactome_enrichment_ordered_acc_NES_decending.csv", row.names = FALSE)


## Upregulated genes 

# Innate immune response pathway
# GSEA neutrophil degranulation
which(reactome_gsea$Description == "Neutrophil degranulation")
gseaplot2(reactome_gsea, 
          geneSetID = 1, 
          title = "Neutrophil degranulation (NES = 2.7250)",
          pvalue_table = TRUE)
ggsave("Reactome_neutrophil_degranulation_GSEA.png", width = 10, height = 8, dpi = 300)

# Interferon alpha/beta signaling
which(reactome_gsea$Description == "Interferon alpha/beta signaling")
gseaplot2(reactome_gsea, 
               geneSetID = 21, 
               title = "Interferon alpha/beta signaling (NES = 2.3662)",
               pvalue_table = TRUE)
ggsave("Interferon alpha_beta signaling_GSEA.png", width = 10, height = 8, dpi = 300)

# Antimicrobial peptides
which(reactome_gsea$Description == "Antimicrobial peptides")
gseaplot2(reactome_gsea, 
          geneSetID = 17, 
          title = "Antimicrobial peptides (NES = 2.4625)",
          pvalue_table = TRUE)
ggsave("Antimicrobial peptides.png", width = 10, height = 8, dpi = 300)

# Hemostasis pathway
# Response to elevated platelet cytosolic Ca2+
which(reactome_gsea$Description == "Response to elevated platelet cytosolic Ca2+")
gseaplot2(reactome_gsea, 
          geneSetID = 36, 
          title = "Response to elevated platelet cytosolic Ca2+ (NES = 2.0445)",
          pvalue_table = TRUE)
ggsave("Response to elevated platelet cytosolic Ca2+.png", width = 10, height = 8, dpi = 300)

# Platelet degranulation
which(reactome_gsea$Description == "Platelet degranulation")
gseaplot2(reactome_gsea, 
          geneSetID = 32, 
          title = "Platelet degranulation (NES = 2.0578)",
          pvalue_table = TRUE)
ggsave("Platelet degranulation.png", width = 10, height = 8, dpi = 300)


# Nervous system pathways
which(reactome_gsea$Description == "GABA receptor activation")
gseaplot2(reactome_gsea, 
          geneSetID = 87, 
          title = "GABA receptor activation (NES = 1.9145)",
          pvalue_table = TRUE)
ggsave("GABA receptor activation.png", width = 10, height = 8, dpi = 300)

## Downregulated
# Metabolism of proteins
# Formation of a pool of free 40S subunits
which(reactome_gsea$Description == "Formation of a pool of free 40S subunits")
gseaplot2(reactome_gsea, 
          geneSetID = 2, 
          title = "Formation of a pool of free 40S subunits (NES = -2.6589)",
          pvalue_table = TRUE)
ggsave("Formation of a pool of free 40S subunits.png", width = 10, height = 8, dpi = 300)

# Eukaryotic Translation Elongation
which(reactome_gsea$Description == "Eukaryotic Translation Elongation")
gseaplot2(reactome_gsea, 
          geneSetID = 3, 
          title = "Eukaryotic Translation Elongation (NES = -2.6429)",
          pvalue_table = TRUE)
ggsave("Eukaryotic Translation Elongation.png", width = 10, height = 8, dpi = 300)

# Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)
which(reactome_gsea$Description == "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)")
gseaplot2(reactome_gsea, 
          geneSetID = 10, 
          title = "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) (NES = -2.5353)",
          pvalue_table = TRUE)
ggsave("Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC).png", width = 15, height = 8, dpi = 300)

## GO-BP
gsea_result <- gseGO(geneList = ranked_genes_entrez,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

dotplot(
  gsea_result,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 35,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "GO-BP",
  #orderBy = "NES",
  label_format = 30
) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16))
ggsave("DESEq2_GO_BP_dot_plot.png", width = 10, height = 15, dpi = 300)

gsea_result <- as.data.frame(gsea_result)
ordered_gsea_result <- gsea_result %>%
  arrange(desc(NES))
write.csv(ordered_gsea_result, "DESEq2_GO_BP_enrichment_ordered_acc_NES_decending.csv", row.names = FALSE)

#================ADAR plots ====================================================

# data
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

# Filter all three ADAR genes
adar_genes <- res_df[res_df$GeneID %in% c("ENSG00000160710|ADAR","ENSG00000197381|ADARB1", "ENSG00000185736|ADARB2"), ]

# Extract gene names/ take out the ensembl_ID
adar_genes$Gene <- sub(".*\\|", "", adar_genes$GeneID)

# plot

# Adding significance information 
adar_genes$Significance <- ifelse(adar_genes$padj < 0.05, "Significant", "Not significant")


# Bar plot instead of lollipop
ggplot(adar_genes, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = Significance)) +
  geom_col(width = 0.2) +  # Bar plot
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("padj=%.4f", padj),
                vjust = ifelse(log2FoldChange < 0, 1.5, -0.5)), 
            size = 4) +
  scale_fill_manual(values = c("Significant" = "steelblue", "Not significant" = "darkred")) +
  labs(x = "", y = "log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

ggsave("ADAR_deseq2_barplot.png", width = 8, height = 8)

#========= volcano plot - all DEGS =============================================

# Data - DESeq2 results
resdata <- read.csv("G:\\Manuscript_2_rerunning_new_SG\\DESeq2\\Aiswarya_DESeq2__-deseq2gene_results-beforeFilters.csv")
resdata <- na.omit(resdata)

# Gene name is set as a column
resdata$gene <- sub(".*\\|", "", resdata$X)

# Adding a column of significance
resdata$significance <- with(resdata, ifelse(padj < 0.05 & abs(log2FoldChange) >= 0.58, "Significant", "Not significant"))

# Order by adjusted p-value to find top 20 genes
resdata_sig<- resdata[order(resdata$padj),]
top20 <- head(resdata_sig, 20)

# volcano plot
g <- ggplot(resdata_sig, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.5) +
  scale_color_manual(values = c("darkcyan", "darkblue")) +
  geom_vline(xintercept=0.58, colour='green', linetype=3) +
  geom_vline(xintercept=-0.58, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  geom_text_repel(data = top20, aes(label = gene), size = 3, max.overlaps = Inf) +
  labs(x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9),
        legend.justification = c("right", "top"))

ggsave("Deseq2_volcano_gene.png", 
       plot = g, width = 8, height = 6, dpi = 300)


#=========================  transcript counts ==================================

# Set the prefix for each output file name
outputPrefix <- "Aiswarya_DESeq2_Transcript_Counts"

# Input count data
data  <- as.matrix(read.csv("transcript_count_matrix.csv", header=T, row.names=1, sep=","))
head(data)

# Create experiment labels (two conditions - Non_critical and Critical)
colData <- read.csv("colData.csv", header=T, row.names=1, sep=",")
head(colData)

# Making sure that the row names in col data match the column names in count data 
all(colnames(data) %in% rownames(colData))

# Making sure that the row names in col data are in the same order as the count data 
all(colnames(data) == rownames(colData))

# Create DESeq input matrix  
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData = colData, 
                              design = ~condition)
dds


# Run DESeq 
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast = c("condition", "Critical","Non_critical"))


### Saving results before applying filters
# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results_before_filters.csv"))

# Applying filters
ressubset= subset(res,  abs(res$log2FoldChange) >= 0.58 & padj<0.05)
ressubset <- ressubset[order(ressubset$padj),]

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(ressubset), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-padj0.05-FC1.5.csv"))


#===================== ADAR_isoform=============================================
#data (manually pulled from the outputs of deseq2 analysis on transcripts)
df_isoform <- read.csv("ADAR_isofrom.csv", header=T, sep = ",")

# adarp110
p <- ggplot(df_isoform, aes(x = Severity, y = ENST00000368471.ADARp110, color = Severity)) +
  geom_boxplot(fill = NA, size = 1) + 
  geom_jitter(width = 0.2, size = 2, aes(fill = Severity), shape = 21) +  
  scale_color_manual(values = c("tomato", "darkblue")) +  
  scale_fill_manual(values = c("tomato", "darkblue")) +  
  labs(y = "ADARp110 expression", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),  
        legend.position = c(0.9, 0.9),  
        legend.justification = c("right", "top"))

print(p)


comparison_result1 <- 
  p + stat_compare_means(method = "wilcox", paired = FALSE, label = "p.format",
                         comparisons = list(c("Non-critical","Critical")))
ggsave("ADARp110_borderonlyfill.png", 
       plot = comparison_result1, width = 8, height = 6, dpi = 300)

# adarp150

p2 <- ggplot(df_isoform, 
             aes(x = Severity, y = ENST00000368474.ADARp150, color = Severity)) +
  geom_boxplot(fill = NA, size = 1) +  
  geom_jitter(width = 0.2, size = 2, aes(fill = Severity), shape = 21) + 
  scale_color_manual(values = c("tomato", "darkblue")) +  
  scale_fill_manual(values = c("tomato", "darkblue")) +
  labs(y = "ADARp150 expression", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),  
        legend.position = c(0.9, 0.9),  
        legend.justification = c("right", "top"))



comparison_result2 <- 
  p2 + stat_compare_means(method = "wilcox", paired = FALSE, label = "p.format",
                          comparisons = list(c("Non-critical","Critical")))


ggsave("ADARp150.png", 
       plot = comparison_result2, width = 8, height = 6, dpi = 300)

e = 1) +  
  geom_jitter(width = 0.2, size = 2, aes(fill = Severity), shape = 21) + 
  scale_color_manual(values = c("tomato", "darkblue")) +  
  scale_fill_manual(values = c("tomato", "darkblue")) +
  labs(y = "ADARp150 expression", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),  
        legend.position = c(0.9, 0.9),  
        legend.justification = c("right", "top"))



comparison_result2 <- 
  p2 + stat_compare_means(method = "wilcox", paired = FALSE, label = "p.format",
                          comparisons = list(c("Non-critical","Critical")))


ggsave("ADARp150.png", 
       plot = comparison_result2, width = 8, height = 6, dpi = 300)

