# R script
# Running DESeq2 for septic male vs female
# Completed version
# MICB 405 101 (2025W1)
# October 23, 2025


# Loading required packages -----------------------------------------------

library(tidyverse)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(clusterProfiler)
library(EnhancedVolcano)

#Analysis of septic male vs female-----------------------------------------
#Path to files
files <- list.files(pattern = "40ReadsPerGene.out.tab$", full.names = TRUE)

# Read and merge counts
read_STAR <- function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  df <- df[, c(1, 4)] # adjust column (2=unstranded, 3 or 4 for stranded)
  colnames(df) <- c("gene_id", basename(file))
  return(df)
}

count_list <- lapply(files, read_STAR)
counts <- Reduce(function(x, y) merge(x, y, by="gene_id"), count_list)
counts <- counts[!counts$gene_id %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped"), ]

head(counts)


# DESeq2 is a widely used tool to assess truly differentially expressed genes.
# Let's prepare the read counts for use in DESeq2.

# Build a new data frame with gene names and read counts.

# If you are using all stranded data, pick the appropriate column.
# If using non-stranded data, use total column.

counts_deseq <- data.frame(row.names = counts$gene_id,
                                f1septic = counts$`f1-40ReadsPerGene.out.tab`,
                                f2septic = counts$`f2-40ReadsPerGene.out.tab`,
                                f3septic = counts$`f3-40ReadsPerGene.out.tab`,
                                m1septic = counts$`m1-40ReadsPerGene.out.tab`,
                           m2septic = counts$`m2-40ReadsPerGene.out.tab`,
                           m3septic = counts$`m3-40ReadsPerGene.out.tab`)

# DESeq2 also requires the read counts to be in matrices, not data frames.
# Matrices also contain data values in a table format, but the entire table's
# values must all be of the same data type.

# We can convert with as.matrix().
counts_matrix <- as.matrix(counts_deseq)

# DESeq2 also requires a table that provides information about the columns.
# These must be in the order in which they appear on the count matrix.
colnames(counts_matrix)

# Confusingly, column metadata is presented as a data frame. The names of each
# sample should be the row name, and be identical and in order of appearance in
# the read counts matrix.
columns_data <- data.frame(gender = c("female",
                                         "female",
                                         "female",
                                         "male","male","male"),
                           row.names = c("f1septic",
                                         "f2septic",
                                         "f3septic",
                                         "m1septic","m2septic","m3septic"))

# Create our DESeq2 object.
dds_matrix <- DESeqDataSetFromMatrix(countData = counts_matrix, # Matrix 
                                     colData = columns_data, # Metadata
                                     design = ~gender)

# Set control condition using the relevel function, to make sure this is the control
dds_matrix$gender <- relevel(dds_matrix$gender, ref = "female")

# Run DESeq2
dds <- DESeq(dds_matrix)

# With this dds object, we can perform a wide variety of interesting analyses.


# Analyses covered in Tutorial 7 ------------------------------------------

# Sanity Check 1: Make a PCA plot with the data ---------------------------

# Perform log transformation on our count data
rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
pcaData <- plotPCA(rld, intgroup = c("gender"),returnData=TRUE) 
plotPCA(rld, intgroup = c("gender"))
# Get percent variance for labeling axes
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot with ggplot2, adding sample labels
PCAplot <- pcaData %>% 
  ggplot(aes(x = PC1, y = PC2, colour = gender)) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-30,30)) +
  scale_x_continuous(limits = c(-50,30)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(linewidth = 2)) +
  labs(colour = "Gender")

# Further analyses of differentially expressed genes ----------------------

# names of the results that DESeq2 calculated
resultsNames(dds)

# Now we will extract the results for our comparison between male and female septic patients
res <- results(dds, name = "gender_male_vs_female") |> as.data.frame()
# we save it as a data frame so we can apply tidyverse dplyr functions better!
head(res)

# Getting rid of NA values
res_no_NA <- res |> drop_na()
head(res_no_NA)
dim(res_no_NA) # Look at how many rows you filtered out!

# Filtering for adjusted p-values <=0.05
res_filtered <- res_no_NA |> filter(padj <= 0.05)
head(res_filtered)
dim(res_filtered) # Look at how many rows you filtered out!

# Pull genes with more than 2x higher/lower expression
res_filtered_final <- res_filtered |>
  filter(log2FoldChange <= -1 | log2FoldChange >= 1)
# The '|' stands for OR here!
head(res_filtered_final)
dim(res_filtered_final)

# Top 10 genes positively differentially expressed.
top10_genes <- res_filtered_final |>
  arrange(desc(log2FoldChange)) |>
  # We use the desc() function to organize the column in descending order.
  head(n = 10)
top10_genes

# How about the 10 genes negatively differentially expressed?
bot10_genes <- res_filtered_final |>
  arrange(log2FoldChange) |>
  # Since we don't use desc(), the column is organized in ascending order.
  head(n = 10)
bot10_genes
# Write CSV
write.csv(res, file = "results_no_filtered.csv")

## Make volcano plot by whether it is significant + large change
vol_plot <- res_filtered_final %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#Analysis for interaction----------------------------------------------------------------
#Path to files
files_interaction <- list.files(path="C:/Users/forss/Documents/R Files/MICB 405 - Tien/MICB 405",pattern = "ReadsPerGene.out.tab$", full.names = TRUE)

# Read and merge counts
read_STAR_int <- function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  df <- df[, c(1, 4)] # adjust column (2=unstranded, 3 or 4? for stranded)
  colnames(df) <- c("gene_id", basename(file))
  return(df)
}

count_list_int <- lapply(files_interaction, read_STAR_int)
counts_int <- Reduce(function(x, y) merge(x, y, by="gene_id"), count_list_int)
counts_int <- counts_int[!counts_int$gene_id %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped"), ]

head(counts_int)


# DESeq2 is a widely used tool to assess truly differentially expressed genes.
# Let's prepare the read counts for use in DESeq2.

# Build a new data frame with gene names and read counts.

# If you are using all stranded data, pick the appropriate column.
# If using non-stranded data, use total column.

counts_deseq_int <- data.frame(row.names = counts_int$gene_id,
                           f1septic = counts_int$`f1-40ReadsPerGene.out.tab`,
                           f2septic = counts_int$`f2-40ReadsPerGene.out.tab`,
                           f3septic = counts_int$`f3-40ReadsPerGene.out.tab`,
                           m1septic = counts_int$`m1-40ReadsPerGene.out.tab`,
                           m2septic = counts_int$`m2-40ReadsPerGene.out.tab`,
                           m3septic = counts_int$`m3-40ReadsPerGene.out.tab`,
                           f1healthy = counts_int$f1_healthy_ReadsPerGene.out.tab,
                           f2healthy=counts_int$f2_healthy_ReadsPerGene.out.tab,
                           f3healthy=counts_int$f3_healthy_ReadsPerGene.out.tab,
                           m1healthy=counts_int$m1_healthy_ReadsPerGene.out.tab,
                           m2healthy=counts_int$m2_healthy_ReadsPerGene.out.tab,
                           m3healthy=counts_int$m3_healthy_ReadsPerGene.out.tab
                      )

# DESeq2 also requires the read counts to be in matrices, not data frames.
# Matrices also contain data values in a table format, but the entire table's
# values must all be of the same data type.

# We can convert with as.matrix().
counts_matrix_int <- as.matrix(counts_deseq_int)

# DESeq2 also requires a table that provides information about the columns.
# These must be in the order in which they appear on the count matrix.
colnames(counts_matrix_int)

# Confusingly, column metadata is presented as a data frame. The names of each
# sample should be the row name, and be identical and in order of appearance in
# the read counts matrix.
columns_data_int <- data.frame(gender = c("female",
                                      "female",
                                      "female",
                                      "male","male","male","female",
                                      "female",
                                      "female",
                                      "male","male","male"), 
                               condition = c("septic","septic","septic","septic","septic","septic",
                                             "healthy","healthy","healthy","healthy","healthy","healthy"),
                           row.names = c("f1septic",
                                         "f2septic",
                                         "f3septic",
                                         "m1septic","m2septic","m3septic",
                                         "f1healthy",
                                         "f2healthy",
                                         "f3healthy",
                                         "m1healthy","m2healthy","m3healthy"))

# Create our DESeq2 object.

dds_matrix_int <- DESeqDataSetFromMatrix(countData = counts_matrix_int, 
                              colData = columns_data_int, 
                              design = ~ gender + condition + gender:condition)
dds_matrix_int$gender <- relevel(dds_matrix_int$gender, ref = "female")

# Run DESeq2
dds_int <- DESeq(dds_matrix_int)

# With this dds object, we can perform a wide variety of interesting analyses.


# Analyses covered in Tutorial 7 ------------------------------------------

# Sanity Check 1: Make a PCA plot with the data ---------------------------

# Perform log transformation on our count data
rld_int <- rlog(dds_int)
vsd_int <- vst(dds_matrix_int, blind=FALSE)

# Generate a PCA plot with DESeq2's plotPCA function
pcaData_int <- plotPCA(rld_int, intgroup = c("gender","condition"),returnData=TRUE) 
plotPCA(rld_int, intgroup = c("gender","condition"))
res_interation <- results(dds_int, name = "gendermale.conditionseptic") #heatmaps

# Generate a PCA plot with DESeq2's plotPCA function

# Get percent variance for labeling axes
percentVar_int <- round(100 * attr(pcaData_int, "percentVar"))

# Plot with ggplot2, adding sample labels
pcaplot_int <- pcaData_int %>% 

    ggplot(aes(x = PC1, y = PC2, colour = gender, shape=condition)) +
  geom_point(size = 4) +
  scale_y_continuous(limits = c(-20,40)) +
  scale_x_continuous(limits = c(-20,40)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(face = "bold",size=18),
    legend.title = element_text(face = "bold",size=14),
    legend.text = element_text(size = 14),
    panel.border = element_rect(linewidth = 2)) +
  labs(colour = "Gender",shape="Condition")
#heatmap
topGenes <- head(order(res_interation$padj), 30)
mat <- assay(vsd_int)[topGenes, ]
heatmap <-pheatmap(mat, annotation_col = columns_data_int)
ggsave("heatmap_septic.png", 
       plot=heatmap, dpi=300, width=10, height=7
)
#volcano
volcano <- EnhancedVolcano(res_interation,
                          lab = rownames(res_interation),
                          title = NULL,
                          subtitle = NULL,
                          legendLabSize = 12,
                          legendIconSize=3,
                          x = 'log2FoldChange',
                          y = 'padj',
                                 pCutoff = 0.05,
                                 FCcutoff = 1)
ggsave("volcano_plot.png",
       plot = volcano,
       width = 10,        # width in inches
       height = 8,        # height in inches
       dpi = 300)
combined <- pcaplot_int + volcano+ plot_layout(ncol = 2)
combined
ggsave("pca_volcano_septic.png", plot=combined,width=12.5, height=7,dpi=300)

sig_int <- subset(res_interation, padj < 0.05)
sig_int$symbol <- mapIds(org.Hs.eg.db,
                         keys = rownames(sig_int),
                         column = "SYMBOL",
                         keytype = "SYMBOL",
                         multiVals = "first")
library(clusterProfiler)
geneList <- res_interation$log2FoldChange
names(geneList) <- rownames(res_interation)
geneList <- sort(geneList, decreasing = TRUE)

gsea_go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 keyType      = "SYMBOL",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05)
mito <- read_excel("C:/Users/forss/Documents/R Files/MICB 405 - Tien/MICB 405/Human.MitoCarta.3.0.xls")
mito_genes <- mito$Symbol

