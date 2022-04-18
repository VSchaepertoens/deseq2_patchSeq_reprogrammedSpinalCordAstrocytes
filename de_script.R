## Gene-level differential expression analysis using DESeq2

## Setup
### Bioconductor and CRAN libraries used
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(here)
library(pheatmap)
library(ggrepel)
library(reshape2)
library(dplyr) 
library(RColorBrewer)

### Load gene counts, meta data, and annotations------------------------------------------------------------------
load(file = "data/counts_spinalCord.Rdata")

# Load the annotation table for GRCm38 (Genome Reference Consortium Mouse Build 38)
tx2gene <- read.delim(here("preprocessing/MouseIndexes/geneInfo.tab"),header=FALSE)

tx2gene <- tx2gene[-1:-5, ]
# Take a look at it 
tx2gene %>% View()

#check countMatrix (GeneCounts) columnnames correspond with rownames of meta file(meta)
all(rownames(meta) == colnames(GeneCounts))

## Create DESeq2Dataset object --------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = GeneCounts,
                              colData = meta,
                              design = ~ tf)
dds
View(counts(dds))

## Count normalization median of ratios method ----------------------------------------------------
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

## Sample-level quality control -------------------------------------------------------------------

### Transform counts for data visualization

# Regularized logarithm
# FALSE argument results in transformation biased to sample condition information
#rld <- rlog(dds, blind=TRUE) # commented out, takes too long

# Variance stabilizing transformations
#FALSE argument results in transformation biased to sample condition information
vsd <- vst(dds, blind=TRUE)


### PCA
### Plot PCA 
plotPCA(vsd, intgroup="tf")
plotPCA(vsd, intgroup="mapping")

# Plot 3rd and 4th principal components
# Input is a matrix of log transformed values
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = tf))


### Plot hierarchical clustering heatmap
### Extract the vsd matrix from the object
vsd_mat <- assay(vsd)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq

### Compute pairwise correlation values
vsd_cor <- cor(vsd_mat)    ## cor() is a base R function

head(vsd_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
df <- as.data.frame(colData(dds)[,c("tf")])
colnames(df) <- c("Transcription factor")
pheatmap(vsd_cor,show_rownames=FALSE,annotation=df)

## Run analysis-------------------------------------------------------------------------------------
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT",reduced = ~ 1, useT=TRUE,minmu=1e-6,minReplicatesForReplace=Inf)

## Check the size factors
sizeFactors(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## Define contrasts, extract results table, and shrink the log2 fold changes
## Cntrl vs Ngn2------------------------------------------------------------------------------------

contrast_ngn2Cntrl <- c("tf","Ngn2","Cntrl")

res_table_ngn2Cntrl <- results(dds, contrast=contrast_ngn2Cntrl, alpha = 0.05)

class(res_table_ngn2Cntrl)

mcols(res_table_ngn2Cntrl, use.names=T)

res_table_ngn2Cntrl %>% data.frame() %>% View()

## Save the unshrunken results to compare
res_table_ngn2Cntrl_unshrunken <- res_table_ngn2Cntrl

# Apply fold change shrinkage
res_table_ngn2Cntrl <- lfcShrink(dds, contrast=contrast_ngn2Cntrl, res=res_table_ngn2Cntrl, type='normal')
res_table_ngn2Cntrl

#MA plot The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes tested. The genes that are significantly DE are colored to be easily identified. 
#plotMA(res_table_ngn2Cntrl_unshrunken)
plotMA(res_table_ngn2Cntrl)

## Summarize results
summary(res_table_ngn2Cntrl, alpha = 0.01)

### Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

res_table_ngn2Cntrl_tb <- res_table_ngn2Cntrl %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#sig_ngn2Cntrl <- res_table_ngn2Cntrl_tb %>%
#   filter(padj < padj.cutoff)

sig_ngn2Cntrl <- res_table_ngn2Cntrl_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_ngn2Cntrl

## Ascl1 vs Cntrl ----------------------------------------------------------------------------------

contrast_ascl1Cntrl <- c("tf","Ascl1","Cntrl")

res_table_ascl1Cntrl <- results(dds, contrast=contrast_ascl1Cntrl, alpha = 0.05)

class(res_table_ascl1Cntrl)

mcols(res_table_ascl1Cntrl, use.names=T)

res_table_ascl1Cntrl %>% data.frame() %>% View()

## Save the unshrunken results to compare
res_table_ascl1Cntrl_unshrunken <- res_table_ascl1Cntrl

# Apply fold change shrinkage
res_table_ascl1Cntrl <- lfcShrink(dds, contrast=contrast_ascl1Cntrl, res=res_table_ascl1Cntrl, type='normal')
res_table_ascl1Cntrl

#MA plot The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes tested. The genes that are significantly DE are colored to be easily identified. 
#plotMA(res_table_ascl1Cntrl_unshrunken)
plotMA(res_table_ascl1Cntrl)

## Summarize results
summary(res_table_ascl1Cntrl, alpha = 0.01)

### Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

res_table_ascl1Cntrl_tb <- res_table_ascl1Cntrl %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#sig_ascl1Cntrl <- res_table_ascl1Cntrl_tb %>%
#   filter(padj < padj.cutoff)

sig_ascl1Cntrl <- res_table_ascl1Cntrl_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_ascl1Cntrl

## Ngn2 vs Ascl1-------------------------------------------------------------------------------------

contrast_tf <- c("tf","Ngn2","Ascl1")

res_table_tf <- results(dds, contrast=contrast_tf, alpha = 0.01, lfcThreshold = 0.58)

class(res_table_tf)

mcols(res_table_tf, use.names=T)

res_table_tf %>% data.frame() %>% View()

## Save the unshrunken results to compare
res_table_tf_unshrunken <- res_table_tf

# Apply fold change shrinkage
res_table_tf <- lfcShrink(dds, contrast=contrast_tf, res=res_table_tf, type='normal')
res_table_tf

#MA plot The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes tested. The genes that are significantly DE are colored to be easily identified. 
#plotMA(res_table_tf_unshrunken)
plotMA(res_table_tf)

## Summarize results
summary(res_table_tf, alpha = 0.01)

### Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

res_table_tf_tb <- res_table_tf %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#sig_tf <- res_table_tf_tb %>%
#   filter(padj < padj.cutoff)

sig_tf <- res_table_tf_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_tf


## Visualization------------------------------------------------------------------------------------

# DESeq2 creates a matrix when you use the counts() function
## Convert normalized_counts to a data frame and transfer the row names to a new column called "gene"
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene")

colnames(normalized_counts) <- sub("X", "", colnames(normalized_counts))

# Merge together (ensembl IDs) the normalized counts data frame with a subset of the annotations in the tx2gene data frame (only the columns for ensembl gene IDs and gene symbols)
grcm38annot <- tx2gene %>% 
  dplyr::select(V1, V2) %>% 
  dplyr::distinct()

## Bring in a column of gene symbols
normalized_counts <- merge(normalized_counts, grcm38annot, by.x="gene", by.y="V1")

# Now create a tibble for the normalized counts
#normalized_counts <- normalized_counts %>%
#  as_tibble()

normalized_counts 

# Data wrangling to name row names with gene IDs and saving the table

normalized_counts_geneIDs <- normalized_counts[ ,-1]
rownames(normalized_counts_geneIDs)=make.names(normalized_counts_geneIDs$V2, unique = TRUE)
normalized_counts_geneIDs <- normalized_counts_geneIDs[ ,-30]

#write.table(normalized_counts_geneIDs, file="data/normalized_counts_geneIDs.txt", sep="\t", quote=T)


## Single gene expression------------------------------------------------------------------------------------------------------------------
# Find the Ensembl ID of Ascl1 and Neurog2
#grcm38annot[grcm38annot$V2 == "Neurog2", "V1"]
#grcm38annot[grcm38annot$V2 == "Ascl1", "V1"]

# Save plotcounts to a data frame object
ngn2 <- plotCounts(dds, gene="ENSMUSG00000027967", intgroup="tf", returnData=TRUE) 
ascl1 <- plotCounts(dds, gene="ENSMUSG00000020052", intgroup="tf", returnData=TRUE) 

# Plot the Ngn2 normalized counts, using the samplenames (rownames(ngn2) as labels)
ggplot(ngn2, aes(x = tf, y = count, color = tf)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(ngn2))) + 
  theme_bw() +
  ggtitle("Ngn2") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot the Ngn2 normalized counts, using the samplenames (rownames(ngn2) as labels)
ggplot(ascl1, aes(x = tf, y = count, color = tf)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(ascl1))) + 
  theme_bw() +
  ggtitle("Ascl1") +
  theme(plot.title = element_text(hjust = 0.5))

  
## Multiple genes using ggplot2-----------------------------------------------------------------------------------------------------------
normCounts <- normalized_counts_geneIDs

normCounts$name <- rownames(normCounts)

colori <- c("blue","red","violet")

#
PanNeuron <- c("Ascl1","Neurog2","Aqp4","Sox9","Snap25", "Syp")

CountPanN <- subset(normCounts, normCounts$name %in% PanNeuron)

mdataPanN <- melt(CountPanN, id=c("name"))

## Rearranging columns
CountPanN2 <- CountPanN %>%
  select(1, 2, 7, 8, 9,13, 17, 18, 19, 20, 21, 27, 28, 29, 30, 3, 4, 5, 6, 11, 12, 22, 23, 24, 25, 26, 10, 14, 15, 16,)
mdataPanN <- melt(CountPanN2, id=c("name"))

meta$newOrder <- c(1,2,15,16,17,18,3,4,5,26,19,20,6,27,28,29,7,8,9,10,11,21,22,23,24,25,12,13,14)
meta$oldOrder <- rownames(meta)
meta2 <- arrange(meta,newOrder)

mdataPanN$type <- meta2$tf[names=mdataPanN$variable]

mdataPanN$valuep <-log2(mdataPanN$value)

# plot panneuronal genes
ggplot(mdataPanN) +
  geom_count(aes(x=variable, y=name, size=valuep, fill= type ), color='black', stroke=0.1, shape=21) + #, fill=type
  # scale_size_continuous(range = c(0,5)) +
  scale_size_area(max_size=4) +
  scale_fill_manual(values = colori) +
  #scale_y_discrete(limits = PanNeuron, position = "bottom") +
  scale_y_discrete(limits = PanNeuron, position = "right")+ #for y axis only right/left possible
  #scale_fill_manual(breaks=meta(pop_colors$Step2), values=pop_colors$Step2) +
  # scale_color_manual(breaks=meta(pop_colors$Step2), values=pop_colors$Step2) +
  xlab("") + ylab("") +
  coord_flip() + 
  #scale_x_discrete(position = "bottom") + 
  #  scale_y_discrete(position = "bottom") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0, size=12))

## Heatmap of all significant genes-------------------------------------------------------------------------------------------------------------
## Ngn2 vs Cntrl
### Extract normalized expression for significant genes from the Ngn2 and Cntrl samples (2:3,8:11,14:22,28:30)
norm_ngn2Cntrlsig <- normalized_counts[,c(1:3,8:11,14:22,28:30)] %>% 
  filter(gene %in% sig_ngn2Cntrl$gene)  

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_ngn2Cntrlsig[2:19], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = df, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Ascl1 vs Cntrl
### Extract normalized expression for significant genes from the Ascl1 and Cntrl samples (4:7,11:13,15:17,23:27)
norm_ascl1Cntrlsig <- normalized_counts[,c(1,4:7,11:13,15:17,23:27)] %>% 
  filter(gene %in% sig_ascl1Cntrl$gene)  

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_ascl1Cntrlsig[2:16], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = df, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Ngn2 vs Ascl1
### Extract normalized expression for significant genes from the Ngn2 and Ascl1 samples (2:10,12:14,18:30)
norm_TFsig <- normalized_counts[,c(1:10,12:14,18:30)] %>% 
  filter(gene %in% sig_tf$gene)  

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_TFsig[2:26], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = df, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Volcano plots (potential further visualization) -------------------------------------------------------------------------------------------------------------

## GO terms analysis (potential further functional analysis) -------------------------------------------------------------------------------------------------------------

## Likelihood ratio test (potential further hypothesis testing) -------------------------------------------------------------------------------------------------------------
