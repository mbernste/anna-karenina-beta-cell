#Import gene expression data to DESeq2 to analyze RNAseq with different visualization tools
#Make sure to edit "GROUP" in metaData based on clustering from dendrogram
#Run this before Venn diagram and volcano plot
library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")

#Use raw data
countData <- read.csv("/Users/matthewbernstein/Development/anna_karenina_beta_cell/data_for_DESeq2/Raw_Matrix_filter_rmDup.csv", header = TRUE, sep = ",")
head(countData)

#metaData without E13.5 samples, separate diabetic group into T1D and T2D
metaData <- read.csv("/Users/matthewbernstein/Development/anna_karenina_beta_cell/data_for_DESeq2/AKP_Metadata_v2.csv", header = TRUE, sep = ",")

head(metaData)

# Round the data
rnd <- function(x) trunc(x+sign(x)*0.5)
rnd_countData <- rnd(countData[,-1])
rnd_countData <- data.frame(GENE_ID = countData[,1],rnd_countData)

#Converting design variables into factors (not necessary)
metaData$COND <- as.factor(metaData$COND)
metaData$BATCH <- as.factor(metaData$BATCH)

run_deseq <- function(counts_mtx, meta_data, cond1, cond2, out_f) {
  dds <- DESeqDataSetFromMatrix(
    countData=counts_mtx, 
    colData=meta_data, 
    design=~COND, 
    tidy = TRUE
  )
  keep <- rowSums(counts(dds) >= 2) > 5
  dds <- dds[keep,]
  dds$COND <- relevel(dds$COND, ref = cond1)
  dds <- DESeq(dds)
  res <- results(dds)
  head(results(dds, tidy=TRUE))
  res <- results(dds, contrast=c("COND",cond1, cond2))
  res <- lfcShrink(dds, contrast=c("COND",cond1, cond2), res = res, type = "ashr")
  write.csv(res, out_f)
  EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "pvalue", xlim = c(-10,10), axisLabSize = 12, title = "", subtitle = "", captionLabSize = 8, pCutoff = 0.05, FCcutoff = 1.5, pointSize = 3.0, labFace = "bold", labSize = 3.0, legendPosition = 'none', legendLabSize = 10)

}

run_deseq(rnd_countData, metaData, "WT", "T2H", "/Users/matthewbernstein/Development/anna_karenina_beta_cell/results/WT_vs_T2H.csv")
run_deseq(rnd_countData, metaData, "T2H", "T2D", "/Users/matthewbernstein/Development/anna_karenina_beta_cell/results/T2H_vs_T2D.csv")
run_deseq(rnd_countData, metaData, "WT", "T2D", "/Users/matthewbernstein/Development/anna_karenina_beta_cell/results/WT_vs_T2D.csv")
