## set wd
setwd("/projectnb2/bf528/users/dreadlocks/project_3/")

## list files
count_file <- list.files("ryan_tmp/result/readcount/", pattern = ".txt$", full.names = T)

## cbind tables
count_table <- read.table(count_file[1], header = T)
count_table <- count_table[,c(1,7)]

for (i in count_file[-1]){
  
  tmp <- read.table(i, header = T)
  tmp <- tmp[,7]
  
  count_table <- cbind(count_table, tmp)
  
}
colnames(count_table) <- c("gene","SRR1177963","SRR1177964",
                           "SRR1177965","SRR1177997",
                           "SRR1177999","SRR1178002",
                           "SRR1178014","SRR1178021","SRR1178047")
head(count_table)

boxplot(count_table[,-1], outline = F, xlab = "Sample",
        ylab = "Count")

## deg analysis
library(DESeq2)
library(apeglm)

my_metadata <- read.csv("ryan_tmp/sample/toxgroup_6_rna_info.csv")
my_metadata

mean(colnames(count_table)[-1] %in% my_metadata$Run) # check all files are matched with metadata

control_matrix <- read.csv("/project/bf528/project_3/samples/control_counts.csv")
all.equal(count_table$gene, control_matrix$Geneid) # check two table are in a same order by gene

control_matrix <- control_matrix[, my_metadata$Run[my_metadata$mode_of_action == "Control"]]
head(control_matrix)

count_table <- cbind(count_table, control_matrix)
head(count_table)

count_table <- subset(count_table,rowSums(count_table==0)==0) # remove all zero count
head(count_table)

rownames(count_table) <- count_table$gene
count_table$gene <- NULL
head(count_table)

## deg object
all.equal(colnames(count_table), my_metadata$Run) # check the order
count_table <- count_table[my_metadata$Run]
all.equal(colnames(count_table), my_metadata$Run)

index_AhR <- with(my_metadata,
                  mode_of_action == "AhR" | 
                    (mode_of_action == "Control" & vehicle == "CMC_.5_%"))
index_CAR_PXR <- with(my_metadata,
                      mode_of_action == "CAR/PXR" | 
                        (mode_of_action == "Control" & vehicle == "CORN_OIL_100_%"))
index_PPARA <- with(my_metadata,
                    mode_of_action == "PPARA" | 
                      (mode_of_action == "Control" & vehicle == "CMC_.5_%"))
meta_data_AhR <- my_metadata[index_AhR,]
meta_data_CAR_PXR <- my_metadata[index_CAR_PXR,]
meta_data_PPARA <- my_metadata[index_PPARA,]

meta_data_AhR$mode_of_action <- factor(meta_data_AhR$mode_of_action)
meta_data_CAR_PXR$mode_of_action <- factor(meta_data_CAR_PXR$mode_of_action)
meta_data_PPARA$mode_of_action <- factor(meta_data_PPARA$mode_of_action)

dds_AhR <- DESeqDataSetFromMatrix(
  countData = count_table[meta_data_AhR$Run],
  colData = meta_data_AhR,
  design= ~ mode_of_action
)
dds_AhR$mode_of_action <- relevel(dds_AhR$mode_of_action, ref='Control')
dds_AhR <- DESeq(dds_AhR)

res_AhR <- results(dds_AhR, alpha = 0.05)
summary(res_AhR)
res_AhR_shrink <- lfcShrink(dds_AhR, coef = 2, type="apeglm")
summary(res_AhR_shrink)
res_AhR_shrink[order(res_AhR_shrink$padj, decreasing = F),][1:10,]

dds_CAR_PXR <- DESeqDataSetFromMatrix(
  countData = count_table[meta_data_CAR_PXR$Run],
  colData = meta_data_CAR_PXR,
  design= ~ mode_of_action
)
dds_CAR_PXR$mode_of_action <- relevel(dds_CAR_PXR$mode_of_action, ref='Control')
dds_CAR_PXR <- DESeq(dds_CAR_PXR)

res_CAR_PXR <- results(dds_CAR_PXR, alpha = 0.05)
summary(res_CAR_PXR)
res_CAR_PXR_shrink <- lfcShrink(dds_CAR_PXR, coef = 2, type="apeglm")
summary(res_CAR_PXR_shrink)
res_CAR_PXR_shrink[order(res_CAR_PXR_shrink$padj, decreasing = F),][1:10,]

dds_PPARA <- DESeqDataSetFromMatrix(
  countData = count_table[meta_data_PPARA$Run],
  colData = meta_data_PPARA,
  design= ~ mode_of_action
)
dds_PPARA$mode_of_action <- relevel(dds_PPARA$mode_of_action, ref='Control')
dds_PPARA <- DESeq(dds_PPARA)

res_PPARA <- results(dds_PPARA, alpha = 0.05)
summary(res_PPARA)
res_PPARA_shrink <- lfcShrink(dds_PPARA, coef = 2, type="apeglm")
summary(res_PPARA_shrink)
res_PPARA_shrink[order(res_PPARA_shrink$padj, decreasing = F),][1:10,]

deg_AhR <- subset(res_AhR_shrink, padj <= 0.05)
deg_AhR <- deg_AhR[order(deg_AhR$padj, decreasing = F),]

deg_CAR_PXR <- subset(res_CAR_PXR_shrink, padj <= 0.05)
deg_CAR_PXR <- deg_CAR_PXR[order(deg_CAR_PXR$padj, decreasing = F),]

deg_PPARA <- subset(res_PPARA_shrink, padj <= 0.05)
deg_PPARA <- deg_PPARA[order(deg_PPARA$padj, decreasing = F),]

## save files
write.csv(res_AhR_shrink, file = "ryan_tmp/result/deseq/res_AhR_shrink.csv",
          row.names = T)
write.csv(deg_AhR, file = "ryan_tmp/result/deseq/deg_AhR.csv",
          row.names = T)

write.csv(res_CAR_PXR_shrink, file = "ryan_tmp/result/deseq/res_CAR_PXR_shrink.csv",
          row.names = T)
write.csv(deg_CAR_PXR, file = "ryan_tmp/result/deseq/deg_CAR_PXR.csv",
          row.names = T)

write.csv(res_PPARA_shrink, file = "ryan_tmp/result/deseq/res_PPARA_shrink.csv",
          row.names = T)
write.csv(deg_PPARA, file = "ryan_tmp/result/deseq/deg_PPARA.csv",
          row.names = T)

write.csv(counts(dds_AhR, normalized=TRUE), 
          file = 'ryan_tmp/result/deseq/AhR_normalized_counts.csv',
          row.names = T)
write.csv(counts(dds_CAR_PXR, normalized=TRUE), 
          file = 'ryan_tmp/result/deseq/CAR_PXR_normalized_counts.csv',
          row.names = T)
write.csv(counts(dds_PPARA, normalized=TRUE), 
          file = 'ryan_tmp/result/deseq/PPARA_normalized_counts.csv',
          row.names = T)



## plot
par(mfrow=c(1,3))
hist(res_AhR_shrink$log2FoldChange[res_AhR_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "AhR_vs_Control")
hist(res_CAR_PXR_shrink$log2FoldChange[res_CAR_PXR_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "CAR_PXR_vs_Control")
hist(res_PPARA_shrink$log2FoldChange[res_PPARA_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "PPARAvs_Control")

par(mfrow=c(1,3))
plot(res_AhR_shrink$log2FoldChange[res_AhR_shrink$padj <= 0.05], 
     -log10(res_AhR_shrink$pvalue[res_AhR_shrink$padj <= 0.05]), xlab = "Log2FC", 
     ylab = "-log10(p_value)", main = "AhR_vs_Control")
plot(res_CAR_PXR_shrink$log2FoldChange[res_CAR_PXR_shrink$padj <= 0.05], 
     -log10(res_CAR_PXR_shrink$pvalue[res_CAR_PXR_shrink$padj <= 0.05]), 
     xlab = "Log2FC", ylab = "-log10(p_value)", main = "CAR_PXR_vs_Control")
plot(res_PPARA_shrink$log2FoldChange[res_PPARA_shrink$padj <= 0.05], 
     -log10(res_PPARA_shrink$pvalue[res_PPARA_shrink$padj <= 0.05]), 
     xlab = "Log2FC", ylab = "-log10(p_value)", main = "PPARAvs_Control")
