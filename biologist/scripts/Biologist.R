## set wd
setwd("/projectnb2/bf528/users/dreadlocks/project_3/")

library(magrittr)
library(pheatmap)

## load the datasets
deg_AhR <- read.csv("ryan_tmp/result/deseq/deg_AhR.csv",
                    row.names = 1)
deg_CAR_PXR <- read.csv("ryan_tmp/result/deseq/deg_CAR_PXR.csv",
                            row.names = 1)
deg_PPARA <- read.csv("ryan_tmp/result/deseq/deg_PPARA.csv",
                            row.names = 1)

count_AhR <- read.csv("ryan_tmp/result/deseq/AhR_normalized_counts.csv",
                      row.names = 1)
count_CAR_PXR <- read.csv("ryan_tmp/result/deseq/CAR_PXR_normalized_counts.csv",
                        row.names = 1)
count_PPARA <- read.csv("ryan_tmp/result/deseq/PPARA_normalized_counts.csv",
                      row.names = 1)

## check the order
all.equal(count_AhR %>% row.names(), 
          count_CAR_PXR %>% row.names(), 
          count_PPARA %>% row.names())

count_table <- cbind(count_AhR, count_CAR_PXR, count_PPARA) # combine tables
count_table <- count_table[, !duplicated(colnames(count_table))]

my_metadata <- read.csv("ryan_tmp/sample/toxgroup_6_rna_info.csv")
mean(colnames(count_table) %in% my_metadata$Run)

annotation_col = data.frame(Vehicle = factor(my_metadata$vehicle[match(colnames(count_table), 
                                                                       my_metadata$Run)]))
sample_name <- my_metadata$mode_of_action[match(colnames(count_table), 
                                                my_metadata$Run)]
colnames(count_table) <- paste0(sample_name, "_", 1:3)
colnames(count_table)[10:12] <- paste0(rep("Control", 3), "_", 4:6)
head(count_table)
rownames(annotation_col) <- colnames(count_table)

## normalizing and filtering gene counts
count_table <- count_table %>% scale() # scale each column in table

cv_index <- apply(count_table, 1, function(x) sd(x)/mean(x))

## heatmap
pheatmap(count_table[cv_index > 0.2,] %>% as.matrix(),
         scale = "row",
         show_rownames = F,
         annotation_col = annotation_col)
