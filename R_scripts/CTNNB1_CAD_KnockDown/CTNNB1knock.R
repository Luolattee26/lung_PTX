rm(list = ls())
setwd('C:/Users/Administrator/Desktop/lung_PTX/R_scripts/CTNNB1_CAD_KnockDown/')
getwd()


library(ggplot2)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)
library(limma)
# cutoff:p.adj<0.05 and abs(lfc)>1
raw_lfc_data <- read.csv('./GSE199835_HCT116_processed_data.txt', header = F, sep = '\t') %>% 
  .[-1, ]
dim(raw_lfc_data)
colnames(raw_lfc_data) <- raw_lfc_data[1, ]
raw_lfc_data <- raw_lfc_data[-1, ]
# filter up|dowm genes
# gene id transfer
genes_up <- raw_lfc_data %>% filter(.$`log2(FC)` > 0)
up_ENTREZID <- bitr(genes_up$`Gene ID`, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
genes_down <- raw_lfc_data %>% filter(.$`log2(FC)` < 0)
down_ENTREZID <- bitr(genes_down$`Gene ID`, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# DEGs
# raw data
df1 <- as.matrix(read.table('raw_featurecount/GSM5988143_featureCounts_HCT116_KO_44_A.txt', sep = '\t', header = F) %>% .[-1, ]) 
df2 <- as.matrix(read.table('raw_featurecount/GSM5988144_featureCounts_HCT116_KO_44_B.txt', sep = '\t', header = F) %>% .[-1, ])
df3 <- as.matrix(read.table('raw_featurecount/GSM5988145_featureCounts_HCT116_KO_47_A.txt', sep = '\t', header = F) %>% .[-1, ])
df4 <- as.matrix(read.table('raw_featurecount/GSM5988146_featureCounts_HCT116_KO_47_B.txt', sep = '\t', header = F) %>% .[-1, ])
df5 <- as.matrix(read.table('raw_featurecount/GSM5988147_featureCounts_HCT116_WT_4_A.txt', sep = '\t', header = F) %>% .[-1, ])
df6 <- as.matrix(read.table('raw_featurecount/GSM5988148_featureCounts_HCT116_WT_4_B.txt', sep = '\t', header = F) %>% .[-1, ])
df7 <- as.matrix(read.table('raw_featurecount/GSM5988149_featureCounts_HCT116_WT_6_A.txt', sep = '\t', header = F) %>% .[-1, ])
df8 <- as.matrix(read.table('raw_featurecount/GSM5988150_featureCounts_HCT116_WT_6_B.txt', sep = '\t', header = F) %>% .[-1, ])

exp_df <- cbind(df1, df2[, 2], df3[, 2], df4[, 2], df5[, 2], df6[, 2], df7[, 2], df8[, 2])
colnames(exp_df) <- c('GeneID', 'KO1', 'KO2', 'KO3', 'KO4', 'WT1', 'WT2', 'WT3', 'WT4')
rownames(exp_df) <- exp_df[, 1]
exp_df <- exp_df[, -1]
exp_num <- matrix(as.numeric(exp_df), nrow = nrow(exp_df))
colnames(exp_num) <- c('KO1', 'KO2', 'KO3', 'KO4', 'WT1', 'WT2', 'WT3', 'WT4')
rownames(exp_num) <- rownames(exp_df)

##### limma DEG
design_MAT <- tribble(
  ~KD, ~WT,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  0, 1,
  0, 1,
  0, 1,
  0, 1
)
design_MAT <- as.data.frame(design_MAT)
rownames(design_MAT) <- colnames(exp_df)

fit1 <- lmFit(exp_num, design_MAT)
df.matrix <- makeContrasts(KD - WT, levels = design_MAT)
fit <- contrasts.fit(fit1, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit, n = Inf, adjust = "BH")
dim(tempOutput %>% filter(abs(logFC) > 1, adj.P.Val < 0.05))

ggplot(tempOutput, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(aes(size = -log10(adj.P.Val), color = -log10(adj.P.Val))) + theme_bw() +
  theme(panel.grid = element_blank()) + scale_color_gradientn(values = seq(0, 1, 0.2), 
                                                              colours = c('#3f56a6', '#3ec2e9', '#feec31', '#f16d68', '#c21d31')) +
  scale_size_continuous(range = c(1, 3.5))

### enrich
## up
# go
go_up <- enrichGO(gene = up_ENTREZID$ENTREZID, OrgDb = 'org.Hs.eg.db',
                  ont = 'all', readable = T, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
up_result_go <- as.data.frame(go_up@result)
# kegg
kegg_up <- enrichKEGG(gene = up_ENTREZID$ENTREZID, organism = 'hsa',
                              keyType = 'kegg', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
up_result_kegg <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
up_result_kegg <- data.frame(up_result_kegg@result)
## down
# go
go_down <- enrichGO(gene = down_ENTREZID$ENTREZID, OrgDb = 'org.Hs.eg.db',
                  ont = 'all', readable = T, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
down_result_go <- as.data.frame(go_down@result)
# kegg
kegg_down <- enrichKEGG(gene = down_ENTREZID$ENTREZID, organism = 'hsa',
                      keyType = 'kegg', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
down_result_kegg <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
down_result_kegg <- data.frame(down_result_kegg@result)


### Plot
## up
# GO
goBP <- subset(up_result_go, subset = (ONTOLOGY == "BP"))[1:10, ]
goCC <- subset(up_result_go, subset = (ONTOLOGY == "CC"))[1:10, ]
goMF <- subset(up_result_go, subset = (ONTOLOGY == "MF"))[1:10, ]
(go_bar <- ggplot(data = goBP, mapping = aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
    scale_fill_gradient(low = '#c82044', high = '#4058a6') + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(go_bar, filename = 'go_barplot_up.pdf', path = './enrichment/', width = 20, units = 'cm')
# KEGG
(kegg_bar <- ggplot(data = up_result_kegg[1:15, ], mapping = aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
    scale_fill_gradient(low = '#c82044', high = '#4058a6') + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(kegg_bar, filename = 'kegg_barplot_up.pdf', path = './enrichment/', width = 20, units = 'cm')
## down
# GO
goBP <- subset(down_result_go, subset = (ONTOLOGY == "BP"))[1:10, ]
goCC <- subset(down_result_go, subset = (ONTOLOGY == "CC"))[1:10, ]
goMF <- subset(down_result_go, subset = (ONTOLOGY == "MF"))[1:10, ]
(go_bar <- ggplot(data = goBP, mapping = aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
    scale_fill_gradient(low = '#c82044', high = '#4058a6') + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(go_bar, filename = 'go_barplot_down.pdf', path = './enrichment/', width = 20, units = 'cm')
# KEGG
(kegg_bar <- ggplot(data = down_result_kegg[1:15, ], mapping = aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
    scale_fill_gradient(low = '#c82044', high = '#4058a6') + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(kegg_bar, filename = 'kegg_barplot_down.pdf', path = './enrichment/', width = 20, units = 'cm')





