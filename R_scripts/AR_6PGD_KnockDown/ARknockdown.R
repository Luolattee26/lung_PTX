rm(list = ls())
setwd('C:/Users/Administrator/Desktop/lung_PTX/R_scripts/AR_6PGD_KnockDown/')
getwd()


library(limma)
library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
raw_cpm <- read.table('./GSE152254_All_data_cpm.txt', header = T, sep = '\t', row.names = 1) %>% .[, c(10:17)]
en_to_symbol <- read.table('./GSE152254_All_data_cpm.txt', header = T, sep = '\t') %>% .[, 1:2]

# limma DEGs
group_list <- factor(c(rep("siCon_48h", 4), rep("siAR_48h", 4)))
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(raw_cpm)

cont.matrix <- makeContrasts(contrasts = paste0('siAR_48h', '-', 'siCon_48h'), levels = design)

# the raw data is normalized,get log2(CPM)
log_cpm <- log2(raw_cpm)
# fit
fit1 <- lmFit(log_cpm, design) # linear fit
fit2 <- contrasts.fit(fit1, cont.matrix) # statistic test
efit <- eBayes(fit2, trend = T)

tempDEGs <- topTable(efit, coef = paste0('siAR_48h', '-', 'siCon_48h'), n = Inf)
DEG_limma_trend <- na.omit(tempDEGs) # final result
PGD_experison <- DEG_limma_trend['ENSG00000142657', ]

DEG_limma_trend$label <- NA
DEG_limma_trend['ENSG00000142657', ]$label <- 'PGD'

# plot
(valcano_plot <- ggplot(DEG_limma_trend, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
  geom_point(aes(size = -log10(adj.P.Val), color = -log10(adj.P.Val))) + theme_bw() +
  theme(panel.grid = element_blank()) + scale_color_gradientn(values = seq(0, 1, 0.2), 
                                                              colours = c('#3f56a6', '#3ec2e9', '#feec31', '#f16d68', '#c21d31')) +
  scale_size_continuous(range = c(1, 2)) + 
  geom_text(aes(label = label), vjust = 1.5))
ggsave(valcano_plot, filename = 'valcano_plot.pdf', path = './', width = 15, height = 20, units = 'cm')  


# enrichment
significant_genes <- DEG_limma_trend %>% filter(abs(logFC) >= 1, adj.P.Val < 0.05)
up_genes <- significant_genes %>% filter(.$logFC > 0)
down_genes <- significant_genes %>% filter(.$logFC < 0)

# get ENTREZID
up_ENTREZID <- bitr(row.names(up_genes), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
down_ENTREZID <- bitr(row.names(down_genes), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

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
                    ont = 'all', readable = T, pAdjustMethod = 'BH', pvalueCutoff = 0.1, qvalueCutoff = 0.25)
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








