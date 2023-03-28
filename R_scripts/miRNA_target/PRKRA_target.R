rm(list = ls())
setwd('C:/Users/Administrator/Desktop/lung_PTX/R_scripts/miRNA_target/')
getwd()

if (require(multiMiR) == F){
  BiocManager::install('multiMiR')
}else{
  library(multiMiR)
}
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)
library(clusterProfiler)

fold_change_data <- read.delim('GSE119845_fold_change.txt', header = T, sep = '\t')
significant_up <- fold_change_data[order(fold_change_data$Fold_Change, decreasing = T), ] %>%
  filter(.$pvalue < 0.05)

#################################################
# heatmap_df <- significant_up[, 5:6]
# pheatmap(heatmap_df, cluster_cols = F, cluster_rows = F, )
#################################################

most_upregurated <- fold_change_data[order(fold_change_data$Fold_Change, decreasing = T), ] %>%
  filter(.$pvalue < 0.05) %>% .[1:11, ]
miRNA_list <- most_upregurated$mirna
# find miRNA targets
# re-Formatting
new_miRNA_list <- vector('list', length = length(miRNA_list))
count = 0
for(miRNA in miRNA_list){
  miRNA_new <- miRNA %>% strsplit('_') %>% .[[1]] %>% .[1:length(.)-1] %>% 
    paste(collapse = '-')
  count <- count + 1
  new_miRNA_list[count] <- miRNA_new
}
new_miRNA_list


target_result <- get_multimir(org = 'hsa',
                              mirna = new_miRNA_list,
                              table = 'validated',
                              summary = T)
all_targets <- target_result@data
all_targets_experiment <- target_result@data[grep('Luciferase|qRT-PCR|Western blot|ChIP-seq|Immunoprecipitaion', 
                                                  target_result@data[, 'experiment']), ]

### Enrichment
## mRNA with experiments support
# GO
go_experiment <- enrichGO(gene = all_targets_experiment$target_entrez, OrgDb = 'org.Hs.eg.db',
               ont = 'all', readable = T, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
go_result_experiment <- as.data.frame(go_experiment@result)
# KEGG
KEGG_experiment <- enrichKEGG(gene = all_targets_experiment$target_entrez, organism = 'hsa',
                          keyType = 'kegg', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
KEGG_result_experiment <- setReadable(KEGG_experiment, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
KEGG_result_experiment <- data.frame(KEGG_result_experiment@result)
## all mRNA
# GO
go_all <- enrichGO(gene = all_targets$target_entrez, OrgDb = 'org.Hs.eg.db',
                          ont = 'all', readable = T, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
go_result_all <- data.frame(go_all)
# KEGG
KEGG_all <- enrichKEGG(gene = all_targets$target_entrez, organism = 'hsa',
                              keyType = 'kegg', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1)
KEGG_result_all <- setReadable(KEGG_all, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
KEGG_result_all <- data.frame(KEGG_result_all@result)

### Plot
## with experiment
# GO
goBP <- subset(go_result_experiment, subset = (ONTOLOGY == "BP"))[1:10, ]
goCC <- subset(go_result_experiment, subset = (ONTOLOGY == "CC"))[1:10, ]
goMF <- subset(go_result_experiment, subset = (ONTOLOGY == "MF"))[1:10, ]
(go_bar <- ggplot(data = goBP, mapping = aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
  scale_fill_gradient(low = 'red', high = 'orange') + 
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(go_bar, filename = 'go_barplot.pdf', path = './enrich_result/', width = 20, units = 'cm')
# KEGG
(kegg_bar <- ggplot(data = KEGG_result_experiment[1:15, ], mapping = aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = 'identity', width = 0.7) + coord_flip() + theme_bw() + 
    scale_fill_gradient(low = 'red', high = 'orange') + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)))
ggsave(kegg_bar, filename = 'kegg_barplot.pdf', path = './enrich_result/', width = 20, units = 'cm')



