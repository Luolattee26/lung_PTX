rm(list = ls())
setwd('C:/Users/Administrator/Desktop/lung_PTX/R_scripts/time_series/')
getwd()


# install packages
library(dplyr)
library(ggplot2)
if (require(MAGeCKFlute) == F){
  BiocManager::install('MAGeCKFlute')
}else{
  library(MAGeCKFlute)
}

# raw data input
raw_count <- read.csv('../../sgRNA_count_test123.count.txt', sep = '\t', header = T)

# QC
sgrna_countsummary <- read.csv('../../sgRNA_count_test123.countsummary.txt', sep = '\t', header = T)
QC <- function(countsummary){
  BarView(countsummary, x = "Label", y = "GiniIndex",ylab = "Gini index", main = "Evenness of sgRNA reads")
  ggsave(path = './count_QC_result', filename = 'GiniIndex.pdf', width = 30, units = 'cm')
  
  countsummary$Missed = log10(countsummary$Zerocounts)
  BarView(countsummary, x = "Label", y = "Missed", fill = "#394E12",ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
  ggsave(path = './count_QC_result', filename = 'Missed_sgrna.pdf', width = 30, units = 'cm')
  
  MapRatesView(countsummary)
  ggsave(path = './count_QC_result', filename = 'MapRates.pdf', width = 30, units = 'cm')
}
QC(sgrna_countsummary)

# cluster and PCA
batch_info_123 <- read.csv('../batch_MAT_123.txt', sep = '\t', header = T)
BatchRemove(mat = raw_count, batchMat = batch_info_123, cluster = T, pca = T, outdir = './Batch_removal_123/')

# translate nonessential gene list into a sgRNA list
sg_gene_list <- read.csv('../Human_GeCKOv2_Library_A_3_mageck.csv', header = F)
nonessential_gene_list <- read.table('../41596_2018_113_MOESM3_ESM.txt', fileEncoding = 'utf-16le', header = F)
# a pick up function
pick_up_nonessential <- function(gene_list, sg_to_gene){
  nonessential_sgrna_list <- data.frame()
  for (i in gene_list$V1) {
    tmp_df <- filter(sg_to_gene, V3 == i)
    nonessential_sgrna_list <- rbind(nonessential_sgrna_list, tmp_df)
  }
  return(nonessential_sgrna_list)
}
nonessential_sgrna_list <- pick_up_nonessential(nonessential_gene_list, sg_gene_list)
write.table(nonessential_sgrna_list$V1, '../nonessential_sgrna_list.txt', quote = F, row.names = F, col.names = T, fileEncoding = 'UTF-8')

# MLE downstream
# input posControl
The_list_of_core_essential_genes <- read.table('../../The_list_of_core_essential_genes.txt', fileEncoding = 'utf-16le', header = T)
FluteMLE('../../test_time.mle.gene_summary.txt', ctrlname = 'D7', treatname = 'D14', pathview.top = 50, cell_lines = 'A549', top = 5, omitEssential = T, proj = 'all_123', 
         enrich_method = 'HGT', scale_cutoff = 2, verbose = T, posControl = The_list_of_core_essential_genes$hgnc_symbol, width = 12, height = 18, norm_method = 'loess')
# select significant genes to rank
raw_MLE_result <- read.csv('../../test_time.mle.gene_summary.txt', sep = '\t', header = T)
selected_result <- read.csv('./MAGeCKFlute_all_123/SelectionData_ScatterView_TreatvsCtrl.txt', header = T, sep = '\t')
genes_significant_list <- (raw_MLE_result %>% filter(.$PTX.fdr < 0.15))$Gene
selected_result_significant <- data.frame()
for(gene in genes_significant_list){
  temp_df <- filter(selected_result, Gene == gene)
  selected_result_significant <- rbind(selected_result_significant, temp_df)
}
# get a rank in significant genes
selected_result_significant <- selected_result_significant %>% mutate(new_rank = rank(.$Diff))

# significant genes PLOT
library(ggrepel)
library(purrr)
# rank plot
rank_plot <- ggplot(selected_result, mapping = aes(x = Rank, y = Diff)) + geom_point(aes(color = group), size = 2.5) + 
  geom_point(data = filter(selected_result_significant, group == 'bottom', new_rank <= 5), size = 4, color = 'purple') + 
  geom_text_repel(data = filter(selected_result_significant, group == 'bottom', new_rank <= 5), aes(label = Gene), size = 4)
# Scatter View
Scatter_View <- ggplot(filter(selected_result, group == 'bottom'), mapping = aes(x = RandomIndex, y = Diff)) + geom_point() +
  geom_point(data = filter(selected_result_significant, group == 'bottom', new_rank <= 5), size = 4, color = 'blue') +
  geom_text_repel(data = filter(selected_result_significant, group == 'bottom', new_rank <= 5), aes(label = Gene), size = 4, color = 'red')
plots <- list(rank_plot, Scatter_View)
filename <- stringr::str_c(c('rank_plot', 'Scatter_View'), '.pdf')
pwalk(list(filename, plots), ggsave, path = './select_significant/', width = 30, height = 25, units = 'cm')





