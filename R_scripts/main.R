rm(list = ls())
setwd('c:/Users/Administrator/Desktop/lung_PTX/R_scripts/')
getwd()


# install package
library(dplyr)
library(ggplot2)
if (require(MAGeCKFlute) == F){
  BiocManager::install('MAGeCKFlute')
}else{
  library(MAGeCKFlute)
}
if (require(sva) == F){
  BiocManager::install('sva')
}else{
  library(sva)
}


# raw data input
raw_count <- read.csv('../sgRNA_count_test123.count.txt', sep = '\t', header = T)
raw_count_only13 <- read.csv('../sgRNA_count_only13.count.txt', sep = '\t', header = T)

# QC
sgrna_countsummary <- read.csv('../sgRNA_count_test123.countsummary.txt', sep = '\t', header = T)
QC <- function(countsummary){
  BarView(countsummary, x = "Label", y = "GiniIndex",ylab = "Gini index", main = "Evenness of sgRNA reads")
  ggsave(path = './count_QC_result', filename = 'GiniIndex.pdf')
  
  countsummary$Missed = log10(countsummary$Zerocounts)
  BarView(countsummary, x = "Label", y = "Missed", fill = "#394E12",ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
  ggsave(path = './count_QC_result', filename = 'Missed_sgrna.pdf')
  
  MapRatesView(countsummary)
  ggsave(path = './count_QC_result', filename = 'MapRates.pdf')
}
QC(sgrna_countsummary)

# cluster and PCA
# all 3 test
batch_info_123 <- read.csv('../batch_MAT_123.txt', sep = '\t', header = T)
removed_result_123 <- BatchRemove(mat = raw_count, batchMat = batch_info_123, cluster = T, pca = T, outdir = './Batch_removal_123/')
raw_count_batch_removed_123 <- removed_result_123$data
write.table(raw_count_batch_removed_123, './Batch_removal_123/raw_count_batch_removed_123.txt', sep = '\t', quote = F, row.names = F, col.names = T)
# test if remove TEST2 samples
batch_info_only13 <- read.csv('../batch_MAT_only13.txt', sep = '\t', header = T)
removed_result_13 <- BatchRemove(mat = raw_count_only13, batchMat = batch_info_only13, cluster = T, pca = T, outdir = './Batch_removal_only13/')
raw_count_batch_removed_13 <- removed_result_13$data
raw_count_batch_removed_13_without_negetive <- raw_count_batch_removed_13[-which(rowSums(raw_count_batch_removed_13 < 0) >= 1), ]
write.table(raw_count_batch_removed_13_without_negetive, './Batch_removal_only13/raw_count_batch_removed_13.txt', sep = '\t', quote = F, row.names = F, col.names = T)
################################################################################
# a negative read count appear after batch removal
# delete all sgRNA with negative read count
################################################################################



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


# downstream functional analysis
library(clusterProfiler)
library(org.Hs.eg.db)
The_list_of_core_essential_genes <- read.table('../The_list_of_core_essential_genes.txt', fileEncoding = 'utf-16le', header = T)
# cell_cycle normalization
FluteMLE('../all_test123.mle.gene_summary.txt', ctrlname = 'DMSO', treatname = 'PTX', pathview.top = 50, cell_lines = 'A549', top = 5, omitEssential = T, proj = 'all_123', 
         enrich_method = 'HGT', scale_cutoff = 2, verbose = T, posControl = The_list_of_core_essential_genes$hgnc_symbol)
# GO&KEGG
# input select result
Selectionall_123squareview_data <- read.csv('./MAGeCKFlute_all_123/Selectionall_123squareview_data.txt', header = T, sep = '\t')
SelectionData_ScatterView_TreatvsCtrl <- read.csv('./MAGeCKFlute_all_123/SelectionData_ScatterView_TreatvsCtrl.txt', header = T, sep = '\t')

# function to generate dif_mat 
dif_mat_generator <- function(mle_result, mle_result_type, position){
  if(mle_result_type == 'tvc'){
    if(position == 'bottom'){
      dif_mat <- filter(mle_result, mle_result$group == 'bottom')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
    if(position == 'top'){
      dif_mat <- filter(mle_result, mle_result$group == 'top')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
  }else if(mle_result_type == 'square'){
    if(position == 'midleft'){
      dif_mat <- filter(mle_result, mle_result$group == 'midleft')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
    if(position == 'topcenter'){
      dif_mat <- filter(mle_result, mle_result$group == 'topcenter')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
    if(position == 'midright'){
      dif_mat <- filter(mle_result, mle_result$group == 'midright')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
    if(position == 'bottomcenter'){
      dif_mat <- filter(mle_result, mle_result$group == 'bottomcenter')
      dif_mat <- as.data.frame(cbind(dif_mat$Symbol, dif_mat$DMSO, dif_mat$PTX, dif_mat$Rank))
      colnames(dif_mat) <- c('Gene', 'DMSO', 'PTX', 'Rank')
      gene_list <- dif_mat$Gene
      rownames(dif_mat) <- gene_list
      dif_mat <- dif_mat[-1]
      return(dif_mat)
    }
  }
}

# function to enrich result visualization
enrich_result_plot <- function(enrich_result, enrich_type, top, plot_type, outdir){
  if(plot_type == 'bar'){
    if(enrich_type == 'go'){
      go_enrich_result <- data.frame(enrich_result)
      write.csv(go_enrich_result, './MAGeCKFlute_all_123/GO/go_enrich_result.csv', quote = F)
      BP_result <- subset(go_enrich_result, subset = (ONTOLOGY == "BP"))[1:top, ]
      CC_result <- subset(go_enrich_result, subset = (ONTOLOGY == "CC"))[1:top, ]
      MF_result <- subset(go_enrich_result, subset = (ONTOLOGY == "MF"))[1:top, ]
      all_result <- rbind(BP_result, CC_result, MF_result)
      all_result <- na.omit(all_result)
      all_result$Description <- factor(all_result$Description,levels = rev(all_result$Description))
      bar_plot <- ggplot(all_result, mapping = aes(x = Description, y = Count, fill = ONTOLOGY)) +
        geom_bar(stat = 'identity', width = 0.3) + coord_flip() + theme_bw() + 
        scale_x_discrete(labels = function(x){
          strwrap(x, width = 45)
        }) + labs(x = 'GO terms', y = 'Number') + 
        theme(axis.title = element_text(size = 13),
              axis.text = element_text(size = 11),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 11),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
      ggsave(bar_plot, path = outdir, filename = 'go_bar_plot.pdf')
    }else if(enrich_type == 'kegg'){
      kegg_enrich_result <- data.frame(enrich_result)
      kegg_enrich_result <- kegg_enrich_result[1:top, ]
      kegg_enrich_result <- na.omit(kegg_enrich_result)
      kegg_enrich_result$Description <- factor(kegg_enrich_result$Description,levels = rev(kegg_enrich_result$Description))
      bar_plot <- ggplot(kegg_enrich_result, mapping = aes(x = Description, y = Count, fill = p.adjust)) +
        geom_bar(stat = 'identity', width = 0.3) + coord_flip() + theme_bw() + 
        scale_x_discrete(labels = function(x){
          strwrap(x, width = 45)
        }) + labs(x = 'KEGG pathways', y = 'Number') + 
        theme(axis.title = element_text(size = 13),
              axis.text = element_text(size = 11),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 11),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
      ggsave(bar_plot, path = outdir, filename = 'kegg_bar_plot.pdf')
    }
  }else if(plot_type == 'point'){
    if(enrich_type == 'kegg'){
      kegg_enrich_result <- data.frame(enrich_result)
      kegg_enrich_result <- kegg_enrich_result[1:top, ]
      kegg_enrich_result <- na.omit(kegg_enrich_result)
      kegg_enrich_result$Description = factor(kegg_enrich_result$Description, levels = rev(kegg_enrich_result$Description))
      point_plot <- ggplot(data = kegg_enrich_result,aes(x = GeneRatio, y = reorder(Description, Count))) +
        geom_point(aes(size = Count, color = -log10(p.adjust)))+
        theme_bw() +
        scale_colour_gradient(low = "blue",high = "red") +
        labs(x = "GeneRatio",y = "", color = expression(-log10(p.adjust)), size = "Count") +
        theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
              legend.title = element_text(size = 11),legend.text = element_text(size = 10))
      point_plot
      ggsave(point_plot, path = outdir, filename = 'kegg_point_plot.pdf', width = 25, height = 25, units = 'cm')
    }else if(enrich_type == 'go'){
      go_enrich_result <- data.frame(enrich_result)
      write.csv(go_enrich_result, './MAGeCKFlute_all_123/GO/go_enrich_result.csv', quote = F)
      BP_result <- subset(go_enrich_result, subset = (ONTOLOGY == "BP"))[1:top, ]
      CC_result <- subset(go_enrich_result, subset = (ONTOLOGY == "CC"))[1:top, ]
      MF_result <- subset(go_enrich_result, subset = (ONTOLOGY == "MF"))[1:top, ]
      all_result <- rbind(BP_result, CC_result, MF_result)
      all_result <- na.omit(all_result)
      go_enrich_result$Description = factor(go_enrich_result$Description, levels = rev(go_enrich_result$Description))
      point_plot <- ggplot(data = go_enrich_result,aes(x = GeneRatio, y = reorder(Description, Count), fill = ONTOLOGY)) +
        geom_point(aes(size = Count, color = -log10(p.adjust)))+
        theme_bw() +
        scale_colour_gradient(low = "blue",high = "red") +
        scale_y_discrete(labels = function(x){
          strwrap(x,width = 40)
        }) +
        labs(x = "GeneRatio",y = "", color = expression(-log10(p.adjust)), size = "Count") +
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10) ,
              legend.title = element_text(size = 11), legend.text = element_text(size = 10))
      ggsave(point_plot, path = outdir, filename = 'go_point_plot.pdf', width = 25, height = 25, units = 'cm')
    }
  }
}

# get go dif_mat
dif_mat_go <- dif_mat_generator(SelectionData_ScatterView_TreatvsCtrl, 'tvc', 'bottom')
dif_mat_go$Rank <- as.numeric(dif_mat_go$Rank)

id_list <- TransGeneID(rownames(dif_mat_go), 'Symbol', 'Entrez')
id_list <- na.omit(id_list)

# perform go enrichment
go_enrich <- enrichGO(
  gene = id_list,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID',
  ont = 'ALL',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = T
)

# get kegg enrichment result from FluteMLE result
kegg_enrich <- read.csv('./MAGeCKFlute_all_123/Enrichment/GroupB_kegg_cell_cycle.txt', header = T, sep = '\t')
kegg_enrich <- filter(kegg_enrich, p.adjust < 0.1, pvalue < 0.05)


# GO&KEGG analysis visualization
enrich_result_plot(kegg_enrich, 'kegg', 10, 'point', './MAGeCKFlute_all_123/KEGG/')
enrich_result_plot(kegg_enrich, 'kegg', 10, 'bar', './MAGeCKFlute_all_123/KEGG/')
enrich_result_plot(go_enrich, 'go', 5, 'point', './MAGeCKFlute_all_123/GO/')
enrich_result_plot(go_enrich, 'go', 5, 'bar', './MAGeCKFlute_all_123/GO/')





# dif_mat_kegg <- dif_mat_generator(SelectionData_ScatterView_TreatvsCtrl, 'tvc', 'bottom')
# dif_mat_kegg$Rank <- as.numeric(dif_mat_kegg$Rank)

# id_list <- TransGeneID(rownames(dif_mat_kegg), 'Symbol', 'Entrez')
# id_list <- na.omit(id_list)

# options(clusterProfiler.download.method = "wininet")
# kegg_enrich <- enrichKEGG(
#   gene = id_list,
#   organism = 'hsa',
#   keyType = 'kegg',
#   pAdjustMethod = 'fdr',
#   pvalueCutoff = 0.01,
#   qvalueCutoff = 0.05,
#   minGSSize = 2,
#   maxGSSize = 200,
#   use_internal_data = F
# )

# kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
# kegg_enrich = data.frame(kegg_enrich@result)



















