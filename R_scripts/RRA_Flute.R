rm(list = ls())
setwd('c:/Users/Administrator/Desktop/lung_PTX/R_scripts/')
getwd()


# install package
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
if (require(MAGeCKFlute) == F){
  BiocManager::install('MAGeCKFlute')
}else{
  library(MAGeCKFlute)
}


RRA_result <- read.csv('RRA_only13/only13.RRA.gene_summary.txt', sep = '\t', header = T)
dif_MAT_enrich <- RRA_result[, 1:6][, -2]
# filter significant genes to ENRICH
dif_MAT_enrich <- filter(dif_MAT_enrich, neg.fdr < 0.05, neg.p.value < 0.05)


# ENRICH
# GO
go_enrich <- enrichGO(
  gene = dif_MAT_enrich$id,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = 'ALL',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = T
)
go_result <- data.frame(go_enrich)

# KEGG
dif_MAT_enrich_kegg <- filter(RRA_result, neg.rank <= 16)

options(clusterProfiler.download.method = "wininet")
id_list <- TransGeneID(dif_MAT_enrich_kegg$id, 'Symbol', 'Entrez')
id_list <- na.omit(id_list)

kegg_enrich <- enrichKEGG(
  gene = id_list,
  organism = 'hsa',
  keyType = 'kegg',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  minGSSize = 2,
  maxGSSize = 200,
  use_internal_data = F
)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
kegg_enrich <- data.frame(kegg_enrich@result)


# function to enrich result visualization
enrich_result_plot <- function(enrich_result, enrich_type, top, plot_type, outdir){
  if(plot_type == 'bar'){
    if(enrich_type == 'go'){
      go_enrich_result <- data.frame(enrich_result)
      write.csv(go_enrich_result, './RRA_only13/GO/go_enrich_result.csv', quote = F)
      BP_result <- subset(go_enrich_result, subset = (ONTOLOGY == "BP"))[1:top, ]
      CC_result <- subset(go_enrich_result, subset = (ONTOLOGY == "CC"))[1:top, ]
      MF_result <- subset(go_enrich_result, subset = (ONTOLOGY == "MF"))[1:top, ]
      all_result <- rbind(BP_result, CC_result, MF_result)
      all_result <- na.omit(all_result)
      all_result$Description <- factor(all_result$Description,levels = rev(all_result$Description))
      bar_plot <- ggplot(all_result, mapping = aes(x = Description, y = Count, fill = ONTOLOGY)) +
        geom_bar(stat = 'identity', width = 0.3) + coord_flip() + theme_bw() + 
        labs(x = 'GO terms', y = 'Number') + 
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
      write.csv(kegg_enrich_result, './RRA_only13/KEGG/kegg_enrich_result.csv', quote = F)
      kegg_enrich_result$Description <- factor(kegg_enrich_result$Description, levels = rev(kegg_enrich_result$Description))
      bar_plot <- ggplot(kegg_enrich_result, mapping = aes(x = Description, y = Count, fill = p.adjust)) +
        geom_bar(stat = 'identity', width = 0.3) + coord_flip() + theme_bw() + 
        labs(x = 'KEGG pathways', y = 'Number') + 
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
      write.csv(kegg_enrich_result, './RRA_only13/KEGG/kegg_enrich_result.csv', quote = F)
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
      all_result$Description = factor(all_result$Description, levels = rev(all_result$Description))
      point_plot <- ggplot(data = all_result,aes(x = GeneRatio, y = reorder(Description, Count), fill = ONTOLOGY)) +
        geom_point(aes(size = Count, color = -log10(p.adjust)))+
        theme_bw() +
        scale_colour_gradient(low = "blue",high = "red") +
        labs(x = "GeneRatio",y = "", color = expression(-log10(p.adjust)), size = "Count") +
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10) ,
              legend.title = element_text(size = 11), legend.text = element_text(size = 10))
      ggsave(point_plot, path = outdir, filename = 'go_point_plot.pdf', width = 25, height = 25, units = 'cm')
    }
  }
}


enrich_result_plot(kegg_enrich, 'kegg', 10, 'point', './RRA_only13/KEGG/')
enrich_result_plot(kegg_enrich, 'kegg', 10, 'bar', './RRA_only13/KEGG/')
enrich_result_plot(go_enrich, 'go', 5, 'point', './RRA_only13/GO/')
enrich_result_plot(go_enrich, 'go', 5, 'bar', './RRA_only13/GO/')


