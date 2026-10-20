

####4.R可视化####

library(SCENIC)
library(SCopeLoomR)
library(pheatmap) 
library(tidyverse)
library(AUCell)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
packageVersion("SCENIC")

scRNA_harmony

scenicLoomPath='D:/B R code/packages/SCENIC/out_SCENIC.loom'
loom <- open_loom(scenicLoomPath)
loom<-open_loom('D:/B R code/packages/SCENIC/out_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
AUC_df = getAUC(regulonAUC[,])
row.names(AUC_df)
rownames(AUC_df) = sub("...$","",rownames(AUC_df))
rownames(AUC_df)
n=t(scale(t(AUC_df)))
n[n>2]=2 
n[n< -2]= -2
dim(n)
scRNA_harmony=scRNA_harmony[,colnames(AUC_df)]
ac=data.frame(group= as.character(Idents(scRNA_harmony)))
n=n[,colnames(n) %in% colnames(scRNA_harmony)]

rownames(ac)=colnames(n) 
cg=read.table('D:/BaiduSyncdisk/R code/packages/SCENIC/allTFs_hg38.txt')[,1]
cg_AUC_all=n[rownames(n) %in% cg,]
dim(cg_AUC_all) 
sub_regulonAUC <-cg_AUC_all
cellTypes <-  as.data.frame(subset(scRNA_harmony@meta.data,select = 'celltypes'))

####rss####
dir.create("./pyscenic")
rss <- calcRSS(AUC=AUC_df, 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        "celltypes"])
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 3,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
ggsave("rssPlot.pdf",path = "./pyscenic")

###rss2####
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
n1=scale(t(n),scale = T,center = T)
n=t(n1)
colnames(n)=factor(Idents(scRNA_harmony),levels = Idents(scRNA_harmony))
col_fun =colorRamp2::colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
anno_col <- color_qualitative2[1:length(colnames(n))]
names(anno_col) <- as.character(colnames(n))
column_ha = HeatmapAnnotation(cluster = colnames(n), col =list(cluster = anno_col))
dir.create("./pyscenic")
pdf(file="./pyscenic/Heatmap of specific TFs.pdf",width = length(anno_col)/1.5,height =length(anno_col)*1.2)
print(Heatmap(n,name = "Z-score",top_annotation = column_ha,
              cluster_columns = F,cluster_rows = F,
              row_title = "Top 10 Marker genes",
              row_names_gp = gpar(fontface = 'italic',fontsize = 10),
              row_names_side = 'left',
              border = T,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_side = 'top',
              column_names_rot = 45,
              col = col_fun))
dev.off()



####Heatmap Featureplot####
pbmc <- CreateSeuratObject(counts = cg_AUC_all)
pbmc<- AddMetaData(pbmc, metadata = scRNA_harmony@meta.data)
pbmc@reductions=scRNA_harmony@reductions
Idents(pbmc)=pbmc$celltypes
allmarkers=FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- allmarkers%>% group_by(cluster) %>%  filter(avg_log2FC > 0& p_val_adj<0.05) %>%top_n(w =avg_log2FC, n = 10)
write.csv(top10,file="./pyscenic/Top10 TF.csv")
mean_gene_exp <- AverageExpression(pbmc,
                                   features = top10$gene,assays = "RNA",
                                   group.by = 'celltypes') %>%data.frame() %>%as.matrix()
# add colnames
colnames(mean_gene_exp) <-str_sub(colnames(mean_gene_exp),5)
colnames(mean_gene_exp)=gsub("\\.\\.", "_", colnames(mean_gene_exp))
mean_gene_exp=as.matrix(mean_gene_exp)
# Z-score
mean_gene_matrix <- as.matrix(mean_gene_exp)
a=t(mean_gene_matrix)
scaled_matrix <- scale(a, center = TRUE, scale = TRUE)
htdf <- t(scaled_matrix)
# color
col_fun =colorRamp2::colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
# top annotation
anno_col <- color_qualitative2[1:length(colnames(htdf))]
names(anno_col) <- as.character(colnames(htdf))
column_ha = HeatmapAnnotation(cluster = colnames(htdf), col =list(cluster = anno_col))
pdf(file="./pyscenic/Top10_TF_zxj.pdf",width = length(anno_col)/1.5,height =length(anno_col)*1.2)
print(Heatmap(htdf,name = "Z-score",top_annotation = column_ha,
              cluster_columns = F,cluster_rows = F,
              row_title = "Top 10 Marker genes",
              row_names_gp = gpar(fontface = 'italic',fontsize = 10),
              row_names_side = 'left',
              border = T,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_side = 'top',
              column_names_rot = 45,
              col = col_fun))
dev.off()
genename=top10 %>% top_n(w = avg_log2FC, n = 3)
genename=genename$gene
genes=NULL
if (is.null(genes)) {
  genes <- genename
}
for (gene in unique(genes)){
  print(plot_featureplot(pbmc,features =gene))
  ggsave(paste(gene,".pdf"),width = 5,height = 5,path = "./pyscenic/")
}
plot_dotplot(pbmc,features= unique(genes))
ggsave("./pyscenic/Dotplot.pdf")


####reference data####
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl



https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
https://github.com/aertslab/pySCENIC/tree/master/resources/mm_mgi_tfs.txt

