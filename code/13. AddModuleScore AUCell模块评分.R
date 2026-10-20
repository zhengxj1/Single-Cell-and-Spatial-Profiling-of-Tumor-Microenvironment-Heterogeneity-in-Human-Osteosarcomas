####AddModuleScore####
library(Seurat)
library(tidyverse)
library(ggpubr)


genelist=list(geneset1,geneset2,YAP)
names=names(geneSets)

groups=unique(Idents(scRNA_harmony))
comparisons <- combn(groups, 2, simplify = FALSE)

scRNA_harmony@meta.data=dplyr::select(scRNA_harmony@meta.data,-matches("Geneset"))
scRNA_harmony@meta.data=dplyr::select(scRNA_harmony@meta.data,-matches("geneset"))
scRNA_harmony$seurat=Idents(scRNA_harmony)
scRNA_harmony<- AddModuleScore(
  scRNA_harmony,
  features = geneSets,
  assay = "RNA",name =paste0("Geneset_",names))
DoHeatmap(scRNA_harmony,features = names)

cluster_cols <- grep("^Geneset_", names(scRNA_harmony@meta.data))
names(scRNA_harmony@meta.data)[cluster_cols] <- names
df=scRNA_harmony@meta.data 
####AUCell####
source("F:/BaiduSyncdisk/R code/scRNA_scripts/genesets.R")

library(AUCell)
h.human <- msigdbr(species=species,category=category)
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}
geneSets=h.sets["HALLMARK_HYPOXIA"]

aucell=function(scRNA_harmony,geneSets){
  exprMatrix <- scRNA_harmony@assays$RNA$data
  cells_rankings <- AUCell_buildRankings(exprMatrix,splitByBlocks=T) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
if (!dir.exists("./AUCell")) {
  dir.create("./AUCell")
}
save(cells_AUC,file="./AUCell/cells_AUC.rda")

for(name in names(geneSets)){
  aucs <- as.numeric(getAUC(cells_AUC)[name, ])
scRNA_harmony$AUCell=aucs 
plot=scCustomize::FeaturePlot_scCustom(scRNA_harmony, features ="AUCell", slot = "counts",split.by = "group",
                     colors_use = color_gradient4)
plot1=plot[[1]]+theme_CNS_blank()
plot2=plot[[2]]+theme_CNS_blank()
plot=plot1 /plot2
ggsave(paste(name,".pdf"),width = 5,height =8,path = "./AUCell")
##boxplot
df=scRNA_harmony@meta.data
normalize <- function(x){return((x-min(x))/(max(x)-min(x)))}
q_low <- quantile(df$AUCell, 0.01)  # 底部分位数

q_high <- quantile(df$AUCell, 0.99) 
df[, !duplicated(colnames(df))]%>%
  ggplot(aes(x = celltypes, y = normalize(AUCell),
             fill =celltypes)) +geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  scale_fill_manual(values = color_qualitative2)+ylab("Score")+
  theme_CNS_blank()+theme_coord_flip()+coord_flip()
ggsave(paste(name,"_celltype_boxplot.pdf"),path = "./AUCell")

df[, !duplicated(colnames(df))]%>%
  ggplot(aes(x = celltypes, y = normalize(AUCell),
             fill =group)) +geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  scale_fill_manual(values = color_qualitative3)+ ylab("Score")+
  ggsignif::geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = wilcox.test)+
  theme_CNS_blank()+theme_coord_flip()+coord_flip()
ggsave(paste(name,"_group_boxplot.pdf"),path = "./AUCell")
df[, !duplicated(colnames(df))]%>%
  ggplot(aes(x = group, y = normalize(AUCell),
             fill =group)) +geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  scale_fill_manual(values = color_qualitative3)+ ylab("Score")+
  theme_CNS_blank()+theme_coord_flip()+coord_flip()
ggsave(paste(name,"_group_group_boxplot.pdf"),path = "./AUCell")
}
}
aucell(scRNA_harmony,geneSets)
aucell(scRNA_harmony,result)




set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
aucMat <- getAUC(cells_AUC)
aucs_t <- data.frame(t(aucMat))
aucs_t$celltypes <- scRNA_harmony@meta.data[rownames(aucs_t), "celltypes"]
aucs_t <- aucs_t %>% group_by(celltypes) %>% dplyr::summarise_each(funs = mean)


write.csv(aucs_t ,"./outputs/AUC.csv")
test1 <- aucs_t
test1$celltypes<- factor(test1$celltypes)
test1$Mc <- ifelse(grepl("^(0|2|4)$", test1$celltypes),"down",ifelse(grepl("^(1|5)$", test1$celltypes),"up","others"))
col = c('#003399', '#999999')
names(col) = c('CD56brightCD16lo', 'CD56dimCD16hi')

g1 <- test1  %>% ggplot() + geom_bar( aes(x =celltypes ,y=cytotoxicity_genes, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(cytotoxicity_genes)),linetype="dashed") + 
  ylab("Cytotoxicity score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12,angle = 90),axis.title.y = element_text(angle = 90, vjust = 1), 
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = color_qualitative5) +
  scale_colour_manual(values = color_qualitative5) 
g1
g2 <- test1  %>% ggplot() + geom_bar( aes(x =celltypes ,y=inflammatory_genes, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(inflammatory_genes)),linetype="dashed") + 
  ylab("Inflammatory score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12,angle = 90),axis.title.y = element_text(angle = 90, vjust = 1), 
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = color_qualitative5) +
  scale_colour_manual(values = color_qualitative5) 
g2
g3 <- test1  %>% ggplot() + geom_bar( aes(x =celltypes ,y=stress_genes, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(stress_genes)),linetype="dashed") + 
  ylab("Stress score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12,angle = 90),axis.title.y = element_text(angle = 90, vjust = 1), 
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = color_qualitative5) +
  scale_colour_manual(values = color_qualitative5) 
g3
g <- g1 / g2 / g3
g
save(g,file="./outputs/Score_AUCell.rda")
ggsave("./outputs/Score_AUCell.pdf", g, height = 5, width =8)
####可视化####

normalize <- function(x){return((x-min(x))/(max(x)-min(x)))}
q_low <- quantile(df$YAP, 0.01)  # 底部分位数
q_high <- quantile(df$YAP, 0.99) 
df[, !duplicated(colnames(df))]%>%
  ggplot(aes(x = celltypes, y = normalize(YAP),
             fill =orig.ident)) +geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  scale_fill_manual(values = color_qualitative3)+
  theme_pubr(base_size = 15) +
  theme(legend.position = "right", axis.title.x = element_blank())+ ylab("Score")+ coord_cartesian(ylim = c(0.1, 0.6))

  ggsave("YAP-all_cluster.tiff",width = 10,height = 5,dpi = 300,path = "./outputs")
  df[, !duplicated(colnames(df))]%>%
    ggplot(aes(x = orig.ident, y = normalize(YAP),
               fill =orig.ident)) +geom_boxplot(outlier.shape = NA, alpha = 0.75)+
    scale_fill_manual(values = color_qualitative3)+
    theme_pubr(base_size = 15) +
    theme(legend.position = "right", axis.title.x = element_blank())+ ylab("Score")+ coord_cartesian(ylim = c(0.1, 0.4))
  
  ggsave("YAP-all_id.tiff",width = 5,height = 5,dpi = 300,path = "./outputs")
FeaturePlot(scRNA_harmony,features = "YAP",cols =palette_scp(n=40,type = "continious",palette = "Spectral") )



groups <- unique(scRNA_harmony$celltypes)
comparisons <- expand.grid(groups, groups)
comparisons <- comparisons[comparisons$Var1 != comparisons$Var2, ]
comparisons <- split(comparisons, 1:nrow(comparisons))
#箱形图显示编码每个亚群基底膜组分的基因的样品水平表达
scRNA_harmony@meta.data %>%
  ggplot(aes(x = groups, y = YAP,
             fill = groups)) +
  theme_pubr(base_size = 15) +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_compare_means(comparisons = comparisons ,
                     size = 5,
                     hide.ns = T,label.y = c(0.65,0.8)) +
  scale_fill_manual(values = color_qualitative2) +
  rotate_x_text(angle = 45) +
  ylab("Basement Membranes") +
  guides(fill = "none") +  expand_limits(y = 0.9)
scRNA_harmony@meta.data %>%
  ggplot(aes(x = orig.ident, y = geneset1,
             fill = orig.ident)) +
  theme_pubr(base_size = 15) +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_compare_means(comparisons = comparisons ,
                     size = 5,
                     hide.ns = T,label.y = c(0.65,0.8)) +
  scale_fill_manual(values = color_qualitative2) +
  rotate_x_text(angle = 45) +
  ylab("Basement Membranes") +
  guides(fill = "none") +  expand_limits(y = 0.9)

####
comparisons <- list(c("C1", "D3"))
scRNA_harmony@meta.data %>%
  ggplot(aes(x = samples, y = geneset1)) +
  theme_pubr(base_size = 15) +
  facet_wrap(~celltypes,ncol = 4) +
  geom_boxplot(aes(fill = celltypes), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_compare_means(comparisons = comparisons,
                     size = 5,method = "wilcox.test",
                     label.y = 0.9) +
  scale_fill_manual(values = color_qualitative2) +
  rotate_x_text(angle = 45) +
  ylab("Interstitial Collagens") +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") + expand_limits(y = 1.1)

FeaturePlot(scRNA_harmony,features = "geneset1")

####Heatmap####
metadata=scRNA_harmony@meta.data
Idents(scRNA_harmony)=scRNA_harmony$secondtype
exprMatrix <- scRNA_harmony@assays$RNA$data
cells_rankings <- AUCell_buildRankings(exprMatrix,splitByBlocks=T) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
aucMat <- getAUC(cells_AUC)
pbmc <- CreateSeuratObject(counts =aucMat)
pbmc=AddMetaData(pbmc,metadata =metadata)
group.by="secondtype" 
mean_gene_exp <- AverageExpression(pbmc,
                                   group.by = group.by) %>%data.frame() %>%as.matrix()
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
anno_col <- Cat_Col[1:length(colnames(htdf))]
names(anno_col) <- as.character(colnames(htdf))
column_ha = HeatmapAnnotation(cluster = colnames(htdf), col =list(cluster = anno_col))
pdf(file="Heatmap_signature.pdf",width = 6,height = 5)
print(Heatmap(htdf,name = "Z-score",top_annotation = column_ha,
              cluster_columns = F,cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              row_names_side = 'right',
              border = T,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_side = 'top',
              column_names_rot = 45,
              col = col_fun))
dev.off()


















