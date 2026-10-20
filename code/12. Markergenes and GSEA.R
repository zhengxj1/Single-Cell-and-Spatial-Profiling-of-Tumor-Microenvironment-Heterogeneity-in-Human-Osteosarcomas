####All markers and clusters and DEG and GO/KEGG####
####VolcanoPlot####
allmarkers=FindAllMarkers(scRNA_harmony, only.pos = F, min.pct = 0.25)
VolcanoPlot(allmarkers,"Overall")
save(allmarkers,file="./Allmakers and GSEA/VolcanoPlot/allmarkers_PN.rda")
save(scRNA_harmony,file="scRNA_harmony.rda")
library(scRNAtoolVis)
jjVolcano(diffData = allmarkers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          #aesCol = c('purple','orange'), #设置点的颜色
          tile.col = Cat_Col, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
)
ggsave("jjVolcano.pdf",path = "./Allmakers and GSEA/VolcanoPlot")
#同celltype不同group对比
for (cluster in unique(scRNA_harmony$celltypes)) {
  sc=SplitObject(scRNA_harmony,split.by = "celltypes")
  subset=sc[[cluster]]
  Idents(subset)=subset$group
  deg=FindAllMarkers(subset, only.pos = F, min.pct = 0.25)
  VolcanoPlot(deg,project = cluster)
}


####不同亚群不同处理差异基因pseudobulk+DESeq2####
padjcutoff = 0.05 # 自定义
Log2FCcutoff= 0.5 
for (cluster in unique(scRNA_harmony$celltypes)) {
  sc=SplitObject(scRNA_harmony,split.by = "maintype")
  subset=sc[[cluster]]
  Idents(subset)=subset$group
  bs = split(colnames(subset),subset$orig.ident)
  ct = do.call(
    cbind,lapply(names(bs), function(x){ 
      # x=names(bs)[[1]]
      kp =colnames(subset) %in% bs[[x]]
      rowSums( as.matrix(subset@assays$RNA@counts[, kp]  ))
    })
  )
  colnames(ct) <- names(bs)
  phe = unique(subset@meta.data[,c('orig.ident','group')])#样本信息，分组信息，根据自己的metadata进行修改
  phe
  group_list = phe[match(names(bs),phe$orig.ident),'group']
  table(group_list)   
  exprSet = ct
  dim(exprSet)
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 1),]
  dim(exprSet)
  table(group_list)
  #第一步构建DEseq对象
  colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)  
  res <- results(dds,contrast = c('group_list', 'Post-chemotherapy', 'Pre-chemotherapy'))
  res<- as.data.frame(res)###提取差异分析结果（包含所有基因）
  ###或者这行命令res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res=subset(res,padj!="NA") ###去除padj为NA的行
  res=subset(res,log2FoldChange!="NA") ###去除log2FoldChange为NA的行
  res <- res[order(res$padj),]
  res1=res
  res1[which(res1$log2FoldChange > Log2FCcutoff & res1$padj <padjcutoff ),'sig'] <- 'up'
  res1[which(res1$log2FoldChange < (-Log2FCcutoff) & res1$padj <padjcutoff),'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) <= Log2FCcutoff | res1$padj >= padjcutoff),'sig'] <- 'none'
  DEG_DEseq2=res1
  plot_volcano=function(data,label=NA,sig=NA){
    data=na.omit(data)
    data=arrange(data,desc(-log10(padj)))
    top5=filter(data,sig=="up") %>% top_n(w=-log10(padj),n=10)
    bottom5=filter(data,sig=="down") %>% top_n(w=-log10(padj),n=10)
    label=c(row.names(top5),row.names(bottom5),label)
    data$label <- c(rep(NA,(nrow(data))))
    data$label[which(rownames(data)%in%label)]=row.names(data)[which(rownames(data)%in%label)]
    ggplot(data)+
      geom_point(alpha = 1, size = 2, aes(x=log2FoldChange, y=-log10(padj),color=sig,fill=sig))+
      scale_color_manual(values =c("#80B1D3","#bfc0c1","#FB8072"))+ #改变颜色
      scale_size_continuous(range = c(1, 3))+ #指定散点大小渐变
      geom_vline(xintercept = 0, linetype = "dashed") + #横向水平参考线
      geom_hline(yintercept =  -log10(0.05), linetype = "dashed")+ #纵向垂直参考线
      theme_bw()+ theme(panel.grid = element_blank())+ #改变主题
      geom_text_repel(aes(x=log2FoldChange, y=-log10(padj),label = label), size = 3.5,force = 20,
                      min.segment.length = 0,segment.colour = "grey40",max.overlaps = 100, key_glyph = draw_key_point)+
      labs(x = "log2FC",y = "-log10(P-value)",title = cluster)+
      theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),# 坐标轴标签和标题
            legend.text = element_text(size = 11),legend.title = element_text(size = 13),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold") ) # 图例标签和标题
  }
  plot_volcano(DEG_DEseq2)
  ggsave(file=paste0(cluster,'.pdf'),width = 4,height = 4)
}


####findmarkers####
dir.create('Allmakers and GSEA')
allmarkers=FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, test.use = "MAST",logfc.threshold = 0.25)
allmarkers=FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(allmarkers,file="./Allmakers and GSEA/allmarkers_Positive.rda")
top10 <- allmarkers%>% group_by(cluster) %>%  filter(avg_log2FC > 0& p_val_adj<0.05) %>%top_n(w =avg_log2FC, n =5)
write.csv(top10,file="./Allmakers and GSEA/Allmakers.csv")

####Heatmap####
library(ComplexHeatmap)
group.by="celltypes"
heatmap_plot_markers=function(scRNA_harmony,markers,group.by=group.by){
  dir.create("./Allmakers and GSEA")
  mean_gene_exp <- AverageExpression(scRNA_harmony,assays = 'SCT',slot = 'data',
                                   features = markers,
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
col_fun =colorRamp2::colorRamp2(c(min(htdf), median(htdf), max(htdf)), c("#0099CC", "white", "#CC0033"))
anno_col <-Cat_Col[1:length(colnames(htdf))]
names(anno_col) <- as.character(colnames(htdf))
column_ha = HeatmapAnnotation(cluster = colnames(htdf), col =list(cluster = anno_col))
pdf(file="./Allmakers and GSEA/Heatmap of averge allmarkers.pdf",width = length(anno_col),height =length(anno_col)*0.9)
print(Heatmap(htdf,name = "Z-score",top_annotation = column_ha,
              cluster_columns = F,cluster_rows = F,
              row_title = "Marker genes",
              row_names_gp = gpar(fontface = 'italic',fontsize = 10),
              row_names_side = 'left',
              border = T,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_side = 'top',
              column_names_rot = 45,
              col = col_fun))
dev.off()
}
top=top10$gene
heatmap_plot_tops=function(scRNA_harmony,group.by=group.by){
  mean_gene_exp <- AverageExpression(scRNA_harmony,assays = 'SCT',
                                     features = top,
                                     group.by = group.by,
                                     slot = 'data') %>%data.frame() %>%as.matrix()
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
  col_fun =colorRamp2::colorRamp2(c(min(htdf), median(htdf), max(htdf)), c("#0099CC", "white", "#CC0033"))
  anno_col <-Cat_Col[1:length(colnames(htdf))]
  names(anno_col) <- as.character(colnames(htdf))
  column_ha = HeatmapAnnotation(cluster = colnames(htdf), col =list(cluster = anno_col))
  pdf(file="./Allmakers and GSEA/Heatmap of tops.pdf",width = 5,height =length(anno_col))
  print(Heatmap(htdf,name = "Z-score",top_annotation = column_ha,
                cluster_columns = F,cluster_rows = F,
                row_title = "TOP DEGs",
                row_names_gp = gpar(fontface = 'italic',fontsize = 6),
                row_names_side = 'left',
                border = T,
                rect_gp = gpar(col = "white", lwd = 1),
                column_names_side = 'top',
                column_names_rot = 45,
                col = col_fun))
  dev.off()
}
heatmap_plot_markers(scRNA_harmony,markers,group.by)
heatmap_plot_tops(scRNA_harmony,group.by)
####GO KEGG####
perform_diff_analysis_and_plot <- function(allmarkers, OrgDb = "org.Hs.eg.db") {
  dir.create("./Allmakers and GSEA/")
  library(msigdbr)
  library(ComplexHeatmap)
  library(clusterProfiler) 
  top10 <- allmarkers%>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
  DEGlist <- split(top10, top10$cluster)
  for (i in 1:length(DEGlist)) {
    DEG <- DEGlist[[i]]
    GSEAgenelist <- bitr(DEG$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=OrgDb)
    genewithFC <- merge(DEG, GSEAgenelist, by.y="SYMBOL", by.x="gene")
    genewithFC1 <- genewithFC[,"avg_log2FC"]
    names(genewithFC1) <- genewithFC[,"ENTREZID"]
    
    genewithFC1 <- sort(genewithFC1, decreasing=TRUE)
    
    egoBP1 <- enrichGO(GSEAgenelist$ENTREZID, OrgDb = OrgDb, ont="BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable=TRUE)
    egoBP <- simplify(egoBP1, cutoff=0.7, by="p.adjust", select_fun=min)
    
    pdf(file=paste0("./Allmakers and GSEA/barplotBP_", names(DEGlist[i]), ".pdf"))
    print(barplot(egoBP, showCategory=20, title="EnrichmentGO_BP"))
    dev.off()
    egokegg <- enrichKEGG(GSEAgenelist$ENTREZID,organism = "hsa", keyType = "kegg",pvalueCutoff = 0.1,
)
    pdf(file=paste0("./Allmakers and GSEA/barplotkegg_", names(DEGlist[i]), ".pdf"))
    print(barplot(egokegg, showCategory=20, title="KEGG"))
    dev.off()
  }
}
# 用法示例
perform_diff_analysis_and_plot(allmarkers, OrgDb="org.Hs.eg.db")

####GSEA of Halmalk####

FindAllMarkers_heatmap=function(scRNA_harmony,species,category){
  dir.create("./Allmakers and GSEA/")
  library(msigdbr)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(singleseqgset)
  GSEA=function(scRNA_harmony,species=species,category=category){
    h.human <- msigdbr(species=species,category=category)
    save(h.human,file='h.human.rda')
    h.names <- unique(h.human$gs_name)
    h.sets <- vector("list",length=length(h.names))
    names(h.sets) <- h.names
    for (i in names(h.sets)) {
      h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
    }
    logfc.data <- logFC(cluster.ids=Idents(scRNA_harmony),expr.mat=scRNA_harmony@assays$RNA@data)
    names(logfc.data)
    gse.res <- wmw_gsea(expr.mat=scRNA_harmony@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
    res.stats <- gse.res[["GSEA_statistics"]]
    res.stats=as.matrix(res.stats)
    row.names(res.stats)=str_sub(row.names(res.stats),10)
    ann_colors = list("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","black","grey")
    return(res.stats) 
  }
  print("starting GSEA of Hallmark")
  res.stats=GSEA(scRNA_harmony,species=species,category=category)
  save(res.stats,file="./Allmakers and GSEA/res.stats.rda")
  anno_col <-Cat_Col[1:length(colnames(res.stats))]
  names(anno_col) <- as.character(colnames(res.stats))
  column_ha = HeatmapAnnotation(cluster = colnames(res.stats), col =list(cluster = anno_col))
  col = colorRamp2::colorRamp2(c(min(res.stats), median(res.stats), max(res.stats)), c("#0099CC", "white", "#CC0033"))
  pdf("./Allmakers and GSEA/GSEA_hallmalk.pdf",width = length(anno_col),height =length(anno_col)*2)
  pdf("./Allmakers and GSEA/GSEA_hallmalk.pdf")
  print(Heatmap(res.stats,column_title="GSEA",top_annotation = column_ha,
                cluster_columns = F,cluster_rows = T,show_column_names = FALSE,show_heatmap_legend=F,
                row_names_side = NULL,show_row_dend = F,
                border = T,
                column_names_side = 'top',col = col,
                column_names_rot = 45,
                rect_gp = gpar(col = "white", lwd = 1),
                heatmap_legend_param = list(
                  title="Pathway activity",
                  border='black')))
  dev.off()
  }
#species="Mus musculus"  "Homo sapiens" ###Homo sapiens or Mus musculus
#category="H"
species="Mus musculus"
species="Homo sapiens"
category="H"
FindAllMarkers_heatmap(scRNA_harmony,species=species,category=category)

####GSVA BP####
species=species ###Homo sapiens or Mus musculus
library(GSVA)
library(limma)
library(SeuratData)
DefaultAssay(scRNA_harmony)="SCT"
case=c("0","2","4")
con=c("1","5")
logFC=0
P=0.05
gsva_average=function(scRNA_harmony){
  matrix=AverageExpression(scRNA_harmony,slot = 'data',
                           group.by = 'celltypes') [[1]]
  matrix= matrix[rowSums(matrix)>0,]
  matrix=as.matrix(matrix)
  GO_bp = msigdbr(species =species,
                  category = "C5", #GO在C5
                  subcategory = "GO:BP") %>% 
    dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol，也可以选择ID，根据自己数据需求来，主要为了方便
  GO_bp_Set = GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
  gsva_scroe_averge <- gsva(expr = matrix, method="gsva",
                            gset.idx.list = GO_bp_Set,
                            kcdf="Gaussian", #查看帮助函数选择合适的kcdf方法 counts选 Poisson，data选Gaussian
                            parallel.sz = 5)
  write.table(gsva_scroe_averge, './Allmakers and GSEA/gsva_scroe_averge_BP.csv', row.names=T, col.names=NA, sep=",")
  return(gsva_scroe_averge)
}
gsva_H=function(scRNA_harmony,geneSets){
  matrix=AverageExpression(scRNA_harmony,slot = 'data',
                           group.by = 'celltypes') [[2]]
  matrix= matrix[rowSums(matrix)>0,]
  matrix=as.matrix(matrix)
  gsva_scroe_averge <- gsva(expr = matrix, method="gsva",
                            gset.idx.list = geneSets,
                            kcdf="Gaussian", #查看帮助函数选择合适的kcdf方法 counts选 Poisson，data选Gaussian
                            parallel.sz = 5)
  write.table(gsva_scroe_averge, './Allmakers and GSEA/gsva_scroe_averge_custome.csv', row.names=T, col.names=NA, sep=",")
  return(gsva_scroe_averge)
}
gsva_sc=function(scRNA_harmony){
  matrix=scRNA_harmony@assays$RNA@data
  matrix=as.matrix(matrix)
  GO_bp = msigdbr(species =species,
                  category = "C5", #GO在C5
                  subcategory = "GO:BP") %>% 
    dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol，也可以选择ID，根据自己数据需求来，主要为了方便
  GO_bp_Set = GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)
  gsva_scroe_sc <- gsva(expr = matrix, method="gsva",
                        gset.idx.list = GO_bp_Set,
                        kcdf="Gaussian", #查看帮助函数选择合适的kcdf方法 counts选 Poisson，data选Gaussian
                        parallel.sz = 5)
  write.table(gsva_scroe_sc, './Allmakers and GSEA/gsva_scroe_sc_BP.csv', row.names=T, col.names=NA, sep=",")
  return(gsva_scroe_averge_sc)
}
gsva_scroe_averge=gsva_average(scRNA_harmony)
gsva_scroe_averge=gsva_custome(scRNA_harmony,geneSets=geneSets)
gsva_sc(scRNA_harmony)
# color
gsva_scroe_averge=t(scale(t(gsva_scroe_averge)))
col_fun =colorRamp2::colorRamp2(c(min(gsva_scroe_averge), median(gsva_scroe_averge), max(gsva_scroe_averge)), c("#0099CC", "white", "#CC0033"))
# top annotation
anno_col <- Cat_Col[1:length(colnames(gsva_scroe_averge))]
names(anno_col) <- as.character(colnames(gsva_scroe_averge))
column_ha = HeatmapAnnotation(cluster = colnames(gsva_scroe_averge), col =list(cluster = anno_col))
pdf(file="./Allmakers and GSEA/GSVA.pdf",width = 10)
print(Heatmap(gsva_scroe_averge,column_title="Cell Cluster",top_annotation = column_ha,
              cluster_columns = F,cluster_rows = T,
              row_names_side = NULL,
              border = T,
              column_names_side = 'top',
              column_names_rot = 45,
              col = col_fun ),
              rect_gp = gpar(col = "white", lwd = 1)
)
dev.off()
pheatmap(gsva_scroe_averge, show_colnames = T, 
         scale = "row",angle_col = "90",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


group_info=scRNA_harmony@meta.data[,group_by]
group_list <- data.frame(celltype = colnames(gsva_scroe), group = group_info)
group_list_use=dplyr::filter(group_list,group%in%c(con,case))
group_list_use$group <- ifelse(group_list_use$group %in% con, "con", ifelse(group_list_use$group %in% case, "case",NA ))
group_list_use=dplyr::arrange(group_list_use,group)
design <- model.matrix(~0+factor(group_list_use$group))
colnames(design)=levels(factor(c("case","con")))
rownames(design)=group_list_use$celltype

contrast.matrix<-makeContrasts(contrasts = "case-con",levels = design)
gsva_scroe2=gsva_scroe[,row.names(design)]
fit <- lmFit(gsva_scroe2,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
diff$group <- ifelse( diff$logFC > logFC & diff$P.Value < P ,"up" ,
                      ifelse(diff$logFC < -logFC & diff$P.Value < P ,"down","noSig")
)
write.csv(diff,file = "./Allmakers and GSEA/gsva diff.csv")
diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>% #注释掉该行即可
  arrange(t) %>% 
  rownames_to_column("ID")
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))
ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x =0,nudge_y =0,hjust =diff2$hjust,
            size = 3)+ #设置字体大小
  labs(x = "GO BP", #设置标题 和 坐标
       y=paste0("t value of GSVA score\n","Cluster(0,2,4)-Cluster(1,5)"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none" #去掉legend
  )

ggsave("./Allmakers and GSEA/gsva results.tiff",width = 12,height = 5)


####GSVA REACTOME####

h_gene_sets = msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
Reactome_pathwaysH = split(x = h_gene_sets$human_gene_symbol, f = h_gene_sets$gs_name)
matrix=as.matrix(scRNA_harmony@assays$RNA@data)
row.names(matrix)=str_to_upper(row.names(matrix))
REACTOME_gsva_res <- gsva(matrix, gset.idx.list = Reactome_pathwaysH, method="gsva",
                          min.sz=10, max.sz=500)
pheatmap(REACTOME_gsva_res[1:20,], show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
write.table(gsva_scroe, './GSVAresults/gsva_scroe_REACTOME.csv', row.names=T, col.names=NA, sep=",")
f=factor(scRNA_harmony$celltypes)
design <- model.matrix(~0+f)
names=str_sub(colnames(design),2)
colnames(design)=str_c("Cluster",1:length(colnames(design)))
REACTOME_gsva.fit <- lmFit(REACTOME_gsva_res, design)
REACTOME_gsva.fit <- eBayes(REACTOME_gsva.fit)
groups=colnames(design)
comparisons <- expand.grid(groups = groups, groups2 = groups)
contrasts <- apply(comparisons, 1, function(x) paste(x[2], "vs", x[1]))

contrast.matrix <- makeContrasts(Cluster1.v.ALL = Cluster1-(Cluster2+Cluster3)/2,
                                 Cluster2.v.ALL = Cluster2-(Cluster1+Cluster3)/2,
                                 Cluster3.v.ALL = Cluster3-(Cluster1+Cluster2)/2,
                                 Cluster1.v.Cluster2 = Cluster1 - Cluster2,
                                 Cluster1.v.Cluster3 = Cluster1-Cluster3,
                                 Cluster2.v.Cluster3 = Cluster1-Cluster3,
                                 levels=design)
REACTOME_gsva.fit2 <- contrasts.fit(REACTOME_gsva.fit, contrast.matrix)
REACTOME_gsva.fit2 <- eBayes(REACTOME_gsva.fit2)
res2 <- decideTests(REACTOME_gsva.fit, p.value=0.01, method = "global")
summary(res2)
REACTOME_tt.res <- list(
  Cluster1 = topTable(REACTOME_gsva.fit2, coef = "Cluster1.v.ALL", n=1270),
  Cluster2 = topTable(REACTOME_gsva.fit2, coef = "Cluster2.v.ALL", n=1270),
  Cluster3 = topTable(REACTOME_gsva.fit2, coef = "Cluster3.v.ALL", n=1270)
)
REACTOME_tt.res_df <- do.call(rbind, REACTOME_tt.res)
REACTOME_tt.res_df$Cluster <- do.call(rbind, strsplit(rownames(REACTOME_tt.res_df), ".", 2))[,1]
REACTOME_tt.res_df$Pathway <- do.call(rbind, strsplit(rownames(REACTOME_tt.res_df), ".", 2))[,2]
REACTOME_tt.res_df$Pathway.label <- gsub("REACTOME_", "", REACTOME_tt.res_df$Pathway)
REACTOME_tt.res_df$Pathway.label <- gsub("_", " ", REACTOME_tt.res_df$Pathway.label)
REACTOME_tt.res_subset <- REACTOME_tt.res_df[REACTOME_tt.res_df$logFC > 0, ] %>%
  group_by(Cluster) %>%
  top_n(wt = -adj.P.Val, n = 5)
REACTOME_tt.res_subset <- REACTOME_tt.res_subset[-10, ]
REACTOME_tt.res_subset <- rbind(REACTOME_tt.res_subset, REACTOME_tt.res_df["Alveolar.REACTOME_TRP_CHANNELS", ])
REACTOME_tt.res_subset %>%
  ggplot(aes(x = reorder(Pathway.label, logFC), y = logFC, label =  paste0("adj.P=", signif(adj.P.Val, 3)),
             fill = Cluster)) +
  theme_pubr(base_size = 15)+
  theme(axis.title.y = element_blank(),legend.title = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "top", strip.text = element_blank(), legend.key.size = unit(10, "pt"))+
  geom_bar(stat = "identity", alpha = 0.75)+
  facet_grid(Cluster~., scales = "free", space = "free") +
  geom_text(hjust = 1.1, colour = "white", size = 3) +
  scale_fill_manual(values = color_qualitative2) +
  coord_flip()


####heatmap of marker gene 每个细胞####
scale_data=ScaleData(scRNA_harmony,features = top10$gene,assays = "RNA")
DoHeatmap(scale_data, features = rev(top10$gene),
          assay = "RNA",
          angle = 45,
          
          
          label = T, disp.max = 2) &
  theme(axis.text.y = element_text(face = "italic", size = 10),legend.position = "none") 

####GSVA####
scRNAseq_data=scRNA_harmony
gene_sets=CAF
output_filename
library(GSVA)
library(dplyr)
library(pheatmap)

calculate_gsva_and_plot_with_metadata <- function(scRNAseq_data, gene_sets, output_filename) {
  # 1. 检查元数据中的celltypes列是否存在
  if (!"celltypes" %in% colnames(scRNAseq_data@meta.data)) {
    stop("celltypes not found in the metadata of the scRNAseq_data object.")
  }
  
  # 2. 计算每个细胞的基因集GSVA评分
  gsva_scores <- gsva(
    as.matrix(scRNAseq_data@assays$RNA@data),
    gene_sets,
    min.sz = 5,parallel.sz = 8,
    kcdf = "Gaussian",
    method = "gsva",
    mx.diff = TRUE,
    verbose = FALSE
  )
  
  # 3. 转换为数据框并添加亚群信息
  gsva_scores_df <- as.data.frame(t(gsva_scores))
  cell_type_vector <- scRNAseq_data@meta.data$celltypes
  gsva_scores_df$cell_type <- cell_type_vector
  
  # 4. 计算亚群内的GSVA评分平均值
  gsva_avg_by_subtype <- gsva_scores_df %>%
    group_by(cell_type) %>%summarise(across(everything(), mean, na.rm = TRUE))
  gsva_avg_by_subtype <- as.data.frame(gsva_avg_by_subtype)
  row.names(gsva_avg_by_subtype) <- gsva_avg_by_subtype$cell_type
  gsva_avg_by_subtype <- gsva_avg_by_subtype[,-1]
  gsva_avg_by_subtype =scale(gsva_avg_by_subtype)
  gsva=t(gsva_avg_by_subtype)
  # 5. 热图展示亚群间的差异
  color_mapping <- colorRamp2::colorRamp2(
    c(min(gsva), 
      median(apply(gsva, 1, median, na.rm = TRUE)), 
      max(gsva)
    ), 
    c("#0099CC", "white", "#CC0033")
  )
  anno_col <-Cat_Col[1:ncol(gsva)]
  names(anno_col) <- as.character(colnames(gsva))
  column_ha = HeatmapAnnotation(cluster = colnames(gsva), col =list(cluster = anno_col))
  pdf(file="./Allmakers and GSEA/Heatmap of tops.pdf",width = length(anno_col)/2,height =length(anno_col)*0.7)
  print(Heatmap(gsva,name = "Z-score",top_annotation = column_ha,
                cluster_columns = F,cluster_rows = F,
                row_title = "TOP DEGs",
                row_names_gp = gpar(fontface = 'italic',fontsize = 6),
                row_names_side = 'left',
                border = T,
                rect_gp = gpar(col = "white", lwd = 1),
                column_names_side = 'top',
                column_names_rot = 45,
                col = color_mapping))
  dev.off()
}

expr <- AverageExpression(scRNA_harmony, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
head(expr)
gsva.kegg <- gsva(expr, CAF, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=8)
head(gsva.kegg)
pheatmap(gsva.kegg, show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))























