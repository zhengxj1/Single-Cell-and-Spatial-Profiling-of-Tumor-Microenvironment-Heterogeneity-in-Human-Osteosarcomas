##cell proportion
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(ggbeeswarm)
####各样本比例图####
scRNA_harmony$celltypes=Idents(scRNA_harmony)
scedata=scRNA_harmony

cell_ratio_wide <- Cellratio %>%
  pivot_wider(
    names_from = Var2,  # 样本名称作为列名
    values_from = Freq   # 比例为值
  ) %>%
  column_to_rownames("Var1")  # 将细胞类型设为行名
# 提取目标亚群的比例数据
cellA <- as.numeric(cell_ratio_wide["CellTypeA", ])
cellB <- as.numeric(cell_ratio_wide["CellTypeB", ])

# 计算皮尔森相关系数及p值
cor_result <- cor.test(cellA, cellB, method = "pearson")
cat("Pearson r =", cor_result$estimate, "p-value =", cor_result$p.value)
ggplot(plot_data, aes(x = CellTypeA, y = CellTypeB)) +
  geom_point(color = "darkgreen", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # 添加回归线
  labs(
    title = paste0("Pearson Correlation: r = ", round(cor_result$estimate, 2)),
    x = "CellTypeA Proportion",
    y = "CellTypeB Proportion"
  ) +
  theme_classic()

proportion=function(scedata){
  dir.create("./outputs")
Cellratio <- prop.table(table(scedata$celltypes, scedata$orig.ident), margin = 2) %>% as.data.frame()#计算各组样本不同细胞群比例
colourCount = length(unique( as.character(Cellratio$Var1)))
ggplot(Cellratio,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7)+ scale_y_continuous(expand=c(0,0))+labs(y="Percentage(%)")+
  scale_x_discrete()+
  scale_fill_manual(values=Cat_Col,name="Celltype")+theme_CNS_blank()+xlab(NULL)+ylab(NULL)+
  coord_flip()
ggsave("proportion-ori.pdf",width=4,dpi=400,path = "./outputs/")
tab=table(scedata$celltypes)
tabl=data.frame(count=tab,
               prop=round(tab/sum(tab)*100))
data=tabl
data$ang=seq(from=(360/nrow(data))/1.5,
             to = (1.5* (360/nrow(data))) - 360,
             length.out = nrow(data)) + 80
data$ang[which(data$ang < -90)] <- (180+ data$ang)[which(data$ang < -90)]
data$x=seq(1,160,length.out=nrow(data))
data$xmin=data$x-0.4
data$xmax=data$x+0.4
ggplot(data, aes(x = 3,
               y = prop.Freq, 
               fill = prop.Var1)) +
  geom_col(width = 1.5,
           color = 'white') +
  coord_polar(theta = "y") +
  scale_fill_manual(values = Cat_Col) +xlim(c(0.2, 3.8)) +
  geom_text(aes(label = paste(prop.Freq,'%')),
          position = position_stack(vjust = 0.5),
          size = 4)+ labs(title = paste("Total = ",sum(data$count.Freq)),fill="Celltype")+theme_CNS_blank()+xlab(NULL)+ylab(NULL)+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  ggsave("celltype prop.pdf",width = 5,height = 5,path = "./outputs/")
write.csv(tabl,file = "./outputs/cell number.csv")
Cellratio2 <- prop.table(table(scedata$celltypes, scedata$group), margin = 2) %>% as.data.frame()
ggplot(Cellratio2,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5)+ 
  theme_classic() +
  labs(y = 'Ratio')+
  scale_fill_manual(name="Celltype",values=Cat_Col)+theme_CNS_blank()+xlab(NULL)+theme_angle()
ggsave("proportion-Group_CELLTYPE_Location.pdf",width = 5,height=5,dpi=400,path = "outputs/")

Cellratio3 = prop.table(table(scedata$orig.ident,scedata$celltypes), margin = 2) %>% as.data.frame()
ggplot(Cellratio3,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5)+ 
  labs(y = 'Ratio',x=NULL)+
  scale_fill_manual(name=NULL,values=color_patient) +theme_CNS_blank()+theme_angle()
ggsave("proportion-patient.pdf",width = length(unique(scedata$celltypes)),height=5,dpi=400,path = "outputs/")

Cellratio4 <- prop.table(table(scedata$group,scedata$celltypes), margin = 2) %>% as.data.frame()
  ggplot(Cellratio4,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5)+  theme_CNS_blank()+xlab(NULL)+
  scale_fill_manual(name="Group",values=color_group2)   + labs(y = 'Ratio',X=NULL)+theme_angle()
ggsave("proportion-Group2.pdf",width = length(unique(scedata$celltypes)),height=5,path = "outputs/")

Cellratio5 <- prop.table(table(scedata$celltypes, scedata$Location), margin = 2) %>% as.data.frame()
ggplot(Cellratio5,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5)+ 
  theme_classic() +
  labs(y = 'Ratio')+
  scale_fill_manual(name="Celltype",values=Cat_Col)+theme_CNS_blank()+xlab(NULL)+theme_angle()
ggsave("proportion-LOcat_CELLTYPE.pdf",width = 5,height=5,dpi=400,path = "outputs/")
Cellratio6 <- prop.table(table(scedata$Location,scedata$celltypes), margin = 2) %>% as.data.frame()
ggplot(Cellratio6,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5)+  theme_CNS_blank()+xlab(NULL)+
  scale_fill_manual(name="Group",values=color_group)   + labs(y = 'Ratio',X=NULL)+theme_angle()
ggsave("proportion-LOC_BOX.pdf",width = length(unique(scedata$celltypes)),height=5,path = "outputs/")
####分组比例图#####
compaired <- list(c("Pre-chemotherapy","Post-chemotherapy"))
scedata=scedata
Cellratio <- prop.table(table(scedata$celltypes, scedata$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据

rownames(cellper) <- cellper[,1]

meta<- FetchData(scedata, vars = c("orig.ident", "Treatment",'Location'))
unique_samples <- unique(meta)

merged_data <- merge(cellper, unique(meta), by.y = "orig.ident",by.x = "Var2")
data_long <- merged_data %>%
  gather(key = "Cell_Type", value = "Expression", -Var2, -Treatment,-Location)

p <- ggplot(data_long, aes(x = Cell_Type, y = Expression, fill = Treatment)) +
  geom_boxplot() +geom_point( position=position_dodge(width=0.75), cex=3, size=3) +  # 点，与箱体图对齐

  ylab("Percent") +
  labs(title = NULL) +
  theme_CNS_blank()+theme_angle()+scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                                                 "#355783", "#ef9a80")) +  xlab("") 
p + stat_compare_means(aes(group = Treatment), label="p.signif",method = "t.test")
ggsave("OC自己数据化疗.pdf")
p <- ggplot(data_long, aes(x = Cell_Type, y = Expression, fill = Location)) +
  geom_boxplot() +geom_point( position=position_dodge(width=0.75), cex=3, size=3) +  # 点，与箱体图对齐
  
  ylab("Percent") +
  labs(title = NULL) +
  theme_CNS_blank()+theme_angle()+scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                                                "#355783", "#ef9a80")) +  xlab("") 
p + stat_compare_means(aes(group = Location), label="p.signif",method = "kruskal.test")
ggsave("自己数据部位.pdf")
# 显示图形
print(p)
ggsave("Location.pdf")


samples <- data.frame(scedata@meta.data)
samples=as.data.frame(distinct(samples,orig.ident,.keep_all = T))
cellper=merge(cellper,samples,by.x="Var2",by.y="orig.ident")
full=cellper
clusters=unique(scedata$celltypes)
pplist = list()
for(cluster in unique(clusters)){
  cellper_  = cellper %>% dplyr::select(one_of(c("Var2",'group',cluster)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  colors <- color_group
  pplist[[cluster]] = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_boxplot(aes(color=group))+
    geom_beeswarm(aes(color=group),cex = 3,size=3, priority = "descending") +
    scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                  "#355783", "#ef9a80"))+ xlab("")+ylab("Percent")+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = wilcox.test)+
    labs(title = paste0("",cluster)) +
    theme_CNS_blank()+theme_angle()
  }
for(i in 1:length(unique(clusters))){
plot_grid(pplist[[i]])
ggsave(filename=paste0(paste0("cell proportion_Loc_",unique(clusters)[i]),".pdf"),width = 4,height=5,dpi=300,path = "outputs/")
}
full=full[,c("Var2",levels(Cellratio$Var1),"group")]
fulllong=pivot_longer(full,
  cols = levels(Cellratio$Var1), # 选择要长格式化的列
  names_to = "celltypes",          # 新的列名，用于存储原始列名
  values_to = "Value"              # 新的列名，用于存储原始列的值
) %>% group_by(celltypes)


ggplot(fulllong, aes(x=celltypes, y=Value,color=group)) + 
  geom_boxplot(position=position_dodge(width=0.75)) +  # 箱体图，宽度调整以对齐
  geom_point( position=position_dodge(width=0.75), cex=3, size=3) +  # 点，与箱体图对齐
  scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                "#355783", "#ef9a80")) + 
  xlab("") + 
  ylab("Percent") +
  ggsignif::geom_signif(comparisons = compaired, step_increase = 0.1,
                        map_signif_level = T, test = wilcox.test) +  # 显示P值
  labs(title = NULL) +
  theme_CNS_blank()+theme_angle()
ggsave("proportion-All_grouphis.pdf",height=5,path = "outputs/")
####亚群相关性热图####
library(ggcorrplot)
dir.create("./outputs/correlation")
Cellratio <- prop.table(table(scedata$celltypes, scedata$orig.ident), margin = 2) %>% as.data.frame()
Cellratio=pivot_wider(Cellratio,id_cols = "Var2",names_from = "Var1",values_from = "Freq") %>% as.data.frame()
row.names(Cellratio)=as.character(Cellratio[,1])
colnames(Cellratio)=as.character(colnames(Cellratio))
data=Cellratio[,-1]
new_colnames <- colnames(data)
for (i in seq_along(new_colnames)) {
  if (nchar(new_colnames[i]) == 1) {
    new_colnames[i] <- paste("Cluster_", new_colnames[i], sep = "")
  }
}
colnames(data) <- new_colnames
methods <- c("pearson", "spearman", "kendall")
results <- list()
for (method in methods) {
  cor_matrix <- cor(data, method = method)
  p_value_matrix <- cor_pmat(data, method = method)
  
  results[[method]] <- list(cor_matrix = cor_matrix, p_value_matrix = p_value_matrix)
}
for (method in methods) {
  cor_matrix <- results[[method]]$cor_matrix
  p_value_matrix <- results[[method]]$p_value_matrix
  ggcorrplot(cor_matrix,colors = c("#053061","#F5F6F6","#67001F"), lab_size = 2.5, #p.mat = p_value_matrix,
             hc.order = F,ggtheme=theme_classic(),
             type = "lower",
             lab = TRUE)+ theme_CNS_blank(base_size = 12)+
    labs(x=NULL,y=NULL,title = 'Correlation')+theme_angle()
  ggsave(paste("correlation_",method,".pdf"),width = 5,height=5,dpi=400,path = "outputs/correlation")
}
####亚群相似性图####
dir.create("./outputs/similarity")
av=AverageExpression(scedata,group.by = "celltypes",assays="RNA")
av <- av[[1]]
av=as.data.frame(av)
head(av)
#选出标准差最大的1000个基因
cg=names(tail(sort(apply(av,1,sd)),1000))
#“查看这1000个基因再各细胞群的表达矩阵
for (method in methods) {
  corr_mat <- cor(av[cg,], method = method)
  range_corr <- range(corr_mat)
  corr_mat_scaled <- matrix(nrow = nrow(corr_mat), ncol = ncol(corr_mat))
  colnames(corr_mat_scaled) <- colnames(corr_mat)
  row.names(corr_mat_scaled) <- row.names(corr_mat)
  
  for (i in 1:nrow(corr_mat)) {
    for (j in 1:ncol(corr_mat)) {
      corr_mat_scaled[i, j] <- (corr_mat[i, j] - range_corr[1]) / (range_corr[2] - range_corr[1]) * 2 - 1
    }
  }
  
  results[[method]] <- list(corr_mat_scaled, range_corr)
}
for (method in methods) {
  corr_mat_scaled<- results[[method]][[1]]
  range_corr=results[[method]][[1]]
ggcorrplot(corr_mat_scaled, colors = c("#053061","#F5F6F6","#67001F"),lab_size = 2.5,#p.mat = p_value_matrix,
           hc.order = F,ggtheme=theme_classic(),
           type = "lower",
           lab = TRUE)+ theme_CNS_blank()+
  labs(x=NULL,y=NULL,title = 'Similarity')+theme_angle()
  ggsave(paste("similarity_",method,".pdf"),width = 5,height=5,dpi=400,path = "outputs/similarity")
}
}
proportion(scRNA_harmony)

library(vcd)
library(gplots)
meta=scRNA_harmony@meta.data
balloonplot(table(meta$orig.ident,meta$SCT_snn_res.0.5))
mosaic(table(meta$orig.ident,meta$SCT_snn_res.0.5),shade=T,legend=T)

unique(Idents(scRNA_harmony))
scedata=OS
scedata$celltypes=Idents(scedata)


FeaturePlot(scedata,features = 'nFeature_RNA')

sed






