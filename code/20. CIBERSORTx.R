####CIBERSORTx####
####CIBERSORTx####
library(Seurat)
library(data.table)
library(tidyverse)
CIBERSORTx=function(seuratObj){
  #从Seurat对象中制作单细胞表达矩阵（行为基因，列为细胞）
  set.seed(123456)
  Idents(seuratObj)=seuratObj$celltypes
  seuratObj=subset(seuratObj, downsample = 400)
  print(dim(seuratObj))
  exp_matrix = seuratObj@assays$RNA@counts#表达矩阵
  metadata = seuratObj@meta.data#metadata
  identical(row.names(seuratObj@meta.data),colnames(exp_matrix))#TRUE
  exp_matrix = exp_matrix %>% as.data.frame() %>% 
    rownames_to_column(var = "GeneSymbol") %>%
    dplyr::select(GeneSymbol,everything())
  exp_matrix[1:3,1:12]
  colnames(exp_matrix)[-1] = as.character(seuratObj$celltypes)
  fwrite(exp_matrix, "./Step1.CibersortX_scRNAassay_secondtype.txt",sep = "\t",row.names =F)
}
CIBERSORTx(scRNA_harmony)

####boxplot####
dir.create("Cibersortx")
TME_data <- read.csv("CIBERSORTx_Results.csv",header = T,row.names = 1)
TME_data = TME_data[,c(1:c(ncol(TME_data)-3))]
data=TME_data
data$ratio=data$Osteoblastic-data$Chondroblastic
data$type <- ifelse(data$ratio > 0, "Osteoblastic", "Chondroblastic")
data=data[,c(1:31,33)]
data1=pivot_longer(data,cols = colnames(data)[1:31],values_to = "count")
library(ggpubr)


p <- ggplot(data1, aes(x=name, y=count, fill = type)) +
  geom_boxplot() +geom_point( position=position_dodge(width=0.75), cex=3, size=1) +  # 点，与箱体图对齐
  
  ylab("Percent") +
  labs(title = NULL) +
  theme_CNS_blank()+theme_angle()+scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781","#355783", "#ef9a80")) +  xlab("") 
p
p + stat_compare_means(aes(group = type), label="p.signif",method = "t.test")
ggsave('boxplot.pdf',width = 10,height = 5)


# 5.1 数据粗处理
TME_New = TME_data
TME_New$group <- "T"
TME_New$sample <- row.names(TME_New)
# 5.2 融合数据
TME_New = reshape2::melt(TME_New)
colnames(TME_New)=c("Group","Sample","Celltype","Composition") #设置行名
plot_order = TME_New %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)
plot_order

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)

ggplot(TME_New,aes(x=Celltype,y=Composition)) + #ggplot作图
  geom_boxplot(aes(color=Celltype))+
  geom_jitter(aes(color=Celltype),size=2,width = 0.3,height = 0.1,alpha=0.8,fill="transparent")+
  scale_color_manual(values = color_more)+ xlab("")+ylab("Percent")+labs(title = "TARGET") + theme_CNS_blank(base_size = 8)
ggsave("Cibersortx/boxplot.pdf",width=5,height=5)
####相关性####
library(ggcorrplot)
data=TME_data
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
dir.create( "correlation")
for (method in methods) {
  cor_matrix <- results[[method]]$cor_matrix
  p_value_matrix <- results[[method]]$p_value_matrix
  ggcorrplot(cor_matrix,colors = c("#053061","#F5F6F6","#67001F"), lab_size = 2.5, #p.mat = p_value_matrix,
             hc.order = F,ggtheme=theme_classic(),
             type = "lower",
             lab = TRUE)+ theme_CNS_blank(base_size = 12)+
    labs(x=NULL,y=NULL,title = 'Correlation')+theme_angle()
  ggsave(paste("correlation_",method,".pdf"),width = 20,height=20,dpi=400,path = "correlation")
}

data=pivot_longer(data,values_to = "value",cols=colnames(data)[1:32])
ggplot(data, aes(x=name, y=value,color=type)) + 
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

####生存####
library("survival")
Variables=colnames(TME_data)
rownames(sur)=str_c(rownames(sur),'-01R')
CoxData=merge(TME_data,sur,by = "row.names", all = FALSE)
CoxData=column_to_rownames(CoxData,var = 'Row.names')
OScoxSummary<-lapply(Variables,function(x){
  library("survival")
  print(x)
  OS<-na.omit(CoxData[,c("Time","Status",x)])
  res.cox<-coxph(Surv(Time,Status)~get(x),data=OS)
  coxSummary<-summary(res.cox)
  coxSummary
  OScoxSummary<-data.frame(id=x,
                             N=nrow(OS),
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR_95L=coxSummary$conf.int[,"lower .95"],
                             HR_95H=coxSummary$conf.int[,"upper .95"],
                             C_index=coxSummary$concordance["C"],
                             C_index_se=coxSummary$concordance["se(C)"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  return(OScoxSummary)
  
})
UnivariateOScoxSummary<-do.call(rbind,OScoxSummary)
UnivariateOScoxSummary
UnivariateOScoxSummary%>%
  mutate(Change=case_when(
    pvalue<0.05& HR>1~"Risk factor",
    pvalue<0.05& HR<1~"Protective factor",
    pvalue=0.05~"Not sig"
  ))->UnivariateOScoxSummary

p<-ggplot(UnivariateOScoxSummary, aes(x =HR, y = id)) +
  geom_errorbarh(aes(xmin = HR_95L, xmax = HR_95H,color=Change), height = 0.2,size=1) +
  geom_point(alpha=0.3,aes(size=pvalue)) +
  labs(x = "OR (95%CI)", y = "") +
  geom_vline(xintercept = 1.0,color="grey",linetype=2,linewidth=1)+
  scale_color_manual(values = c("#377eb8","#4daf4a","#e41a1c"))+
  theme_bw()+
  theme(axis.text = element_text(size = 12,color = "black"),
        panel.border = element_rect(linewidth = 1))

p

####corrplot####
library(corrplot)
corda.cell.type <- psych::corr.test(TME_data, TME_data, method = "sp")
r.df = as.matrix(corda.cell.type$r)
p.df = as.matrix(corda.cell.type$p.adj)

pdf(file="./Cibersortx/corrplot.pdf",width =3.5,height =3)
corrplot(as.matrix(r.df),
         method = "color",
         col=rev(COL2("RdBu", 200)),
         type = c("upper"),
         p.mat = p.df,#是否添加了P值，如果是NULL，则sig.level到pch.cex参数无效
         sig.level = 0.05,#P值检验水平设置
         insig = c("label_sig"),
         pch.cex = 1, 
         tl.col = "black",
         tl.cex = 0.8,#文本标签的大小
         na.label = "",#NA值如何表示
         na.label.col = "#bcb9b5",#表示NA值符号的颜色
)
dev.off()




ggplot(UnivariateOScoxSummary, aes(x =OR, y = Variables)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI,color=Change), height = 0.2,size=1) +
  geom_point(alpha=0.3,aes(size=P_Value)) +
  xlim(min(or_ci), max(or_ci)) +
  labs(x = "OR (95%CI)", y = "") +
  geom_vline(xintercept = 1.0,color="grey",linetype=2,linewidth=1)+
  scale_color_manual(values = c("#377eb8","#4daf4a","#e41a1c"))+
  theme_bw()+
  theme(axis.text = element_text(size = 12,color = "black"),
        panel.border = element_rect(linewidth = 1))







