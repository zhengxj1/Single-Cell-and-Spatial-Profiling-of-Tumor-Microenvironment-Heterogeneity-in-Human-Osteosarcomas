####dotplot featureplot vlnplot####

library(tidyverse)
library(MySeuratWrappers)
library(scCustomize)
markers=c("Krt14","Krt13","Krt5","Krt4","Col17a1","Dst","Top2a","Mki67","Stra6","Krt15","Dmkn","Sult2b1","Sprr1a","S100a9",
          "Sox2","Tp63")
markers=osteo
scdata=qread('3 sndata_all_merge.qs')
Cat_Col=color_qualitative3
Seurat_plot=function(scdata,group_by,markers){
  dir.create("./Seurat_plot")
  dir.create("./Seurat_plot/featureplot")
  library(tidyverse)
  library(scCustomize)
  library(scop)
  Idents(scdata)=scdata$orig.ident
  MySeuratWrappers::VlnPlot(scdata, features = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt","pct_counts_rb"), 
                            pt.size=0,ncol = 2,cols = color_patient,x.lab = "")
  ggsave("1. QCvln.ident.pdf",width =12, path = "./Seurat_plot")
  Idents(scdata)=scdata$celltypes
  MySeuratWrappers::VlnPlot(scdata, features =c("nCount_RNA", "nFeature_RNA", "pct_counts_mt","pct_counts_rb"), pt.size=0,ncol = 2,cols = Cat_Col,x.lab = "")
  ggsave("1. QCvln.pdf",path = "./Seurat_plot",height =7,width = 5)
  CellDimPlot(
    srt = scdata,
    group.by = c("orig.ident"),
    reduction = "UMAP",palcolor =color_patient,
    theme_use = "theme_blank"
  )
  ggsave("2. DimPlot_by_ident.pdf",width = 6,height = 5,path = "./Seurat_plot",dpi = 1000) 
  CellDimPlot(
    srt = scdata,
    group.by = c("celltypes"),
    reduction = "UMAP",palcolor =Cat_Col,
    theme_use = "theme_blank"
  )
  ggsave("2.DimPlot_by_celltype.pdf",width = 6,height = 5,path = "./Seurat_plot") 
 
  CellDimPlot(
    srt = scdata,
    group.by = c("celltypes"),label = T,
    reduction = "UMAP",palcolor =Cat_Col,
    theme_use = "theme_blank"
  )
  ggsave("3. DimPlot.pdf",width = 7,height = 5,path = "./Seurat_plot") 
  CellDimPlot(
    srt = scdata,
    group.by = c("celltypes"),split.by = 'group',
    reduction = "UMAP",palcolor =Cat_Col,
    theme_use = "theme_blank"
  )
  ggsave("3.DimPlot_split_Group.pdf",width =10,height = 5,path = "./Seurat_plot") 
  markers=unique(markers)
  FeatureDimPlot(
    srt = scdata,
    features = markers,
    reduction = "UMAP",
    theme_use = "theme_blank"
  )
  ggsave("4.FeaturePlot.pdf",path = "./Seurat_plot")
  MySeuratWrappers::VlnPlot(scdata, features = unique(markers),
                            pt.size=0,stacked=T,
                            cols = Cat_Col,#颜色  
                            direction = "horizontal", #水平作图  
                            x.lab = '', y.lab = '')+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave("4.VlnPlot.pdf",path = "./Seurat_plot",width = length(markers)/4,height = 5) 
  DotPlot(scdata, features = markers,group.by = "celltypes")+coord_flip()+theme_CNS_blank()  +theme_angle()+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
  ggsave("5.1 DotPlot.pdf",path = "./Seurat_plot",height = length(markers)/4) 
  DotPlot(scdata, features = markers)+theme_CNS_blank()  +theme_angle()+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
  ggsave("5.2 DotPlot.pdf",path = "./Seurat_plot",height = length(markers)/4) 
  plots=list()
  all_cells <- colnames(scdata)
  for (i in 1:length(markers)) {
    tryCatch({
      plot <-SCP::FeatureDimPlot(
        srt = scdata,pt.size = 0.001,
        features = markers[i],sample(all_cells, size = round(0.1 * length(all_cells))),
        reduction = "UMAP",
        theme_use = "theme_blank"
      )
      ggsave(filename = paste0("FeaturePlot_", paste0(markers[i], ".pdf")), width = 2.5, height = 2,limitsize = FALSE, path = "./Seurat_plot/featureplot")
      plots[[i]] <- plot
    }, error = function(e) {
      cat("Error encountered for marker:", markers[i], "\n")
    })
  }
  p<-patchwork::wrap_plots(compact(plots), ncol = 3)
  ggsave(p,file="Featureplot_All.pdf",width = 10,height = length(markers)*1,dpi =300,limitsize = FALSE,path = "./Seurat_plot/featureplot" )
}
group_by="group"
Seurat_plot(scdata,"ANP",markers)
Seurat_plot(scdata,"maintype",markers)