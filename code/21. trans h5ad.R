library(SeuratDisk)
library(Seurat)
library(dior)
####R export####
SaveH5Seurat(scRNA_harmony, filename = "scRNA_harmony.h5Seurat")
Convert("scRNA_harmony.h5Seurat", dest = "h5ad")



##读入
Convert("oc_sc.h5ad", dest = "h5seurat", overwrite = TRUE)
scRNA_harmony <- LoadH5Seurat("oc_sc.h5Seurat")


####diopy R Export, Python inport####
library(dior)
library(SingleCellExperiment)
original_counts <- GetAssayData(scRNA_harmony, assay = "RNA", slot = "counts")
original_metadata <- scRNA_harmony@meta.data
new_seurat <- CreateSeuratObject(counts = original_counts, meta.data = original_metadata)
scdata=JoinLayers(scdata)
scdata[["RNA"]] <-as(object=scdata[["RNA"]], Class ="Assay")
dior::write_h5(scdata,file="scRNA_harmony.h5",object.type = "seurat")

adata = diopy.input.read_h5(file = 'scRNA_harmony.h5')

#### diopy Python Export, R inport####
diopy.output.write_h5(adata, 
                      file = 'data_py.h5')
scRNA_harmony <- read_h5(file = "data_py.h5",
                assay.name = 'RNA', 
                target.object = 'seurat')
setlevel <- function(adata, batch_key, min_cells = 10) {
  celltype_counts <- table(adata@meta.data[[batch_key]])
  ordered_celltypes <- names(sort(celltype_counts, decreasing = TRUE))
  adata@meta.data[[batch_key]] <- factor(adata@meta.data[[batch_key]], levels = ordered_celltypes, ordered = TRUE)
  return(adata)
}

scRNA_harmony=setlevel(scRNA_harmony,'celltypes')

scRNA_harmony$percent_mito=scRNA_harmony$pct_counts_mt
scRNA_harmony$batch=scRNA_harmony$orig.ident
scRNA_harmony$nFeature_RNA=scRNA_harmony$n_genes_by_counts
scRNA_harmony$
table(scRNA_harmony$celltypes)
table(scRNA_harmony$maintype)
