remove_small_samples <- function(seurat_obj, meta_key='orig.ident', min_cells = 200) {
  seurat_obj$orig.ident=as.factor(seurat_obj$orig.ident)
  sample_counts <- table(seurat_obj@meta.data[[meta_key]])
  
  # 3. 筛选出细胞数大于等于 min_cells 的样本
  valid_samples <- names(sample_counts[sample_counts >= min_cells])
  
  # 4. 从 Seurat 对象中过滤出符合条件的细胞
  seurat_filtered <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[meta_key]] %in% valid_samples])
  # 5. 返回过滤后的 Seurat 对象
  return(seurat_filtered)
}
scRNA_harmony=remove_small_samples(OS)
table(scRNA_harmony$orig.ident)
