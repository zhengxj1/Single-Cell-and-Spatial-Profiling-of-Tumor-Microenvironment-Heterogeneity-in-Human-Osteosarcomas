#Cellchat
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)

#Step1. 构建cellchat对象
cellchat_obj=scdata
cellchat_obj=createCellChat(object =cellchat_obj, 
                            meta = cellchat_obj@meta.data,
                            group.by = "celltypes")

#Step2. 加载CellChatDB数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat_obj@DB <- CellChatDB.use

##Step3. 对表达数据进行预处理
### 1.5 对表达数据进行预处理，用于细胞间通讯分析
# subset the expression data of signaling genes for saving computation cost
cellchat_obj <- subsetData(cellchat_obj) # This step is necessary even if using the whole database

cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

#Step4. 计算通讯概率，推断细胞通讯网络
cellchat_obj  <- computeCommunProb(cellchat_obj ,population.size = F)
cellchat_obj <- filterCommunication(cellchat_obj , min.cells = 10)
saveRDS(cellchat_obj, file = 'cellchat_res.rds')
cellchat_obj = readRDS(file = 'cellchat_res.rds')
#Step5. 提取预测的细胞通讯网络为data frame
### cellchat取子集
#barcode.use = sample(row.names(cellchat@meta),100)
#cellchat.subset = subsetCellChat(cellchat,cells.use = barcode.use)
#作者提供了函数subsetCommunication去提取用户感兴趣的配受体对强度为data.frame：
#获取所有的配受体对以及其通讯概率
df.net <- subsetCommunication(cellchat_obj)
dir.create("./cellchat/")
write.csv(df.net, "./cellchat/cellchat_cell-cell_communications.all.csv")
#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat_obj,slot.name = "netP")

#Step6. 在信号通路水平推断细胞通讯 
###推断出的每个配体-受体对的细胞间通讯网络和每个信号通路分别存储在' net '和' netP '插槽中。
cellchat_obj <- computeCommunProbPathway(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)
head(cellchat_obj@net)
head(cellchat_obj@netP)

#Step7. 计算加和的cell-cell通讯网络 用户还可以通过设置sources.use 和targets.use来计算细胞亚群之间的加和网络
netVisual_bubble(cellchat_obj, sources.use = c(3,4), 
                 targets.use = c(1,2,3,4), remove.isolate = FALSE)

groupSize <- as.numeric(table(cellchat_obj@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize,
                 weight.scale =TRUE, label.edge= F,
                 title.name = "Interaction weights/strength")

netAnalysis_contribution(cellchat, signaling = 'EGF')

pathways.show = 'EGF'

#由于细胞间通讯网络的复杂性，我们可以对每个细胞亚群发出的信号进行检测。
#这里我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边权值：
mat <- cellchat_obj@net$weight
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}


#Step8. 使用层次图（Hierarchical plot），圆圈图（Circle plot）或和弦图（Chord diagram）可视化每个信号通路
pathways.show <- c("MIF") 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",label.edge= T)

# Circle plot show L-R pairs
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
LR.show
#"MIF_CD74_CXCR4"
# Vis
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = "MIF_CD74_CXCR4", layout = "circle")









