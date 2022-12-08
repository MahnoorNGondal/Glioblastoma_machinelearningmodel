setwd("C:/Users/chris/Desktop/final_590")

library(Seurat)
library(dplyr)
library(patchwork)

#loading data
count_matrix <- Read10X_h5("Glioma_GSE131928_10X_expression.h5", use.names = TRUE, unique.features = TRUE)
glioma_metadata <- read.table("Glioma_GSE131928_10X_CellMetainfo_table.tsv", header = TRUE, sep = "\t")

row.names(glioma_metadata) <- glioma_metadata$Cell
#creating seurat object
count_seurat <- CreateSeuratObject(counts = count_matrix, meta.data = glioma_metadata)

#from seurat tutorial
#count_seurat[["percent.mt"]] <- PercentageFeatureSet(count_seurat, pattern = "^MT-")
#count_seurat <- NormalizeData(count_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
count_seurat <- FindVariableFeatures(count_seurat, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(count_seurat), 10)
all.genes <- rownames(count_seurat)
count_seurat <- ScaleData(count_seurat, features = all.genes)
count_seurat <- RunPCA(count_seurat, features = VariableFeatures(object = count_seurat))
count_seurat <- FindNeighbors(count_seurat, dims = 1:10)
count_seurat <- FindClusters(count_seurat, resolution = 0.5)
#creating UMAP
count_seurat <- RunUMAP(count_seurat, dims = 1:10)
DimPlot(count_seurat, reduction = "umap")
DimPlot(count_seurat, reduction = "umap", group.by = "Celltype..malignancy.")
DimPlot(count_seurat, reduction = "umap", group.by = "Cluster")
DimPlot(count_seurat, reduction = "umap", group.by = "Celltype..major.lineage.")
DimPlot(count_seurat, reduction = "umap", group.by = "Celltype..minor.lineage.")

#using FindAllMarkers
markers_myclusters <- FindAllMarkers(count_seurat) 
Idents(count_seurat = count_seurat) <- count_seurat@meta.data$'Cluster'
markers_ogcluster <- FindAllMarkers(count_seurat)
Idents(count_seurat = count_seurat) <- count_seurat@meta.data$'Celltype..malignancy.'
markers_malignancy2 <- FindAllMarkers(count_seurat)
Idents(count_seurat = count_seurat) <- count_seurat@meta.data$'Celltype..major.lineage.'
markers_majlineage <- FindAllMarkers(count_seurat)
Idents(count_seurat = count_seurat) <- count_seurat@meta.data$'Celltype..minor.lineage.'
markers_minlineage <- FindAllMarkers(count_seurat)

write.table(markers_malignancy, "markers_myclusters.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(markers_malignancy2, "markers_malignancy.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(markers_ogcluster, "markers_originalclusters.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(markers_majlineage, "markers_majlineage.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(markers_minlineage, "markers_minlineage.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#top/bottom 50 malignant/min lineage genes
markers_malignancy2$avg_log2FC_abs <- abs(markers_malignancy2$avg_log2FC)
mal_markers_ordered <- markers_malignancy2[order(markers_malignancy2$avg_log2FC_abs),]
mal_markers_top50 <- mal_markers_ordered[1:50,]
mal_markers_bottom50 <- mal_markers_ordered[5190:5240,]

markers_minlineage$avg_log2FC_abs <- abs(markers_minlineage$avg_log2FC)
min_markers_ordered <- markers_minlineage[order(markers_minlineage$avg_log2FC_abs),]
min_markers_top50 <- min_markers_ordered[1:50,]
min_markers_bottom50 <- min_markers_ordered[5190:5240,]

write.table(mal_markers_top50, "malignancy_top50.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(mal_markers_bottom50, "malignancy_bottom50.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(min_markers_top50, "min_lineage_top50.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(min_markers_bottom50, "min_lineage_bottom50.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#quick stats for report
table(glioma_metadata$Celltype..malignancy.)
#Immune cells Malignant cells          Others 
#4166            8998             394 
#0.30727          0.66366         0.02906
table(glioma_metadata$Celltype..major.lineage.)/13558
# AC-like Malignant             CD8Tex          Malignant MES-like Malignant 
#262                197                367               5203 
#Mono/Macro NPC-like Malignant    Oligodendrocyte OPC-like Malignant 
#3969               1368                394               1798 
# AC-like Malignant             CD8Tex          Malignant MES-like Malignant 
#0.01932438         0.01453017         0.02706889         0.38375867 
#Mono/Macro NPC-like Malignant    Oligodendrocyte OPC-like Malignant 
#0.29274229         0.10089984         0.02906033         0.13261543 
table(glioma_metadata$Celltype..minor.lineage.)/13558
# AC-like Malignant             CD8Tex                 M1          Malignant 
#262                197               3679                367 
#MES-like Malignant           Monocyte NPC-like Malignant    Oligodendrocyte 
#5203                290               1368                394 
#OPC-like Malignant 
#1798
# AC-like Malignant             CD8Tex                 M1          Malignant 
#0.01932438         0.01453017         0.27135271         0.02706889 
#MES-like Malignant           Monocyte NPC-like Malignant    Oligodendrocyte 
#0.38375867         0.02138959         0.10089984         0.02906033 
#OPC-like Malignant 
#0.13261543 

#generating bar chart for report
cell_types <- read.table("cell-types.txt", header = TRUE, sep = "\t")
library(ggplot2)
cell_types$Cell_type <- as.character(cell_types$Cell_type)
cell_types$Cell_type <- factor(cell_types$Cell_type, levels=c("M1", "MES-like Malignant", "Malignant", "AC-like Malignant", "NPC-like Malignant", "Monocyte", "CD8Tex", "OPC-like Malignant", "Oligodendrocyte"))
ggplot(cell_types, aes(Number_of_cells, Cell_type)) +
  geom_bar(aes(fill = Cell_type), stat = "identity", position = "dodge") + 
  ggtitle("Cell types within single cell data") + labs(x="Number of Cells", y="Cell Type") 

#gene expression values
celltype_data <- AverageExpression(count_seurat, group.by = 'Celltype..minor.lineage.', slot = 'data')
celltype_counts <- AverageExpression(count_seurat, group.by = 'Celltype..minor.lineage.', slot = 'counts')
malignancy_data <- AverageExpression(count_seurat, group.by = 'Celltype..malignancy.', slot = 'data')
malignancy_counts <- AverageExpression(count_seurat, group.by = 'Celltype..malignancy.', slot = 'counts')
celltype_data <- data.frame(celltype_data)
malignancy_data <- data.frame(malignancy_data)
celltype_data$gene <- rownames(celltype_data)
malignancy_data$gene <- rownames(malignancy_data)
#(counts and data ended up being the same)
write.table(celltype_data, "celltype_data.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(celltype_counts, "celltype_counts.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(malignancy_data, "malignancy_data.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(malignancy_counts, "malignancy_counts.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#getting top/bottom 50 from average expression, so generating those datasets above first
mal_markers_top50 <- merge(mal_markers_top50,malignancy_data,by="gene")
mal_markers_bottom50 <- merge(mal_markers_bottom50,malignancy_data,by="gene")
min_markers_top50 <- merge(min_markers_top50,celltype_data,by="gene")
min_markers_bottom50 <- merge(min_markers_bottom50,celltype_data,by="gene")

saveRDS(count_seurat, "count_seurat.rds")
