library(Seurat)
library(Azimuth)
library(SeuratData)

query <- RunAzimuth(
  query = "你的10x_filtered_feature_bc_matrix.h5",
  reference = "pbmcref"
)

anno <- data.frame(
  cell_id = colnames(query),
  predicted.celltype.l1 = query$predicted.celltype.l1,
  predicted.celltype.l2 = query$predicted.celltype.l2,
  predicted.celltype.l2.score = query$predicted.celltype.l2.score
)

write.csv(anno, "azimuth_pbmc_annotations.csv", row.names = FALSE)