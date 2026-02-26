############################################################
# PDAC MULTI-DATASET scRNA-seq INTEGRATION PIPELINE
# Datasets:
# - GSE278688
# - GSE212966
# - GSE217845
# - GSE300304
############################################################

############################
# Load required libraries
############################
library(Seurat)        # Core single-cell framework
library(dplyr)         # Data manipulation
library(ggplot2)       # Visualization
library(patchwork)     # Combine plots
library(org.Hs.eg.db)  # Human gene annotation
library(GEOquery)      # GEO access

#____________________________________________________________#
# DATASET 1 — GSE278688
#____________________________________________________________#

setwd("DIRECTION/TO/FILE/GSE278688")     
list.files()  # Check files exist

# Read sparse matrix (genes x cells)
data1 <- ReadMtx(
  mtx = "GSE278688_sc_matrix.mtx.gz",
  features = "GSE278688_sc_features.tsv.gz",
  cells = "GSE278688_sc_barcodes.tsv.gz"
)

# Create Seurat object
sr_1 <- CreateSeuratObject(counts = data1)

# Add metadata
meta <- read.csv("GSE278688_sc_metadata.csv.gz", row.names = 1)
sr_1 <- AddMetaData(sr_1, metadata = meta)

# Calculate mitochondrial percentage (QC metric)
sr_1[["percent.mt"]] <- PercentageFeatureSet(sr_1, pattern = "^MT-")

# Standardize sample naming: PDAC05-Tumor-GEX → Tumor-05
sr_1$orig.ident <- sub(
  "PDAC(\\d+)-(Normal|Tumor|PBMC)-GEX",
  "\\2-\\1",
  sr_1$orig.ident
)

# Fix specific formatting issue
sr_1$orig.ident[sr_1$orig.ident == "PDAC05"] <- "Normal-05"

# Set identities to sample names
Idents(sr_1) <- "orig.ident"

# Order samples biologically: Normal → Tumor → PBMC
sample_type <- sub("-.*", "", sr_1$orig.ident)
type_order <- c("Normal", "Tumor", "PBMC")

ordered_levels <- sr_1$orig.ident[
  order(match(sample_type, type_order),
        as.numeric(sub(".*-", "", sr_1$orig.ident)))
]
ordered_levels <- unique(ordered_levels)

sr_1$orig.ident <- factor(sr_1$orig.ident, levels = ordered_levels)
Idents(sr_1) <- "orig.ident"

# QC visualization
VlnPlot(sr_1, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(sr_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(sr_1, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Filter low-quality cells
sr_1 <- subset(
  sr_1,
  subset = nCount_RNA <= 20000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 10
)

# Remove genes expressed in fewer than 3 cells
counts <- GetAssayData(sr_1, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= 3
sr_1 <- sr_1[keep_genes, ]

#____________________________________________________________#
# DATASET 2 — GSE212966
#____________________________________________________________#

setwd("DIRECTION/TO/FILE/GSE212966/GSE212966_RAW")

files <- list.files(pattern = "matrix.mtx.gz")
samples <- gsub("_matrix.mtx.gz", "", files)

# Load each sample separately
sample_list <- lapply(samples, function(s) {
  mat <- ReadMtx(
    mtx = paste0(s, "_matrix.mtx.gz"),
    features = paste0(s, "_genes.tsv.gz"),
    cells = paste0(s, "_barcodes.tsv.gz")
  )
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3)
  obj$sample <- s
  return(obj)
})

names(sample_list) <- samples

# Merge samples
sr_2 <- merge(sample_list[[1]], y = sample_list[-1], add.cell.ids = samples)

# QC metric
sr_2[["percent.mt"]] <- PercentageFeatureSet(sr_2, pattern = "^MT-")

# Clean sample names → Tumor-XX / Normal-XX
sr_2$orig.ident <- sub(".*_", "", sr_2$orig.ident)
num <- as.numeric(sub(".*?(\\d+)", "\\1", sr_2$orig.ident))
type <- ifelse(grepl("^PDAC", sr_2$orig.ident), "Tumor", "Normal")
sr_2$orig.ident <- paste0(type, "-", sprintf("%02d", num))

Idents(sr_2) <- "orig.ident"

# QC filtering
sr_2 <- subset(
  sr_2,
  subset = nCount_RNA <= 60000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 8000 &
    percent.mt <= 20
)

#____________________________________________________________#
# DATASET 3 — GSE217845
#____________________________________________________________#

setwd("DIRECTION/TO/FILE/GSE217845/GSE217845_RAW")

tumor_samples <- c(
  "GSM6727546_PDAC_47_Tumor",
  "GSM6727547_PDAC_48_Tumor",
  "GSM6727548_PDAC_50_Tumor",
  "GSM6727549_PDAC_51_Tumor",
  "GSM6727550_PDAC_55_Tumor",
  "GSM6727551_PDAC_60_Tumor",
  "GSM6727552_PDAC_47_PB",
  "GSM6727553_PDAC_48_PB",
  "GSM6727554_PDAC_50_PB",
  "GSM6727555_PDAC_51_PB",
  "GSM6727556_PDAC_55_PB",
  "GSM6727557_PDAC_60_PB"
)

sample_list_3 <- lapply(tumor_samples, function(s) {
  mat <- ReadMtx(
    mtx = paste0(s, "_matrix.mtx.gz"),
    features = paste0(s, "_features.tsv.gz"),
    cells = paste0(s, "_barcodes.tsv.gz")
  )
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3)
  obj$sample <- s
  return(obj)
})

names(sample_list_3) <- tumor_samples

sr_3 <- merge(sample_list_3[[1]], y = sample_list_3[-1], add.cell.ids = tumor_samples)

sr_3[["percent.mt"]] <- PercentageFeatureSet(sr_3, pattern = "^MT-")

# Rename → Tumor-XX / PBMC-XX
num <- sub(".*PDAC_(\\d+)_.*", "\\1", sr_3$orig.ident)
type <- ifelse(grepl("_PB$", sr_3$orig.ident), "PBMC", "Tumor")
sr_3$orig.ident <- paste0(type, "-", num)

Idents(sr_3) <- "orig.ident"

# QC filtering
sr_3 <- subset(
  sr_3,
  subset = nCount_RNA <= 30000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 25
)

#____________________________________________________________#
# DATASET 4 — GSE300304
#____________________________________________________________#

setwd("DIRECTION/TO/FILE/GSE300304/GSE300304_RAW")

files <- c(
  "GSM9058075_pt27_Pancreas_CD45Pos_live-GEX_S6.h5",
  "GSM9058085_037_Panc_CD45_filtered_feature_bc_matrix.h5",
  "GSM9058095_036_panc_CD45_filtered_feature_bc_matrix.h5",
  "GSM9058077_pt27_PBMC-GEX_S4.h5",
  "GSM9058087_037_PBMC_filtered_feature_bc_matrix.h5",
  "GSM9058097_036_PBMC_filtered_feature_bc_matrix.h5"
)

# Rename samples based on filename
rename_sample <- function(filename) {
  num <- sub(".*?(pt|_)(\\d+).*", "\\2", filename)
  num <- as.numeric(num)
  type <- ifelse(grepl("PBMC", filename, ignore.case = TRUE), "PBMC", "Tumor")
  paste0(type, "-", num)
}

obj_list <- lapply(files, function(f) {
  counts <- Read10X_h5(f)
  obj <- CreateSeuratObject(counts = counts)
  obj$orig.ident <- rename_sample(f)
  return(obj)
})

sr_4 <- merge(obj_list[[1]], y = obj_list[-1])
Idents(sr_4) <- "orig.ident"

# QC + gene filtering
sr_4[["percent.mt"]] <- PercentageFeatureSet(sr_4, pattern = "^MT-")
sr_4 <- JoinLayers(sr_4)

counts <- GetAssayData(sr_4, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= 3
sr_4 <- subset(sr_4, features = rownames(counts)[keep_genes])

sr_4 <- subset(
  sr_4,
  subset = nCount_RNA <= 25000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 20
)

#____________________________________________________________#
# MERGE ALL DATASETS
#____________________________________________________________#

sr_2 <- JoinLayers(sr_2)
sr_3 <- JoinLayers(sr_3)

combined <- merge(
  sr_1,
  y = list(sr_2, sr_3, sr_4),
  add.cell.ids = c("Dataset1","Dataset2","Dataset3","Dataset4"),
  project = "PDAC_combined"
)

combined <- JoinLayers(combined)

saveRDS(combined, file = "combined_PDAC_seurat.rds")

#____________________________________________________________#
# NORMALIZATION + DIMENSION REDUCTION
#____________________________________________________________#

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, features = VariableFeatures(combined))
combined <- RunPCA(combined, features = VariableFeatures(combined))

ElbowPlot(combined, ndims = 50)

saveRDS(combined, "combined_PDAC_postPCA.rds")

#____________________________________________________________#
# HARMONY BATCH CORRECTION
#____________________________________________________________#

combined <- harmony::RunHarmony(
  object = combined,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:40
)

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:40)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:40)

DimPlot(combined, reduction = "umap", label = TRUE)

saveRDS(combined, "combined_PDAC_postHarmony_40PCs.rds")

#____________________________________________________________#
# CLUSTER MARKERS
#____________________________________________________________#

markers <- FindAllMarkers(
  combined,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.25
)

write.csv(markers, "cluster_markers_level1.csv")

#____________________________________________________________#
# AUTOMATED CELL TYPE ANNOTATION (SingleR)
#____________________________________________________________#

sce <- as.SingleCellExperiment(combined)
ref <- celldex::HumanPrimaryCellAtlasData()

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
combined$SingleR_label <- pred$labels

DimPlot(combined, group.by = "SingleR_label", label = TRUE)

#____________________________________________________________#
# MANUAL ANNOTATION SUPPORT
#____________________________________________________________#

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

nice_table <- top_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(
    genes = paste0(gene, " (", round(avg_log2FC,2), ")", collapse = ", ")
  )

write.csv(nice_table, "cluster_markers.csv")

DotPlot(combined, features = unique(top_markers$gene)) + RotatedAxis()