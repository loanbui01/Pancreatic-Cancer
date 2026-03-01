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
library(presto)        # Faster Wilconxon DEG
library(pheatmap)      # Visualization
library(stringr)
library(scales)
library(lme4)
library(broom.mixed)
library(forcats)

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

# QC visualization
VlnPlot(sr_2, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3, pt.size = 0)

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

# QC visualization
VlnPlot(sr_3, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3, pt.size = 0)

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

# QC visualization
VlnPlot(sr_4, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3, pt.size = 0)

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

# I would recommend making the UMAP based on sample type after you have saved this object: 
# Inspect unique names first (optional)
unique(combined$orig.ident)

# Create a new merged metadata column
combined$Sample_Group <- str_remove(combined$orig.ident, "-.*")

# Plot UMAP with merged grouping
DimPlot(combined,
        reduction = "umap",
        group.by = "Sample_Group")

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

table(Idents(combined), combined$SingleR_label)       #Check agreement between clusters and singleR
tab <- table(Idents(combined), combined$SingleR_label)
tab_df <- as.data.frame.matrix(tab)
write.csv(tab_df, "cluster_vs_SingleR.csv")

#____________________________________________________________#
# MANUAL ANNOTATION SUPPORT
#____________________________________________________________#

canonical_markers <- c(
  "CD79A", "CD3D", "NKG7", "COL1A1",
  "CD68", "CD14", "CLDN5", "KRT18",
  "TPSAB1", "SOX2", "PPBP", "CLEC4C"
)

# Checking 12 canonical markers

plots <- FeaturePlot(
  combined,
  features = canonical_markers,
  ncol = 3,
  combine = FALSE   # IMPORTANT: returns list of plots
)

plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    ggtitle(canonical_markers[i]) +
    theme(
      plot.title = element_text(
        size = 18,        # Larger title
        face = "bold",
        hjust = 0.5
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    )
})

wrap_plots(plots, ncol = 3)

# You would then inspect the canonical markers before doing manual annotation: 

# Make sure identities are clusters
Idents(combined) <- "seurat_clusters"

# Visualize marker distribution across clusters
DotPlot(combined, features = canonical_markers) +
  RotatedAxis() +
  theme_classic()

# Named vector mapping cluster number -> cell type
annotation_map <- c(
  "0" = "T cells", "1" = "T cells", "2" = "B cells", "3" = "NK cells",
  "4" = "Epithelial cells", "5" = "T cells", "6" = "Fibroblasts",
  "7" = "Macrophages/Monocytes", "8" = "T cells",
  "9" = "Macrophages/Monocytes", "10" = "T cells",
  "11" = "Endothelial cells", "12" = "Fibroblasts",
  "13" = "B cells", "14" = "T cells",
  "15" = "Epithelial cells", "16" = "B cells",
  "17" = "Epithelial cells", "18" = "Mast cells",
  "19" = "T cells", "20" = "Platelets",
  "21" = "T cells", "22" = "T cells",
  "23" = "T cells", "24" = "T cells",
  "25" = "Epithelial cells", "26" = "T cells",
  "27" = "Neuronal-like cells", "28" = "Specialized immune cells",
  "29" = "T cells", "30" = "Epithelial cells",
  "31" = "Epithelial cells", "32" = "Epithelial cells",
  "33" = "B cells", "34" = "B cells",
  "35" = "Epithelial cells", "36" = "B cells"
)

# Convert clusters to character
clusters_char <- as.character(combined$seurat_clusters)

# Map clusters to labels
manual_labels <- annotation_map[clusters_char]

# IMPORTANT: attach cell names
names(manual_labels) <- colnames(combined)

# Add metadata
combined <- AddMetaData(
  combined,
  metadata = manual_labels,
  col.name = "ManualAnnotation"
)

# Set identities
Idents(combined) <- "ManualAnnotation"

# Plot
DimPlot(combined, label = TRUE, repel = TRUE)

#____________________________________________________________#
# FIGURE 7 - CELL TYPE COMPOSITION
#____________________________________________________________#

meta <- as.data.frame(combined@meta.data)

meta <- meta %>%
  filter(!is.na(ManualAnnotation),
         tissue %in% c("Adjacent_normal", "Tumor", "PBMC"))

cell_order <- meta$ManualAnnotation %>%
  table() %>%
  sort(increasing = TRUE) %>%
  names()

colnames(combined@meta.data)

meta$ManualAnnotation <- factor(meta$ManualAnnotation, levels = cell_order)
combined$ManualAnnotation <- meta$ManualAnnotation

# % Cell type by Condition
df_condition <- meta %>%
  filter(tissue %in% c("Adjacent_normal", "Tumor")) %>%
  group_by(ManualAnnotation, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ManualAnnotation) %>%
  mutate(percent = n / sum(n) * 100)

p1 <- ggplot(df_condition,
             aes(x = ManualAnnotation, y = percent, fill = tissue)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100)) +
  scale_fill_manual(
    values = c("Adjacent_normal" = "blue",
               "Tumor" = "red"),
    name = "Tissue Type",
    labels = c("Normal", "Tumor")
  ) +
  labs(title = "% Cell Type by Condition",
       x = "Cell Type",
       y = "Percentage (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

p1

# % Cell type by Tissue
df_tissue <- meta %>%
  group_by(ManualAnnotation, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ManualAnnotation) %>%
  mutate(percent = n / sum(n) * 100)

p2 <- ggplot(df_tissue,
             aes(x = ManualAnnotation, y = percent, fill = tissue)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100)) +
  scale_fill_manual(
    values = c("Adjacent_normal" = "blue",
               "Tumor" = "red",
               "PBMC" = "yellow"),
    name = "Tissue Type",
    labels = c("Normal", "PBMC", "Tumor")
  ) +
  labs(title = "% Cell Type by Tumor Locations",
       x = "Cell Type",
       y = "Percentage (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

p2

# Total Cell Numbers
meta$ManualAnnotation <- as.character(unlist(meta$ManualAnnotation))
str(meta$ManualAnnotation)

df_total <- as.data.frame(table(meta$ManualAnnotation))
colnames(df_total) <- c("ManualAnnotation", "total_cells")

p3 <- ggplot(df_total,
             aes(x = total_cells, y = ManualAnnotation)) +
  labs(title = "Total Number of Cell Types",
       x = "Cell Counts",
       y = "Cell Type") + 
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic()

p3

# Combine
p1 | p2 | p3

#____________________________________________________________#
# FIGURE 9 - GLMM 
#____________________________________________________________#

