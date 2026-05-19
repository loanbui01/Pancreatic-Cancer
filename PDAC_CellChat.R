library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)
library(GEOquery)
library(presto)
library(pheatmap)
library(stringr)
library(scales)
library(lme4)
library(broom.mixed)
library(forcats)

#______________________________________________________________________________#
#GSE278688

setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC/GSE278688")     
list.files()

data1 <- ReadMtx(
  mtx = "GSE278688_sc_matrix.mtx.gz",
  features = "GSE278688_sc_features.tsv.gz",
  cells = "GSE278688_sc_barcodes.tsv.gz"
)

sr_1 <- CreateSeuratObject(counts = data1)

meta <- read.csv("GSE278688_sc_metadata.csv.gz", row.names = 1)
sr_1 <- AddMetaData(sr_1, metadata = meta)

sr_1
head(sr_1)

sr_1[["percent.mt"]] <- PercentageFeatureSet(
  sr_1,
  pattern = "^MT-"
)

sr_1$orig.ident
unique(sr_1$orig.ident)

#_____________Trial lb_______________________________________________#
# Convert factor to character
# sr_1$orig.ident <- as.character(sr_1$orig.ident)

# Replace the incorrect sample name
# sr_1$orig.ident[sr_1$orig.ident == "PDAC05"] <- "PDAC05-PBMC-GEX"

# Convert back to factor
# sr_1$orig.ident <- as.factor(sr_1$orig.ident)

# Check again
# unique(sr_1$orig.ident)

#____________________________________________________________#

sr_1$orig.ident <- sub(
  "PDAC(\\d+)-(Normal|Tumor|PBMC)-GEX",
  "\\2-\\1",
  sr_1$orig.ident
)

sr_1$orig.ident[sr_1$orig.ident == "PDAC05"] <- "PBMC-05"

Idents(sr_1) <- "orig.ident"

###Group sample types
# Extract sample type
sample_type <- sub("-.*", "", sr_1$orig.ident)

# Define desired order
type_order <- c("Normal", "Tumor", "PBMC")

# Order by type first, then number
ordered_levels <- sr_1$orig.ident[
  order(match(sample_type, type_order),
        as.numeric(sub(".*-", "", sr_1$orig.ident)))
]

# Make unique (important)
ordered_levels <- unique(ordered_levels)

# Set factor order
sr_1$orig.ident <- factor(sr_1$orig.ident, levels = ordered_levels)

# Update Seurat identities
Idents(sr_1) <- "orig.ident"

# Extract sample type from each level
sample_levels <- levels(sr_1$orig.ident)
sample_type <- sub("-.*", "", sample_levels)

# Assign colors based on group
sample_colors <- ifelse(sample_type == "Tumor", "red",
                        ifelse(sample_type == "Normal", "blue",
                               ifelse(sample_type == "PBMC", "yellow", "gray")))

# Name the vector so Seurat matches correctly
names(sample_colors) <- sample_levels


VlnPlot(
  sr_1,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "orig.ident",
  cols = sample_colors,
  ncol = 3,
  pt.size = 0
)

FeatureScatter(sr_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(sr_1, feature1 = "nCount_RNA", feature2 = "percent.mt")

sr_1 <- subset(
  sr_1,
  subset =
    nCount_RNA <= 20000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 10
)

counts <- GetAssayData(
  sr_1,
  assay = "RNA",
  layer = "counts"
)

keep_genes <- Matrix::rowSums(counts > 0) >= 3

sr_1 <- sr_1[keep_genes, ]

VlnPlot(
  sr_1,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

sr_1
#______________________________________________________________________________#
#GSE212966

setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC/GSE212966/GSE212966_RAW")     
list.files()

files <- list.files(pattern = "matrix.mtx.gz")
files

samples <- gsub("_matrix.mtx.gz", "", files)

sample_list <- lapply(samples, function(s) {
  
  mat <- ReadMtx(
    mtx = paste0(s, "_matrix.mtx.gz"),
    features = paste0(s, "_genes.tsv.gz"),
    cells = paste0(s, "_barcodes.tsv.gz")
  )
  
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3)   #Filter for 3 cells
  
  obj$sample <- s   # add sample label
  
  return(obj)
})

sample_list

names(sample_list) <- samples

sr_2 <- merge(
  sample_list[[1]],
  y = sample_list[-1],
  add.cell.ids = samples
)

sr_2

sr_2[["percent.mt"]] <- PercentageFeatureSet(sr_2, pattern = "^MT-")

sr_2$orig.ident
unique(sr_2$orig.ident)

sr_2$orig.ident <- sub(".*_", "", sr_2$orig.ident)

# Extract number
num <- as.numeric(sub(".*?(\\d+)", "\\1", sr_2$orig.ident))

# Convert type
type <- ifelse(grepl("^PDAC", sr_2$orig.ident), "Tumor", "Normal")

# Rebuild name
sr_2$orig.ident <- paste0(type, "-", sprintf("%02d", num))

Idents(sr_2) <- "orig.ident"

sr_2 <- JoinLayers(sr_2)

sr_2$orig.ident <- factor(sr_2$orig.ident)

# Extract sample type from each level
sample_levels <- levels(sr_2$orig.ident)
sample_type <- sub("-.*", "", sample_levels)

# Assign colors based on group
sample_colors <- ifelse(sample_type == "Tumor", "red",
                        ifelse(sample_type == "Normal", "blue",
                               "gray"))

# Name the vector so Seurat matches correctly
names(sample_colors) <- sample_levels

VlnPlot(
  sr_2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "orig.ident",
  cols = sample_colors,
  ncol = 3,
  pt.size = 0
)

sr_2 <- subset(
  sr_2,
  subset =
    nCount_RNA <= 60000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 8000 &
    percent.mt <= 20
)

sr_2
#______________________________________________________________________________#
#GSE217845

setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC/GSE217845/GSE217845_RAW")     
list.files()

tumor_samples <- c(
  "GSM6727545_LPDAC_30_Tumor",
  "GSM6727548_PDAC_50_Tumor",
  "GSM6727550_PDAC_55_Tumor",
  "GSM6727551_PDAC_60_Tumor",
  "GSM6727554_PDAC_50_PB",
  "GSM6727556_PDAC_55_PB",
  "GSM6727557_PDAC_60_PB"
)

sample_list_3 <- lapply(tumor_samples, function(s) {
  
  mat <- ReadMtx(
    mtx = paste0(s, "_matrix.mtx.gz"),
    features = paste0(s, "_features.tsv.gz"),
    cells = paste0(s, "_barcodes.tsv.gz")
  )
  
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3)   #Filter for 3 cells
  
  obj$sample <- s   # add sample label
  
  return(obj)
})

sample_list_3

names(sample_list_3) <- tumor_samples

sr_3 <- merge(
  sample_list_3[[1]],
  y = sample_list_3[-1],
  add.cell.ids = tumor_samples
)

sr_3

sr_3[["percent.mt"]] <- PercentageFeatureSet(sr_3, pattern = "^MT-")

sr_3$orig.ident
unique(sr_3$orig.ident)

# Extract patient number
num <- sub(".*PDAC_(\\d+)_.*", "\\1", sr_3$orig.ident)

# Determine sample type
type <- ifelse(grepl("_PB$", sr_3$orig.ident),
               "PBMC",
               "Tumor")

# Rebuild names
sr_3$orig.ident <- paste0(type, "-", num)

Idents(sr_3) <- "orig.ident"

sr_3 <- JoinLayers(sr_3)

# Get unique sample names
sample_levels <- unique(sr_3$orig.ident)

# Extract type and number from levels
sample_type <- sub("-.*", "", sample_levels)
sample_number <- as.numeric(sub(".*-", "", sample_levels))

# Define desired biological order
type_order <- c("Tumor", "PBMC")

# Order levels
ordered_levels <- sample_levels[
  order(
    match(sample_type, type_order),
    sample_number
  )
]

# Apply ordered factor
sr_3$orig.ident <- factor(sr_3$orig.ident, levels = ordered_levels)
Idents(sr_3) <- "orig.ident"


# -----------------------------
# Create color vector
# -----------------------------

sample_colors <- ifelse(sample_type == "Tumor", "red",
                        ifelse(sample_type == "PBMC", "yellow",
                               "gray"))

names(sample_colors) <- sample_levels

# Reorder color vector to match new levels
sample_colors <- sample_colors[ordered_levels]


# -----------------------------
# Plot
# -----------------------------

VlnPlot(
  sr_3,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "orig.ident",
  cols = sample_colors,
  ncol = 3,
  pt.size = 0
)

# QC filtering
sr_3 <- subset(
  sr_3,
  subset = nCount_RNA <= 30000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 25
)

sr_3

################################################################################
setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC/GSE300304/GSE300304_RAW")     
list.files()

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

# Get unique sample names
sample_levels <- unique(sr_4$orig.ident)

# Extract type and number
sample_type   <- sub("-.*", "", sample_levels)
sample_number <- as.numeric(sub(".*-", "", sample_levels))

# Define desired order
type_order <- c("Tumor", "PBMC")

# Order levels biologically
ordered_levels <- sample_levels[
  order(
    match(sample_type, type_order),
    sample_number
  )
]

# Apply ordered factor
sr_4$orig.ident <- factor(sr_4$orig.ident, levels = ordered_levels)
Idents(sr_4) <- "orig.ident"


sample_colors <- ifelse(sample_type == "Tumor", "red",
                        ifelse(sample_type == "PBMC", "yellow", "gray"))

names(sample_colors) <- sample_levels

# Reorder colors to match factor levels
sample_colors <- sample_colors[ordered_levels]


# QC + gene filtering
sr_4[["percent.mt"]] <- PercentageFeatureSet(sr_4, pattern = "^MT-")
sr_4 <- JoinLayers(sr_4)

counts <- GetAssayData(sr_4, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= 3
sr_4 <- subset(sr_4, features = rownames(counts)[keep_genes])

VlnPlot(
  sr_4,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "orig.ident",
  cols = sample_colors,
  ncol = 3,
  pt.size = 0
)

sr_4 <- subset(
  sr_4,
  subset = nCount_RNA <= 25000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 20
)
#______________________________________________________________________________#

sr_1
sr_2
sr_3
sr_4

#You should join layers for everything before merging. 

sr_2 <- JoinLayers(sr_2)
sr_3 <- JoinLayers(sr_3)

combined <- merge(
  sr_1,
  y = list(sr_2, sr_3, sr_4),
  add.cell.ids = c("Dataset1","Dataset2","Dataset3","Dataset4"),
  project = "PDAC_combined"
)

combined <- JoinLayers(combined)

combined

saveRDS(combined, file = "C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC/Updated_combined_PDAC_seurat.rds")

#______________________________________________________________________________#
# setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC")
# combined <- readRDS("combined_PDAC_seurat.rds")

meta <- combined@meta.data
head(meta)

write.csv(meta, "updated_PDAC_cell_metadata.csv")

#______________________________________________________________________________#
library(SingleR)
library(celldex)
library(Seurat)
library(harmony)
library(SingleCellExperiment)

# setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC")
# combined <- readRDS("combined_PDAC_seurat.rds")
dim(combined)

#Normalize
combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


#ggplot(hvf_info, aes(x = mean, y = variance.standardized, color = selected)) +
#geom_point(alpha = 0.5) +
#theme_classic() +
#labs(
#title = "Highly Variable Gene Selection",
#x = "Mean Expression",
#y = "Standardized Variance"
#) +
#theme(
#plot.title = element_text(size = 20, face = "bold"),
#axis.title.x = element_text(size = 18),
#axis.title.y = element_text(size = 18),
#axis.text.x  = element_text(size = 16),
#axis.text.y  = element_text(size = 16)
#)

#Scale data
#combined <- ScaleData(
#combined,
#features = rownames(combined)
#)     #Error: cannot allocate vector of size 89.4 Gb
#So we should only scale the top 2000 HVGs

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined)
)

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined)
)

# Extract PCA standard deviations
pca_stdev <- combined[["pca"]]@stdev

# Calculate percent variance explained
percent_var <- (pca_stdev^2) / sum(pca_stdev^2) * 100

# Calculate cumulative variance
cumulative_var <- cumsum(percent_var)

plot(
  cumulative_var,
  type = "b",
  xlab = "Principal Components",
  ylab = "Cumulative Variance Explained (%)",
  main = "Cumulative Variance Explained by PCs"
)
abline(h = 95, col = "red", lty = 2)

VariableFeaturePlot(combined)

ElbowPlot(combined, ndims = 50)

#saveRDS(combined, file = "combined_PDAC_postPCA.rds")

#combined <- readRDS("combined_PDAC_postPCA.rds")

combined <- harmony::RunHarmony(
  object = combined,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:40
)

Reductions(combined)
DimPlot(combined, reduction = "harmony")


combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:40
)

combined <- FindClusters(
  combined,
  resolution = 0.5
)

length(unique(Idents(combined)))

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:40
)

DimPlot(combined, reduction = "umap", label = TRUE) #Plot cluster
DimPlot(combined, reduction = "umap", group.by = "orig.ident") #Plot batch mixing

# Inspect unique names first (optional)
unique(combined$orig.ident)

# Create a new merged metadata column
combined$Sample_Group <- str_remove(combined$orig.ident, "-.*")

# Plot UMAP with merged grouping
DimPlot(combined,
        reduction = "umap",
        group.by = "Sample_Group") +
  scale_color_manual(values = c(
    "Tumor" = "red",
    "Normal" = "blue",
    "PBMC" = "yellow"
  ))

saveRDS(combined, "updated_combined_PDAC_postHarmony_40PCs_res05.rds")

markers <- FindAllMarkers(
  combined,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.25
)

write.csv(markers, "updated_cluster_markers_corrected_level1_res05.csv")


sce <- as.SingleCellExperiment(combined)

ref <- celldex::HumanPrimaryCellAtlasData()

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

combined$SingleR_label <- pred$labels

DimPlot(combined, group.by = "SingleR_label", label = TRUE)

saveRDS(combined, file = "updated_combined_post_SingleR.rds")

# setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC")
# combined <- readRDS("combined_post_SingleR.rds")

table(Idents(combined), combined$SingleR_label)       #Check agreement between clusters and singleR
tab <- table(Idents(combined), combined$SingleR_label)
tab_df <- as.data.frame.matrix(tab)
write.csv(tab_df, "updated_cluster_vs_SingleR_res05.csv")

################################################################################
#manual annotation
library(dplyr)

head(markers)
dim(markers)

top_markers <- markers %>%                  #Extract top markers per cluster
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

nice_table <- top_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(
    genes = paste0(gene, " (", round(avg_log2FC,2), ")", collapse = ", ")
  )

View(nice_table)

write.csv(nice_table, "updated_cluster_markers_res05.csv")

DotPlot(combined, features = unique(top_markers$gene)) + RotatedAxis()

VlnPlot(combined, features = c("nFeature_RNA","percent.mt"), group.by = "seurat_clusters")
###############################################################################
setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC")
combined <- readRDS("combined_post_SingleR.rds")

canonical_markers <- c(
  "CD79A", "CD3D", "NKG7", "COL1A1",
  "CD68", "CD14", "CLDN5", "KRT18",
  "TPSAB1", "SOX2", "PPBP", "CLEC4C"
)

canon_name <- c(
  "B cells", "T cells", "NK cells", "Fibroblasts",
  "Macrophage", "Monocytes", "Endothelial cells", "Epithelial cells",
  "Mast cells", "Neuronal-like cells", "Platelets", "Specialized immune cells"
)

plots <- FeaturePlot(
  combined,
  features = canonical_markers,
  ncol = 3,
  combine = FALSE   # IMPORTANT: returns list of plots
)

plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    ggtitle(paste(canonical_markers[i], "-", canon_name[i])) +
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

#Make sure identities are clusters
Idents(combined) <- "seurat_clusters"

#Visualize marker distribution across clusters
# DotPlot(combined, features = canonical_markers) +
  # RotatedAxis() +
  # theme_classic()

# Named vector mapping cluster number -> cell type
annotation_map <- c(
  "0" = "T cells", "1" = "T cells", "2" = "B cells", "3" = "NK cells",
  "4" = "T cells", "5" = "Monocytes", "6" = "Epithelial cells",
  "7" = "T cells", "8" = "Fibroblasts",
  "9" = "T cells", "10" = "Endothelial cells",
  "11" = "B cells", "12" = "Macrophages",
  "13" = "Fibroblasts", "14" = "Epithelial cells",
  "15" = "B cells", "16" = "Mast cells",
  "17" = "Monocytes", "18" = "T cells",
  "19" = "T cells", "20" = "Platelets",
  "21" = "T cells", "22" = "T cells",
  "23" = "Neuronal-like cells", "24" = "T cells",
  "25" = "T cells", "26" = "B cells",
  "27" = "Neuronal-like cells", "28" = "Epithelial cells",
  "29" = "Epithelial cells", "30" = "B cells",
  "31" = "Epithelial cells"
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

saveRDS(combined, file = "combined_post_ManualAnnotation.rds")

################################################################################
combined = readRDS("combined_post_ManualAnnotation.rds")

meta <- combined@meta.data %>%
  filter(!is.na(ManualAnnotation),
         tissue %in% c("Tumor", "Adjacent_normal", "PBMC"))

meta$ManualAnnotation <- as.character(unlist(meta$ManualAnnotation))
meta$tissue <- as.character(unlist(meta$tissue))

cell_order <- meta %>%
  dplyr::count(ManualAnnotation) %>%
  dplyr::arrange(n) %>%
  dplyr::pull(ManualAnnotation)

meta$ManualAnnotation <- factor(meta$ManualAnnotation,
                                levels = cell_order)

df_total <- meta %>%
  dplyr::count(ManualAnnotation)

p1 <- ggplot(df_total,
             aes(x = n, y = ManualAnnotation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total Number of Cell Types",
       x = "Cell Counts",
       y = "Cell Type") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

p1

df_locations <- meta %>%
  group_by(ManualAnnotation, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ManualAnnotation) %>%
  mutate(percent = n / sum(n) * 100)

p2 <- ggplot(df_locations,
             aes(x = percent,
                 y = ManualAnnotation,
                 fill = tissue)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("Tumor" = "red",
               "Adjacent_normal" = "blue",
               "PBMC" = "yellow"),
    labels = c("Normal", "PBMC", "Tumor"),  
    name = "Tissue Type"
  ) +
  labs(title = "% Cell Types by Condition",
       x = "Percentage (%)",
       y = "Cell Type") +
  theme_classic()  +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

p2

df_conditions <- meta %>%
  filter(tissue %in% c("Tumor", "Adjacent_normal")) %>%
  group_by(ManualAnnotation, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ManualAnnotation) %>%
  mutate(percent = n / sum(n) * 100)

p3 <- ggplot(df_conditions,
             aes(x = percent,
                 y = ManualAnnotation,
                 fill = tissue)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Tumor" = "red",
                               "Adjacent_normal" = "blue"),
    labels = c("Normal", "Tumor"),  
    name = "Tissue Type")+
  labs(title = "% Cell Types by Condition",
       x = "Percentage (%)",
       y = "Cell Type") +
  theme_classic()  +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

p3

################################################################################

meta <- as.data.frame(combined@meta.data)

meta <- meta %>%
  filter(tissue %in% c("Tumor", "Adjacent_normal"),
         !is.na(ManualAnnotation)) %>%
  mutate(
    tumor_status = ifelse(tissue == "Tumor", 1, 0),
    patient = factor(patients)
  )

cell_types <- unique(meta$ManualAnnotation)

results_list <- list()

#----------------------------------------
# Fit GLMM per cell type
#----------------------------------------
for (ct in cell_types) {
  
  df_ct <- meta %>%
    mutate(success = ifelse(ManualAnnotation == ct, 1, 0)) %>%
    group_by(patient, tumor_status) %>%
    summarise(
      successes = sum(success),
      total = n(),
      failures = total - successes,
      .groups = "drop"
    )
  
  # Skip if no variation
  if (length(unique(df_ct$successes)) <= 1) next
  
  model <- glmer(
    cbind(successes, failures) ~ tumor_status + (1 | patient),
    data = df_ct,
    family = binomial
  )
  
  tidy_res <- tidy(model, effects = "fixed", conf.int = TRUE)
  
  tumor_row <- tidy_res %>%
    filter(term == "tumor_status") %>%
    mutate(
      cell_type = ct,
      OR = exp(estimate),
      CI_low = exp(conf.low),
      CI_high = exp(conf.high)
    )
  
  results_list[[ct]] <- tumor_row
}

results_df <- bind_rows(results_list)

#----------------------------------------
# Multiple testing correction
#----------------------------------------
results_df <- results_df %>%
  mutate(
    FDR = p.adjust(p.value, method = "fdr")
  )
manual_order <- c(
  "Specialized immune cells",
  "Neuronal-like cells",
  "Mast cells",
  "Platelets",
  "Endothelial cells",
  "Monocytes",
  "Macrophages",
  "Fibroblasts",
  "Epithelial cells",
  "B cells",
  "NK cells",
  "T cells"
)

results_df$cell_type <- factor(
  results_df$cell_type,
  levels = manual_order
)

results_df <- results_df %>%
  mutate(
    FDR = p.adjust(p.value, method = "fdr"),
    
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    direction = ifelse(OR > 1,
                       "Enriched in Tumor",
                       "Depleted in Tumor")
  )

#----------------------------------------
# Forest Plot
#----------------------------------------
ggplot(results_df,
       aes(x = OR,
           y = cell_type)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_low,
                     xmax = CI_high),
                 height = 0.2) +
  scale_x_log10() +
  labs(x = "Odds Ratio (Tumor vs Normal)",
       y = "Cell Type") +
  theme_classic()  +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )


ggplot(results_df,
       aes(x = OR,
           y = cell_type,
           color = direction)) +
  
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  
  geom_point(size = 3) +
  
  geom_errorbarh(aes(xmin = CI_low,
                     xmax = CI_high),
                 height = 0.2) +
  
  geom_text(aes(label = significance),
            hjust = -0.4,
            size = 6,
            color = "black") +
  
  scale_color_manual(values = c(
    "Enriched in Tumor" = "red",
    "Depleted in Tumor" = "blue"
  )) +
  
  scale_x_log10() +
  
  labs(x = "Odds Ratio (Tumor vs Normal)",
       y = "Cell Type",
       color = "Direction") +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

################################################
#CellChat
################################################

combined
colnames(combined@meta.data)
head(combined@meta.data)
levels(Idents(combined))

# Relabel clusters 23 and 27 in ManualAnnotation
combined$ManualAnnotation <- as.character(combined$ManualAnnotation)

combined$ManualAnnotation[combined$seurat_clusters == 23] <- "Schwann cells"
combined$ManualAnnotation[combined$seurat_clusters == 27] <- "Neuroendocrine cells"

combined$ManualAnnotation <- as.factor(combined$ManualAnnotation)

# Update Idents to match
Idents(combined) <- "ManualAnnotation"

# Confirm
table(combined$ManualAnnotation)

# Check cell type distribution per condition
table(combined$Sample_Group, combined$ManualAnnotation)

###
saveRDS(combined, file = "updated_pre_CellChat.rds")
setwd("C:/Users/huybu/OneDrive/Desktop/Dayton/Research/PDAC")
combined <- readRDS("updated_pre_CellChat.rds")

# Packages
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jinworks/CellChat")

library(CellChat)
library(Seurat)
library(tidyverse)

# Subset to Tumor and Normal conditions only, dropping PBMC
combined_sub <- subset(combined, Sample_Group %in% c("Tumor", "Normal"))

# Check how many NA cells exist in tumor
sum(is.na(combined_sub$ManualAnnotation))

# Remove cells with NA annotation
combined_sub <- subset(combined_sub, 
                       cells = colnames(combined_sub)[!is.na(combined_sub$ManualAnnotation)])

# Confirm no NAs remain
sum(is.na(combined_sub$ManualAnnotation))

# Split into two separate Seurat objects by condition
seurat_tumor  <- subset(combined_sub, Sample_Group == "Tumor")
seurat_normal <- subset(combined_sub, Sample_Group == "Normal")

# Extract normalized data and cell type labels for each condition
# Tumor
cellchat_tumor <- createCellChat(object = seurat_tumor,
                                 group.by = "ManualAnnotation",
                                 assay = "RNA")

# Normal
cellchat_normal <- createCellChat(object = seurat_normal,
                                  group.by = "ManualAnnotation",
                                  assay = "RNA")

# Confirm
cellchat_tumor
cellchat_normal

# Load the human CellChat database
CellChatDB <- CellChatDB.human

# Show what's available in the database
showDatabaseCategory(CellChatDB)

# Subset the database to focus on Secreted Signaling only
# This is the most relevant for TME interactions
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

# Assign the database to both objects
cellchat_tumor@DB  <- CellChatDB.use
cellchat_normal@DB <- CellChatDB.use

# Preprocess for each condition

# Tumor
cellchat_tumor <- subsetData(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedGenes(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedInteractions(cellchat_tumor)

# Normal
cellchat_normal <- subsetData(cellchat_normal)
cellchat_normal <- identifyOverExpressedGenes(cellchat_normal)
cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

# Drop unused factor levels for tumor
cellchat_tumor@idents <- droplevels(cellchat_tumor@idents, 
                                    exclude = setdiff(levels(cellchat_tumor@idents), 
                                                      unique(cellchat_tumor@idents)))

# Drop unused factor levels for normal
cellchat_normal@idents <- droplevels(cellchat_normal@idents, 
                                     exclude = setdiff(levels(cellchat_normal@idents), 
                                                       unique(cellchat_normal@idents)))

# Confirm levels look clean
levels(cellchat_tumor@idents)
levels(cellchat_normal@idents)

# Compute communication probabilities for each condition

# Tumor
cellchat_tumor <- computeCommunProb(cellchat_tumor,
                                    type = "triMean")

# Normal
cellchat_normal <- computeCommunProb(cellchat_normal,
                                     type = "triMean")

# Filter out communications with very few cells
# (minimum 10 cells per cell group — CellChat default)

cellchat_tumor  <- filterCommunication(cellchat_tumor, min.cells = 10)
cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)

# Compute pathway-level communication probabilities
# This aggregates individual L-R interactions into signaling pathways
cellchat_tumor  <- computeCommunProbPathway(cellchat_tumor)
cellchat_normal <- computeCommunProbPathway(cellchat_normal)

# Calculate aggregated interaction weights per cell type
cellchat_tumor  <- aggregateNet(cellchat_tumor)
cellchat_normal <- aggregateNet(cellchat_normal)

# Merge tumor and normal CellChat objects for comparative analysis
object.list <- list(Normal = cellchat_normal, Tumor = cellchat_tumor)

cellchat_merged <- mergeCellChat(object.list, 
                                 add.names = names(object.list))

# Confirm
cellchat_merged

# Compare total number and strength of interactions
compareInteractions(cellchat_merged, 
                    show.legend = FALSE,
                    group = c(1, 2))

# Compare interaction strength (weight) instead of count
compareInteractions(cellchat_merged, 
                    show.legend = FALSE,
                    group = c(1, 2),
                    measure = "weight")

# Define your populations of interest
pops_of_interest <- c("Schwann cells", "Neuroendocrine cells")

# Differential number of interactions — Tumor vs Normal
# Positive = gained in Tumor, Negative = lost in Tumor
netVisual_diffInteraction(cellchat_merged,
                          weight.scale = TRUE,
                          label.edge = FALSE,
                          sources.use = pops_of_interest,
                          targets.use = pops_of_interest)

# Differential interaction strength — Tumor vs Normal
netVisual_diffInteraction(cellchat_merged,
                          weight.scale = TRUE,
                          label.edge = FALSE,
                          measure = "weight",
                          sources.use = pops_of_interest,
                          targets.use = pops_of_interest)


# Compute network centrality scores for each condition
cellchat_tumor  <- netAnalysis_computeCentrality(cellchat_tumor)
cellchat_normal <- netAnalysis_computeCentrality(cellchat_normal)

# Update the merged object
object.list <- list(Normal = cellchat_normal, Tumor = cellchat_tumor)
cellchat_merged <- mergeCellChat(object.list, 
                                 add.names = names(object.list))

# Extract all significant interactions where Schwann cells are senders
# and T cells are receivers, comparing Tumor vs Normal

netAnalysis_signalingRole_scatter(cellchat_merged)

# Extract ligand-receptor pairs between Schwann and T cells
# for each condition separately

# Tumor
df_tumor_sw_t <- subsetCommunication(cellchat_tumor,
                                     sources.use = "Schwann cells",
                                     targets.use = "T cells")

# Normal  
df_normal_sw_t <- subsetCommunication(cellchat_normal,
                                      sources.use = "Schwann cells",
                                      targets.use = "T cells")

# View results
head(df_tumor_sw_t)
head(df_normal_sw_t)