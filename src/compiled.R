# setup environment ----------------------------------------------------------

library(tidyverse)
library(Seurat)
# source('~/westerlund/R/functions.R')
# source('~/westerlund/R/seurat_wrappers.R')
source('~/Projects/westerlund/R/utilities.R')
source('~/Projects/westerlund/R/plotting.R')
# source('~/westerlund/R/granuloma_specific.R')
set.seed(1)
gc()

# data preprocessing -----------------------------------------------------------

# load in barcodes from lane A
bctbl <- readRDS("data/raw/A_barcodetable.rds")
colSums(bctbl)
a_bc_counts <- bctbl[, c(14:19, 26:31)]
a_bc_counts <- Matrix::Matrix(as.matrix(a_bc_counts), sparse = TRUE)
colnames(a_bc_counts) <- paste0("BC", 1:12)
a_bc_counts[1:5, 1:5]
rownames(a_bc_counts) <- paste0(rownames(a_bc_counts), "-1")

# load in barcodes from lane B
bctbl <- readRDS("data/raw/B_barcodetable.rds")
b_bc_counts <- bctbl[,  c(14:19, 26:31)]
b_bc_counts <- Matrix::Matrix(as.matrix(b_bc_counts), sparse = TRUE)
colnames(b_bc_counts) <- paste0("BC", 1:12)
b_bc_counts[1:5, 1:5]
rownames(b_bc_counts) <- paste0(rownames(b_bc_counts), "-1")

# merge barcodes
bc_counts <- rbind(a_bc_counts, b_bc_counts)

# load in counts
a_counts <- Read10X_h5("data/raw/A_filtered_feature_bc_matrix.h5")
table(rownames(a_bc_counts) %in% colnames(a_counts))
b_counts <- Read10X_h5("data/raw/B_filtered_feature_bc_matrix.h5")
table(rownames(b_bc_counts) %in% colnames(b_counts))

# take intersection of LMO BCs and count BCs
a_overlap <- intersect(rownames(a_bc_counts), colnames(a_counts))
b_overlap <- intersect(rownames(b_bc_counts), colnames(b_counts))

# subset all matrices
a_counts <- a_counts[, a_overlap]
a_bc_counts <- a_bc_counts[a_overlap, ]
b_counts <- b_counts[, b_overlap]
b_bc_counts <- b_bc_counts[b_overlap, ]

# rename barcodes
colnames(a_counts) <- substr(colnames(a_counts), 1, 16)
colnames(b_counts) <- substr(colnames(b_counts), 1, 16)
rownames(a_bc_counts) <- substr(rownames(a_bc_counts), 1, 16)
rownames(b_bc_counts) <- substr(rownames(b_bc_counts), 1, 16)

colnames(a_counts) <- paste0(colnames(a_counts), "_", "A")
colnames(b_counts) <- paste0(colnames(b_counts), "_", "B")
rownames(a_bc_counts) <- paste0(rownames(a_bc_counts), "_", "A")
rownames(b_bc_counts) <- paste0(rownames(b_bc_counts), "_", "B")

# merge matrices
bc_counts <- rbind(a_bc_counts, b_bc_counts)
gex_counts <- RowMergeSparseMatrices(a_counts, b_counts)

# object preprocessing -----------------------------------------------------------

# create object on LMO counts first
object <- CreateSeuratObject(counts = Matrix::t(bc_counts), project = "AB", assay = "LMO",
                             min.cells = 0, min.features = 0, names.delim = "_", names.field = 2)

# add GEX counts
object[["RNA"]] <- CreateAssayObject(gex_counts, min.cells = 0, min.features = 0)
object$cell_barcode <- Cells(object)
object$lane <- object$orig.ident

# demux LMO barcodes -------------------------------------------------------------------
Idents(object) <- "lane"
a <- subset(object, idents = "A")
a <- NormalizeData(a, assay = "LMO", normalization.method = "CLR")
a <- HTODemux(a, assay = "LMO", positive.quantile = 0.999, init = 13)
table(a$LMO_classification.global)

b <- subset(object, idents = "B")
b <- NormalizeData(b, assay = "LMO", normalization.method = "CLR")
b <- HTODemux(b, assay = "LMO", positive.quantile = 0.999, init = 13)
table(b$LMO_classification.global)

object$sample_barcode <- "Negative"
object$sample_barcode[Cells(a)] <- as.character(a$hash.ID)
object$sample_barcode[Cells(b)] <- as.character(b$hash.ID)
table(object$sample_barcode)

rm(a, b)

# add sample metadata ----------------------------------------------------------------

meta <- read_csv("data/raw/sample_metadata.csv")
meta <- merge(object[[]] %>% select(cell_barcode, sample_barcode), meta, by.x = "sample_barcode", by.y = "barcode")
meta <- meta %>% select(-sample_barcode)
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_barcode
meta
object <- AddMetaData(object, meta[, -1])
object[[]]

object
object <- object[, !is.na(object$unique_sample_id) & object$mtb_condition != "PDIMKO"]
table(object$unique_sample_id)
saveRDS(object, "data/processed/base_object.rds")

# ambient gene identification --------------------------------------------------------------

counts <- list(a_counts, b_counts)
gene_sums <- data.frame(total_umis = rowSums(gex_counts), gene = rownames(gex_counts))

amb <- pbapply::pblapply(counts, function(x) {
  amb <- DropletUtils::estimateAmbience(x)
  amb <- qdapTools::vect2df(amb)
  colnames(amb) <- c("gene", "est_amb")
  return(amb)
})

amb[[1]]$lane <- "A"
amb[[2]]$lane <- "B"
amb <- bind_rows(amb)
hist(log10(amb$est_amb))

threshold <- quantile(amb$est_amb, 0.997)
gene_stats <- amb %>%
  group_by(gene) %>%
  summarize(mean_est_amb = mean(est_amb), n = sum(est_amb >= threshold))
gene_stats <- merge(gene_stats, gene_sums, by = "gene")
gene_stats
a <- gene_stats %>%
  ggplot(aes(x = mean_est_amb, y = total_umis)) +
  geom_point(size = 0.5, alpha = 0.5) +
  #geom_hline(yintercept = 8) +
  geom_vline(xintercept = threshold) +
  #scale_x_sqrt(breaks = c(0.0001, 0.001, seq(0.002, 0.008, 0.002))) +
  scale_x_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Mean estimated ambient fraction", y = "Total detected UMIs per gene") +
  GeneralTheme(18) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles()
a
SavePlot(a, filename = "ambient_genes", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

write_csv(gene_stats, "data/reported/ambient_gene_stats.csv")
saveRDS(amb, "data/processed/amb_results.rds")

ambient_genes <- as.character(amb$gene[amb$est_amb >= threshold])
saveRDS(ambient_genes, "data/processed/ambient_genes.rds")

# filtering and doublet detection -----------------------------------------------------------------

DefaultAssay(object) <- "RNA"
Idents(object) <- "unique_sample_id"

table(Idents(object))
<- obj <- object
obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent_ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp(s|l)")
VlnPlot(obj, features = c("percent_mito", "percent_ribo"))
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

b <- ggplot(obj[[]], aes(x = percent_mito, y = percent_ribo, color = lane)) +
  geom_point(size = 0.5, alpha = 0.5) +
  ggthemes::scale_color_ptol(name = "Lane") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black", size = 1) +
  labs(x = "% mitochondrial UMIs", y = "% ribosomal UMIs") +
  GeneralTheme(18) +
  SpaceAxisTitles() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.9, 0.9))
b
SavePlot(b, filename = "percent_tech_genes", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

coef(lm(nFeature_RNA ~ nCount_RNA, data = obj[[]]))
c <- ggplot(obj[[]], aes(x = nCount_RNA, y = nFeature_RNA, color = lane)) +
  geom_point(size = 0.5, alpha = 0.5) +
  ggthemes::scale_color_ptol(name = "Lane") +
  geom_hline(yintercept = 300, linetype = "dashed", color = "black", size = 1) +
  #geom_abline(slope = 0.17, intercept = 200, linetype = "dashed", color = "black", size = 1) +
  labs(x = "# UMIs", y = "# genes detected") +
  scale_x_log10() +
  scale_y_log10() +
  GeneralTheme(18) +
  SpaceAxisTitles() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = "none")
c
SavePlot(c, filename = "library_sizes", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

obj <- subset(obj, subset = nFeature_RNA >= 300 & percent_mito <= 10)
obj

scrubbed <- Scrub(obj, batch.var = "lane")
scrubbed_scores <- scrubbed[[]] %>% select(scrublet_score, scrublet_label)
obj <- AddMetaData(obj, scrubbed_scores)
obj$scrublet_score <- as.numeric(obj$scrublet_score)
obj$scrublet_label <- as.logical(obj$scrublet_label)

saveRDS(obj, "data/processed/filtered_object.rds")

head(obj[[]])
sample_meta <- obj[[]] %>%
  mutate(unique_sample_id = paste0(unique_sample_id, lane)) %>%
  group_by(unique_sample_id) %>%
  summarize(
    `Lane` = unique(lane),
    `Mouse strain` = unique(mouse_condition),
    `Infection condition` = unique(mtb_condition),
    `Number of cells` = n(),
    `Mean UMIs` = mean(nCount_RNA),
    `Mean genes detected` = mean(nFeature_RNA),
    `Mean % mito. UMIs` = mean(percent_mito),
    `Mean % ribo. UMIs` = mean(percent_ribo))
write_csv(sample_meta, file = "data/reported/sample_meta.csv")

# dimensionality reduction and clustering ---------------------------------

DefaultAssay(obj) <- "RNA"
obj <- ReduceDims(obj, n.var.features = 3000, n.pcs = 50, remove.genes = ambient_genes)
obj <- Embed(seurat.obj = obj, reduction = "pca", dims.use = 30, knn.use = 20)
obj <- Walktrap(obj, col.name = "pca_walktrap", graph.name = "pca_snn_20")
obj <- SetIdent(obj, value = "pca_walktrap")
obj <- BuildClusterTree(obj, reorder = TRUE, reorder.numeric = TRUE)
obj[["pca_walktrap"]] <- Idents(obj)

d <- DimPlot(object = obj, group.by = "pca_walktrap", label = TRUE, raster = F) +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + NoLegend() +
  colorspace::scale_color_discrete_qualitative("Dark 3") +
  labs(title = "Walktrap clusters") +
  theme(
    axis.text = element_blank(),
    axis.title.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank())
d
SavePlot(d, filename = "walktrap_umap", root = "", h = 5, w = 5, s = 1, save.data = FALSE)


# obj <- AutoCluster(seurat.obj = obj, graph.name = "pca_snn_20", col.name = "pca_leiden", mod.percentile = 1)
# obj <- SetIdent(obj, value = "pca_leiden")
# obj <- BuildClusterTree(obj, reorder = TRUE, reorder.numeric = TRUE)
# obj[["pca_leiden"]] <- Idents(obj)

DimPlot(obj, group.by = c("pca_walktrap", "pca_leiden"), label = TRUE) + NoLegend()

Idents(obj) <- "pca_walktrap"
markers <- presto::wilcoxauc(obj, group_by = "pca_walktrap", seurat_assay = "RNA", assay = "data")
top_markers <- markers %>%
  group_by(feature) %>%
  top_n(1, avgExpr) %>%
  group_by(group) %>%
  top_n(10, auc) %>%
  arrange(as.numeric(group), desc(auc))
VlnPlot(obj, features = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"),
        sort = TRUE, pt.size = 0, ncol = 2)

table(obj$pca_walktrap)

plots <- lapply(levels(obj), function(x) {
  p <- DimPlot(object = obj, group.by = "pca_walktrap", label = FALSE, raster = F,
               sizes.highlight = 1, cols.highlight = c("black"),
               cells.highlight = obj$cell_barcode[obj$pca_walktrap == x]) +
    GeneralTheme(12) + RemoveBackgrounds() + SpaceAxisTitles() + NoLegend() +
    labs(title = x) +
    theme(
      axis.text = element_blank(),
      axis.title.x  = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank())
  return(p)
})
e <- wrap_plots(plots, ncol = 7, nrow = 6)
e
SavePlot(e, filename = "walktrap_clusterspecific_umap", root = "", h = 12, w = 12, s = 1, save.data = FALSE)

# transfer tabula muris data ----------------------------------------------

tabulamuris <- readRDS("data/external/tabulamuris.rds")
DimPlot(tabulamuris, group.by = "cell_ontology_class", split = "method", label = TRUE) + NoLegend()
tabulamuris <- tabulamuris[, Cells(tabulamuris)[tabulamuris$method == "droplet" & !is.na(tabulamuris$cell_ontology_class)]]
anchors <- FindTransferAnchors(reference = tabulamuris,
                               query = obj,
                               dims = 1:30,
                               n.trees = 20)
scores <- MappingScore(anchors, ndim = 20)
predictions <- TransferData(anchorset = anchors,
                            refdata = tabulamuris$cell_ontology_class,
                            n.trees = 20,
                            dims = 1:30)
predictions$mapping_score <- scores
obj <- AddMetaData(obj, metadata = predictions)
View(as.data.frame(table(obj$predicted.id, obj$pca_walktrap)))

DimPlot(obj, group.by = c("pca_walktrap", "predicted.id"), label = TRUE, repel = TRUE, ncol = 1) + NoLegend()

xvar = "predicted.id"
fillvar = "pca_walktrap"
res_props <- obj[[]] %>%
  dplyr::group_by(!!sym(xvar), !!sym(fillvar)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::group_by(!!sym(fillvar), !!sym(xvar))
head(res_props)

e2 <- ggplot(res_props, aes(x = prop, y = predicted.id)) +
  geom_point() +
  facet_wrap(~ pca_walktrap) +
  GeneralTheme(6) +
  theme(panel.grid.major.y = element_line(color = "gray90", size = 0.5))
e2
SavePlot(e2, filename = "tabulamuris_prediction_props", root = "", h = 12, w = 12, s = 1, save.data = FALSE)

# annotation of cell types --------------------------------------------------------------

annotations <- readxl::read_xlsx("data/processed/celltype_assignments.xlsx")
annotations <- merge(annotations, obj[[]][, c("pca_walktrap", "cell_barcode")])
rownames(annotations) <- annotations$cell_barcode
obj <- AddMetaData(obj, metadata = annotations[, -4])
DimPlot(obj, group.by = "cell_type_l2", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()

# run decontX for poor cell identification -------------------------------------------------------------

library(celda)
decon_results <- decontX(x = obj@assays$RNA@counts, z = obj$cell_type_l2)
obj$contam_fraction <- decon_results$contamination
f <- FeaturePlot(obj, features = "contam_fraction", reduction = "pca_umap", cols = c("gray80", "black")) +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + RemoveAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "DecontX contamination fraction")
f
SavePlot(f, filename = "decontx_contamination_umap", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

VlnPlot(obj, features = "contam_fraction", group.by = "pca_walktrap")
obj[[]] %>% group_by(pca_walktrap) %>% summarize(mean = mean(contam_fraction)) %>% arrange(desc(mean))

quantile(obj$contam_fraction, seq(0, 1, 0.1))
table(obj$contam_fraction >= 0.25)

# add decontX counts
obj[["DECON"]] <- CreateAssayObject(counts = decon_results$decontXcounts)
saveRDS(decon_results, file = "data/processed/decontx.rds")

# add decontX counts
obj[["DECON"]] <- CreateAssayObject(counts = decon_results$decontXcounts)
DefaultAssay(obj) <- "DECON"
obj <- ReduceDims(obj, n.var.features = 3000, n.pcs = 50, pca.name = "decon_pca", remove.genes = ambient_genes)
obj <- Embed(seurat.obj = obj, reduction = "decon_pca", dims.use = 30, knn.use = 20)
obj <- Walktrap(obj, col.name = "decon_walktrap", graph.name = "decon_pca_snn_20")
obj <- SetIdent(obj, value = "decon_walktrap")
obj <- BuildClusterTree(obj, reorder = TRUE, reorder.numeric = TRUE)
obj[["decon_walktrap"]] <- Idents(obj)
DimPlot(obj, group.by = "decon_walktrap", reduction = "decon_pca_umap", label = TRUE, repel = FALSE) + NoLegend()


?# reannotate --------------------------------------------------------------

table(obj$cell_type_l2)
sub <- obj[, obj$cell_type_l2 %in% c("T cell 1", "T cell 2", "T cell 3", "T cell 4",
                                     "T cell 5", "NK cell 1", "NK cell 2", "NK cell 3",
                                     "Proliferating lymphocyte", "B cell 1", "B cell 2")]
DefaultAssay(sub) <- "RNA"
sub <- ReduceDims(sub, n.var.features = 3000, n.pcs = 50, remove.genes = ambient_genes)
sub <- Embed(seurat.obj = sub, reduction = "pca", dims.use = 30, knn.use = 20)
sub <- Walktrap(sub, col.name = "lymphocyte_walktrap", graph.name = "pca_snn_20")
sub <- SetIdent(sub, value = "lymphocyte_walktrap")
sub <- BuildClusterTree(sub, reorder = TRUE, reorder.numeric = TRUE)
sub[["lymphocyte_walktrap"]] <- Idents(sub)

DimPlot(sub, group.by = c("cell_type_l2", "lymphocyte_walktrap"), label = TRUE) + NoLegend()
VlnPlot(sub, features = c("Cd8a", "Cd3g", "Cd4", "Klrb1", "Klrc1"), stack = TRUE, flip = TRUE)
markers <- presto::wilcoxauc(sub, group_by = "lymphocyte_walktrap", seurat_assay = "RNA", assay = "data")
top_markers <- markers %>%
  group_by(feature) %>%
  top_n(1, avgExpr) %>%
  group_by(group) %>%
  top_n(10, auc) %>%
  arrange(as.numeric(group), desc(auc))

annotations <- readxl::read_xlsx("data/processed/celltype_assignments.xlsx", sheet = 2)
annotations <- merge(annotations, sub[[]][, c("lymphocyte_walktrap", "cell_barcode")])
rownames(annotations) <- annotations$cell_barcode
sub <- AddMetaData(sub, metadata = annotations[, -4])
DimPlot(sub, group.by = "cell_type_l3", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()

obj$cell_type_l3 <- obj$cell_type_l2
obj$cell_type_l3[sub$cell_barcode] <- sub$cell_type_l3
DimPlot(obj, group.by = "cell_type_l3", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()

table(obj$cell_type_l2)
sub <- obj[, obj$cell_type_l2 %in% c("Alveolar macrophage", "Proliferating alveolar macrophage", "Classical monocyte",
                                     "Non-classical monocyte 1", "Non-classical monocyte 2", "Non-classical monocyte 3",
                                     "Myeloid cell 1", "Myeloid cell 2", "Neutrophil", "Dendritic cell", "Recruited macrophage")]
DefaultAssay(sub) <- "RNA"
sub <- ReduceDims(sub, n.var.features = 3000, n.pcs = 50, remove.genes = ambient_genes)
sub <- Embed(seurat.obj = sub, reduction = "pca", dims.use = 30, knn.use = 20)
sub <- Walktrap(sub, col.name = "myeloid_walktrap", graph.name = "pca_snn_20")
sub <- SetIdent(sub, value = "myeloid_walktrap")
sub <- BuildClusterTree(sub, reorder = TRUE, reorder.numeric = TRUE)
sub[["myeloid_walktrap"]] <- Idents(sub)

DimPlot(sub, group.by = c("cell_type_l2", "myeloid_walktrap"), label = TRUE) + NoLegend()
VlnPlot(sub, features = c("Csf3r", "Mrc1", "Ly6c", "C1qa"), stack = TRUE, flip = TRUE)
VlnPlot(sub, features = c("nCount_RNA"), stack = F, flip = F)
markers <- presto::wilcoxauc(sub, group_by = "myeloid_walktrap", seurat_assay = "RNA", assay = "data")
top_markers <- markers %>%
  group_by(feature) %>%
  top_n(1, avgExpr) %>%
  group_by(group) %>%
  top_n(10, auc) %>%
  arrange(as.numeric(group), desc(auc))

annotations <- readxl::read_xlsx("data/processed/celltype_assignments.xlsx", sheet = 3)
annotations <- merge(annotations, sub[[]][, c("myeloid_walktrap", "cell_barcode")], by = "myeloid_walktrap")
rownames(annotations) <- annotations$cell_barcode
sub <- AddMetaData(sub, metadata = annotations[, -3])
DimPlot(sub, group.by = "cell_type_l3", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()

obj$cell_type_l3[sub$cell_barcode] <- sub$cell_type_l3
DimPlot(obj, group.by = "cell_type_l3", label = TRUE, repel = FALSE, ncol = 1) + NoLegend()

# umap plotting -----------------------------------------------------------

g <- DimPlot(object = obj, group.by = "cell_type_l1", reduction = "pca_umap",
             label = F, raster = F, cells = obj$cell_barcode[obj$cell_type_l1 != "Nondescript"]) +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + RemoveAxes() +
  colorspace::scale_color_discrete_qualitative("Dark 3") +
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  theme(legend.text = element_text(size = 10))
g
SavePlot(g, filename = "cell_type_l1_umap", root = "", h = 5, w = 7, s = 1, save.data = FALSE)

h <- DimPlot(object = obj, group.by = "cell_type_l2", reduction = "pca_umap",
             label = F, raster = F, cells = obj$cell_barcode[obj$cell_type_l2 != "Nondescript"]) +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + RemoveAxes() +
  colorspace::scale_color_discrete_qualitative("Dark 3") +
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  theme(legend.text = element_text(size = 10))
h
SavePlot(h, filename = "cell_type_l2_umap", root = "", h = 5, w = 10, s = 1, save.data = FALSE)

h2 <- DimPlot(object = obj, group.by = "cell_type_l3", reduction = "pca_umap",
             label = F, raster = F, cells = obj$cell_barcode[obj$cell_type_l3 != "Nondescript"]) +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + RemoveAxes() +
  colorspace::scale_color_discrete_qualitative("Dark 3") +
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  theme(legend.text = element_text(size = 10))
h2
SavePlot(h2, filename = "cell_type_l3_umap", root = "", h = 5, w = 10, s = 1, save.data = FALSE)

mouse = "TLR2KO"
mtb = "H37Rv"
sub_obj <- obj[, obj$cell_barcode[obj$cell_type_l3 != "Nondescript"]]
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$mouse_condition == mouse & sub_obj$mtb_condition == mtb]
umap_df <- as.data.frame(sub_obj@reductions$pca_umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$pUMAP_1, umap_df_sub$pUMAP_2, n = 200)
i <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = pUMAP_1, y = pUMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density") +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = glue("{mouse} {mtb}"))
SavePlot(plot = i, filename = glue("density_umap_{mouse}_{mtb}"), w = 6, h = 5, s = 1, save.data = FALSE, root = "")

mouse = "TLR2KO"
mtb = "Uninfected"
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$mouse_condition == mouse & sub_obj$mtb_condition == mtb]
umap_df <- as.data.frame(sub_obj@reductions$pca_umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$pUMAP_1, umap_df_sub$pUMAP_2, n = 200)
j <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = pUMAP_1, y = pUMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density") +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = glue("{mouse} {mtb}"))
SavePlot(plot = j, filename = glue("density_umap_{mouse}_{mtb}"), w = 6, h = 5, s = 1, save.data = FALSE, root = "")

mouse = "Null"
mtb = "H37Rv"
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$mouse_condition == mouse & sub_obj$mtb_condition == mtb]
umap_df <- as.data.frame(sub_obj@reductions$pca_umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$pUMAP_1, umap_df_sub$pUMAP_2, n = 200)
k <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = pUMAP_1, y = pUMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density") +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = glue("{mouse} {mtb}"))
SavePlot(plot = k, filename = glue("density_umap_{mouse}_{mtb}"), w = 6, h = 5, s = 1, save.data = FALSE, root = "")

mouse = "Null"
mtb = "Uninfected"
cells_to_highlight <- sub_obj$cell_barcode[sub_obj$mouse_condition == mouse & sub_obj$mtb_condition == mtb]
umap_df <- as.data.frame(sub_obj@reductions$pca_umap@cell.embeddings)
umap_df_sub <- umap_df[cells_to_highlight, ]
umap_df_sub$density <- get_density(umap_df_sub$pUMAP_1, umap_df_sub$pUMAP_2, n = 200)
l <- ggplot() +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "black", size = 2) +
  ggrastr::geom_point_rast(data = umap_df, mapping = aes(x = pUMAP_1, y = pUMAP_2), color = "gray95", size = 1) +
  ggrastr::geom_point_rast(data = umap_df_sub, mapping = aes(x = pUMAP_1, y = pUMAP_2, color = density), size = 1) +
  GeneralTheme(18) +
  scale_color_viridis_c(option = "G", direction = -1, name = "Density") +
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", title = glue("{mouse} {mtb}"))
SavePlot(plot = l, filename = glue("density_umap_{mouse}_{mtb}"), w = 6, h = 5, s = 1, save.data = FALSE, root = "")

# find and save markers ---------------------------------------------------

obj <- readRDS("data/processed/annotated_object.rds")
markers <- presto::wilcoxauc(obj, group_by = "cell_type_l3", seurat_assay = "RNA", assay = "data")
top_markers <- markers %>%
  filter(padj <= 0.05 & auc >= 0.6)
write_csv(top_markers, file = "data/reported/celltype_l3_markers.csv")

# plot markers
use_obj <- obj[, obj$cell_barcode[obj$cell_type_l3 != "Nondescript"]]
obj$cell_type_l3[obj$cell_type_l3 == "Lung endothelial cell"] <- "Nondescript"
use_obj <- obj[, obj$cell_barcode[obj$cell_type_l3 != "Nondescript"]]

Idents(use_obj) <- "cell_type_l3"
use_obj <- BuildClusterTree(use_obj, reorder = TRUE, reorder.numeric = F)
use_obj[["cell_type_l3"]] <- Idents(use_obj)

cell_type_order <- levels(Idents(use_obj))
markers <- presto::wilcoxauc(use_obj, group_by = "cell_type_l3", seurat_assay = "RNA", assay = "data")
subset_markers <- markers %>%
  filter(!grepl("^Hsp|^Dna|^mt-", feature) & logFC >= 0.1) %>%
  group_by(feature) %>%
  top_n(1, avgExpr) %>%
  group_by(group) %>%
  top_n(1, logFC) %>%
  arrange(group, desc(logFC))
markers_df <- markers[markers$feature %in% subset_markers$feature, ]
a <- markers_df %>% mutate(feature = factor(feature, levels = subset_markers$feature)) %>%
  ggplot() +
  geom_tile(aes(x = feature, y = group, fill = auc), color = "black", size = 0.25) +
  geom_point(aes(x = feature, y = group, size = pct_in), shape = 21, fill = "transparent", color = "black") +
  colorspace::scale_fill_continuous_diverging("Blue-Red 3", name = "AUROC", mid = 0.5) +
  #colorspace::scale_fill_continuous_sequential("Grays", name = "logFC") +
  scale_size_continuous(range = c(0.1, 4), name = "% detected") +
  #colorspace::scale_fill_continuous_sequential("Grays", name = "Average\nExpr.") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "") +
  GeneralTheme(14) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, angle = 90),
        axis.ticks.y = element_blank())
a
SavePlot(plot = a, filename = "cell_type_l3_markers", w = 12, h = 6, s = 1, save.data = FALSE, root = "")

# save object -------------------------------------------------------------

saveRDS(obj, "data/processed/annotated_object.rds")
obj <- readRDS("data/processed/annotated_object.rds")

# proportion plotting ---------------------------------------------------------------------

obj$condition <- paste(obj$mouse_condition, obj$mtb_condition, sep = " ")
table(obj$condition)
obj$condition <- as.factor(obj$condition)
levels(obj$condition) <- c("WT H37Rv", "WT Uninfected", "Tlr2-/- H37Rv", "Tlr2-/- Uninfected")

counts <- obj[[]] %>% group_by(lane, unique_sample_id, condition) %>% summarize(n = n())
m <- counts %>% unite("sample_id", lane, unique_sample_id) %>%
  ggplot() +
  geom_boxplot(aes(x = n, y = condition), width = 0.5) +
  geom_jitter(aes(x = n, y = condition, fill = sample_id), shape = 21, color = "white", size = 4,
              position = position_jitter(height = 0.1)) +
  colorspace::scale_fill_discrete_qualitative("Dark 3", guide = FALSE) +
  labs(x = "Number of cells", y = "") +
  GeneralTheme(base_size = 18) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles()
m
SavePlot(plot = m, filename = "cell_counts", w = 6, h = 5, s = 1, save.data = FALSE, root = "")

xvar = "condition"
fillvar = "cell_type_l3"
res_props <- obj[[]] %>%
  dplyr::group_by(!!sym(xvar), !!sym(fillvar)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::group_by(!!sym(fillvar), !!sym(xvar))
head(res_props)

n <- ggplot() +
  geom_bar(aes_string(x = "prop", y = xvar, fill = fillvar), width = 0.8,
           data = res_props %>% filter(!!sym(fillvar) != "Nondescript"), stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  ggthemes::scale_fill_tableau("Tableau 20", name = "Cell type") +
  labs(x = "Fractional abundance", y = "") +
  GeneralTheme(base_size = 18) +
  RemoveBackgrounds(outline = F) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 0.5),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0, unit = "in"),
    panel.grid.major = element_line(color = "gray80", size = 1, linetype = "solid")
  ) + SpaceAxisTitles()
n
SavePlot(plot = n, filename = "cell_type_l1_props_barplot", w = 7, h = 4, s = 1, save.data = FALSE, root = "")

xvar <- "cell_type_l3"
yvar <- "condition"

res_props <- obj[[]] %>%
  dplyr::group_by(!!sym(xvar), !!sym(yvar)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(prop = n / sum(n))
res_props

sample_count <- obj[[]] %>%
  dplyr::select(condition, unique_sample_id) %>%
  group_by(unique_sample_id) %>%
  mutate(batch_frequency = sum(n())) %>%
  distinct()
sample_count

count_table <- obj[[]] %>%
  dplyr::select(unique_sample_id, !!sym(xvar)) %>%
  group_by(unique_sample_id, !!sym(xvar)) %>%
  dplyr::count(name = "count", .drop = FALSE)
count_table
count_table <- merge(count_table, sample_count, by = "unique_sample_id", sort = FALSE)
count_table <- count_table %>% mutate(other = batch_frequency - count)

formula = cbind(count, other) ~ cell_type_l3 * condition
model1 <- glm(formula = formula, family = 'binomial', data = count_table)
emm1 <- emmeans::emmeans(model1, specs = revpairwise ~ condition | cell_type_l3)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

emm2 <- emmeans::emmeans(model1, specs = ~ cell_type_l3)
emm2 %>%
  summary(type = 'response') %>%
  dplyr::select(cell_type_l3, prob) -> mean_probs
c_results %>% left_join(mean_probs) -> m_results

m_results$adj_p <- p.adjust(m_results$p.value, method = "fdr")
m_results$sig <- ifelse(m_results$adj_p <= 0.05, TRUE, FALSE)

write_csv(m_results, file = "data/reported/cell_abundance_lm_results.csv")

unique(m_results$contrast)
contrast_use <- unique(m_results$contrast)[2]
contrast_use
o <- m_results %>%
  mutate(logp = -log10(p.value)) %>%
  mutate(logp = ifelse(logp > 20, 20, logp)) %>%
  filter(contrast == contrast_use) %>%
  ggplot(aes(x = odds.ratio, y = logp, fill = prob, label = cell_type_l3)) +
  geom_point(shape = 21, size = 4, alpha = 0.8) +
  #geom_crossbar(aes(xmin = odds.ratio - SE, xmax = odds.ratio + SE)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray80") +
  ggrepel::geom_text_repel(size = 2, box.padding = 0.5, force = 5, max.overlaps = Inf) +
  #ylim(c(-10, 60)) +
  scale_x_continuous(trans = scales::log10_trans(), breaks = c(0, 1, 5, 10)) +
  colorspace::scale_fill_continuous_sequential("Blues 3", name = "Average\nabundance") +
  labs(x = glue::glue("Odds Ratio {contrast_use}"), y = "-log10(adj. P value)") +
  GeneralTheme(14)
plot(o)
SavePlot(plot = o, filename = "oddsratio_tlr2wt", w = 6, h = 5, s = 1, save.data = FALSE, root = "")

obj$batch <- paste(obj$unique_sample_id, obj$lane, sep = "_")
res_props <- obj[[]] %>%
  dplyr::group_by(batch, condition, cell_type_l3) %>%
  dplyr::summarise(n = n()) %>%
  group_by(batch) %>%
  dplyr::mutate(prop = n / sum(n))
res_props

o2 <- res_props %>%
  ggplot(aes(x = prop, y = condition)) +
  geom_boxplot() +
  geom_point(size = 0.5, alpha = 0.8) +
  labs(x = "Proportion of cells within each sample", y = "Condition") +
  GeneralTheme(14) +
  scale_x_continuous(trans = scales::log10_trans()) +
  facet_wrap(~ cell_type_l3, ncol = 7)
o2
SavePlot(plot = o2, filename = "cell_type_l3_props_boxplot_split", w = 1, h = 10, s = 1, save.data = FALSE, root = "")

# tlr2 plotting -----------------------------------------------------------

p <- FeaturePlot(obj, features = "Tlr2", reduction = "pca_umap", cols = c("gray80", "black"),
                 order = TRUE, min.cutoff = "q10", max.cutoff = "q90") +
  GeneralTheme(18) + RemoveBackgrounds() + SpaceAxisTitles() + RemoveAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "Tlr2 expression in WT mice")
p
SavePlot(p, filename = "tlr2_umap", root = "", h = 5, w = 6, s = 1, save.data = FALSE)

# markers between conditions ----------------------------------------------

celltypes <- unique(obj$cell_type_l3)
markers <- pbapply::pblapply(celltypes, function(x) {
  boolean_vec <- obj$cell_type_l3 == x & obj$mtb_condition == "H37Rv"
  if(sum(boolean_vec) > 99) {
    markers <- presto::wilcoxauc(obj@assays$RNA@data[, obj$cell_barcode[boolean_vec]], obj$mouse_condition[boolean_vec])
    #top_markers <- markers %>% group_by(feature) %>% filter(padj <= 0.05)
    #top_markers$celltype <- x
    markers$celltype <- x
    return(markers)
  }
})
markers <- bind_rows(markers)
write_csv(markers, "data/reported/cell_type_l3_specific_tlr2vsnull_infectionmarkers.csv")

markers <- markers %>% filter(padj <= 0.05 & !celltype %in% c("Nondescript", "Smooth muscle cell 2"))
markers_df <- markers %>% group_by(celltype) %>% summarize(n = n())
q <- ggplot(markers_df, aes(y = fct_reorder(celltype, n), x = n)) +
  geom_col(color = "black", fill = "gray80") +
  labs(x = "Number of DE genes between Tlr2-/- and WT", y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 700)) +
  GeneralTheme(18) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles()
SavePlot(q, filename = "number_degenes_percelltypel3", root = "", h = 5, w = 6, s = 1, save.data = FALSE)

# enrichment --------------------------------------------------------------

enrich_genes <- markers %>%
  filter(celltype == "Recruited macrophage" & group == "Null" & padj <= 0.05) %>%
  top_n(100, logFC) %>%
  pull(feature)

clipr::write_clip(enrich_genes)

# library(enrichR)
# dbs <- listEnrichrDbs()
# dbs_use <- dbs$libraryName[c(122, 96, 168, 149, 146)]
# enriched <- enrichr(genes = enrich_genes, databases = dbs[1])

markers <- markers %>% filter(celltype == "Recruited macrophage")
markers_df <- markers %>% filter(group == "Null") %>% mutate(neglogp = -log10(padj), sig = ifelse(padj <= 0.05, T, F))
write_csv(markers_df, "data/reported/cell_type_l3_recruitedmaconly_specific_tlr2vsnull_infectionmarkers.csv")
# markers_df <- markers_df %>% mutate(label = ifelse(feature %in% c("Sod2", "Isg15", "Irgm1", "C3ar1", "Ifi47", "Zbp1",
#                                                                   "Irf7", "Tnfaip2", "Cebpb", "Prdx5",
#                                                                   "Trem2", "Cd9", "Lgals3", "Fabp5", "Abcg1"), T, F))
markers_df <- markers_df %>% mutate(label = ifelse(feature %in% c("Sod2", "Ptgs2", "Cebpb", "Clec4e", "Tnfaip2", "Prdx5"), T, F))
r <- ggplot(markers_df, aes(x = logFC, y = neglogp, color = sig)) +
  geom_point() +
  ggrepel::geom_text_repel(data = markers_df %>% filter(label == T), min.segment.length = 0.05, box.padding = 0.5, nudge_y = 1, nudge_x = 0.4,
                           mapping = aes(x = logFC, y = neglogp, color = sig, label = feature)) +
  labs(x = "Log(Fold Change)\nTlr2-/-     <>     WT", y = "-Log10(adj. P)") +
  GeneralTheme(18) +
  scale_color_manual(values = c("gray80", "black")) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles() +
  NoLegend()
r
SavePlot(r, filename = "degenes_recruitedmacs", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

Idents(obj) <- "mouse_condition"
StartFuture()
markers <- FindMarkers(obj[, obj$cell_type_l3 == "Recruited macrophage"], ident.1 = "TLR2KO", slot = "data", return.thresh = 1,
                       logfc.threshold = 0, test.use = "wilcox", min.pct = 0, only.pos = FALSE, base = 2)
StopFuture()
markers <- markers %>% rownames_to_column("gene")
markers$padj <- p.adjust(markers$p_val, method = "fdr")
markers_df <- markers %>% mutate(neglogp = -log10(padj), sig = ifelse(padj <= 0.05, T, F))
write_csv(markers_df, "data/reported/cell_type_l3_recruitedmaconly_specific_tlr2vsnull_infectionmarkers_seurat.csv")
# markers_df <- markers_df %>% mutate(label = ifelse(feature %in% c("Sod2", "Isg15", "Irgm1", "C3ar1", "Ifi47", "Zbp1",
#                                                                   "Irf7", "Tnfaip2", "Cebpb", "Prdx5",
#                                                                   "Trem2", "Cd9", "Lgals3", "Fabp5", "Abcg1"), T, F))
markers_df <- markers_df %>% mutate(label = ifelse(gene %in% c("Sod2", "Ptgs2", "Cebpb", "Clec4e", "Tnfaip2", "Prdx5"), T, F))
r <- ggplot(markers_df, aes(x = avg_log2FC, y = neglogp, color = sig)) +
  geom_point() +
  ggrepel::geom_text_repel(data = markers_df %>% filter(label == T), min.segment.length = 0.05, box.padding = 0.5, nudge_y = 1, nudge_x = -0.5,
                           mapping = aes(x = avg_log2FC, y = neglogp, color = sig, label = gene)) +
  labs(x = "Log2(Fold Change)\nWT     <>     Tlr2-/-", y = "-Log10(adj. P)") +
  GeneralTheme(18) +
  scale_color_manual(values = c("gray80", "black")) +
  RemoveBackgrounds(outline = TRUE) +
  SpaceAxisTitles() +
  NoLegend()
r
SavePlot(r, filename = "degenes_recruitedmacs", root = "", h = 5, w = 5, s = 1, save.data = FALSE)

tsv1 <- read_tsv(file = "data/reported/GO_Biological_Process_2021_table.txt")
tsv1$db <- "GO Biological Process 2021"
tsv2 <- read_tsv(file = "data/reported/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_table.txt")
tsv2$db <- "ENCODE"
tsv3 <- read_tsv(file = "data/reported/ChEA_2016_table.txt")
tsv3$db <- "ChEA 2016"
tsv4 <- read_tsv(file = "data/reported/MSigDB_Hallmark_2020_table.txt")
tsv4$db <- "MSigDB Hallmark"
enr <- bind_rows(tsv1, tsv2, tsv3, tsv4)
enr <- janitor::clean_names(enr)
enr <- enr %>% filter(adjusted_p_value <= 0.05)
enr <- enr %>% filter(db != "ENCODE") %>% group_by(db) %>% top_n(-3, adjusted_p_value)

s <- ggplot(enr %>% mutate(adjusted_p_value = -log10(adjusted_p_value)),
            aes(x = adjusted_p_value, y = fct_reorder(term, adjusted_p_value),
                size = odds_ratio, fill = db)) +
  geom_point(shape = 21, color = "black") +
  labs(y = "Term", x = "-Log10(adj. P)") +
  GeneralTheme(18) +
  scale_size_continuous(name = "Odds Ratio") +
  ggthemes::scale_fill_ptol(name = "Database") +
  RemoveBackgrounds(outline = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  SpaceAxisTitles()
s
SavePlot(s, filename = "tl2degenes_recruitedmac_enrichment", root = "", h = 4, w = 10, s = 1, save.data = FALSE)


# loadings ----------------------------------------------------------------

ls <- read_csv("data/external/mtb_WTvsTLR2KO.csv", sheet = 1)
ls <- janitor::clean_names(ls)
ls <- ls %>% gather(key = "pc", value = "loading", -column1)
ls <- ls %>% group_by(pc) %>% top_n(50, loading)
sigs <- split(ls$column1, ls$pc)

library(UCell)
obj <- AddModuleScore_UCell(obj, features = sigs, seed = 1)

df <- obj[[]] %>%
  select(batch, cell_type_l3, condition, pc_1_UCell, pc_2_UCell) %>%
  group_by(condition, batch, cell_type_l3) %>%
  summarize(pc1 = mean(pc_1_UCell),
         pc2 = mean(pc_2_UCell)) %>%
  gather(key = "pc", value = "score", -cell_type_l3, -batch, -condition)
df

wilcox.test(df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "WT H37Rv" & df$pc == "pc1"],
            df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "Tlr2-/- H37Rv" & df$pc == "pc1"])
wilcox.test(df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "WT H37Rv" & df$pc == "pc2"],
            df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "Tlr2-/- H37Rv" & df$pc == "pc2"])

t <- df %>%
  filter(cell_type_l3 == "Recruited macrophage" & condition %in% c("WT H37Rv", "Tlr2-/- H37Rv")) %>%
  ggplot(mapping = aes(x = condition, y = score, color = condition)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "", y = "Score") +
  ggthemes::scale_color_ptol(name = "Condition", guide = FALSE) +
  GeneralTheme(14) +
  SpaceAxisTitles() +
  facet_wrap(~ pc, scales = "free_y")
plot(t)
SavePlot(t, filename = "pcscores_recruitedmac", root = "", h = 5, w = 6, s = 1, save.data = FALSE)


# infected_genes ----------------------------------------------------------

obj <- readRDS("data/processed/annotated_object.rds")
ls <- read_csv("data/external/mtb_WTvsTLR2KO.csv")
ls <- janitor::clean_names(ls)
ls <- ls %>% mutate(group = ifelse(avg_log2fc < 0, "TLR2KO", "WT")) %>% mutate(avg_log2fc = abs(avg_log2fc))
ls <- ls %>% group_by(group) %>% top_n(50, avg_log2fc)
sigs <- split(ls$x1, ls$group)

library(UCell)
obj <- UCell::AddModuleScore_UCell(obj, features = sigs)

obj[[]]
df <- obj[[]] %>%
  select(batch, cell_type_l3, condition, TLR2KO_UCell, WT_UCell) %>%
  group_by(condition, batch, cell_type_l3) %>%
  summarize(tlr2ko = mean(TLR2KO_UCell),
            wt = mean(WT_UCell)) %>%
  gather(key = "score_type", value = "score", -cell_type_l3, -batch, -condition)
df

wilcox.test(df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "WT H37Rv" & df$score_type == "tlr2ko"],
            df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "Tlr2-/- H37Rv" & df$score_type == "tlr2ko"])
wilcox.test(df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "WT H37Rv" & df$score_type == "wt"],
            df$score[df$cell_type_l3 == "Recruited macrophage" & df$condition == "Tlr2-/- H37Rv" & df$score_type == "wt"])

t <- df %>%
  filter(cell_type_l3 == "Recruited macrophage" & condition %in% c("WT H37Rv", "Tlr2-/- H37Rv")) %>%
  separate(batch, c("genetics", "mtb", "mouse", "lane"), sep = "_", remove = FALSE) %>%
  ggplot(mapping = aes(x = condition, y = score, group = condition, color = lane)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3) +
  labs(x = "", y = "Score") +
  ggthemes::scale_color_ptol(name = "Mouse") +
  GeneralTheme(14) +
  SpaceAxisTitles() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~ score_type, scales = "free_y")
plot(t)
SavePlot(t, filename = "compgenesbwinfectedinvitro_recruitedmac", root = "", h = 4, w = 7, s = 1, save.data = FALSE)


# hallmark scores ---------------------------------------------------------

tnfa <- read_csv(file = "data/external/hallmark_tnf.txt", skip = 2, col_names = FALSE)
ifng <- read_csv(file = "data/external/hallmark_ifng.txt", skip = 2, col_names = FALSE)
ifna <- read_csv(file = "data/external/hallmark_ifna.txt", skip = 2, col_names = FALSE)

library(nichenetr)
tnfa <- nichenetr::convert_human_to_mouse_symbols(tnfa$X1)
ifng <- nichenetr::convert_human_to_mouse_symbols(ifng$X1)
ifna <- nichenetr::convert_human_to_mouse_symbols(ifna$X1)

sigs <- list(tnfa = tnfa, ifng = ifng, ifna = ifna)
obj <- UCell::AddModuleScore_UCell(obj, features = sigs)

obj[[]]
df <- obj[[]] %>%
  select(batch, cell_type_l3, condition, tnfa_UCell, ifng_UCell, ifna_UCell) %>%
  group_by(condition, batch, cell_type_l3) %>%
  summarize(TNFA = mean(tnfa_UCell),
            IFNG = mean(ifng_UCell),
            IFNA = mean(ifna_UCell)) %>%
  gather(key = "score_type", value = "score", -cell_type_l3, -batch, -condition)
df

res <- df %>%
  filter(cell_type_l3 == "Recruited macrophage" & condition %in% c("WT H37Rv", "Tlr2-/- H37Rv")) %>%
  group_by(score_type) %>%
  do(w = wilcox.test(score ~ condition, data = ., paired = FALSE)) %>%
  summarise(score_type, wilcox = w$p.value)
res

u <- df %>%
  filter(cell_type_l3 == "Recruited macrophage" & condition %in% c("WT H37Rv", "Tlr2-/- H37Rv")) %>%
  separate(batch, c("genetics", "mtb", "mouse", "lane"), sep = "_", remove = FALSE) %>%
  ggplot(mapping = aes(x = condition, y = score, group = condition, color = lane)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "", y = "Score") +
  ggthemes::scale_color_ptol(name = "Mouse") +
  GeneralTheme(14) +
  SpaceAxisTitles() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~ score_type, scales = "free_y")
plot(u)

grid <- t + u + patchwork::plot_layout(widths = c(2, 3))
grid
SavePlot(grid, filename = "scores_recruitedmac", root = "", h = 5, w = 10, s = 1, save.data = FALSE)


# compare to pisu et al ---------------------------------------------------

obj <- readRDS("data/processed/annotated_object.rds")
obj_f <- obj[, obj$cell_type_l1 %in% c("Granulocyte", "Myeloid cell")]
pisu <- readRDS("data/external/GSE167232_mtb_integrated.rds")
pisu[["ann_level_5"]] <- Idents(pisu)
DimPlot(pisu, label = TRUE)

DefaultAssay(pisu) <- "RNA"
pisu <- ScaleData(pisu, features = rownames(pisu))

vargenes <- VariableFeatures(obj_f)
pisu_vargenes <- VariableFeatures(pisu)
total_vargenes <- intersect(pisu_vargenes, vargenes)

anchors <- FindTransferAnchors(reference = pisu, query = obj_f,
                               features = total_vargenes,
                               dims = 1:30,
                               n.trees = 20,
                               reduction = "pcaproject")
mapping <- MappingScore(anchors, ndim = 30)
predictions <- TransferData(anchorset = anchors,
                            refdata = pisu$ann_level_5,
                            n.trees = 50,
                            dims = 1:30)
predictions$mapping_score <- mapping
predictions$mapping_score
predictions$cell_barcode <- rownames(predictions)
write_csv(predictions, "data/processed/pisu_predictions.csv")
obj_f <- AddMetaData(obj, metadata = predictions)
obj_f[[]]

g <- DimPlot(object = obj_f, group.by = c("predicted.id"), reduction = "pca_umap",
             label = T, raster = F, cells = obj_f$cell_barcode[obj_f$cell_type_l1 != "Nondescript"]) +
  GeneralTheme(10)
g
SavePlot(g, filename = "predicted_umap", root = "", h = 5, w = 7, s = 1, save.data = FALSE)

obj_f[[]]
df <- obj_f[[]] %>%
  select(batch, cell_type_l1, cell_type_l3, predicted.id, 67:80) %>%
  filter(cell_type_l1 %in%c("Granulocyte", "Myeloid cell")) %>%
  filter(cell_type_l3 != "Nondescript") %>%
  pivot_longer(cols = 5:18, names_to = "pisu_subset", names_pattern = "prediction.score.(.*)", values_to = "score")
pal <- circlize::colorRamp2(
  seq(-4, 4, length.out = 11),
  rev(RColorBrewer::brewer.pal(11, "RdBu")))
df %>%
  group_by(cell_type_l3, pisu_subset) %>%
  summarize(score = mean(score)) %>%
  ungroup() %>%
  heatmap(.row = cell_type_l3, .column = pisu_subset, .value = score, scale = "column",
          palette_value = pal,
          #row_km = 4,
          #column_km = 4,
          column_title = "Pisu subsets projected",
          row_title = "Cell type (L3) scored",
          border_gp = grid::gpar(col = "black", lty = 1),
          rect_gp = grid::gpar(col = "black", lty = 1),
          heatmap_legend_param = list(title = "Scaled\nscore",
                                      border = "black",
                                      title_position = "topleft")) %>%
  save_pdf("plots/pisu_predicted_heatmap.pdf", width = 8, height = 5, units = c("in"))

library(tidyHeatmap)
pal <- circlize::colorRamp2(
  seq(-2, 3, length.out = 11),
  rev(RColorBrewer::brewer.pal(11, "RdBu")))
tbl %>%
  ungroup() %>%
  heatmap(Var1, Var2, prop, .scale = "column",
          palette_value = pal,
          border = TRUE,
          column_title = "Cynomolgus subsets",
          row_title = "Esaulova et al. 2021 Subsets",
          border_gp = grid::gpar(col = "black", lty = 1),
          rect_gp = grid::gpar(col = "black", lty = 1),
          column_names_rot = 45,
          heatmap_legend_param = list(title = "Column-scaled\npredicted\nproportion",
                                      border = "black",
                                      title_position = "topleft")) %>%
  save_pdf("plots/esakha2022_predicted_heatmap.pdf", width = 8, height = 5, units = c("in"))
