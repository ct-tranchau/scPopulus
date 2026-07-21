## ===========================================================================
## 00 -- CCA integration of the poplar snRNA-seq samples + omg cell-type annotation.
##   * reads the six filtered 10x matrices (3 Control, 3 MeJA)
##   * QC floor (min.cells = 3, min.features = 200), 5000 HVGs per sample
##   * CCA anchors -> IntegrateData -> PCA(50) -> UMAP/clusters on dims 1:34, res 0.25
##   * FindAllMarkers, then omg (OrthoMarkerGeneGroups) for cell-type prediction
## Produces the integrated object used by scripts 10-40 (final cell_type labels are
## curated from the omg predictions and added to the object before running 10).
## ===========================================================================
library(Seurat)
library(ggplot2)

## input / output folders: from `Rscript script.R <input_dir> <output_dir>`, else defaults
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output"

## six filtered 10x sample directories (place under input/)
sample_dirs <- c(
  ctr1 = file.path(input_dir, "filtered-ctr", "filtered_ctr1"),
  ctr2 = file.path(input_dir, "filtered-ctr", "filtered_ctr2"),
  ctr3 = file.path(input_dir, "filtered-ctr", "filtered_ctr3"),
  mej1 = file.path(input_dir, "filtered-mej", "filtered_mej1"),
  mej2 = file.path(input_dir, "filtered-mej", "filtered_mej2"),
  mej3 = file.path(input_dir, "filtered-mej", "filtered_mej3"))

## ---- read each sample, normalize, find variable features ----
objs <- lapply(names(sample_dirs), function(nm) {
  obj <- CreateSeuratObject(Read10X(sample_dirs[[nm]]), project = nm,
                            min.cells = 3, min.features = 200)
  obj$sample     <- nm
  obj$condition  <- ifelse(grepl("^ctr", nm), "ctr", "mej")
  obj$orig.ident <- nm
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # nuclear ref -> ~0
  obj <- NormalizeData(obj)
  FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
})
names(objs) <- names(sample_dirs)

## ---- CCA integration ----
features <- SelectIntegrationFeatures(objs, nfeatures = 5000)
anchors  <- FindIntegrationAnchors(objs, anchor.features = features)
dat.integrated <- IntegrateData(anchorset = anchors)

## ---- dimensionality reduction + clustering (dims 1:34, res 0.25) ----
DefaultAssay(dat.integrated) <- "integrated"
dat.integrated <- ScaleData(dat.integrated, verbose = FALSE)
dat.integrated <- RunPCA(dat.integrated, npcs = 50, verbose = FALSE)
dat.integrated <- RunUMAP(dat.integrated, dims = 1:34)
dat.integrated <- FindNeighbors(dat.integrated, dims = 1:34)
dat.integrated <- FindClusters(dat.integrated, resolution = 0.25)

## ---- save object + UMAPs ----
saveRDS(dat.integrated, file.path(output_dir, "scRNA_integrated_all_samples_PC34_res0.25_annotated.rds"))
ggsave(file.path(output_dir, "UMAP_clusters.pdf"),     DimPlot(dat.integrated, reduction = "umap", label = TRUE),                 width = 8, height = 6)
ggsave(file.path(output_dir, "UMAP_by_sample.pdf"),    DimPlot(dat.integrated, reduction = "umap", group.by = "sample"),          width = 8, height = 6)
ggsave(file.path(output_dir, "UMAP_by_condition.pdf"), DimPlot(dat.integrated, reduction = "umap", group.by = "condition"),       width = 8, height = 6)

## ---- cluster marker genes (input for annotation) ----
DefaultAssay(dat.integrated) <- "RNA"
dat.integrated[["RNA"]] <- JoinLayers(dat.integrated[["RNA"]])
dat.integrated <- NormalizeData(dat.integrated, verbose = FALSE)
Idents(dat.integrated) <- "seurat_clusters"
markers <- FindAllMarkers(dat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(output_dir, "marker_genes_all_clusters.csv"), row.names = FALSE)

## ---- cell-type annotation with omg (OrthoMarkerGeneGroups) ----
if (!requireNamespace("omg", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("LiLabAtVT/OrthoMarkerGeneGroups", subdir = "omg", upgrade = "never")
}
library(omg)

## omg reference stores genes as "Potri_001G..." (underscore); convert the dot form
markers$gene <- sub("\\.", "_", markers$gene)
res <- omg(markers, fdr = 0.05, top_n = 200,
           outdir = file.path(output_dir, "omg_output"), write_files = TRUE, verbose = TRUE)
print(res$predictions)

sessionInfo()
