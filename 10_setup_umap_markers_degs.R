## ===========================================================================
## 10 -- Setup, UMAP, cluster marker genes, and MeJA-vs-Control DEGs.
##   * loads libraries and the integrated + annotated Seurat object
##   * makes the UMAP (cell types + condition overlay)
##   * finds marker genes for each cluster              -> output/marker_genes_all_clusters.csv
##   * finds DEGs (MeJA vs Control) within each cluster  -> output/DEG_MeJA_vs_Control_per_cluster.csv
## Run from the repo root; scripts 20-40 read the tables it writes.
## Replace with your input and output directory.
## ===========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)

## input / output folders: from `Rscript script.R <input_dir> <output_dir>`, else defaults
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output"

## load & prepare: expression figures use the RNA assay (log-normalized)
obj <- readRDS(file.path(input_dir, "scRNA_integrated_all_samples_PC34_res0.25_annotated.rds"))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])        # Seurat v5: collapse per-sample layers
obj <- NormalizeData(obj, verbose = FALSE)
obj$Condition <- factor(ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control"),
                        levels = c("Control", "MeJA"))

## ---- UMAP: annotated cell types + Control/MeJA overlay ----
p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
              label = TRUE, repel = TRUE, raster = FALSE, pt.size = 0.2) + ggtitle(NULL)
ggsave(file.path(output_dir, "UMAP_annotated.pdf"), p1, width = 10, height = 8)

p2 <- DimPlot(obj, reduction = "umap", group.by = "Condition", shuffle = TRUE, pt.size = 0.35,
              raster = FALSE, cols = c(Control = "#0072B2", MeJA = "#C0392B")) + ggtitle(NULL) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  theme(axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16))
p2$layers[[1]]$aes_params$alpha <- 0.5
ggsave(file.path(output_dir, "UMAP_condition_overlay_MeJA_vs_Control.pdf"), p2, width = 9, height = 7.5)
ggsave(file.path(output_dir, "UMAP_condition_overlay_MeJA_vs_Control.png"), p2, width = 9, height = 7.5, dpi = 300)

## ---- marker genes for each cluster (FindAllMarkers, positive) ----
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(output_dir, "marker_genes_all_clusters.csv"), row.names = FALSE)

## ---- DEGs: MeJA vs Control within each cluster ----
md    <- obj@meta.data
map   <- md %>% distinct(seurat_clusters, cell_type) %>% mutate(cluster = as.character(seurat_clusters))
clord <- as.character(sort(unique(as.numeric(map$cluster))))
deg_list <- list()
for (cl in clord) {
  sub <- subset(obj, cells = rownames(md)[as.character(md$seurat_clusters) == cl])
  Idents(sub) <- sub$Condition
  if (sum(sub$Condition == "MeJA") >= 3 && sum(sub$Condition == "Control") >= 3) {
    dm <- FindMarkers(sub, ident.1 = "MeJA", ident.2 = "Control",
                      logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
    dm$gene <- rownames(dm); dm$cluster <- cl; deg_list[[cl]] <- dm
  }
  rm(sub); gc(verbose = FALSE)
}
all_deg <- do.call(rbind, deg_list)
write.csv(all_deg, file.path(output_dir, "DEG_MeJA_vs_Control_per_cluster.csv"), row.names = FALSE)

## per-cluster up-DEG counts + labels (used by script 20)
up <- all_deg %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>% count(cluster, name = "n_up")
counts <- data.frame(cluster = clord) %>% left_join(up, by = "cluster") %>%
  mutate(n_up = ifelse(is.na(n_up), 0L, n_up), label = map$cell_type[match(cluster, map$cluster)])
write.csv(counts, file.path(output_dir, "up_DEG_MeJA_counts_per_cluster.csv"), row.names = FALSE)

sessionInfo()
