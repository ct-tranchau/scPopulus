## ===========================================================================
## 40 -- All cell wall / secondary-metabolism figures in one file:
##   (A) dot plot of cell wall genes across cell type x condition
##   (B) module score (AddModuleScore) bar chart, mean +/- SE, Control vs MeJA
##   (C) clustered heatmap (row-scaled) across cell types
##   (D) FeaturePlot grid (MeJA vs Control) for EXPR/PAL1/LAC14/CAD9
## ===========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(pheatmap)
library(RColorBrewer)
library(scattermore)

## input / output folders: from `Rscript script.R <input_dir> <output_dir>`, else defaults
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output"

## cell wall / secondary-metabolism genes; id -> label ("" = show ID only)
cw_lab <- c("Potri.006G069600"="","Potri.003G083200"="EXPR/EXLB1","Potri.003G072800"="PME3",
  "Potri.002G018300"="CAD9","Potri.008G106400"="DIOX1/DOX1/PADOX-1","Potri.001G105200"="GPX6/PHGPX",
  "Potri.016G132800"="","Potri.002G065300"="","Potri.007G106400"="GONST1","Potri.019G088600"="LAC14",
  "Potri.006G126800"="PAL1","Potri.019G101900"="EXPB3","Potri.018G094200"="AT4CL2","Potri.015G003600"="RCI3A",
  "Potri.006G094100"="LAC7","Potri.019G130700"="C4H/REF3","Potri.019G088500"="LAC14","Potri.004G144600"="",
  "Potri.013G019800"="CESA9","Potri.016G091100"="PAL1","Potri.005G135300"="","Potri.019G049500"="4CL3",
  "Potri.003G214500"="","Potri.001G304800"="CCoAOMT1","Potri.016G132732"="")
cw <- names(cw_lab)

obj <- readRDS(file.path(input_dir, "scRNA_integrated_all_samples_PC34_res0.25_annotated.rds"))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
cw  <- cw[cw %in% rownames(obj)]

md  <- obj@meta.data
ct  <- sub("^[0-9]+:\\s*", "", as.character(md$cell_type))
cnd <- ifelse(grepl("mej", md$condition, ignore.case = TRUE), "mej", "ctr")
ord <- data.frame(cl = as.numeric(as.character(md$seurat_clusters)), lin = ct) %>%
  group_by(lin) %>% summarise(cl = min(cl), .groups = "drop") %>% arrange(cl)
lins <- ord$lin

## ---- (A) dot plot: genes on y, cell type x condition on x ----
obj$ct_cond <- factor(paste0(ct, "_", cnd), levels = c(paste0(lins, "_ctr"), paste0(lins, "_mej")))
dd <- DotPlot(obj, features = cw, group.by = "ct_cond", assay = "RNA")$data
dd$features.plot <- factor(dd$features.plot, levels = rev(cw))
pA <- ggplot(dd, aes(id, features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "lightgrey", high = "#67001f", name = "Average Expression") +
  scale_size(range = c(0, 8), name = "Percent Expressed") +
  labs(title = "Secondary metabolism and cell wall gene expression across root cell types", x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14), legend.text = element_text(size = 13),
        panel.grid.major = element_line(color = "grey92"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "Cellwall_secondary_genes_dotplot_by_celltype_condition.pdf"), pA, width = 13, height = 8)
ggsave(file.path(output_dir, "Cellwall_secondary_genes_dotplot_by_celltype_condition.png"), pA, width = 13, height = 8, dpi = 300)

## ---- (B) module score (AddModuleScore) bar chart, mean +/- SE ----
set.seed(42)
obj <- AddModuleScore(obj, features = list(cw), name = "CellWall_Lignin_Module", assay = "RNA")
dfB <- data.frame(lineage = factor(ct, levels = lins),
                  Condition = factor(ifelse(cnd == "mej", "MeJA", "Control"), levels = c("Control", "MeJA")),
                  score = obj$CellWall_Lignin_Module1)
summ <- dfB %>% group_by(lineage, Condition) %>%
  summarise(mean = mean(score), se = sd(score)/sqrt(n()), .groups = "drop")
pB <- ggplot(summ, aes(lineage, mean, fill = Condition)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.8), width = 0.25, linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(Control = "#0072B2", MeJA = "#C0392B")) +
  labs(title = "MeJA enhances cell wall and secondary metabolism module activity",
       x = "Cell type", y = "Mean cell wall/lignin module score") +
  theme_bw(base_size = 17) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17), axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 21, face = "bold"), plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
        legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "Cellwall_module_score_barplot.pdf"), pB, width = 13.5, height = 7.6)
ggsave(file.path(output_dir, "Cellwall_module_score_barplot.png"), pB, width = 13.5, height = 7.6, dpi = 300)

## ---- (C) clustered heatmap (per-gene z-score of mean expression per cell type) ----
dat <- GetAssayData(obj, assay = "RNA", layer = "data")[cw, , drop = FALSE]
mat <- sapply(lins, function(l) Matrix::rowMeans(dat[, ct == l, drop = FALSE]))
rownames(mat) <- cw; colnames(mat) <- lins
z <- t(scale(t(mat))); z[is.na(z)] <- 0; z <- pmax(pmin(z, 3), -3)
rlab <- ifelse(cw_lab[cw] == "", cw, paste0(cw_lab[cw], " / ", cw))
ph <- pheatmap(z, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
               breaks = seq(-3, 3, length.out = 102), cluster_rows = TRUE, cluster_cols = FALSE,
               labels_row = rlab, fontsize_row = 10, fontsize_col = 11, angle_col = 45,
               border_color = "grey80", treeheight_row = 45, treeheight_col = 0,
               main = "Cell wall / secondary metabolism gene expression across root cell types", silent = TRUE)
ggsave(file.path(output_dir, "Heatmap_cellwall_secondary_by_celltype.pdf"), ph$gtable, width = 8.8, height = 7)
ggsave(file.path(output_dir, "Heatmap_cellwall_secondary_by_celltype.png"), ph$gtable, width = 8.8, height = 7, dpi = 300)

## ---- (D) FeaturePlot grid (MeJA top / Control bottom), red ----
genesF <- c("Potri.003G083200"="EXPR / EXLB1\nPotri.003G083200","Potri.006G126800"="PAL1\nPotri.006G126800",
            "Potri.019G088600"="LAC14\nPotri.019G088600","Potri.002G018300"="CAD9\nPotri.002G018300")
idsF <- names(genesF); idsF <- idsF[idsF %in% rownames(obj)]
emb   <- Embeddings(obj, "umap")
condF <- ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control")
exprF <- FetchData(obj, vars = idsF, layer = "data")
long <- do.call(rbind, lapply(idsF, function(g) data.frame(
  UMAP_1 = emb[,1], UMAP_2 = emb[,2], condition = condF, gene = genesF[[g]], expr = exprF[[g]], row.names = NULL)))
long$gene <- factor(long$gene, levels = unname(genesF[idsF]))
long$condition <- factor(long$condition, levels = c("MeJA", "Control"))
cap <- ceiling(as.numeric(quantile(long$expr[long$expr > 0], 0.99)) * 2) / 2
long$expr_c <- pmin(long$expr, cap); long <- long[order(long$expr_c), ]
pD <- ggplot(long, aes(UMAP_1, UMAP_2, color = expr_c)) +
  geom_scattermore(pointsize = 10, pixels = c(4096, 4096)) + coord_fixed(1) +
  facet_grid(condition ~ gene, switch = "y") +
  scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, cap),
    breaks = seq(0, cap, by = if (cap <= 3) 0.5 else 1), name = "Expression\n(log-norm)",
    guide = guide_colorbar(barwidth = 12, barheight = 0.7, title.position = "top",
                           title.hjust = 0.5, ticks.colour = "grey30", frame.colour = "grey30")) +
  theme_bw(base_size = 19) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(2, "pt"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.x = element_text(size = 15, lineheight = 0.9, face = "bold"),
        strip.text.y.left = element_text(size = 20, face = "bold", angle = 90),
        legend.position = "bottom", legend.title = element_text(size = 17), legend.text = element_text(size = 16))
ggsave(file.path(output_dir, "FeaturePlot_paper_grid_EXPR_PAL1_LAC14_CAD9.pdf"), pD, width = 13, height = 7.5)
ggsave(file.path(output_dir, "FeaturePlot_paper_grid_EXPR_PAL1_LAC14_CAD9.png"), pD, width = 13, height = 7.5, dpi = 400, bg = "white")

sessionInfo()
