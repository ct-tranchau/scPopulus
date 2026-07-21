## ===========================================================================
## 30 -- JA-signaling / defense figures:
##   * dot plot of JA/defense genes across cell type x condition
##   * JA/defense module score (AddModuleScore) boxplot by cell type x condition,
##     with per-cell-type Wilcoxon Control-vs-MeJA significance stars
## ===========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)

## input / output folders: from `Rscript script.R <input_dir> <output_dir>`, else defaults
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output"

## JA-signaling / defense genes (dot-plot y-axis order)
ja <- c("Potri.009G142300","Potri.003G083200","Potri.002G039100","Potri.005G032400",
        "Potri.001G288600","Potri.003G067600","Potri.T011200","Potri.013G153400",
        "Potri.001G154200","Potri.004G007500","Potri.T131500","Potri.018G008500",
        "Potri.013G102800")

obj <- readRDS(file.path(input_dir, "scRNA_integrated_all_samples_PC34_res0.25_annotated.rds"))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
ja <- ja[ja %in% rownames(obj)]

md  <- obj@meta.data
ct  <- sub("^[0-9]+:\\s*", "", as.character(md$cell_type))
cnd <- ifelse(grepl("mej", md$condition, ignore.case = TRUE), "mej", "ctr")
ord <- data.frame(cl = as.numeric(as.character(md$seurat_clusters)), lin = ct) %>%
  group_by(lin) %>% summarise(cl = min(cl), .groups = "drop") %>% arrange(cl)

## ---- (A) JA/defense dot plot: genes on y, cell type x condition on x ----
obj$ct_cond <- factor(paste0(ct, "_", cnd), levels = c(paste0(ord$lin, "_ctr"), paste0(ord$lin, "_mej")))
dd <- DotPlot(obj, features = ja, group.by = "ct_cond", assay = "RNA")$data
dd$features.plot <- factor(dd$features.plot, levels = rev(ja))
pA <- ggplot(dd, aes(id, features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "lightgrey", high = "red", name = "Average Expression") +
  scale_size(range = c(0, 8), name = "Percent Expressed") +
  labs(title = "JA signaling and defense gene expression across root cell types", x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14), legend.text = element_text(size = 13),
        panel.grid.major = element_line(color = "grey92"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "JA_defense_genes_dotplot_by_celltype_condition.pdf"), pA, width = 18, height = 6.5)
ggsave(file.path(output_dir, "JA_defense_genes_dotplot_by_celltype_condition.png"), pA, width = 18, height = 6.5, dpi = 300)

## ---- (B) module score (AddModuleScore) boxplot + Wilcoxon stats ----
set.seed(42)
obj <- AddModuleScore(obj, features = list(ja), name = "JAdefense", assay = "RNA")
df <- data.frame(lineage = factor(ct, levels = ord$lin),
                 Condition = factor(ifelse(cnd == "mej", "MeJA", "Control"), levels = c("Control", "MeJA")),
                 score = obj$JAdefense1)

stats <- do.call(rbind, lapply(levels(df$lineage), function(ln) {
  s <- df[df$lineage == ln, ]; c1 <- s$score[s$Condition == "Control"]; c2 <- s$score[s$Condition == "MeJA"]
  if (length(c1) < 3 || length(c2) < 3) return(NULL)
  sp <- sqrt(((length(c1)-1)*var(c1) + (length(c2)-1)*var(c2)) / (length(c1)+length(c2)-2))
  data.frame(lineage = ln, wilcox_p = wilcox.test(c2, c1)$p.value, cohens_d = (mean(c2)-mean(c1))/sp) }))
stats$wilcox_padj <- p.adjust(stats$wilcox_p, "BH")
stats$sig <- ifelse(stats$wilcox_padj < 1e-4, "****", ifelse(stats$wilcox_padj < 1e-3, "***",
             ifelse(stats$wilcox_padj < 1e-2, "**", ifelse(stats$wilcox_padj < 0.05, "*", "ns"))))

ylim <- quantile(df$score, c(0.005, 0.995))
ann  <- data.frame(lineage = factor(levels(df$lineage), levels = levels(df$lineage)),
                   y = ylim[2]*0.96, sig = stats$sig[match(levels(df$lineage), stats$lineage)])
pB <- ggplot(df, aes(lineage, score, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7, linewidth = 0.3) +
  geom_text(data = ann, aes(lineage, y, label = sig), inherit.aes = FALSE, size = 6, fontface = "bold") +
  scale_fill_manual(values = c(Control = "#0072B2", MeJA = "#C0392B")) +
  coord_cartesian(ylim = ylim) + geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(title = "JA/defense module scores across cell types and conditions",
       subtitle = "Wilcoxon rank-sum, Control vs MeJA per cell type (BH-adjusted): **** p<1e-4",
       x = "Cell type", y = "JA/defense module score") +
  theme_classic(base_size = 17) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17), axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 21, face = "bold"), plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5),
        legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 17))
ggsave(file.path(output_dir, "JA_defense_module_score_boxplot_stats.pdf"), pB, width = 13.5, height = 7.8)
ggsave(file.path(output_dir, "JA_defense_module_score_boxplot_stats.png"), pB, width = 13.5, height = 7.8, dpi = 300)

sessionInfo()
