## ===========================================================================
## 20 -- Downstream of the MeJA-vs-Control DEGs (from 10):
##   * up-regulated DEG bar chart (Cortex clusters merged)
##   * up/down diverging bar chart per cluster
##   * UpSet plot of up-DEG intersections across lineages
##   * cell-type-resolved GO Biological Process enrichment dot plot (Control vs MeJA)
## Bar charts/UpSet read only the DEG table; GO uses biomaRt + org.At.tair.db.
## ===========================================================================
library(ggplot2)
library(dplyr)
library(UpSetR)
library(clusterProfiler)
library(org.At.tair.db)
library(biomaRt)
options(timeout = 600)

## input / output folders: from `Rscript script.R <input_dir> <output_dir>`, else defaults
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output"

deg <- read.csv(file.path(output_dir, "DEG_MeJA_vs_Control_per_cluster.csv"), stringsAsFactors = FALSE)
lab <- read.csv(file.path(output_dir, "up_DEG_MeJA_counts_per_cluster.csv"), stringsAsFactors = FALSE)
deg$cluster <- as.character(deg$cluster); lab$cluster <- as.character(lab$cluster)
deg$lineage <- sub("^[0-9]+:\\s*", "", lab$label[match(deg$cluster, lab$cluster)])
ord <- lab %>% mutate(cl = as.numeric(cluster), lin = sub("^[0-9]+:\\s*", "", label)) %>%
  group_by(lin) %>% summarise(cl = min(cl), .groups = "drop") %>% arrange(cl)
lin_levels <- ord$lin

## ---- (A) up-regulated DEGs per lineage (Cortex merged) ----
up <- deg %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  group_by(lineage) %>% summarise(n_up = n_distinct(gene), .groups = "drop")
up$lineage <- factor(up$lineage, levels = lin_levels)
pA <- ggplot(up, aes(lineage, n_up)) +
  geom_col(fill = "#C0392B", width = 0.75) +
  geom_text(aes(label = n_up), vjust = -0.3, size = 6.5) +
  labs(x = "Cell type", y = "Number of up-regulated DEGs (MeJA)") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 21, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggsave(file.path(output_dir, "barchart_up_DEG_MeJA_per_cluster.pdf"), pA, width = 12, height = 7.5)
ggsave(file.path(output_dir, "barchart_up_DEG_MeJA_per_cluster.png"), pA, width = 12, height = 7.5, dpi = 300)

## ---- (B) up/down diverging bar chart ----
cnt <- deg %>% mutate(sig = p_val_adj < 0.05) %>% group_by(cluster) %>%
  summarise(Up = sum(sig & avg_log2FC > 0.25), Down = sum(sig & avg_log2FC < -0.25), .groups = "drop")
plotdf <- bind_rows(data.frame(cluster = cnt$cluster, dir = "Up",   n =  cnt$Up),
                    data.frame(cluster = cnt$cluster, dir = "Down", n = -cnt$Down))
plotdf$label <- factor(sub("^[0-9]+:\\s*", "", lab$label[match(plotdf$cluster, lab$cluster)]), levels = lin_levels)
pB <- ggplot(plotdf, aes(label, n, fill = dir)) +
  geom_col(width = 0.75) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = abs(n), vjust = ifelse(n >= 0, -0.3, 1.3)), size = 4.2) +
  scale_fill_manual(values = c(Up = "#C0392B", Down = "#0072B2"), name = NULL,
                    labels = c("Up-regulated (MeJA)", "Down-regulated (MeJA)")) +
  labs(x = "Cluster", y = "Number of DEGs  (down  |  up)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.title = element_text(size = 16, face = "bold"), legend.position = "top") +
  scale_y_continuous(labels = abs)
ggsave(file.path(output_dir, "barchart_up_down_DEG_MeJA_per_cluster.pdf"), pB, width = 12, height = 8)
ggsave(file.path(output_dir, "barchart_up_down_DEG_MeJA_per_cluster.png"), pB, width = 12, height = 8, dpi = 300)

## ---- (C) UpSet of up-DEG sets across lineages ----
upg  <- deg %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)
sets <- lapply(split(upg$gene, upg$lineage), unique)
sets <- sets[order(sapply(sets, length), decreasing = TRUE)]
draw <- function() print(upset(fromList(sets), sets = rev(names(sets)), keep.order = TRUE,
  nintersects = 40, order.by = "freq", mb.ratio = c(0.6, 0.4),
  main.bar.color = "#C0392B", sets.bar.color = "#0072B2", matrix.color = "#2C3E50",
  point.size = 2.6, line.size = 0.9,
  mainbar.y.label = "Intersection size (shared up-DEGs)", sets.x.label = "Upregulated genes per lineage",
  text.scale = c(2.2, 2.0, 2.0, 1.8, 2.0, 1.6)))
pdf(file.path(output_dir, "UpSet_MeJA_up_DEGs_lineages.pdf"), width = 14, height = 8); draw(); dev.off()
png(file.path(output_dir, "UpSet_MeJA_up_DEGs_lineages.png"), width = 14, height = 8, units = "in", res = 300); draw(); dev.off()

## ---- (D) GO BP enrichment dot plot (Control-up vs MeJA-up, Arabidopsis orthologs) ----
## poplar -> Arabidopsis (TAIR) ortholog map, fetched once from Ensembl Plants then cached
if (file.exists(file.path(output_dir, "poplar_ath_ortho.rds"))) {
  o <- readRDS(file.path(output_dir, "poplar_ath_ortho.rds"))
} else {
  m <- useEnsemblGenomes(biomart = "plants_mart", dataset = "ptrichocarpa_eg_gene")
  o <- getBM(attributes = c("ensembl_gene_id", "athaliana_eg_homolog_ensembl_gene"), mart = m)
  o$gene <- sub("\\.v[0-9]+\\.[0-9]+$", "", o$ensembl_gene_id)   # strip .v4.1 suffix
  o$ath  <- toupper(o$athaliana_eg_homolog_ensembl_gene)
  o <- unique(o[!is.na(o$ath) & o$ath != "", c("gene", "ath")])
  saveRDS(o, file.path(output_dir, "poplar_ath_ortho.rds"))
}
p2a <- split(o$ath, o$gene); to_ath <- function(g) unique(unlist(p2a[g], use.names = FALSE))
universe_ath <- unique(o$ath)
up_sets   <- lapply(split(deg$gene[deg$p_val_adj<0.05 & deg$avg_log2FC> 0.25],
                          deg$lineage[deg$p_val_adj<0.05 & deg$avg_log2FC> 0.25]), unique)
down_sets <- lapply(split(deg$gene[deg$p_val_adj<0.05 & deg$avg_log2FC< -0.25],
                          deg$lineage[deg$p_val_adj<0.05 & deg$avg_log2FC< -0.25]), unique)
terms <- c(
  "GO:0042254"="ribosome biogenesis","GO:0002181"="cytoplasmic translation",
  "GO:0042546"="cell wall biogenesis","GO:0009664"="plant-type cell wall organization",
  "GO:0005976"="polysaccharide metabolic process","GO:0009611"="response to wounding",
  "GO:0006979"="response to oxidative stress","GO:0009620"="response to fungus",
  "GO:0009414"="response to water deprivation","GO:0009753"="response to jasmonic acid",
  "GO:0019748"="secondary metabolic process","GO:0044550"="secondary metabolite biosynthetic process",
  "GO:0009698"="phenylpropanoid metabolic process","GO:0009699"="phenylpropanoid biosynthetic process",
  "GO:0009808"="lignin metabolic process","GO:0009809"="lignin biosynthetic process")
enr <- function(genes) {
  ath <- to_ath(genes); if (length(ath) < 5) return(NULL)
  er <- enrichGO(ath, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP",
                 universe = universe_ath, pvalueCutoff = 1, qvalueCutoff = 1,
                 minGSSize = 3, maxGSSize = 5000)
  r <- as.data.frame(er); r[r$ID %in% names(terms), c("ID", "Count", "p.adjust")]
}
collect <- function(sets, panel) do.call(rbind, lapply(names(sets), function(ln) {
  r <- enr(sets[[ln]]); if (is.null(r) || !nrow(r)) return(NULL)
  data.frame(panel = panel, lineage = ln, term = terms[r$ID], Count = r$Count, padj = r$p.adjust) }))
df <- rbind(collect(down_sets, "Control"), collect(up_sets, "MeJA"))
df <- df[df$padj < 0.05, ]; df$neglog10 <- -log10(df$padj)
df$lineage <- factor(df$lineage, levels = lin_levels)
df$panel   <- factor(df$panel, levels = c("Control", "MeJA"))
df$term    <- factor(df$term, levels = rev(unname(terms)))
pD <- ggplot(df, aes(lineage, term, size = Count, color = neglog10)) +
  geom_point() + facet_grid(. ~ panel) +
  scale_color_viridis_c(option = "viridis", name = expression(-log[10]*"(adjusted p-value)")) +
  scale_size_continuous(range = c(2, 9), name = "Gene count") +
  labs(x = "Cell type", y = "GO biological process") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 19, face = "bold"), strip.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "GO_BP_celltype_dotplot_MeJA_vs_Control.pdf"), pD, width = 23, height = 8.5)
ggsave(file.path(output_dir, "GO_BP_celltype_dotplot_MeJA_vs_Control.png"), pD, width = 23, height = 8.5, dpi = 300)

sessionInfo()
