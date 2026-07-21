## Render one combined PDF per gene (Control | MeJA on top, violin below) for download.
## Display stays PNG (browsers can't show PDF in <img>); these are the downloadables.
library(Seurat)
library(ggplot2)
library(scattermore)
library(patchwork)

obj   <- readRDS("../shiny_app/webapp_data.rds")
genes <- read.csv("../shiny_app/genes.csv", stringsAsFactors = FALSE)
genes <- genes[genes$gene_id %in% rownames(obj), ]
dir.create("pdf", showWarnings = FALSE)

obj$Condition <- factor(ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control"),
                        levels = c("Control", "MeJA"))
emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$Condition <- obj$Condition
ct_levels <- unique(obj$cell_type)[order(as.numeric(sub(":.*", "", unique(obj$cell_type))))]
obj$cell_type <- factor(obj$cell_type, levels = ct_levels)

feat <- function(gene, cond, cap) {
  d <- emb[emb$Condition == cond, ]
  d$expr <- pmin(FetchData(obj, vars = gene, layer = "data")[rownames(d), 1], cap)
  d <- d[order(d$expr), ]
  ggplot(d, aes(UMAP_1, UMAP_2, color = expr)) +
    geom_scattermore(pointsize = 3) +
    scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, cap), name = "log-norm") +
    coord_fixed() + ggtitle(cond) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18), legend.position = "right")
}
vln <- function(gene) {
  d <- data.frame(expr = FetchData(obj, vars = gene, layer = "data")[, 1],
                  cell_type = obj$cell_type, Condition = obj$Condition)
  ggplot(d, aes(cell_type, expr, fill = Condition)) +
    geom_violin(scale = "width", linewidth = 0.2, position = position_dodge(0.85)) +
    scale_fill_manual(values = c(Control = "#0072B2", MeJA = "#C0392B")) +
    labs(x = NULL, y = "log-norm expression") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.title = element_blank())
}

ug <- unique(genes$gene_id)
for (i in seq_along(ug)) {
  g <- ug[i]; e <- FetchData(obj, vars = g, layer = "data")[, 1]
  cap <- max(0.5, ceiling(as.numeric(quantile(e[e > 0], 0.99, na.rm = TRUE)) * 2) / 2)
  p <- (feat(g, "Control", cap) | feat(g, "MeJA", cap)) / vln(g) +
    plot_annotation(title = g, theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))) +
    plot_layout(heights = c(1, 0.8))
  ggsave(file.path("pdf", paste0(g, ".pdf")), p, width = 12, height = 11)
  if (i %% 25 == 0) cat(i, "/", length(ug), "\n")
}
cat("rendered", length(ug), "PDFs\n")
