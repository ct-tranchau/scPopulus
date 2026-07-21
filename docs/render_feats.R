## Re-render only the Control / MeJA feature PNGs with the color legend at the BOTTOM
## (more horizontal space for the UMAP scatter). Reads ../shiny_app/webapp_data.rds.
library(Seurat)
library(ggplot2)
library(scattermore)

obj   <- readRDS("../shiny_app/webapp_data.rds")
genes <- read.csv("../shiny_app/genes.csv", stringsAsFactors = FALSE)
genes <- genes[genes$gene_id %in% rownames(obj), ]

obj$Condition <- factor(ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control"),
                        levels = c("Control", "MeJA"))
emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$Condition <- obj$Condition

feat <- function(gene, cond, cap) {
  d <- emb[emb$Condition == cond, ]
  d$expr <- pmin(FetchData(obj, vars = gene, layer = "data")[rownames(d), 1], cap)
  d <- d[order(d$expr), ]
  ggplot(d, aes(UMAP_1, UMAP_2, color = expr)) +
    geom_scattermore(pointsize = 2.6) +
    scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, cap), name = "log-norm") +
    coord_fixed() + ggtitle(cond) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "bottom",
          legend.title = element_text(size = 13), legend.text = element_text(size = 12),
          legend.key.width = unit(1.3, "cm"), legend.key.height = unit(0.35, "cm"))
}

ug <- unique(genes$gene_id)
for (i in seq_along(ug)) {
  g <- ug[i]; e <- FetchData(obj, vars = g, layer = "data")[, 1]
  cap <- max(0.5, ceiling(as.numeric(quantile(e[e > 0], 0.99, na.rm = TRUE)) * 2) / 2)
  ggsave(file.path("images", paste0(g, "_Control.png")), feat(g, "Control", cap), width = 5, height = 5, dpi = 90, bg = "white")
  ggsave(file.path("images", paste0(g, "_MeJA.png")),    feat(g, "MeJA",    cap), width = 5, height = 5, dpi = 90, bg = "white")
  if (i %% 25 == 0) cat(i, "/", length(ug), "\n")
}
cat("re-rendered", length(ug), "feature PNG pairs (bottom legend)\n")
