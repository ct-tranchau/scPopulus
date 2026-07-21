## Render everything for the static GitHub Pages site (docs/):
##   - interactive plotly UMAP widget (umap_widget.html, hover keeps working)
##   - per-gene Control / MeJA / violin PNGs (images/)
##   - genes.js (gene list for the dropdown)
## Reads ../shiny_app/webapp_data.rds and ../shiny_app/genes.csv.
library(Seurat)
library(ggplot2)
library(scattermore)
library(plotly)
library(htmlwidgets)

obj   <- readRDS("../shiny_app/webapp_data.rds")
genes <- read.csv("../shiny_app/genes.csv", stringsAsFactors = FALSE)
genes <- genes[genes$gene_id %in% rownames(obj), ]
dir.create("images", showWarnings = FALSE)

obj$Condition <- factor(ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control"),
                        levels = c("Control", "MeJA"))
emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$Condition <- obj$Condition
ct_levels <- unique(obj$cell_type)[order(as.numeric(sub(":.*", "", unique(obj$cell_type))))]
obj$cell_type <- factor(obj$cell_type, levels = ct_levels)

## ---- interactive UMAP (plotly widget) ----
umap_df <- data.frame(UMAP_1 = emb$UMAP_1, UMAP_2 = emb$UMAP_2, cell_type = obj$cell_type)
cent <- do.call(rbind, lapply(ct_levels, function(l) {
  d <- umap_df[umap_df$cell_type == l, ]
  data.frame(cell_type = l, UMAP_1 = median(d$UMAP_1), UMAP_2 = median(d$UMAP_2))
}))
umap_plotly <- plot_ly(umap_df, x = ~UMAP_1, y = ~UMAP_2, color = ~cell_type,
                       colors = scales::hue_pal()(length(ct_levels)),
                       type = "scattergl", mode = "markers", marker = list(size = 4),
                       hoverinfo = "text", text = ~cell_type, showlegend = FALSE) %>%
  add_annotations(data = cent, x = ~UMAP_1, y = ~UMAP_2, text = ~cell_type, showarrow = FALSE,
                  font = list(size = 12, color = "black"), bgcolor = "rgba(255,255,255,0.7)") %>%
  layout(paper_bgcolor = "white", plot_bgcolor = "white",
         xaxis = list(title = "UMAP 1", zeroline = FALSE),
         yaxis = list(title = "UMAP 2", zeroline = FALSE, scaleanchor = "x", scaleratio = 1))
saveWidget(umap_plotly, "umap_widget.html", selfcontained = FALSE, libdir = "umap_lib",
           title = "Cell types")
cat("saved umap_widget.html\n")

## ---- per-gene Control / MeJA / violin images ----
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
vln <- function(gene) {
  d <- data.frame(expr = FetchData(obj, vars = gene, layer = "data")[, 1],
                  cell_type = obj$cell_type, Condition = obj$Condition)
  ggplot(d, aes(cell_type, expr, fill = Condition)) +
    geom_violin(scale = "width", linewidth = 0.2, position = position_dodge(0.85)) +
    scale_fill_manual(values = c(Control = "#0072B2", MeJA = "#C0392B")) +
    labs(x = NULL, y = "log-norm expression") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.text.y = element_text(size = 12),
          legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 14))
}

ug <- unique(genes$gene_id)
for (i in seq_along(ug)) {
  g <- ug[i]; e <- FetchData(obj, vars = g, layer = "data")[, 1]
  cap <- max(0.5, ceiling(as.numeric(quantile(e[e > 0], 0.99, na.rm = TRUE)) * 2) / 2)
  ggsave(file.path("images", paste0(g, "_Control.png")), feat(g, "Control", cap), width = 5, height = 5, dpi = 90, bg = "white")
  ggsave(file.path("images", paste0(g, "_MeJA.png")),    feat(g, "MeJA",    cap), width = 5, height = 5, dpi = 90, bg = "white")
  ggsave(file.path("images", paste0(g, "_violin.png")),  vln(g),                  width = 11, height = 3.8, dpi = 90, bg = "white")
  if (i %% 25 == 0) cat(i, "/", length(ug), "\n")
}
cat("rendered", length(ug), "genes\n")

## ---- genes.js: [{id,label,cat}] ordered JA -> cell wall -> clusters 0..14 ----
g <- genes; g$label <- ifelse(g$symbol == "", g$gene_id, paste0(g$symbol, "  (", g$gene_id, ")"))
clus <- unique(g$category[grepl("^[0-9]+:", g$category)]); clus <- clus[order(as.numeric(sub(":.*", "", clus)))]
g$category <- factor(g$category, levels = c("JA / defense", "Cell wall / secondary metabolism", clus))
g <- g[order(g$category), ]
rows <- apply(g, 1, function(r) sprintf("{id:%s,label:%s,cat:%s}",
             shQuote(r["gene_id"]), shQuote(trimws(r["label"])), shQuote(as.character(r["category"]))))
cat("var GENES=[", paste(rows, collapse = ","), "];", sep = "", file = "genes.js")
cat("wrote genes.js (", nrow(g), "rows )\nDone.\n")
