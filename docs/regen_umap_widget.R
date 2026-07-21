## Regenerate only the UMAP plotly widget with tighter margins (less top whitespace).
library(Seurat)
library(plotly)
library(htmlwidgets)

obj <- readRDS("../shiny_app/webapp_data.rds")
ct_levels <- unique(obj$cell_type)[order(as.numeric(sub(":.*", "", unique(obj$cell_type))))]
obj$cell_type <- factor(obj$cell_type, levels = ct_levels)

## exact publication palette from plot_umap_final.R (cluster 0..14 order),
## mapped by cluster number so the hover UMAP matches UMAP_annotated.pdf
pal_hex <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#7B3FA0", "#9C755F", "#FF9DA7",
             "#B07AA1", "#EDC948", "#76B7B2", "#00A0B0", "#2E5A87", "#8CD17D", "#8C1D40", "#B6992D")
pal <- setNames(pal_hex[order(as.numeric(sub(":.*", "", ct_levels)))], ct_levels)
emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb) <- c("UMAP_1", "UMAP_2")
umap_df <- data.frame(UMAP_1 = emb$UMAP_1, UMAP_2 = emb$UMAP_2, cell_type = obj$cell_type)
cent <- do.call(rbind, lapply(ct_levels, function(l) {
  d <- umap_df[umap_df$cell_type == l, ]
  data.frame(cell_type = l, UMAP_1 = median(d$UMAP_1), UMAP_2 = median(d$UMAP_2))
}))
cent$lab <- sub(":.*", "", cent$cell_type)   # short cluster number only, to avoid label overlap

umap_plotly <- plot_ly(umap_df, x = ~UMAP_1, y = ~UMAP_2, color = ~cell_type,
                       colors = pal,
                       type = "scattergl", mode = "markers", marker = list(size = 4),
                       hoverinfo = "text", text = ~cell_type, showlegend = FALSE) %>%
  add_annotations(data = cent, x = ~UMAP_1, y = ~UMAP_2, text = ~lab, showarrow = FALSE,
                  font = list(size = 13, color = "black"), bgcolor = "rgba(255,255,255,0.75)") %>%
  layout(paper_bgcolor = "white", plot_bgcolor = "white",
         margin = list(t = 6, b = 40, l = 50, r = 12),          # tight top margin
         xaxis = list(title = "UMAP 1", zeroline = FALSE),
         yaxis = list(title = "UMAP 2", zeroline = FALSE, scaleanchor = "x", scaleratio = 1))
saveWidget(umap_plotly, "umap_widget.html", selfcontained = FALSE, libdir = "umap_lib", title = "Cell types")
cat("regenerated umap_widget.html\n")
