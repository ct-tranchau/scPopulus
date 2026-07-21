## Re-render only the violin images at higher resolution (sharp at full page width).
library(Seurat)
library(ggplot2)

obj   <- readRDS("../shiny_app/webapp_data.rds")
genes <- read.csv("../shiny_app/genes.csv", stringsAsFactors = FALSE)
genes <- genes[genes$gene_id %in% rownames(obj), ]

obj$Condition <- factor(ifelse(grepl("mej", obj$condition, ignore.case = TRUE), "MeJA", "Control"),
                        levels = c("Control", "MeJA"))
ct_levels <- unique(obj$cell_type)[order(as.numeric(sub(":.*", "", unique(obj$cell_type))))]
obj$cell_type <- factor(obj$cell_type, levels = ct_levels)

vln <- function(gene) {
  d <- data.frame(expr = FetchData(obj, vars = gene, layer = "data")[, 1],
                  cell_type = obj$cell_type, Condition = obj$Condition)
  ggplot(d, aes(cell_type, expr, fill = Condition)) +
    geom_violin(scale = "width", linewidth = 0.2, position = position_dodge(0.85)) +
    scale_fill_manual(values = c(Control = "#0072B2", MeJA = "#C0392B")) +
    labs(x = NULL, y = "log-norm expression") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 15),
          legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 15))
}

ug <- unique(genes$gene_id)
for (i in seq_along(ug)) {
  ggsave(file.path("images", paste0(ug[i], "_violin.png")), vln(ug[i]),
         width = 16, height = 5.8, dpi = 200, bg = "white")
  if (i %% 25 == 0) cat(i, "/", length(ug), "\n")
}
cat("re-rendered", length(ug), "violin images (hi-res)\n")
