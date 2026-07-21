# A Single-Cell Transcriptomic Atlas of the *Populus* Root Uncovers Cell Type-Specific Methyl Jasmonate Response

Single-nucleus RNA-seq analysis of how the defense phytohormone methyl jasmonate
(MeJA) reprograms the *Populus* (poplar) root. Comparing **control (ctr)** and
**MeJA-treated (mej)** nuclei within each annotated cell type resolves the
cell-type-specific induction of jasmonate signaling, defense, and cell wall /
secondary-metabolism (phenylpropanoid, lignin) pathways.

## Layout

```
.
├── 00_integration_annotation.R     # CCA integration of the 6 samples + omg cell-type annotation
├── 10_setup_umap_markers_degs.R    # object, UMAP, cluster markers, MeJA-vs-Control DEGs
├── 20_deg_barchart_upset_go.R      # DEG bar charts, UpSet, GO BP enrichment dot plot
├── 30_ja_dotplot_and_modulescore.R # JA/defense dot plot + module score boxplot
├── 40_cellwall_all_plots.R         # all cell wall / secondary-metabolism figures
│
├── marker_genes_all_clusters_PC34_res0.25.csv   # cluster marker genes (FindAllMarkers)
├── up_down_DEG_MeJA_counts_per_cluster.csv       # up / down DEG counts per cluster
├── omg_output_PC34_res0.25_top200/               # omg cell-type annotation output
│
├── shiny_app/                       # interactive gene-expression web app (R Shiny)
│
├── input/    # raw inputs (not tracked) — filtered 10x matrices / the Seurat object
└── output/   # where the scripts write figures & tables by default (not tracked)
```

## How to run

Run the scripts from the repo root, in order:

- **`00`** integrates the raw matrices into one object. Its omg predictions are
  used to curate the `cell_type` labels, which are added to the object.
- **`10`–`40`** start from that annotated object. **Run `10` first** — it writes
  the marker and DEG tables that `20` reads.

Each script uses `input/` and `output/` by default, or takes them as arguments:

```sh
# defaults: input/ and output/
Rscript 10_setup_umap_markers_degs.R

# custom paths (input first, output second)
Rscript 10_setup_umap_markers_degs.R /path/to/input /path/to/output
```

You can also just open a script in RStudio and run it — no arguments needed.

## Inputs (place in `input/`)

- For `00`: the six filtered 10x matrices under `filtered-ctr/filtered_ctr{1,2,3}`
  and `filtered-mej/filtered_mej{1,2,3}`.
- For `10-40`: `scRNA_integrated_all_samples_PC34_res0.25_annotated.rds` — the
  integrated + annotated Seurat object (RNA + integrated assays, UMAP,
  `cell_type`, `condition`, `seurat_clusters`) produced by `00`.

The raw matrices and the integrated Seurat object are large and are **not**
committed here — they are deposited in a public database (TBD).

## Included results

The key derived tables are committed so they can be browsed without re-running
the pipeline (figures and `.rds` objects are not tracked):

- **`marker_genes_all_clusters_PC34_res0.25.csv`** — marker genes for each
  cluster (`Seurat::FindAllMarkers`, positive markers). Columns: `p_val`,
  `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`, `cluster`, `gene`.
- **`up_down_DEG_MeJA_counts_per_cluster.csv`** — number of up- and
  down-regulated MeJA-vs-Control DEGs (`p_val_adj < 0.05`, `|avg_log2FC| > 0.25`)
  per cluster. Columns: `cluster`, `Up`, `Down`.
- **`omg_output_PC34_res0.25_top200/`** — cell-type annotation from
  [omg](https://github.com/LiLabAtVT/OrthoMarkerGeneGroups) (top-200 markers/cluster):
  - `cell_type_predictions.csv` — per-cluster prediction, with the
    `consolidated_cell_type` call and a `prediction_confidence`.
  - `compare_15species_all.csv` / `compare_15species_heatmap.pdf` — orthogroup
    marker matches across 15 reference species.
  - `extract_table_significant.csv` — significant orthogroup–cell-type matches.
  - `pairwise/` — per reference species × tissue significance tables and heatmaps.

Running the scripts regenerates equivalent tables under `output/`.

## Requirements

R (≥ 4.4) with `Seurat` (v5), `ggplot2`, `dplyr`, `Matrix`, `pheatmap`,
`RColorBrewer`, `UpSetR`, `scattermore`, and for GO enrichment
`clusterProfiler`, `org.At.tair.db`, `biomaRt`. Script `00` also uses
`omg` ([LiLabAtVT/OrthoMarkerGeneGroups](https://github.com/LiLabAtVT/OrthoMarkerGeneGroups),
installed from GitHub on first run).

## Contact

Questions: Miaomiao Li (lim7@ornl.gov), Tran Chau (tnchau@vt.edu).
