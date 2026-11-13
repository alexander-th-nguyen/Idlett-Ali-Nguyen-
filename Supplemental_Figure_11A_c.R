# load packages
library(conflicted)
library(SummarizedExperiment)
library(forcats)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(PhosR)
library(imputeLCMD)
conflicts_prefer(base::setdiff)

# function for enrichment plots
enrich_curve <- function(
    gene.set.id, gsea_res, dge,
    plot.stats = TRUE,
    # Positioning of stats
    just = "left",
    x = 1.02,
    y = 0.95,
    y.gap = 0.07
) {
  entry <- gsea_res@result |>
    dplyr::filter(ID == gene.set.id)
  
  nes <- entry |>
    dplyr::pull(NES) %>%
    scales::label_number(0.0001)(.) %>%
    paste0("NES = ", .)
  es <- entry |>
    dplyr::pull(enrichmentScore) %>%
    scales::label_number(0.0001)(.) %>%
    paste0("ES = ", .)
  qval <- entry |>
    dplyr::pull(qvalue) %>%
    scales::label_number(0.0001)(.) %>%
    paste0("q = ", .)
  
  enrichplot::gseaplot(
    x = gsea_res,
    geneSetID = gene.set.id,
    by = "runningScore",
    title = gene.set.id
  ) +
    ggplot2::annotation_custom(
      grid::textGrob(nes,  x, y,               just = just)
    ) +
    ggplot2::annotation_custom(
      grid::textGrob(es,   x, y - y.gap,       just = just)
    ) +
    ggplot2::annotation_custom(
      grid::textGrob(qval, x, y - (2 * y.gap), just = just)
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "lightgrey")
    )
}

## import summarized experiment ####
summarized_experiment <- readRDS("Replacing_with_your_file_path")

## import dge ####
twopointfive_vs_baseline_dge <- readRDS("Replacing_with_your_file_path")
row.names(twopointfive_vs_baseline_dge) <- twopointfive_vs_baseline_dge$protein_id
stopifnot(!(any(duplicated(twopointfive_vs_baseline_dge$protein_id))))

# store organism
organism <- "Homo sapiens"

## GSEA for 1 hour versus 0 hour ####

# check duplicated or NA counts
sum(is.na(twopointfive_vs_baseline_dge$gene))
sum(duplicated(twopointfive_vs_baseline_dge$gene))

# filter for non-NA gene rows
twopointfive_vs_baseline_dge <- twopointfive_vs_baseline_dge |>
  dplyr::filter(!is.na(gene))

# check number of duplicates
twopointfive_vs_baseline_dge$gene |>
  duplicated() |>
  sum()

# find duplicate gene rows
twopointfive_vs_baseline_dge[duplicated(twopointfive_vs_baseline_dge$gene) | duplicated(twopointfive_vs_baseline_dge$gene, fromLast = TRUE), ]

# get protein ids of duplicate gene rows
dupe_ids <- twopointfive_vs_baseline_dge[duplicated(twopointfive_vs_baseline_dge$gene) | duplicated(twopointfive_vs_baseline_dge$gene, fromLast = TRUE), ] |>
  dplyr::pull(protein_id)

# find top ave expression for a gene group
twopointfive_vs_baseline_dge |>
  dplyr::filter(protein_id %in% dupe_ids) |>
  dplyr::slice_max(order_by = AveExpr, n = 1, by = gene)

# filter duplicate gene rows
nrow(twopointfive_vs_baseline_dge)
twopointfive_vs_baseline_dge <- twopointfive_vs_baseline_dge |>
  dplyr::slice_max(order_by = AveExpr, n = 1, by = gene)
nrow(twopointfive_vs_baseline_dge)

# check to make sure there are no NA or duplicate gene rows
stopifnot(!all(is.na(twopointfive_vs_baseline_dge$gene)))
stopifnot(!all(duplicated(twopointfive_vs_baseline_dge$gene)))

stopifnot(exprs = {
  "gene" %in% colnames(twopointfive_vs_baseline_dge)
  "logFC" %in% colnames(twopointfive_vs_baseline_dge)
})

paths_hallmark <- msigdbr::msigdbr(
  species = organism,
  category = "H"
) |>
  dplyr::select(gs_name, gene_symbol)

twopointfive_vs_baseline_dge <- twopointfive_vs_baseline_dge |>
  dplyr::arrange(dplyr::desc(logFC))

fc_symbol <- twopointfive_vs_baseline_dge |>
  dplyr::pull(logFC, gene)

stopifnot(exprs = {
  !any(is.na(fc_symbol))
})

set.seed(381272)
gsea_res <- clusterProfiler::GSEA(
  geneList = fc_symbol,
  TERM2GENE = paths_hallmark,
  minGSSize = 10,
  maxGSSize = Inf,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  by = "fgsea",
  nPermSimple = 100000
)

# generate enrichment curve
enrich_curve("HALLMARK_COAGULATION", gsea_res, twopointfive_vs_baseline_dge)