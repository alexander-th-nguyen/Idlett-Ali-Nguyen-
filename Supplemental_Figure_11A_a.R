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

# create paths for protein data and metadata
protein_path <- "Replacing_with_your_file_path"
metadata_path <- "Replacing_with_your_file_path"

## import metadata ####
metadata <- read_xlsx(path = metadata_path)

# set factors for groups
metadata$group <- factor(x = metadata$group, levels = c(
  "Baseline", 
  "1 Hour", 
  "2.5 Hours" 
  ))

## import protein data ####
protein_data <- read_xlsx(path = protein_path, .name_repair = "minimal")

# select necessary columns and match with metadata
protein_data <- protein_data |>
  select(
    "protein_id" = "Protein.Names",
    "gene"       = "Genes",
    "desc"       = "First.Protein.Description",
    all_of(metadata$data_colname)
  )

## convert data into a SummarizedExperiment object

stopifnot(exprs = {
  sum(is.na(protein_data$protein_id)) == 0
  sum(duplicated(protein_data$protein_id)) == 0
  sum(duplicated(metadata$sample)) == 0
  sum(duplicated(metadata$data_colname)) == 0
})

temp_matrix <- protein_data |>
  select(protein_id, all_of(metadata$data_colname)) |>
  column_to_rownames(var = "protein_id") |>
  as.matrix()
stopifnot(identical(colnames(temp_matrix), metadata$data_colname))
colnames(temp_matrix) <- metadata$sample

temp_rowdata <- protein_data |> select(protein_id, gene, desc) |> as.data.frame()
row.names(temp_rowdata) <- temp_rowdata$protein_id
stopifnot(identical(row.names(temp_rowdata), row.names(temp_matrix)))

temp_coldata <- metadata |> as.data.frame()
row.names(temp_coldata) <- temp_coldata$sample
stopifnot(identical(row.names(temp_coldata), colnames(temp_matrix)))

protein <- SummarizedExperiment(
  assays = list(intensity = temp_matrix),
  rowData = temp_rowdata,
  colData = temp_coldata
)

stopifnot(exprs = {
  identical(row.names(protein), row.names(rowData(protein)))
  identical(row.names(protein), rowData(protein)$protein_id)
  identical(colnames(protein), row.names(colData(protein)))
  identical(colnames(protein), colData(protein)$sample)
})

## log2 transformation ####

# check zero and NA expression
sum(assays(protein)[["intensity"]] %in% 0)
sum(is.na(assays(protein)[["intensity"]]))

# convert zeroes to NA
assays(protein)[["intensity"]] <- assays(protein)[["intensity"]] |>
  apply(2, \(x) {
    x[x == 0] <- NA
    x
  })

# log2 transform in a new assay slot
assays(protein)[["log2"]] <- assays(protein)[["intensity"]] |>
  apply(2, \(x) {
    x <- log2(x)
    x
  })

## filter lowly expressed proteins ####
rows_before_num <- nrow(protein)

# keep rows where at least 50% of replicates in each group are non-missing in all 3 of the experimental groups
rows_keep <- PhosR::selectGrps(
  mat = assays(protein)[["log2"]],
  grps = as.character(colData(protein)$group),
  percent = 0.5,
  n = 3
) |>
  row.names()

protein <- protein[rows_keep, ]

rows_after_num <- nrow(protein)

rows_before_num
rows_after_num

## VSN normalization ####

assays(protein)[["normalized"]] <- assays(protein)[["intensity"]]
fit <- vsn::vsnMatrix(x = assays(protein)[["normalized"]])
assays(protein)[["normalized"]] <- vsn::predict(
  object = fit,
  newdata = assays(protein)[["normalized"]],
  log2scale = TRUE,
  useDataInFit = FALSE
)

## minProb imputation ####
set.seed(391283)
assays(protein)[["minProb"]] <- imputeLCMD::impute.MinProb(
  assays(protein)[["normalized"]],
  q = 0.01,
  tune.sigma = 1
)

## run DGE on minProb imputed data set ####

experiment <- assays(protein)[["minProb"]]
group <- colData(protein)$group |>
  fct_recode(
    "Baseline" = "Baseline",
    "Hour_1" = "1 Hour",
    "Hour_2.5" = "2.5 Hours",
  )

design <- stats::model.matrix(~0 + group)
colnames(design) <- sub("^group", "", colnames(design))

fit <- limma::lmFit(
  object = experiment,
  design = design
)

matrix <- limma::makeContrasts(
  "1_Hour_vs_Baseline" = Hour_1 - Baseline,
  "2.5_Hours_vs_Baseline" = Hour_2.5 - Baseline,
  levels = design
)

fit <- limma::contrasts.fit(fit = fit, contrasts = matrix)
fit <- limma::eBayes(fit = fit)

## DGEs ####

# 1 hour vs baseline comparison
dge <- limma::topTable(
  fit = fit,
  adjust = "BH",
  coef = "1_Hour_vs_Baseline",
  number = Inf
) |>
  rownames_to_column(var = "protein_id") |>
  left_join(y = as.data.frame(rowData(protein)), by = "protein_id") |>
  relocate(protein_id, gene, desc)

# 2.5 hours vs baseline comparison
dge2 <- limma::topTable(
  fit = fit,
  adjust = "BH",
  coef = "2.5_Hours_vs_Baseline",
  number = Inf
) |>
  rownames_to_column(var = "protein_id") |>
  left_join(y = as.data.frame(rowData(protein)), by = "protein_id") |>
  relocate(protein_id, gene, desc)