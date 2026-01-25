
############################################
## scCONGAS Core 1 – R script version
## Converted from .Rmd to .R
############################################

## --- Environment and libraries ----
reticulate::use_condaenv("CONGAS", required = TRUE)
library(reticulate)

library(Rcongas)
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(data.table)
library(cowplot)

## --- Read multiome data ----
multi <- Read10X(
  "/group/sottoriva/00-PROCESSED_DATA/2023-ROVIGO/pre_demultiplexing/BTHROHT01_BTHROHT02_BTHROHT04/multiome/LAZ_315/filtered_feature_bc_matrix"
)

myrna_counts <- multi$`Gene Expression`
myatac_counts <- multi$`Peaks`

# Replace ":" with "-" in ATAC peak names
rownames(myatac_counts) <- gsub(":", "-", rownames(myatac_counts))

## --- Read barcode list and subset ----
mybarcode_list <- read.table(
  "/home/fatemeh.nasrabadi/ROVIGO_congas/BTROHT01/Core/CONGAS/BTROHT01_core_barcodes_all.tsv",
  sep = "\t",
  header = TRUE
)

mybarcode_values <- mybarcode_list$X

myrna_counts  <- myrna_counts[, colnames(myrna_counts) %in% mybarcode_values]
myatac_counts <- myatac_counts[, colnames(myatac_counts) %in% mybarcode_values]

## --- Read features ----
myfeatures <- fread(
  "gunzip -c /group/sottoriva/00-PROCESSED_DATA/2023-ROVIGO/pre_demultiplexing/BTHROHT01_BTHROHT02_BTHROHT04/multiome/LAZ_315/filtered_feature_bc_matrix/features.tsv.gz"
)

colnames(myfeatures) <- c("gene_id", "gene", "modality", "chr", "from", "to")

## --- ATAC feature dataframe ----
myatac_featureDF <- data.frame(id = rownames(myatac_counts)) %>%
  tidyr::separate(id, c("chr", "from", "to"), "-")

## --- Create CONGAS tibbles ----
myatac <- create_congas_tibble(
  counts   = myatac_counts,
  modality = "ATAC",
  save_dir = NULL,
  features = myatac_featureDF
)

myatac <- na.omit(myatac)

## --- Metadata construction ----
filtered_barcodes <- mybarcode_list %>%
  filter(X %in% myatac$cell)

metadata_atac <- filtered_barcodes %>%
  mutate(cell = paste0(X, "-ATAC"), type = celltype) %>%
  select(cell, type)

metadata_rna <- filtered_barcodes %>%
  mutate(cell = paste0(X, "-RNA"), type = celltype) %>%
  select(cell, type)

mymetadata <- bind_rows(metadata_atac, metadata_rna)

## --- Normalisation factors ----
mynorm_atac <- Rcongas:::auto_normalisation_factor(myatac) %>%
  mutate(modality = "ATAC")

myrna <- create_congas_tibble(
  counts   = myrna_counts,
  modality = "RNA",
  save_dir = NULL,
  features = features
)

mynorm_rna <- Rcongas:::auto_normalisation_factor(myrna) %>%
  mutate(modality = "RNA")

myatac <- myatac %>% mutate(value = as.integer(value))
myrna  <- myrna  %>% mutate(value = as.integer(value))

## --- Filter RNA genes ----
myrna <- myrna %>% filter_known_genes(what = "r")

all_genes <- unique(myrna$gene)
mito_genes <- all_genes[str_starts(all_genes, "MT-")]
myrna <- myrna %>% filter(!gene %in% mito_genes)

## --- Read CNV segments ----
mysegs <- read.csv(
  "/home/fatemeh.nasrabadi/ROVIGO_congas/BTROHT01/Core/CONGAS/BTROHT01_copre_segs_for_congas.csv",
  sep = "\t",
  row.names = 1
)

mysegments <- mysegs %>%
  filter(chr != "chrX", chr != "chrY")

## --- Initialize CONGAS object ----
myx <- init(
  rna  = myrna,
  atac = myatac,
  segmentation = mysegments,
  rna_normalisation_factors  = mynorm_rna,
  atac_normalisation_factors = mynorm_atac,
  rna_likelihood  = "NB",
  atac_likelihood = "NB",
  description = "Bimodal sample",
  multiome = TRUE
)

## --- QC plots ----
pdf("plot1.pdf", width = 12, height = 8)
ggplot(myx$input$segmentation,
       aes(x = ATAC_nonzerocells, y = ATAC_peaks)) +
  geom_point() +
  theme_bw()
dev.off()

pdf("plot2.pdf", width = 12, height = 8)
ggplot(myx$input$segmentation,
       aes(x = RNA_nonzerocells, y = RNA_genes)) +
  geom_point() +
  theme_bw()
dev.off()

## --- Filter segments ----
myx <- filter_segments(
  myx,
  RNA_nonzerocells  = 300,
  ATAC_nonzerocells = 300
)

## --- Count distributions ----
pdf("countdistribution_RNA_ATAC.pdf", width = 12, height = 8)
plot_data(
  myx,
  what = "histogram",
  segments = get_input(myx, "segmentation") %>%
    ungroup() %>%
    mutate(L = to - from) %>%
    arrange(desc(L)) %>%
    top_n(30) %>%
    pull(segment_id)
)
dev.off()

## --- Annotated distributions ----
myx$input$metadata <- mymetadata

pdf("countdistribution_annotated.pdf", width = 12, height = 8)
plot_data(
  myx,
  to_plot  = "type",
  position = "stack",
  segments = get_input(myx, "segmentation") %>%
    ungroup() %>%
    mutate(L = to - from) %>%
    arrange(desc(L)) %>%
    top_n(30) %>%
    pull(segment_id)
)
dev.off()

## --- Segment selection ----
filt <- segments_selector_congas(myx)

## --- Model configuration ----
k <- 1:4
binom_limits <- c(40, 1000)

hyperparams_filt <- auto_config_run(
  filt,
  k,
  prior_cn = c(0.2, 0.6, 0.1, 0.05, 0.05),
  init_importance = 0.6,
  CUDA = FALSE,
  normal_cells = TRUE
)

hyperparams_filt$binom_prior_limits <- binom_limits

fit_filt <- Rcongas:::fit_congas(
  filt,
  K = k,
  lambdas = 0.5,
  learning_rate = 0.01,
  steps = 5000,
  model_parameters = hyperparams_filt,
  model_selection = "BIC",
  latent_variables = "G",
  CUDA = FALSE,
  temperature = 20,
  same_mixing = TRUE,
  threshold = 0.001
)

## --- Save results ----
saveRDS(
  fit_filt,
  "/home/fatemeh.nasrabadi/ROVIGO_congas/BTROHT01/Core/CONGAS/rcongas_obj_filtered_BTROHT01_core_all_cells.rds"
)

pdf("fitplotb1core.pdf", width = 12, height = 8)
plot_fit(fit_filt, what = "scores")
dev.off()

pdf("CONGAS_clusterb1core.pdf", width = 12, height = 8)
plot_fit(fit_filt, "posterior_CNA")
dev.off()

pdf("density_fitplotb1core.pdf", width = 20, height = 16)
cowplot::plot_grid(
  plotlist = plot_fit(fit_filt, what = "density", highlights = TRUE),
  ncol = 4
)
dev.off()

pdf("heatmap_b1core.pdf", width = 20, height = 16)
plot_fit(
  fit_filt,
  what = "heatmap",
  scale = TRUE,
  scale_min_lim = -2,
  scale_max_lim = 2
)
dev.off()

## --- Export tables ----
df <- data.frame(
  cluster    = fit_filt$best_fit$CNA_real$cluster,
  value      = fit_filt$best_fit$CNA_real$value,
  segment_id = fit_filt$best_fit$CNA_real$segment_id
)

write.table(
  df,
  file = "/home/fatemeh.nasrabadi/ROVIGO_congas/BTROHT01/Core/CONGAS/ploidy_by_segments_BTROHT01_core.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  fit_filt$best_fit$cluster_assignments,
  file = "/home/fatemeh.nasrabadi/ROVIGO_congas/BTROHT01/Core/CONGAS/BTROHT01_CONGAS_Clusters.csv",
  row.names = FALSE
)

############################################
## End of script
############################################
