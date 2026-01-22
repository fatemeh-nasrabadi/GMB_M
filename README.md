### Project Description:

This repository contains the analysis pipeline and results of a study focused on the molecular heterogeneity of glioblastoma (GBM), comparing tumor cells from the core and perimarginal zones. Using single-cell RNA sequencing (scRNA-seq) and single-cell ATAC sequencing (scATAC-seq), the project investigates gene expression and chromatin accessibility profiles of GBM cells.

#### **Folder Breakdown**:

1. **1__10xMultiome_snATAC_snRNA**:

   * Notebooks for data preprocessing, quality control, and integration of gene expression and chromatin accessibility data using 10x Genomics multiome technology (`10x_snATAC_snRNA analysis_Core1.ipynb`, `scCONGAS_Core1.Rmd`).

2. **2_GEX**:

   * Gene expression analysis:

     * `FINAL - Data aggregation - Gex.ipynb`: Aggregates gene expression data.
     * `GEX_VISUALIZATION_SIGNATURES.ipynb`: Visualizes gene expression signatures.
     * `DGE_PATHWAYENRICHMDNT.ipynb`: Performs differential gene expression analysis and pathway enrichment.
     * `HOTSPOT.ipynb` and `LIANA.ipynb`: Conducts hotspot and ligand-receptor interaction analysis.

3. **3_ATAC**:

   * Chromatin accessibility analysis:

     * `AGG_scATAC_10Samples.ipynb`: Aggregates ATAC-seq data.
     * `Peakcalling_motifenrichment_TFsvisualization.ipynb` and `Peaks_annotation_ROVIGO.ipynb`: Perform peak calling and motif enrichment.
     * `Track_plots_SCAAs.Rmd`: Generates chromatin accessibility track plots.

4. **4_WGS**:

   * Whole-genome sequencing analysis:

     * `Cohort_WGS_Oncoplot.ipynb`: Visualizes genomic alterations.
     * `ROVIGO_WGS_samplecheck.ipynb`: Checks WGS sample integrity.
     * `SparseSignatures_ex.Rmd`: Analyzes sparse signatures in WGS data.

5. **All_figures.pdf**:

   * Compilation of all figures summarizing the key findings.

### **Project Objective**:

This project aims to compare the molecular and epigenetic features of glioblastoma tumor cells in the core and perimarginal regions. By using multi-omics techniques, it identifies genes, transcriptional programs, and chromatin accessibility profiles linked to the invasive behavior of perimarginal tumor cells, enhancing our understanding of glioblastoma heterogeneity and informing potential therapeutic strategies.
