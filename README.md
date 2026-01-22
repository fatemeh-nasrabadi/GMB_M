# Project Directory Structure

This repository contains code and documentation for multi-omics data analysis. Below is the directory structure of the project:


├── 1__10xMultiome_snATAC_snRNA
│   ├── 10x_snATAC_snRNA analysis_Core1.ipynb
│   └── scCONGAS_Core1.Rmd
│
├── 2_GEX
│   ├── 2_GEX_Aggrigation
│   │   └── FINAL - Data aggregation - Gex.ipynb
│   ├── 3_GEX_AGGREGATED_VISUALIZATION_SIGNATURES
│   │   └── GEX_VISUALIZATION_SIGNATURES.ipynb
│   ├── 4_GEX_DGE_PathwayENRICHMENT
│   │   └── DGE_PATHWAYENRICHMDNT.ipynb
│   └── 5_GEX_HOTSPOT_LIANA
│   ├── HOTSPOT.ipynb
│   └── LIANA.ipynb
│
├── 3_ATAC
│   ├── 1_ATAC_Aggrigation
│   │   └── AGG_scATAC_10Samples.ipynb
│   ├── 2_Peakcalling_mutif_enrichment
│   │   ├── Peakcalling_motifenrichment_TFsvisualization.ipynb
│   │   └── Peaks_annotation_ROVIGO.ipynb
│   └── Track_plots_SCAAs.Rmd
│
├── 4_WGS
│   ├── Cohort_WGS_Oncoplot.ipynb
│   ├── ROVIGO_WGS_samplecheck.ipynb
│   └── SparseSignatures_ex.Rmd
│
└── All_figures.pdf



#### Project Description:

This repository contains the analysis pipeline and results of a comprehensive study investigating the molecular heterogeneity of glioblastoma (GBM), with a specific focus on comparing tumor cells from the core and perimarginal zones. The project employs an integrated approach using single-cell RNA sequencing (scRNA-seq) and single-cell ATAC sequencing (scATAC-seq) to uncover gene expression and chromatin accessibility profiles of GBM cells.

#### Folder Breakdown:

1. **1__10xMultiome_snATAC_snRNA**:

   * This folder contains the analysis notebooks for single-cell RNA and ATAC sequencing data. The `10x_snATAC_snRNA analysis_Core1.ipynb` and `scCONGAS_Core1.Rmd` scripts focus on data preprocessing, quality control, and integration of the gene expression (GEX) and chromatin accessibility (ATAC) modalities using the 10x Genomics multiome technology.

2. **2_GEX**:

   * This folder organizes various steps in gene expression analysis.

     * **2_GEX_Aggrigation**: Contains the `FINAL - Data aggregation - Gex.ipynb` script, which performs data aggregation for gene expression across multiple samples.
     * **3_GEX_AGGREGATED_VISUALIZATION_SIGNATURES**: Includes `GEX_VISUALIZATION_SIGNATURES.ipynb`, focusing on visualizing gene expression signatures.
     * **4_GEX_DGE_PathwayENRICHMENT**: Houses the `DGE_PATHWAYENRICHMDNT.ipynb` script for differential gene expression analysis and pathway enrichment.
     * **5_GEX_HOTSPOT_LIANA**: Contains scripts for hotspot analysis and ligand-receptor interaction (LIANA) analysis, with `HOTSPOT.ipynb` and `LIANA.ipynb`.

3. **3_ATAC**:

   * This section deals with the analysis of chromatin accessibility using ATAC sequencing.

     * **1_ATAC_Aggrigation**: The `AGG_scATAC_10Samples.ipynb` script aggregates ATAC-seq data from multiple samples.
     * **2_Peakcalling_mutif_enrichment**: Includes scripts for peak calling and motif enrichment (`Peakcalling_motifenrichment_TFsvisualization.ipynb` and `Peaks_annotation_ROVIGO.ipynb`).
     * **Track_plots_SCAAs.Rmd**: Generates track plots for visualizing chromatin accessibility.

4. **4_WGS**:

   * This folder contains the analysis scripts for whole-genome sequencing data.

     * **Cohort_WGS_Oncoplot.ipynb**: Produces oncoplots to visualize the genomic alterations.
     * **ROVIGO_WGS_samplecheck.ipynb**: Validates and checks the integrity of the WGS samples.
     * **SparseSignatures_ex.Rmd**: Explores sparse signatures from the whole-genome sequencing data.

5. **All_figures.pdf**:

   * A compilation of all figures generated from the analysis, summarizing the key findings of the study.

The primary goal of this project is to investigate the molecular and epigenetic differences between glioblastoma tumor cells in the core and perimarginal regions, which are known for their distinct behaviors in tumor progression and recurrence. By applying cutting-edge multi-omics techniques, the study identifies key genes, transcriptional programs, and chromatin accessibility profiles associated with the invasive properties of perimarginal tumor cells. This will help improve our understanding of glioblastoma heterogeneity and potentially lead to more effective therapeutic strategies targeting the perimarginal zone.
