### **Project Description**:

This repository contains the analysis pipeline and results of a study focused on the molecular heterogeneity of glioblastoma (GBM), comparing tumor cells from the core and perimarginal zones. The project utilizes multi-omics techniques, including **single-cell RNA sequencing (scRNA-seq)** and **single-cell ATAC sequencing (scATAC-seq)**, to investigate gene expression and chromatin accessibility profiles of GBM cells. By comparing the core and perimarginal tumor regions, the project aims to identify key genes, transcriptional programs, and chromatin accessibility profiles associated with the invasive behavior of perimarginal tumor cells, advancing our understanding of glioblastoma heterogeneity and informing potential therapeutic strategies.

---

### **Folder Breakdown**:

#### **1__10xMultiome_snATAC_snRNA**:

* **10x_snATAC_snRNA analysis_Core1.ipynb**:

  * **Quality Control (QC) for RNA Data**
  * **RNA Data Filtration**
  * **Dimensionality Reduction and Clustering**
  * **Cell Type Annotation with CellTypist**
  * **Gene Expression Analysis**
  * **Copy Number Variation (CNV) Analysis**
  * **ATAC-Seq Analysis**
  * **ATAC Data Normalization and Dimensionality Reduction**
  * **Integration of RNA and ATAC Data**
  * **Data Export**
  * **Final Outputs**
* **scCONGAS_Core1.Rmd**: Data processing and analysis for multiome integration.

#### **2_GEX** (Gene Expression Analysis):

* **FINAL - Data aggregation - Gex.ipynb**: Data integration for GEX modality and batch effect correction.
* **GEX_VISUALIZATION_SIGNATURES.ipynb**: Visualizes gene expression signatures.
* **DGE_PATHWAYENRICHMDNT.ipynb**: Differential gene expression analysis and pathway enrichment.
* **HOTSPOT.ipynb and LIANA.ipynb**: Hotspot and ligand-receptor interaction analysis.

#### **3_ATAC** (Chromatin Accessibility Analysis):

* **AGG_scATAC_10Samples.ipynb**: Aggregates ATAC-seq data.
* **Peakcalling_motifenrichment_TFsvisualization.ipynb and Peaks_annotation_ROVIGO.ipynb**: Peak calling and motif enrichment.
* **Track_plots_SCAAs.Rmd**: Chromatin accessibility track plots.

#### **4_WGS** (Whole-Genome Sequencing Analysis):

* **Cohort_WGS_Oncoplot.ipynb**: Visualizes genomic alterations.
* **ROVIGO_WGS_samplecheck.ipynb**: Validates WGS samples.
* **SparseSignatures_ex.Rmd**: Sparse signature analysis for WGS data.

#### **All_figures.pdf**:

* Compilation of all figures summarizing the key findings.

---

### **Project Objective**:

The primary goal of this project is to explore and compare the molecular and epigenetic features of glioblastoma tumor cells from the core and perimarginal regions. By leveraging multi-omics technologies, the study identifies genes, transcriptional programs, and chromatin accessibility profiles linked to the invasive behavior of perimarginal tumor cells. The insights gained will help improve our understanding of glioblastoma heterogeneity and potentially lead to more effective therapeutic strategies targeting the perimarginal zone.

