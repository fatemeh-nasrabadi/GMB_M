### **10x Single Cell Multiome Analysis: RNA + ATAC and Genomic Landscape of GBM Cohort**

* **Data Sources**:

  * **10x multiome snRNA-seq** and **snATAC-seq** raw multiome data.
  * **100x WGS analysis**: Genomic landscape of the GBM cohort.

* **Analysis Flow**:

  * Sample analysis of gene expression and chromatin accessibility.
  * Data integration of RNA and ATAC modalities.

* **Downstream Analysis**:

  * **Gene Expression**: Differential expression, pathway enrichment analysis, **Hotspot** (identifies informative genes and gene modules), and **LIANA** (ligand-receptor interaction analysis).
  * **Chromatin Accessibility**: ATAC-seq analysis for regulatory regions.

* **Whole-Genome Sequencing (WGS)**:

  * Analyzes **genomic landscape** of the GBM cohort, identifying genomic alterations and mutational signatures.

### **Folder Breakdown**:

#### **1__10xMultiome_snATAC_snRNA**:

* **10x_snATAC_snRNA analysis_Core1.ipynb**:

  * **Single-Cell Quality Control (QC) for RNA Data**
  * **Single-Cell RNA Data Filtration**
  * **Single-Cell Dimensionality Reduction and Clustering**
  * **Cell Type Annotation with CellTypist for Single-Cell Data**
  * **Single-Cell Gene Expression Analysis**
  * **Single-Cell Copy Number Variation (CNV) Analysis**
  * **Single-Cell ATAC-Seq Analysis**
  * **ATAC Data Normalization and Dimensionality Reduction for Single-Cell Data**
  * **Integration of Single-Cell RNA and ATAC Data**
  * **Data Export**
  * **Final Outputs**

* **scCONGAS_Core1.Rmd**: Unsupervised single-cell clustering based on CNVs and multiome data.

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
