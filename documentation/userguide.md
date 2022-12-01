<!-- omit in toc -->
# GENE Platform Quickstart Guide

<!-- omit in toc -->
## Table of Contents
- [Preparing input files](#preparing-input-files)
  - [Examples](#examples)
    - [**DGE Analysis**](#dge-analysis)
    - [**Batch correction**](#batch-correction)
    - [**Correlation**](#correlation)
    - [**Normalization**](#normalization)
  - [Gene count matrix (counts.tsv)](#gene-count-matrix-countstsv)
    - [Required columns for count matrix](#required-columns-for-count-matrix)
    - [**Example**: RNA-Seq count matrix (counts.tsv)](#example-rna-seq-count-matrix-countstsv)
    - [**Example**: Microarray count matrix (counts.tsv)](#example-microarray-count-matrix-countstsv)
  - [Metadata for samples (coldata.tsv)](#metadata-for-samples-coldatatsv)
    - [Required columns for coldata.tsv](#required-columns-for-coldatatsv)
    - [**Example**: coldata.tsv file](#example-coldatatsv-file)
  - [Filter file (filter.txt)](#filter-file-filtertxt)
  - [Configuration file (config.yml)](#configuration-file-configyml)
- [Differential expression analysis](#differential-expression-analysis)
- [Data preparation](#data-preparation)
  - [Batch correction](#batch-correction-1)
  - [Preprocessing](#preprocessing)
  - [Normalization](#normalization-1)
- [RNA-Seq correlation](#rna-seq-correlation)

<div style="page-break-after: always;"></div>

## Preparing input files
### Examples
#### **DGE Analysis**
RNA-Seq: [counts.tsv](assets/example_data/analysis/rnaseq/counts.tsv)&nbsp; | &nbsp;[coldata.tsv](assets/example_data/analysis/rnaseq/coldata.tsv)&nbsp; | &nbsp;[filter.txt](assets/example_data/analysis/rnaseq/filter.txt)&nbsp; | &nbsp;[config.yml](assets/example_data/analysis/rnaseq/config.yml)  
Microarray: [counts.tsv](assets/example_data/analysis/microarray/counts.tsv)&nbsp; | &nbsp;[coldata.tsv](assets/example_data/analysis/microarray/coldata.tsv)&nbsp; | &nbsp;[filter.txt](assets/example_data/analysis/microarray/filter.txt)&nbsp; | &nbsp;[config.yml](assets/example_data/analysis/microarray/config.yml)

#### **Batch correction**
RNA-Seq: [counts.tsv](assets/example_data/batch_correction/rnaseq/counts.tsv)&nbsp; | &nbsp;[coldata.tsv](assets/example_data/batch_correction/rnaseq/coldata.tsv)

#### **Correlation**
[ohsu-cnl_counts.tsv](assets/example_data/correlation/ohsu-cnl_counts.tsv)

#### **Normalization**
[ohsu-cnl_counts.tsv](assets/example_data/normalization/ohsu-cnl_counts.tsv)&nbsp; | &nbsp;[ohsu-cnl_coldata.tsv](assets/example_data/normalization/ohsu-cnl_coldata.tsv)

### Gene count matrix (counts.tsv)

*Used for*: [Analysis](#differential-expression-analysis)&nbsp; | &nbsp;[Batch correction](#batch-correction)&nbsp; | &nbsp;[RNA-seq correlation](#rna-seq-correlation)&nbsp; | &nbsp;[Normalization](#normalization)

#### Required columns for count matrix
* *symbol* or *Hugo_symbol* column containing gene symbols
* A column for each sample
  * Column name sould be the sample name (if a coldata.tsv is provided, the sample names must match)
  * Column values should be the counts associated with each gene for that sample

Extra columns are okay (they will be ignored).

#### **Example**: RNA-Seq count matrix (counts.tsv)

Hugo_Symbol | sampleName1 | sampleName2 | sampleName3 | sampleName4 | ...
--- | --- | --- | --- | --- | ---
FAM208A | 2523.95	| 2703.47	| 2036.34	| 2204.85 | ...
RADIL | 320.0 | 183.0 | 209.0 | 577.0 | ...
AP1M2 | 34.0 | 8.0 | 24.0 | 2.0 | ...
... | ... | ... | ... | ... | ...

<div style="page-break-after: always;"></div>

#### **Example**: Microarray count matrix (counts.tsv)
symbol | sampleName1 | sampleName2 | sampleName3 | sampleName4 | ...
--- | --- | --- | --- | --- | ---
A1BG | 5.68557 | 6.30914 | 6.5588 | 5.67676 | ...
A1BG-AS1 | 6.35464 | 7.02546 | 7.01993 | 6.82097 | ...
A1CF | 6.0172 | 6.80624 | 5.87632 | 6.0752 | ...
... | ... | ... | ... | ... | ...

### Metadata for samples (coldata.tsv)

*Used for*: [Analysis](#differential-expression-analysis)&nbsp; | &nbsp;[Batch correction](#batch-correction)&nbsp; | &nbsp;[Normalization](#normalization)

#### Required columns for coldata.tsv
* *sample_name* containing the sample names (should match the sample names in the count matrix)
* *condition* containing the condition for each sample (e.g. healthy vs. disease)
* *batch* (only required for batch correction)

Extra columns are okay (they will be ignored).

#### **Example**: coldata.tsv file

sample_name | condition | batch *<span style="font-weight: normal">(only required for batch correction)</span>*
--- | --- | ---
sampleName1 | tumor | 1
sampleName2 | healthy | 2
sampleName3 | healthy | 1
sampleName4 | tumor | 2
... | ... | ...

### Filter file (filter.txt)
**Used for**: [Analysis](#differential-expression-analysis) (optional)
You can provide a text file with a list of genes to include in the "filtered output" of the analysis. The file should contain one gene symbol per line, for example:
  
```
FAM208A
RADIL
AP1M2
NEB
TEPP
CD151
AC009403.2
CAMP
```

### Configuration file (config.yml)
**Used for**: [Analysis](#differential-expression-analysis) (optional)  
A configuration file can be provided to specify the parameters for the analysis. The file should be in YAML format. If you do not provide a config.yml file, you can fill them in on the form after uploading your data files. See [here](#examples) to download an example config.yml file. The following parameters can be specified in config.yml:  
**min_expr**: Minimum expression threshold for a gene to be eligible for identification as differentially expressed  
**min_prop**: The minimum proportion of samples (for a gene) in which the expression of that gene must exceed min_expr. For example, if min_prop=1/3, and min_expr=1, then the expression of a gene must exceed 1 in at least 1/3 of samples to be eligible for identification as differentially expressed  
**padj_thresh**: The significance of threshold for the adjusted p-value of a gene’s expression — the expression change significance between conditions (e.g., disease vs. normal). Used in combination with log_fc to determine differential expression  
**log_fc**: Log fold-change. The fold-change (change in expression) in the expression of a gene between conditions (e.g., disease vs. normal). Used in combination with padj_thresh to determine differential expression  
**adj_method**: Method for adjusting expression differential expression p-values. See R documentation for a list of available methods.  
**condition**: Name of the column name for the factor of interest (e.g., “disease status”)  
**contrast_level**: The level (value) of condition for which effect size is being studied (e.g., “tumor”)  
**reference_level**: The level (value) of condition against which effect size is being compared (e.g., “healthy”)  
**use_qual_weights**: Used for microarray analysis. Whether to use quality weights for estimating differential expression. Recommended. See [documentation](https://rdrr.io/bioc/limma/man/arrayWeights.html).  

<div style="page-break-after: always;"></div>

## Differential expression analysis
This section is in not written yet...

<div style="page-break-after: always;"></div>

## Data preparation
### Batch correction
This section is in not written yet...

### Preprocessing
This section is in not written yet...

### Normalization
This section is in not written yet...

<div style="page-break-after: always;"></div>

## RNA-Seq correlation
This section is in not written yet...
