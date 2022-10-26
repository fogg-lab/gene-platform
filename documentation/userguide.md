<!-- omit in toc -->
# GENE Platform Quickstart Guide

<!-- omit in toc -->
## Table of Contents
- [Preparing input files](#preparing-input-files)
  - [Gene count matrix (counts.tsv)](#gene-count-matrix-countstsv)
    - [Required columns for count matrix](#required-columns-for-count-matrix)
    - [**Example**: RNA-seq count matrix (counts.tsv)](#example-rna-seq-count-matrix-countstsv)
    - [**Example**: Microarray count matrix (counts.tsv)](#example-microarray-count-matrix-countstsv)
  - [Metadata for samples (coldata.tsv)](#metadata-for-samples-coldatatsv)
    - [**Example**: coldata.tsv file](#example-coldatatsv-file)
    - [Required columns for coldata.tsv](#required-columns-for-coldatatsv)
  - [Filter file (filter.txt)](#filter-file-filtertxt)
  - [Configuration file (config.yml)](#configuration-file-configyml)
- [Differential expression analysis](#differential-expression-analysis)
- [Data prep tools](#data-prep-tools)
  - [Batch correction](#batch-correction)
  - [Preprocessing](#preprocessing)
  - [Normalization](#normalization)
- [RNA-seq correlation](#rna-seq-correlation)

<div style="page-break-after: always;"></div>

## Preparing input files
### Gene count matrix (counts.tsv)

*Used for*: [Analysis](#differential-expression-analysis), [Batch correction](#batch-correction), [RNA-seq correlation](#rna-seq-correlation), and [Normalization](#normalization)

#### Required columns for count matrix
* *symbol* or *Hugo_symbol* column containing gene symbols
* A column for each sample
  * Column name sould be the sample name (if a coldata.tsv is provided, the sample names must match)
  * Column values should be the counts associated with each gene for that sample

Extra columns are okay (they will be ignored).

#### **Example**: RNA-seq count matrix (counts.tsv)

Hugo_Symbol | sampleName1 | sampleName2 | sampleName3 | sampleName4 | ...
--- | --- | --- | --- | --- | ---
FAM208A | 2523.95	| 2703.47	| 2036.34	| 2204.85 | ...
RADIL | 320.0 | 183.0 | 209.0 | 577.0 | ...
AP1M2 | 34.0 | 8.0 | 24.0 | 2.0 | ...
... | ... | ... | ... | ... | ...

#### **Example**: Microarray count matrix (counts.tsv)
symbol | sampleName1 | sampleName2 | sampleName3 | sampleName4 | ...
--- | --- | --- | --- | --- | ---
A1BG | 5.68557 | 6.30914 | 6.5588 | 5.67676 | ...
A1BG-AS1 | 6.35464 | 7.02546 | 7.01993 | 6.82097 | ...
A1CF | 6.0172 | 6.80624 | 5.87632 | 6.0752 | ...
... | ... | ... | ... | ... | ...

### Metadata for samples (coldata.tsv)

*Used for*: [Analysis](#differential-expression-analysis), [Batch correction](#batch-correction), and [Normalization](#normalization)

#### **Example**: coldata.tsv file

#### Required columns for coldata.tsv
* *sample_name* containing the sample names (should match the sample names in the count matrix)
* *condition* containing the condition for each sample (e.g. healthy vs. disease)
* *batch* (only required for batch correction)

Extra columns are okay (they will be ignored).

sample_name | condition | batch *<span style="font-weight: normal">(only required for batch correction)</span>*
--- | --- | ---
sampleName1 | tumor | 1
sampleName2 | healthy | 2
sampleName3 | healthy | 1
sampleName4 | tumor | 2
... | ... | ...

### Filter file (filter.txt)
**Used for**: [Analysis](#differential-expression-analysis) (optional)
To-do

### Configuration file (config.yml)
**Used for**: [Analysis](#differential-expression-analysis) (optional)
To-do

<div style="page-break-after: always;"></div>

## Differential expression analysis

<div style="page-break-after: always;"></div>

## Data prep tools
### Batch correction
To-do

### Preprocessing
To-do

### Normalization
To-do

<div style="page-break-after: always;"></div>

## RNA-seq correlation
To-do
