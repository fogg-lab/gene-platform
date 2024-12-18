# GENE Platform

A desktop application for gene expression analysis of bulk RNA-seq data, built with Electron and React. This project is a rewrite of the [original Flask-based web application](https://github.com/fogg-lab/gene-platform-archive).

Download and install the latest version on the [releases page](https://github.com/fogg-lab/gene-platform/releases).

## Features

- Exploratory Data Analysis
- Differential expression analysis for bulk RNA-seq data with a `limma-voom` pipeline
- Enrichment analysis via limma's `camera` method
- Interactive data visualization
- Load [harmonized GEO and GDC datasets](https://github.com/fogg-lab/curated-bulk-rnaseq-gene-expression)
- Local data processing via WebAssembly
