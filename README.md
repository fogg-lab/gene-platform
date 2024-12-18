# GENE Platform

A desktop application for gene expression analysis of bulk RNA-seq data, built with Electron and React. This project is a rewrite of the [original Flask-based web application](https://github.com/fogg-lab/gene-platform-archive).

## Features

- Exploratory Data Analysis
- Differential expression analysis for bulk RNA-seq data with a `limma-voom` pipeline
- Enrichment analysis via limma's `camera` method
- Interactive data visualization
- Load [harmonized GEO and GDC datasets](https://github.com/fogg-lab/curated-bulk-rnaseq-gene-expression)
- Local data processing via WebAssembly

## Setup

Prerequisite: [Install Node.js and npm](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm)

### Run from the command line

1. Clone the repository and navigate to the project root directory.
   ```bash
   git clone https://github.com/fogg-lab/gene-platform.git
   cd gene-platform
   ```
2. Install dependencies using `npm`.
   ```bash
   npm install
   ```
3. Run the app in development mode.
   ```bash
   npm run dev
   ```
