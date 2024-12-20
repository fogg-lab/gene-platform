# GENE Platform

**Gene Expression Explorer**

A desktop application for gene expression analysis of bulk RNA-seq data, built with Electron and React. This project is a rewrite of the [original Flask-based web application](https://github.com/fogg-lab/gene-platform-archive).

## Features

- Exploratory Data Analysis
- Differential expression analysis for bulk RNA-seq data with a `limma-voom` pipeline
- Enrichment analysis via limma's `camera` method
- Interactive data visualization
- Load [harmonized GEO and GDC datasets](https://github.com/fogg-lab/curated-bulk-rnaseq-gene-expression)
- Local data processing via WebAssembly

## Setup

Set up and launch GENE on your machine.

### Prerequisites
- Install Git. macOS users can run `git -v` in Terminal for automatic installation. On Windows, download the installer from [git-scm.com](https://git-scm.com) and perform a standard installation.
- Install Node.js. Download and run the installer from [nodejs.org](https://nodejs.org). You may perform a regular installation using the default options. No extra dependencies are necessary.

### Install and Launch GENE
Set up and launch the app from a new command line session. Launch Command Prompt on Windows or Terminal on macOS, and run the four commands below in sequence.

1. Clone this repository.
   ```bash
   git clone https://github.com/fogg-lab/gene-platform.git
   ```
2. Navigate to the project root directory.
   ```bash
   cd gene-platform
   ```
3. Install project dependencies.
   ```bash
   npm install
   ```
4. Launch the app in development mode.
   ```bash
   npm run dev
   ```

## Quickstart demo

Follow [this link](https://cdnapisec.kaltura.com/p/391241/embedPlaykitJs/uiconf_id/44855082?iframeembed=true&entry_id=1_xk6xve4k&config%5Bprovider%5D=%7B%22widgetId%22%3A%221_7iju7bzc%22%7D&config%5Bplayback%5D=%7B%22startTime%22%3A0%7D) to watch a 5-minute demo.
