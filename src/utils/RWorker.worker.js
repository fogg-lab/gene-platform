let webR;

async function initializeWebR() {
  try {
    const { WebR } = await import(/* webpackIgnore: true */ 'https://webr.r-wasm.org/latest/webr.mjs');
    webR = new WebR({
      baseURL: 'https://webr.r-wasm.org/latest/'
    });

    await webR.init();
    console.log('WebR initialized successfully');

    await webR.installPackages(['limma'], {
      repos: ['https://bioc.r-universe.dev', 'https://repo.r-wasm.org']
    });
    console.log('Limma package installed successfully');
  } catch (error) {
    console.error('Failed to initialize WebR:', error);
    throw error;
  }
}

self.onmessage = async function(event) {
  if (!webR) {
    await initializeWebR();
  }

  const { action, data } = event.data;
  for (const key in data) {
    self[key] = data[key];
  }

  try {
    let result;
    switch (action) {
      case 'run_de_analysis':
        // Bind variables to the R environment
        const varNames = [
          'geneSymbols',
          'contrastGroup',
          'referenceGroup',
          'libSizes',
          'numSamples',
          'numGenes'
        ];
        for (const varName of varNames) {
          await webR.objs.globalEnv.bind(varName, self[varName]);
        }

        // coldata is an object with keys cols, rows, and data
        const coldataDataframe = Object.fromEntries(
          coldata.cols.map(colName => [colName, coldata.data.map(row => row[coldata.cols.indexOf(colName)])])
        );
        await webR.objs.globalEnv.bind('coldata', coldataDataframe);
        await webR.objs.globalEnv.bind('coldata_cols', coldata.cols);

        // Convert counts to Uint8Array
        const countsUint8Array = new Uint8Array(counts.buffer);
        await webR.objs.globalEnv.bind('counts', countsUint8Array);

        const _result = await (await webR.evalR(`
          library(limma)
          library(statmod)

          # Convert counts to matrix and set names
          counts <- readBin(counts, what = "integer", n = numGenes * numSamples)
          counts <- t(matrix(counts, nrow = numSamples, ncol = numGenes))
          rownames(counts) <- geneSymbols
          colnames(counts) <- coldata$sample_id

          # Filter lowly expressed genes
          cpm_counts <- sweep(counts, 2, libSizes, FUN="/") * 1e6
          min_count <- 10
          min_total_count <- 15
          group_sizes <- table(ifelse(coldata$sample_id %in% contrastGroup, "contrast", "reference"))
          min_samples <- min(group_sizes)
          keep <- rowSums(cpm_counts >= min_count) >= min_samples
          keep <- keep & rowSums(counts) >= min_total_count
          counts <- counts[keep, ]

          # Create design matrix
          group <- ifelse(coldata$sample_id %in% contrastGroup, "contrastGrp", "referenceGrp")
          group <- factor(group, levels = c("referenceGrp", "contrastGrp"))
          design <- model.matrix(~0+group)
          colnames(design) <- levels(group)

          # Differential expression analysis
          y <- voom(counts, design, lib.size = libSizes)
          fit <- lmFit(y, design)
          contr <- makeContrasts(contrastGrp - referenceGrp, levels = design)
          contr_fit <- contrasts.fit(fit, contr)
          contr_fit <- eBayes(contr_fit)

          # Get topTable results
          tbl <- topTable(contr_fit, sort.by = "P", adjust.method = "BH", n = Inf)
          final_res <- list(res = tbl, row_names = rownames(tbl))
          final_res
        `)).toJs();

        const res = _result.values[0];

        // Get column indexes directly from res.names
        const logFCIndex = res.names.indexOf('logFC');
        const tIndex = res.names.indexOf('t');
        const pIndex = res.names.indexOf('P.Value');
        const padjIndex = res.names.indexOf('adj.P.Val');
        const BIndex = res.names.indexOf('B');

        result = {
          data: Float64Array.from(res.values.flatMap(col => col.values)),
          logFC: res.values[logFCIndex].values,
          t: res.values[tIndex].values,
          p_value: res.values[pIndex].values,
          p_value_adj: res.values[padjIndex].values,
          B: res.values[BIndex].values,
          row_names: _result.values[1].values,
          column_names: res.names
        };
        break;

      case 'run_camera':
        // Bind variables to the R environment
        await webR.objs.globalEnv.bind('numSamples', design.length);
        await webR.objs.globalEnv.bind('numGenes', numGenes);
        await webR.objs.globalEnv.bind('geneSymbols', geneSymbols);
        await webR.objs.globalEnv.bind('geneSetNames', geneSetNames);
        await webR.objs.globalEnv.bind('design', design);
        await webR.objs.globalEnv.bind('contrast', contrast);
        await webR.objs.globalEnv.bind('minSetSize', data.minSetSize || 15);

        // Convert counts to Uint8Array first (quirk of webR)
        const countsUint8Arr = new Uint8Array(counts.buffer);
        await webR.objs.globalEnv.bind('counts', countsUint8Arr);

        // Flatten geneSetSymbols into a single array and record set sizes
        const geneSetSizes = geneSetSymbols.map(set => set.length);
        const flattenedSymbols = geneSetSymbols.flat();

        await webR.objs.globalEnv.bind('flattenedSymbols', flattenedSymbols);
        await webR.objs.globalEnv.bind('geneSetSizes', geneSetSizes);
        await webR.objs.globalEnv.bind('geneSetNames', geneSetNames);

        const _res = await (await webR.evalR(`
          library(limma)

          # Convert counts to matrix and set names
          counts <- readBin(counts, what = "integer", n = numGenes * numSamples)
          counts <- t(matrix(counts, nrow = numSamples, ncol = numGenes))
          rownames(counts) <- geneSymbols

          # Reconstruct the gene sets from flattenedSymbols and geneSetSizes
          index_list <- vector("list", length(geneSetSizes))
          pos <- 1
          for (i in seq_along(geneSetSizes)) {
            set_length <- geneSetSizes[i]
            set_symbols <- flattenedSymbols[pos:(pos + set_length - 1)]
            pos <- pos + set_length

            # Match the set symbols to our geneSymbols
            indices <- match(set_symbols, geneSymbols)
            original_size <- length(set_symbols)
            valid_indices <- indices[!is.na(indices)]
            
            # Only include sets that meet the minimum size requirement
            if (length(valid_indices) >= minSetSize) {
              index_list[[i]] <- structure(valid_indices, original_size = original_size)
            } else {
              index_list[[i]] <- NULL
            }
          }
          names(index_list) <- geneSetNames

          # Filter out empty sets and NULL entries
          valid_sets <- !sapply(index_list, is.null)
          index_list <- index_list[valid_sets]

          # Run camera
          camera_result <- camera(
            y = counts,
            index = index_list,
            design = design,
            contrast = contrast,
            inter.gene.cor = 0.01,
            use.ranks = TRUE
          )

          camera_df <- as.data.frame(camera_result)
          camera_df$Name <- rownames(camera_df)
          camera_df$NGenes <- sapply(index_list, function(x) attr(x, "original_size"))
          camera_df <- camera_df[order(camera_df$FDR), ]

          list(
            NGenes = as.vector(camera_df$NGenes),
            Direction = as.vector(camera_df$Direction),
            PValue = as.vector(camera_df$PValue),
            FDR = as.vector(camera_df$FDR),
            Name = as.vector(camera_df$Name)
          )
        `)).toJs();

        result = {
          NGenes: Array.from(_res.values[0].values),
          Direction: Array.from(_res.values[1].values),
          PValue: Array.from(_res.values[2].values),
          FDR: Array.from(_res.values[3].values),
          Name: Array.from(_res.values[4].values)
        };
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};

export default null;

