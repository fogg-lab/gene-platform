import { WebR } from 'webr';

let webR;

async function initializeWebR() {
  webR = new WebR();
  await webR.init();
  await webR.installPackages(['limma'], {
    repos: ['https://bioc.r-universe.dev', 'https://repo.r-wasm.org']
  });
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
          'geneSymbols',      // Array<string>
          'contrastGroup',   // Array<string>
          'referenceGroup',  // Array<string>
          'libSizes',        // Array<number>
          'numSamples',      // number
          'numGenes'         // number
        ];
        for (const varName of varNames) {
          await webR.objs.globalEnv.bind(varName, self[varName]);
        }

        // coldata is an object with keys cols, rows, and data
        // coldata.cols is an array of strings (column names) starting with 'sample_id'
        // coldata.data is an array of arrays containing [sample_id, ...values] for each sample
        const coldataDataframe = Object.fromEntries(
          coldata.cols.map(colName => [colName, coldata.data.map(row => row[coldata.cols.indexOf(colName)])])
        );
        await webR.objs.globalEnv.bind('coldata', coldataDataframe);
        await webR.objs.globalEnv.bind('coldata_cols', coldata.cols);

        // counts is an Int32Array
        // I found out the hard way that WebR specifically needs a Uint8 view when it comes to TypedArrays
        // (if given a Int32Array it will truncate every 32 bit int down to a single byte)
        const countsUint8Array = new Uint8Array(counts.buffer);
        await webR.objs.globalEnv.bind('counts', countsUint8Array);

        const _result = await (await webR.evalR(`
          library(limma)
          library(statmod)

          # Convert counts (initially a flat raw vector) to a numGenes x numSamples matrix
          counts <- readBin(counts, what = "integer", n = numGenes * numSamples)
          counts <- matrix(counts, nrow = numGenes, ncol = numSamples)
          rownames(counts) <- geneSymbols
          colnames(counts) <- coldata$sample_id

          # Filter lowly expressed genes
          cpm_counts <- counts / rep(libSizes, each = numGenes) * 1e6
          cutoff <- 1
          drop <- which(apply(cpm_counts, 1, max) < cutoff)
          if (length(drop) > 0) {
            counts <- counts[-drop, ]
          }
          # Create design matrix
          group <- ifelse(coldata$sample_id %in% contrastGroup, "contrastGrp", "referenceGrp")
          group <- factor(group, levels = c("referenceGrp", "contrastGrp"))
          design <- model.matrix(~0+group)
          colnames(design) <- levels(group)

          # Limma-voom
          y <- voom(counts, design, lib.size = libSizes)
          fit <- lmFit(y, design)
          contr <- makeContrasts(contrastGrp-referenceGrp, levels = colnames(coef(fit)))
          contr_fit <- contrasts.fit(fit, contr)
          contr_fit <- eBayes(contr_fit)
          list(res = topTable(contr_fit, sort.by = "P", adjust.method = "BH", n = Inf), row_names = rownames(counts))
        `)).toJs();
        const res = _result.values[0];
        result = {
          data: Float64Array.from(res.values.flatMap(col => col.values)),
          logFC: res.values[_result.values[0].names.indexOf('logFC')].values,
          t: res.values[_result.values[0].names.indexOf('t')].values,
          p_value: res.values[_result.values[0].names.indexOf('P.Value')].values,
          p_value_adj: res.values[_result.values[0].names.indexOf('adj.P.Val')].values,
          B: res.values[_result.values[0].names.indexOf('B')].values,
          row_names: _result.values[1].values,
          column_names: res.names
        }
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
