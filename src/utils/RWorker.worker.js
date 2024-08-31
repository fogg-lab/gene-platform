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
        const varNames = ['counts', 'contrastGroup', 'referenceGroup', 'numSamples', 'numGenes'];
        for (const varName of varNames) {
          // counts is a TypedArray, specifically an Int32Array
          // coldata is an object with keys col, rows, and data
          // coldata.col is an array of strings (column names) starting with 'sample_id'
          // coldata.data is an array of arrays containing [sample_id, ...values] for each sample
          // contrastGroup is an array of strings (sample ids)
          // referenceGroup is an array of strings (sample ids)
          // numSamples is an integer
          // numGenes is an integer
          await webR.objs.globalEnv.bind(varName, self[varName]);
        }
        const coldataDataframe = Object.fromEntries(
          coldata.cols.map(colName => [colName, coldata.data.map(row => row[coldata.cols.indexOf(colName)])])
        );
        await webR.objs.globalEnv.bind('coldata', coldataDataframe);
        await webR.objs.globalEnv.bind('coldata_cols', coldata.cols);
        // await webR.objs.globalEnv.bind('coldata_data', coldata.data); // OLD. Need to rewrite below to not use coldata_data

        result = await webR.evalR(`
          library(limma)
          library(statmod)

          # Convert counts (initially a 1D raw vector) to a numGenes x numSamples matrix
          counts <- matrix(counts, nrow = numGenes, ncol = numSamples)
          rownames(counts) <- paste0("gene", 1:numGenes)
          colnames(counts) <- coldata$sample_id

          # Create design matrix
          group <- factor(c(rep("contrast", length(contrastGroup)), rep("reference", length(referenceGroup))))
          design <- model.matrix(~0+group)
          colnames(design) <- levels(group)

          # Add covariates to the design matrix
          # TODO: Add ability to define covariates in UI that will be added here

          # Normalize and transform data
          dge <- DGEList(counts = counts)
          keep <- filterByExpr(dge, design)
          dge <- dge[keep, , keep.lib.sizes=FALSE]
          dge <- calcNormFactors(dge)
          v <- voom(dge, design, plot=FALSE)

          # Fit the model
          fit <- lmFit(v, design)

          # Define contrasts
          contrasts <- makeContrasts(
            contrast = contrastGroup - referenceGroup,
            levels = colnames(design)
          )

          # Compute contrasts
          fit2 <- contrasts.fit(fit, contrasts)
          fit2 <- eBayes(fit2)

          # Get results
          results <- topTable(fit2, coef=1, number=Inf, sort.by="P")
          results_list <- lapply(1:nrow(results), function(i) as.list(results[i,]))
          results_list
        `);
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
