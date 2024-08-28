importScripts('https://webr.r-wasm.org/latest/webr.mjs');

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

  try {
    let result;
    switch (action) {
      case 'run_de_analysis':
        result = await webR.evalR(`
          library(limma)
          library(statmod)

          # Prepare data
          counts <- matrix(c(${data.counts}), nrow = ${data.geneCount}, ncol = ${data.sampleCount}, byrow = TRUE)
          rownames(counts) <- c(${data.geneNames.map(name => `"${name}"`).join(', ')}))
          colnames(counts) <- c(${data.sampleNames.map(name => `"${name}"`).join(', ')})

          # Create design matrix
          group <- factor(c(${data.groups.map(g => `"${g}"`).join(', ')}))
          design <- model.matrix(~0+group)
          colnames(design) <- levels(group)

          # Add covariates to the design matrix
          covariates <- data.frame(${data.covariates.map(cov => `${cov} = c(${data.coldataTable[cov].join(', ')})`).join(', ')}))
          design <- cbind(design, covariates)

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
            ${data.contrastName} = ${data.contrastGroup} - ${data.referenceGroup},
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
      // other R functions can be added here
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
