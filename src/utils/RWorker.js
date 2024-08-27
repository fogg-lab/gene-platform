// RWorker.js
importScripts('https://webr.r-wasm.org/latest/webr.mjs');

let webR;

async function initializeWebR() {
  webR = new WebR();
  await webR.init();
  await webR.installPackages(['limma'], {
    repos: 'https://bioc.r-universe.dev'
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
          # Implement DE analysis using limma
          # Return results as a data frame
        `);
        break;
      // Add more cases for other R functions
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
