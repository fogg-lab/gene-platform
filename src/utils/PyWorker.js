importScripts('https://cdn.jsdelivr.net/pyodide/v0.22.1/full/pyodide.js');

let pyodide;

async function initializePyodide() {
  pyodide = await loadPyodide({
    indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.22.1/full/'
  });
  await pyodide.loadPackage(['numpy', 'pandas', 'scipy', 'sklearn']);
  await pyodide.runPythonAsync(`
    import sys
    sys.path.append('gene-platform-utils/python')
    from gene_platform_utils import plot_eda, plot_de
  `);
}

self.onmessage = async function(event) {
  if (!pyodide) {
    await initializePyodide();
  }

  const { action, data } = event.data;

  try {
    let result;
    switch (action) {
      case 'transform_log':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          counts = np.array(${JSON.stringify(data.counts)})
          transformed = np.log2(counts + 1)
          transformed.tolist()
        `);
        break;
      case 'transform_vst':
        // Implement VST transformation
        break;
      case 'create_pca':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_eda import create_pca_plot
          create_pca_plot(${JSON.stringify(data.counts)}, ${JSON.stringify(data.sample_ids)})
        `);
        break;
      // Add more cases for other EDA functions
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
