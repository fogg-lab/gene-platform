import { getPublicUrl } from '../utils/environment';
import { loadPyodide } from 'pyodide';

let pyodide;

async function initializePyodide() {
  pyodide = await loadPyodide({ indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.2/full/' });
  await pyodide.loadPackage(['numpy', 'scipy', 'scikit-learn', 'micropip']);

  // Install gene-platform-utils
  await pyodide.runPythonAsync(`
    import micropip
    await micropip.install('networkx')
    await micropip.install('plotly')
    await micropip.install('${getPublicUrl()}/py-wheels/gene_platform_utils-0.0.1-py3-none-any.whl')
  `);
}

self.onmessage = async function(event) {
  if (!pyodide) {
    await initializePyodide();
  }

  const { action, data } = event.data;
  for (const key in data) {
    self[key] = data[key];
  }

  try {
    let result;
    switch (action) {
      case 'transform_vst':
        result = (await pyodide.runPythonAsync(`
          from js import expression, numSamples, numGenes
          import numpy as np
          from gene_platform_utils.transformation import vst
          counts = np.asarray(expression).reshape(numSamples, numGenes)
          vst(counts)
        `)).toJs();
        break;
      case 'transform_log2':
        result = (await pyodide.runPythonAsync(`
          from js import counts, numGenes, numSamples
          import numpy as np
          from gene_platform_utils.transformation import log2_1p
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          log2_1p(counts_2d)
        `)).toJs();
        break;
      case 'transform_ln':
        result = (await pyodide.runPythonAsync(`
          from js import counts, numGenes, numSamples
          import numpy as np
          from gene_platform_utils.transformation import ln_1p
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          ln_1p(counts_2d)
        `)).toJs();
        break;
      case 'transform_log10':
        result = (await pyodide.runPythonAsync(`
          from js import counts, numGenes, numSamples
          import numpy as np
          from gene_platform_utils.transformation import log10_1p
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          log10_1p(counts_2d)
        `)).toJs();
        break;
      case 'create_heatmap':
        result = await pyodide.runPythonAsync(`
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          from gene_platform_utils.plot_eda import create_correlation_heatmap
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          create_correlation_heatmap(counts_2d, sample_ids)
        `);
        break;
      case 'create_pca':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_eda import create_pca_plot
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          create_pca_plot(counts_2d, sample_ids)
        `);
        break;
      case 'create_tsne':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_eda import create_tsne_plot
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          counts_2d = np.asarray(counts).reshape(numGenes, numSamples)
          create_tsne_plot(counts_2d, sample_ids)
        `);
        break;
      case 'create_volcano_plot':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_de import create_volcano_plot
          from js import data, fdr, lfc_thresh, pval_thresh, cohort_name, row_names, column_names
          create_volcano_plot(data, row_names, column_names, lfc_thresh, pval_thresh, cohort_name)
        `);
        break;
      case 'create_mean_difference_plot':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_de import create_mean_difference_plot
          from js import data, cohort_name, fdr, row_names, column_names
          create_mean_difference_plot(data, row_names, column_names, fdr, cohort_name)
        `);
        break;
      case 'create_gene_concept_network':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_gsea import gene_concept_network_plot
          from js import gsea_res, de_res, ensembl_to_symbol, color_metric, pvalue_threshold, layout_seed, color_seed
          gene_concept_network_plot(gsea_res, de_res, ensembl_to_symbol, color_metric, pvalue_threshold, layout_seed, color_seed)
        `);
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};

export default null;
