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
    await micropip.install([
      'https://github.com/fogg-lab/gene-platform-utils/releases/download/latest/gene_platform_utils-0.0.1-py3-none-any.whl'
    ])
    from gene_platform_utils import transformation, plot_de, plot_eda, plot_gsea
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
        case 'transform_vst':
          result = await pyodide.runPythonAsync(`
            import numpy as np
            from gene_platform_utils.transformation import vst
            counts = ${data.counts}
            transformed = vst(counts)
            transformed.tolist()
          `);
          break;
      case 'transform_log2':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.transformation import log2_1p
          counts = ${data.counts}
          transformed = log2_1p(counts)
          transformed.tolist()
        `);
        break;
      case 'transform_ln':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.transformation import ln_1p
          counts = ${data.counts}
          transformed = ln_1p(counts)
          transformed.tolist()
        `);
        break;
      case 'transform_log10':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.transformation import log10_1p
          counts = ${data.counts}
          transformed = log10_1p(counts)
          transformed.tolist()
        `);
        break;
      case 'create_heatmap':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.plot_eda import create_correlation_heatmap
          counts = ${data.counts}
          sample_ids = ${JSON.stringify(data.sample_ids)}
          create_correlation_heatmap(counts, sample_ids)
        `);
        break;
      case 'create_pca':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.plot_eda import create_pca_plot
          counts = ${data.counts}
          sample_ids = ${JSON.stringify(data.sample_ids)}
          create_pca_plot(counts, sample_ids)
        `);
        break;
      case 'create_tsne':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.plot_eda import create_tsne_plot
          counts = ${data.counts}
          sample_ids = ${JSON.stringify(data.sample_ids)}
          create_tsne_plot(counts, sample_ids)
        `);
        break;
      case 'create_volcano_plot':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.plot_de import create_volcano_plot
          data = ${data.data}
          row_names = ${JSON.stringify(data.row_names)}
          column_names = ${JSON.stringify(data.column_names)}
          lfc_thresh = ${data.lfc_thresh}
          pval_thresh = ${data.pval_thresh}
          cohort_name = ${JSON.stringify(data.cohort_name)}
          create_volcano_plot(data, row_names, column_names, lfc_thresh, pval_thresh, cohort_name)
        `);
        break;
      case 'create_mean_difference_plot':
        result = await pyodide.runPythonAsync(`
          import numpy as np
          from gene_platform_utils.plot_de import create_mean_difference_plot
          data = ${data.data}
          row_names = ${JSON.stringify(data.row_names)}
          column_names = ${JSON.stringify(data.column_names)}
          fdr = ${data.fdr}
          cohort_name = ${JSON.stringify(data.cohort_name)}
          create_mean_difference_plot(data, row_names, column_names, fdr, cohort_name)
        `);
        break;
      case 'create_gene_concept_network':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_gsea import gene_concept_network_plot
          gsea_res = ${JSON.stringify(data.gsea_res)}
          de_res = ${JSON.stringify(data.de_res)}
          ensembl_to_symbol = ${JSON.stringify(data.ensembl_to_symbol)}
          color_metric = ${JSON.stringify(data.color_metric)}
          pvalue_threshold = ${data.pvalue_threshold}
          layout_seed = ${data.layout_seed}
          color_seed = ${data.color_seed}
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
