import { isElectron } from '../utils/environment';
import { loadPyodide } from 'pyodide';

let pyodide;

async function initializePyodide() {
  pyodide = await loadPyodide({
    indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.4/full/'
  });
  await pyodide.loadPackage(['numpy', 'scipy', 'scikit-learn', 'micropip']);

  // Get the base URL for the wheel
  let publicUrl = isElectron() ? '' : `${self.location.origin}${getPublicUrl()}/static`;
  const wheelPath = `${publicUrl}/py-wheels/gene_platform_utils-0.0.1-py3-none-any.whl`;

  // Install gene-platform-utils
  await pyodide.runPythonAsync(`
    import micropip
    await micropip.install('networkx')
    await micropip.install('plotly')
    await micropip.install('${wheelPath}')
  `);
}

self.onmessage = async function (event) {
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
      case 'check_genes':
        const res = (await pyodide.runPythonAsync(`
          from js import genesQuery, humanGeneReference, mouseGeneReference
          import numpy as np
          from gene_platform_utils.genes import get_gene_info, get_duplicated
          hs_gene_result = get_gene_info(genesQuery, humanGeneReference)
          mm_gene_result = get_gene_info(genesQuery, mouseGeneReference)
          if mm_gene_result["num_matches"] > hs_gene_result["num_matches"]:
            species = "mouse"
            res = mm_gene_result
          else:
            species = "human"
            res = hs_gene_result
          gene_reference_idx = res["reference_indices"]
          not_found_gene_mask = gene_reference_idx == -1
          duplicate_genes_indices = get_duplicated(genesQuery)
          duplicate_gene_mask = np.zeros(len(gene_reference_idx), dtype=bool)
          for dup_indices in duplicate_genes_indices:
            duplicate_gene_mask[dup_indices] = True
          # Gather result codes:
          # 0: reference was found for the gene & the # of occurences of the gene in genesQuery == 1
          # 1: reference not found for the gene & the # of occurences of the gene in genesQuery == 1
          # 2: reference was found for the gene & the # of occurences of the gene in genesQuery >= 2
          # 3: reference not found for the gene & the # of occurences of the gene in genesQuery >= 2
          query_result_code = np.zeros(len(genesQuery), dtype=np.int32)
          query_result_code[not_found_gene_mask] += 1
          query_result_code[duplicate_gene_mask] += 2
          res["queryResultCode"] = query_result_code.tolist()
          res["duplicateGenesIndices"] = duplicate_genes_indices
          res["detectedSpecies"] = species
          res
        `)).toJs();
        result = {};
        result.queryResultGeneIndices = res.reference_indices;
        result.queryResultCode = res.queryResultCode;
        result.duplicateGenesIndices = res.duplicateGenesIndices;
        result.detectedSpecies = res.detectedSpecies;
        result.numFound = res.num_matches;
        break;
      case 'transform_vst':
        result = (await pyodide.runPythonAsync(`
          from js import expression, numSamples, numGenes
          import numpy as np
          from gene_platform_utils.transformation import vst
          counts_2d = np.asarray(expression).reshape(numSamples, numGenes)
          vst(counts_2d).astype(np.float32).flatten()
        `)).toJs();
        break;
      case 'transform_log2':
        result = (await pyodide.runPythonAsync(`
          from js import expression
          import numpy as np
          from gene_platform_utils.transformation import log2_1p
          log2_1p(np.asarray(expression)).astype(np.float32)
        `)).toJs();
        break;
      case 'transform_ln':
        result = (await pyodide.runPythonAsync(`
          from js import expression
          import numpy as np
          from gene_platform_utils.transformation import ln_1p
          ln_1p(np.asarray(expression)).astype(np.float32)
        `)).toJs();
        break;
      case 'transform_log10':
        result = (await pyodide.runPythonAsync(`
          from js import expression, numGenes, numSamples
          import numpy as np
          from gene_platform_utils.transformation import log10_1p
          log10_1p(np.asarray(expression)).astype(np.float32)
        `)).toJs();
        break;
      case 'compute_tmm_effective_library_sizes':
        result = (await pyodide.runPythonAsync(`
          from js import expression, numGenes, numSamples
          import numpy as np
          from gene_platform_utils.between_sample_norm import compute_tmm_effective_library_sizes
          counts_2d = np.asarray(expression).reshape(numSamples, numGenes)
          compute_tmm_effective_library_sizes(counts_2d)
        `)).toJs();
        break;
      case 'create_heatmap':
        result = await pyodide.runPythonAsync(`
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          from gene_platform_utils.plot_eda import create_correlation_heatmap
          counts_2d = np.asarray(counts).reshape(numSamples, numGenes).T
          create_correlation_heatmap(counts_2d, sample_ids)
        `);
        break;
      case 'create_pca':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_eda import create_pca_plot
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          counts_2d = np.asarray(counts).reshape(numSamples, numGenes).T
          create_pca_plot(counts_2d, sample_ids)
        `);
        break;
      case 'create_tsne':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_eda import create_tsne_plot
          from js import counts, numGenes, numSamples, sample_ids
          import numpy as np
          counts_2d = np.asarray(counts).reshape(numSamples, numGenes).T
          create_tsne_plot(counts_2d, sample_ids)
        `);
        break;
      case 'create_volcano_plot':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_de import create_volcano_plot
          from js import data, lfc_thresh, pval_thresh, cohort_name, row_names, column_names
          data_2d = np.asarray(data).reshape(len(column_names), len(row_names)).T
          create_volcano_plot(data_2d, row_names, column_names, lfc_thresh, pval_thresh, cohort_name)
        `);
        break;
      case 'create_mean_difference_plot':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_de import create_mean_difference_plot
          from js import data, cohort_name, fdr, row_names, column_names
          data_2d = np.asarray(data).reshape(len(column_names), len(row_names)).T
          create_mean_difference_plot(data_2d, row_names, column_names, fdr, cohort_name)
        `);
        break;
      case 'create_gene_concept_network':
        result = await pyodide.runPythonAsync(`
          from gene_platform_utils.plot_gsea import gene_concept_network_plot
          from js import gsea_res, de_res, color_metric, pvalue_threshold, layout_seed, color_seed
          de_res_2d = np.asarray(de_res).reshape(len(column_names), len(row_names)).T
          gene_concept_network_plot(gsea_res, de_res_2d, color_metric, pvalue_threshold, layout_seed, color_seed)
        `);
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};

export default null;
