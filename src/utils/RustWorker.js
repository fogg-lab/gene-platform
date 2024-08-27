importScripts('./pkg/gsea_rs.js');

let wasm;

async function initializeWasm() {
  wasm = await import('./pkg/gsea_rs.js');
  await wasm.default();
}

self.onmessage = async function(event) {
  if (!wasm) {
    await initializeWasm();
  }

  const { action, data } = event.data;

  try {
    let result;
    switch (action) {
      case 'run_gsea':
        result = wasm.prerank_rs(
          data.genes,
          data.metric,
          data.geneSets,
          data.weight,
          data.minSize,
          data.maxSize,
          data.nperm,
          BigInt(data.seed)
        );
        break;
    }
    self.postMessage({ status: 'success', result });
  } catch (error) {
    self.postMessage({ status: 'error', error: error.toString() });
  }
};
