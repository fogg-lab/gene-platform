let wasm;

async function initializeWasm() {
  const wasmUrl = 'https://github.com/wigginno/gsea-rs/releases/download/latest/gsea_rs_bg.wasm';
  const jsUrl = 'https://github.com/wigginno/gsea-rs/releases/download/latest/gsea_rs.js';

  const wasmModule = await WebAssembly.instantiateStreaming(fetch(wasmUrl));
  const js = await import(jsUrl);
  wasm = js.default(wasmModule.instance.exports);
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
