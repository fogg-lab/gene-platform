import init, { prerank_rs } from "gsea-rs";

//let wasm;
//async function initializeWasm() {
//  const publicUrl = getPublicUrl();
//  const wasmUrl = `${publicUrl}/wasm/gsea_rs_bg.wasm`;
//  const jsUrl = `${publicUrl}/wasm/gsea_rs.js`;
//  // Load the JS file as a module
//  const module = await import(/* webpackIgnore: true */ jsUrl);
//  // Initialize the WASM module
//  const wasmResponse = await fetch(wasmUrl);
//  const wasmBuffer = await wasmResponse.arrayBuffer();
//  wasm = await module.default(wasmBuffer);
//}

self.onmessage = async function (event) {
  // if (!wasm) {
  //   await init();
  // }
  await init();

  const { action, data } = event.data;

  try {
    let result;
    switch (action) {
      case "run_gsea":
        result = prerank_rs(
            data.genes,
            data.metric,
            data.geneSets,
            data.weight,
            data.minSize,
            data.maxSize,
            data.nperm,
            data.seed.toString()
          );
        break;
    }
    self.postMessage({ status: "success", result });
  } catch (error) {
    self.postMessage({ status: "error", error: error.toString() });
  }
};
