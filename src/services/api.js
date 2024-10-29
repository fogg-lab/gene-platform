import pako from "pako";
import Papa from "papaparse";
const { loadPyodide } = require("pyodide");
import { strReplaceAll } from "../utils/helper";

/**
 * Retrieve a dataset from an external database (GDC, GEO-human, or GEO-mouse).
 * @param {string} dataSrc - The data source ('GDC', 'GEO-human', or 'GEO-mouse').
 * @param {string} datasetID - The identifier of the dataset (GDC project ID or GEO accession number).
 * @returns {Promise<Object>} - A promise that resolves to an object containing the dataset with the following structure:
 *   {
 *     expression: Int32Array, // 1D array representing the expression matrix (samples x genes)
 *     counts: Int32Array, // 1D array representing the transposed expression matrix (genes x samples)
 *     countsTable: {
 *       cols: string[], // Array of column names, including 'Ensembl gene', 'Symbol', and sample IDs
 *       rows: string[], // Array of Ensembl gene IDs
 *       data: number[][] // 2D array of counts values
 *     },
 *     coldataTable: {
 *       cols: string[], // Array of column names from the coldata CSV (bio)
 *       rows: string[], // Array of sample IDs
 *       data: string[][] // 2D array of metadata values
 *     },
 *     genesTable: {
 *       cols: string[], // Array of column names from the genes CSV
 *       rows: string[], // Array of gene IDs (ensembl_id)
 *       data: string[][] // 2D array of gene information
 *     }
 *   }
 */
export async function getExternalDataset(dataSrc, datasetID) {
    const baseUrl = `/datasets/curated-bulk-rnaseq-gene-expression/${dataSrc}`;

    // Fetch and parse coldata
    const coldataResponse = await fetch(`${baseUrl}/coldata/${datasetID}.csv.gz`);
    const coldataArrayBuffer = await coldataResponse.arrayBuffer();
    const coldataUnzipped = pako.ungzip(new Uint8Array(coldataArrayBuffer));
    const coldataText = new TextDecoder().decode(coldataUnzipped);
    let coldataData = Papa.parse(coldataText, { header: true }).data;
    if (coldataData.at(-1).sample_id === "") {
        coldataData = coldataData.slice(0, -1);
    }

    // Fetch and parse genes
    const genesResponse = await fetch(`${baseUrl}/genes/${datasetID}.csv.gz`);
    const genesArrayBuffer = await genesResponse.arrayBuffer();
    const genesUnzipped = pako.ungzip(new Uint8Array(genesArrayBuffer));
    const genesText = new TextDecoder().decode(genesUnzipped);
    let genesData = Papa.parse(genesText, { header: true }).data;
    if (genesData.at(-1).ensembl_gene === "") {
        genesData = genesData.slice(0, -1);
    }

    // Fetch expression data
    const expressionResponse = await fetch(
        `${baseUrl}/expression/${datasetID}.npy.gz`
    );
    const expressionArrayBuffer = await expressionResponse.arrayBuffer();
    const expressionUnzipped = pako.ungzip(new Uint8Array(expressionArrayBuffer));

    // Use pyodide to load and process the NPY data
    let pyodide = await loadPyodide({
        indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/",
    });
    await pyodide.loadPackage("numpy");

    const np = pyodide.pyimport("numpy");
    const io = pyodide.pyimport("io");

    const expressionBuffer = io.BytesIO(pyodide.toPy(expressionUnzipped));
    let expressionArray = np.load(expressionBuffer);
    const expression = new Int32Array(expressionArray.flatten().toJs());
    const counts = new Int32Array(np.transpose(expressionArray).flatten().toJs());

    // Create coldataTable
    const coldataCols = Object.keys(coldataData[0]);
    const coldataTable = {
        cols: coldataCols,
        rows: coldataData.map((row) => row.sample_id),
        data: coldataData.map((row) => coldataCols.map((col) => row[col])),
    };

    // Create genesTable
    const genesCols = Object.keys(genesData[0]);
    const genesTable = {
        cols: genesCols,
        rows: genesData.map((row) => row.ensembl_gene),
        data: genesData.map((row) => genesCols.map((col) => row[col])),
    };

    // Create countsTable
    const countsTable = {
        cols: ["symbol", ...coldataTable.rows],
        rows: genesTable.rows,
        data: (() => {
            const numSamples = coldataTable.rows.length;
            const numGenes = genesTable.rows.length;
            return Array.from({ length: numGenes }, (_, geneIndex) => {
                const hgncSymbol = genesTable.data[geneIndex][2];
                const geneCountsStart = geneIndex * numSamples;
                const geneCounts = counts.slice(
                    geneCountsStart,
                    geneCountsStart + numSamples
                );
                return [hgncSymbol, ...Array.from(geneCounts)];
            });
        })(),
    };

    return {
        expression,
        counts,
        countsTable,
        coldataTable,
        genesTable,
    };
}

/**
 * Retrieve a collection of gene sets from MSigDB.
 * @param {string} species - The species of the gene sets to retrieve ('human' or 'mouse').
 * @param {string} collection - The collection of gene sets to retrieve (e.g. 'H', 'MH', 'C2:CP:BIOCARTA', etc.).
 * @returns {Promise<Object>} - Promise for the gene sets in the collection and their symbols, for example:
 *   {
 *     geneSets: ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", ...],
 *     geneSetSymbols: [["ABCA1", "ACKR3", "AREG", ...], ["ACKR3", "ADM", "ADORA2B", ...], ...]
 *   }
 */
export async function getExternalGeneSetCollection(species, collection) {
    const baseUrl = "/broadinstitute/gsea-msigdb/msigdb/release";
    const speciesCode = species === "mouse" ? "Mm" : "Hs";
    const collectionCode = strReplaceAll(collection, ":", ".").toLowerCase();
    const url = `${baseUrl}/2024.1.${speciesCode}/${collectionCode}.v2024.1.${speciesCode}.json`;
    const response = await fetch(url);
    const data = await response.json();
    // data looks like this:
    // {
    //   "HALLMARK_MITOTIC_SPINDLE":{"collection":"H","systematicName":"M5893",...,"geneSymbols":["ABI1","ABL1",...]},
    //   "HALLMARK_WNT_BETA_CATENIN_SIGNALING":{"collection":"H","systematicName":"M5895",...,"geneSymbols":["ADAM17","AXIN1",...]},
    //   ...
    // }
    // we want to return an object with the following structure:
    // {
    //   geneSets: ["HALLMARK_MITOTIC_SPINDLE", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", ...],
    //   geneSetSymbols: [["ABI1","ABL1",...], ["ADAM17","AXIN1",...], ...]
    // }
    const geneSets = Object.keys(data);
    return {
        geneSets: geneSets,
        geneSetSymbols: geneSets.map((geneSet) => data[geneSet].geneSymbols),
    };
}

/**
 * Retrieve a gene set from MSigDB.
 * @param {string} species - The species of the gene set to retrieve ('human' or 'mouse').
 * @param {string} geneSetName - The name of the gene set to retrieve (e.g. 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', etc.).
 * @returns {Promise<Object>} - Promise for the gene set with the following structure:
 *   {
 *     geneSymbols: ["ABI1", "ABL1", ...]
 *   }
 */
export async function getExternalGeneSet(species, geneSetName) {
    const baseUrl = "/msigdb/gsea/msigdb";
    const url = `${baseUrl}/${species}/download_geneset.jsp?geneSetName=${geneSetName}&fileType=json`;
    const response = await fetch(url);
    const data = await response.json();
    return {
        geneSymbols: data[geneSetName].geneSymbols,
    };
}
