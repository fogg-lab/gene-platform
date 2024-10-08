import pako from "pako";
import Papa from "papaparse";
const { loadPyodide } = require("pyodide");

/**
 * Retrieve a dataset from an external database (GEO or GDC).
 * @param {string} dataSrc - The data source ('GEO' or 'GDC').
 * @param {string} datasetID - The identifier of the dataset (GEO accession number or GDC project ID).
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
    const baseUrl = `/api/curated-bulk-rnaseq-human-gene-expression/${dataSrc}`;

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
    const expressionResponse = await fetch(`${baseUrl}/expression/${datasetID}.npy.gz`);
    const expressionArrayBuffer = await expressionResponse.arrayBuffer();
    const expressionUnzipped = pako.ungzip(new Uint8Array(expressionArrayBuffer));

    // Use pyodide to load and process the NPY data
    let pyodide = await loadPyodide({ indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.2/full/' });
    await pyodide.loadPackage('numpy');

    const np = pyodide.pyimport('numpy');
    const io = pyodide.pyimport('io');

    const expressionBuffer = io.BytesIO(pyodide.toPy(expressionUnzipped));
    let expressionArray = np.load(expressionBuffer);
    const expression = new Int32Array(expressionArray.flatten().toJs());
    const counts = new Int32Array(np.transpose(expressionArray).flatten().toJs());

    // Create coldataTable
    const coldataCols = Object.keys(coldataData[0]);
    const coldataTable = {
        cols: coldataCols,
        rows: coldataData.map(row => row.sample_id),
        data: coldataData.map(row => coldataCols.map(col => row[col]))
    };

    // Create genesTable
    const genesCols = Object.keys(genesData[0]);
    const genesTable = {
        cols: genesCols,
        rows: genesData.map(row => row.ensembl_gene),
        data: genesData.map(row => genesCols.map(col => row[col]))
    };

    // Create countsTable
    const countsTable = {
        cols: ['Ensembl gene', 'Symbol', ...coldataTable.rows],
        rows: genesTable.rows,
        data: (() => {
            const numSamples = coldataTable.rows.length;
            const numGenes = genesTable.rows.length;
            return Array.from({ length: numGenes }, (_, geneIndex) => {
                const ensemblGene = genesTable.rows[geneIndex];
                const symbol = genesTable.data[geneIndex][1];
                const geneCountsStart = geneIndex * numSamples;
                const geneCounts = counts.slice(geneCountsStart, geneCountsStart + numSamples);
                return [ensemblGene, symbol, ...Array.from(geneCounts)];
            });
        })()
    };

    return {
        expression,
        counts,
        countsTable,
        coldataTable,
        genesTable,
    };
}
