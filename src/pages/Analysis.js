import React, { useState, useEffect, useCallback, useRef } from 'react';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
import GeneSetCollectionsPopup from '../components/ui/GeneSetCollectionsPopup';
import WorkerManager from '../utils/Workers';
import { getExternalDataset, getExternalGeneSetCollection } from '../services/api';
import humanGenes from '../assets/genes/human_genes.json';
import mouseGenes from '../assets/genes/mouse_genes.json';
import DifferentialExpressionContent from '../components/ui/DifferentialExpressionContent';
import ExplorationContent from '../components/ui/ExplorationContent';
import GSEAContent from '../components/ui/GSEAContent';
import GSEAInputForm from '../components/form/GSEAInputForm';
import DEAInputForm from '../components/form/DEAInputForm';
import EDAInputForm from '../components/form/EDAInputForm';

const Analysis = () => {
    const [activeTab, setActiveTab] = useState('table');
    const [error, setError] = useState(null);
    const [isVisible, setIsVisible] = useState(false);
    const [dataset, setDataset] = useState(null);
    const [edaData, setEdaData] = useState(null);
    const [deData, setDeData] = useState(null);
    const [gseaData, setGseaData] = useState(null);
    const [currentTable, setCurrentTable] = useState(null);
    const [currentPlot, setCurrentPlot] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [tableColumns, setTableColumns] = useState([]);

    const [contrastGroup, setContrastGroup] = useState({ samples: [] });
    const [referenceGroup, setReferenceGroup] = useState({ samples: [] });

    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);

    const [geneSetCollections, setGeneSetCollections] = useState([]);
    const [gseaParams, setGseaParams] = useState({
        weight: 1,
        minSize: 15,
        maxSize: 500,
        nperm: 1000,
    });

    const [currentStage, setCurrentStage] = useState('exploration');
    const [lockedStates, setLockedStates] = useState({
        exploration: false,
        differential: true,
        enrichment: true
    });

    const updateLockedState = (stage, isLocked) => {
        setLockedStates(prevStates => ({
            ...prevStates,
            [stage]: isLocked
        }));
    };

    const handleAddSamplesToGroup = useCallback((isContrast, samplesToAdd) => {
        setContrastGroup(prevContrastGroup => {
            const updatedContrastGroup = isContrast
                ? [...new Set([...prevContrastGroup.samples, ...samplesToAdd.filter(sample =>
                    !prevContrastGroup.samples.some(s => s.id === sample.id)
                )])]
                : prevContrastGroup.samples.filter(sample => !samplesToAdd.some(s => s.id === sample.id));
            return { ...prevContrastGroup, samples: updatedContrastGroup };
        });
        setReferenceGroup(prevReferenceGroup => {
            const updatedReferenceGroup = !isContrast
                ? [...new Set([...prevReferenceGroup.samples, ...samplesToAdd.filter(sample =>
                    !prevReferenceGroup.samples.some(s => s.id === sample.id)
                )])]
                : prevReferenceGroup.samples.filter(sample => !samplesToAdd.some(s => s.id === sample.id));
            return { ...prevReferenceGroup, samples: updatedReferenceGroup };
        });
    }, []);

    const handleRemoveSamplesFromGroup = useCallback((isContrast, sampleIdsToRemove) => {
        if (isContrast) {
            setContrastGroup(prevGroup => ({
                ...prevGroup,
                samples: prevGroup.samples.filter(sample => !sampleIdsToRemove.includes(sample.id))
            }));
        } else {
            setReferenceGroup(prevGroup => ({
                ...prevGroup,
                samples: prevGroup.samples.filter(sample => !sampleIdsToRemove.includes(sample.id))
            }));
        }

        // Update row colors
        setDataset(prevDataset => {
            if (!prevDataset || !prevDataset.coldataTable) return prevDataset;

            const updatedColdataTable = {
                ...prevDataset.coldataTable,
                data: prevDataset.coldataTable.data.map(row => {
                    const sampleId = row[prevDataset.coldataTable.cols.indexOf('sample_id')];
                    if (sampleIdsToRemove.includes(sampleId)) {
                        return row.map((cell, index) =>
                            index === prevDataset.coldataTable.cols.indexOf('group') ? '' : cell
                        );
                    }
                    return row;
                })
            };

            return {
                ...prevDataset,
                coldataTable: updatedColdataTable
            };
        });
    }, []);

    const handleClearGroup = useCallback((isContrast) => {
        if (isContrast) {
            setContrastGroup({ samples: [] });
        } else {
            setReferenceGroup({ samples: [] });
        }
    }, []);

    useEffect(() => {
        if (dataset) {
            setEdaData({
                tables: {
                    coldata: dataset.coldataTable,
                    counts: dataset.countsTable
                },
                plots: {}
            });
            setCurrentTable('coldata');
        }
    }, [dataset]);

    const handleDatasetSelect = async (type, data) => {
        let dset;
        if (type === 'external' || type === 'upload') {
            setDataset(data);
            dset = data;
        } else if (type === 'example') {
            setIsLoading(true);
            try {
                const exampleData = await getExternalDataset('GDC', 'CDDP_EAGLE-1');
                dset = exampleData;
                setDataset(exampleData);

                handleClearGroup(true);
                handleClearGroup(false);

                setReferenceGroup(
                    {
                        samples: exampleData.coldataTable.data
                            .filter(row => row[exampleData.coldataTable.cols.indexOf('pack_years_smoked')] === '0.0')
                            .map(row => ({
                                id: row[exampleData.coldataTable.cols.indexOf('sample_id')],
                                [exampleData.coldataTable.cols[0]]: row[0]
                            }))
                    }
                );

                setContrastGroup(
                    {
                        samples: exampleData.coldataTable.data
                            .filter(row => {
                                const packYears = parseFloat(row[exampleData.coldataTable.cols.indexOf('pack_years_smoked')]);
                                return !isNaN(packYears) && packYears >= 50;
                            })
                            .map(row => ({
                                id: row[exampleData.coldataTable.cols.indexOf('sample_id')],
                                [exampleData.coldataTable.cols[0]]: row[0]
                            }))
                    }
                );
                const collectionSpecies = "human";
                const collectionId = "C5:GO:BP";
                //const collectionId = "C2:CP:KEGG_LEGACY";
                const newCollection = await getExternalGeneSetCollection(collectionSpecies, collectionId);
                setGeneSetCollections(prev => [...prev, { name: collectionId, data: newCollection }]);
            } catch (error) {
                console.error('Error loading example dataset:', error);
                setError('Failed to load example dataset');
            } finally {
                setIsLoading(false);
            }
        }

        if (dset && dset.coldataTable) {
            const tableData = dset.coldataTable.data.map((row, index) => {
                const obj = {};
                dset.coldataTable.cols.forEach((col, colIndex) => {
                    obj[col] = row[colIndex];
                });
                obj.id = index;
                return obj;
            });

            const columns = dset.coldataTable.cols.map(col => ({ key: col, name: col }));

            setTableData(tableData);
            setTableColumns(columns);
        }
    };

    const handleStageChange = (stage) => {
        setCurrentStage(stage);
        switch (stage) {
            case 'exploration':
                setCurrentTable('coldata');
                setCurrentPlot(null);
                break;
            case 'differential':
                if (deData) {
                    setCurrentTable('de_results');
                    setCurrentPlot('volcano_plot');
                }
                break;
            case 'enrichment':
                if (gseaData) {
                    setCurrentTable('gsea_results');
                    setCurrentPlot('gene_concept_network');
                }
                break;
        }
    };

    const runAnalysis = async (processedData) => {
        setIsLoading(true);
        setProgress(0);

        // If processedData is provided (from EDAInputForm), use it instead of the current dataset
        if (processedData) {
            setDataset(processedData);
        }

        const numInputGenes = (processedData?.countsTable || dataset.countsTable).rows.length;
        const numSamples = (processedData?.coldataTable || dataset.coldataTable).rows.length;
        const geneInfo = await WorkerManager.runTask('py', 'check_genes', {
            genesQuery: (processedData?.countsTable || dataset.countsTable).rows,
            humanGeneReference: humanGenes.genes,
            mouseGeneReference: mouseGenes.genes
        });
        const geneReference = geneInfo.detectedSpecies === "mouse" ? mouseGenes.genes : humanGenes.genes;

        // Create filtered expression and countsTable that excludes duplicate and unmatched genes and uses HGNC symbol
        // Genes that are unique and were mapped to known HGNC symbols using the gene reference have result code 0
        // (1-3 are failure codes. 1: reference for gene not found, 2: gene not unique, 3: not found and not unique)
        const geneIndex = Array.from({ length: numInputGenes }, (_, i) => i);
        const geneMask = geneIndex.map((i) => geneInfo.queryResultCode[i] == 0);
        const mappedGeneIdx = geneIndex.filter((i) => geneMask[i]);
        const mappedGeneSymbols = mappedGeneIdx.map((i) => geneReference[geneInfo.queryResultGeneIndices[i]][2]);
        const numGenes = mappedGeneSymbols.length;
        const mappedGenesExpression = new (processedData?.expression || dataset.expression).constructor(numGenes * numSamples);
        let mappedGenesExpressionIndex = 0;
        let colIndex = 0;
        for (let i = 0; i < (processedData?.expression || dataset.expression).length; i++) {
            if (geneMask[colIndex++]) {
                mappedGenesExpression[mappedGenesExpressionIndex++] = (processedData?.expression || dataset.expression)[i];
            }
            if (colIndex === numInputGenes) colIndex = 0;
        }
        const mappedGenesCountsTable = {
            cols: ["Symbol", ...(processedData?.countsTable || dataset.countsTable).cols],
            rows: mappedGeneSymbols,
            data: mappedGeneIdx.map((i) => [mappedGeneSymbols[i], ...(processedData?.countsTable || dataset.countsTable).data[i]])
        }
        setProgress(20);

        // EDA
        if (currentStage === 'exploration') {
            try {
                const transformedCounts = await WorkerManager.runTask('py', 'transform_vst', {
                    expression: mappedGenesExpression,
                    numSamples: numSamples,
                    numGenes: numGenes,
                });
                setProgress(40);
                const pcaPlot = await WorkerManager.runTask('py', 'create_pca', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: mappedGenesCountsTable.cols.slice(1)
                });
                setProgress(60);
                const tsnePlot = await WorkerManager.runTask('py', 'create_tsne', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: mappedGenesCountsTable.cols.slice(1)
                });
                setProgress(80);
                const heatmap = await WorkerManager.runTask('py', 'create_heatmap', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: mappedGenesCountsTable.cols.slice(1)
                });

                setEdaData({
                    tables: {
                        coldata: (processedData?.coldataTable || dataset.coldataTable),
                        counts: mappedGenesCountsTable
                    },
                    plots: { pca: pcaPlot, tsne: tsnePlot, heatmap: heatmap }
                });
                setCurrentTable('coldata');
                updateLockedState('differential', false);
            } catch (error) {
                console.error("Error in exploration analysis:", error);
                setError("An error occurred during the exploration analysis. Please try again.");
            }
        }

        // DE Analysis
        if (currentStage === 'differential') {
            try {
                // Calculate effective library sizes
                const effectiveLibSizes = Array.from(await WorkerManager.runTask('py', 'compute_tmm_effective_library_sizes', {
                    expression: mappedGenesExpression,
                    numSamples: numSamples,
                    numGenes: numGenes
                }));

                // Get the sample IDs from contrast and reference groups
                const selectedSampleIds = new Set([
                    ...contrastGroup.samples.map(s => s.id),
                    ...referenceGroup.samples.map(s => s.id)
                ]);

                const sampleIds = mappedGenesCountsTable.cols.slice(1).filter(col => selectedSampleIds.has(col));
                const filteredColdata = {
                    cols: dataset.coldataTable.cols,
                    data: sampleIds.map(sampleId => {
                        const rowIndex = dataset.coldataTable.data.findIndex(row => row[0] === sampleId);
                        return dataset.coldataTable.data[rowIndex];
                    })
                }

                const countsColMask = new Array(mappedGenesCountsTable.cols.length).fill(false);
                countsColMask[0] = true;
                countsColMask[1] = true;
                sampleIds.forEach(sampleId => {
                    const index = mappedGenesCountsTable.cols.indexOf(sampleId);
                    countsColMask[index] = true;
                });
                const filteredCountsTable = {
                    cols: [mappedGenesCountsTable.cols[0], mappedGenesCountsTable.cols[1], ...sampleIds],
                    data: mappedGenesCountsTable.data.map(row => row.filter((_, index) => countsColMask[index]))
                }
                const filteredCounts = Int32Array.from(filteredCountsTable.data.flatMap(row => row.slice(1)));
                const filteredEffectiveLibSizes = effectiveLibSizes.filter((_, index) => countsColMask[index + 1]);
                const filteredNumSamples = filteredColdata.data.length;

                setProgress(40);
                const deResults = await WorkerManager.runTask('r', 'run_de_analysis', {
                    counts: filteredCounts,
                    geneSymbols: mappedGeneSymbols,
                    coldata: filteredColdata,
                    contrastGroup: contrastGroup.samples.map(sample => sample.id),
                    referenceGroup: referenceGroup.samples.map(sample => sample.id),
                    libSizes: filteredEffectiveLibSizes,
                    numSamples: filteredNumSamples,
                    numGenes: numGenes
                });
                setDataset({ ...dataset, deResults: deResults });
                const deTable = {
                    cols: ['symbol', 'logFC', 't', 'p_value', 'p_value_adj'],
                    data: [...Array(deResults.row_names.length)].map((_, i) => [deResults.row_names[i], deResults.logFC[i], deResults.t[i], deResults.p_value[i], deResults.p_value_adj[i]])
                };

                setProgress(60);
                const volcanoPlot = await WorkerManager.runTask('py', 'create_volcano_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    pval_thresh: 0.05,   // todo: make this dynamic or allow the user to set it in the analysis input form
                    lfc_thresh: 1.5,    // todo: make this dynamic or allow the user to set it in the analysis input form
                    cohort_name: 'cohort name goes here 🚧'  // todo: allow this to be set by the user
                });

                setProgress(80);
                const meanDifferencePlot = await WorkerManager.runTask('py', 'create_mean_difference_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    fdr: 0.05,   // todo: allow this to be set by the user
                    cohort_name: 'cohort name goes here 🚧'  // todo: allow this to be set by the user
                });

                setDeData({ table: deTable, plots: { meanDifference: meanDifferencePlot, volcano: volcanoPlot } });
                updateLockedState('enrichment', false);
            } catch (error) {
                console.error("Error in differential expression analysis:", error);
                setError("An error occurred during the differential expression analysis. Please try again.");
            }
        }

        // GSEA
        if (currentStage === 'enrichment') {
            try {
                const combinedGeneSets = geneSetCollections.reduce((acc, collection) => {
                    return {
                        ...acc, ...collection.data.geneSets.reduce((setAcc, setName, index) => {
                            setAcc[setName] = collection.data.geneSetSymbols[index];
                            return setAcc;
                        }, {})
                    };
                }, {});
                setProgress(40);
                const gseaResults = await WorkerManager.runTask('rust', 'run_gsea', {
                    genes: mappedGeneSymbols,
                    metric: dataset.deResults.t,
                    geneSets: combinedGeneSets,
                    weight: gseaParams.weight,
                    minSize: gseaParams.minSize,
                    maxSize: gseaParams.maxSize,
                    nperm: gseaParams.nperm,
                    seed: Date.now()
                });
                setProgress(60);
                const geneConceptNetwork = await WorkerManager.runTask('py', 'create_gene_concept_network', {
                    gsea_res: gseaResults,
                    de_res: dataset.deResults,
                    color_metric: 'P.Value',
                    pvalue_threshold: 0.05,
                    layout_seed: Date.now(),
                    color_seed: Date.now()
                });
                setProgress(80);
                setGseaData({ table: gseaResults, plots: { geneConceptNetwork: geneConceptNetwork } });
            } catch (error) {
                console.error("Error in gene set enrichment analysis:", error);
                setError("An error occurred during the gene set enrichment analysis. Please try again.");
            }
        }
        setProgress(100);
        setIsLoading(false);
    };

    const renderTable = () => {
        if ((!tableData || !tableColumns.length) && currentStage !== 'differential') {
            return (
                <div className='analysisContentGuide'>
                    <h1>Analysis has not yet been run</h1>
                </div>
            );
        }

        let currentTableData = tableData;
        let currentTableColumns = tableColumns;

        if (currentStage === 'differential' && deData && deData.table) {
            const rawData = deData.table.data;
            const cols = deData.table.cols;
            currentTableData = rawData.map((row, index) => {
                const obj = { id: index };
                cols.forEach((col, colIndex) => {
                    obj[col] = row[colIndex];
                });
                return obj;
            });
            currentTableColumns = cols.map(col => ({ key: col, name: col }));
            return (
                <div>
                    <h1>DE Results</h1>
                    <DataTable
                        data={currentTableData}
                        columns={currentTableColumns}
                        contrastGroup={contrastGroup}
                        referenceGroup={referenceGroup}
                        onAddSamplesToGroup={handleAddSamplesToGroup}
                        onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                        requiresToolbar={false}
                    />
                </div>
            );
        }

        const requiresToolbar = currentStage !== 'exploration';
        console.log("Current stage: ", currentStage);

        return (
            <DataTable
                data={currentTableData}
                columns={currentTableColumns}
                contrastGroup={contrastGroup}
                referenceGroup={referenceGroup}
                onAddSamplesToGroup={handleAddSamplesToGroup}
                onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                requiresToolbar={requiresToolbar}
            />
        );
    };

    const handleRemoveGeneSetCollection = (index) => {
        setGeneSetCollections(prev => prev.filter((_, i) => i !== index));
    };

    if (error) {
        return <div>Error: {error}</div>;
    }

    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                {currentStage === 'exploration' && (
                    <EDAInputForm
                        setIsVisible={setIsVisible}
                        onDatasetSelect={handleDatasetSelect}
                        onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                        runAnalysis={runAnalysis}
                        isLoading={isLoading}
                        handleStageChange={handleStageChange}
                        currentStage={currentStage}
                        edaData={edaData}
                    />
                )}
                {currentStage === 'differential' && (
                    <DEAInputForm
                        contrastGroup={contrastGroup}
                        referenceGroup={referenceGroup}
                        onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                        runAnalysis={runAnalysis}
                        handleStageChange={handleStageChange}
                        currentStage={currentStage}
                        deData={deData}
                    />
                )}
                {currentStage === 'enrichment' && (
                    <GSEAInputForm
                        setIsVisible={setIsVisible}
                        onDatasetSelect={handleDatasetSelect}
                        runAnalysis={runAnalysis}
                        isLoading={isLoading}
                        onAddGeneSetCollection={(species, collectionId) => getExternalGeneSetCollection(species, collectionId)
                            .then(newCollection => setGeneSetCollections(prev => [...prev, { name: collectionId, data: newCollection }]))}
                        geneSetCollections={geneSetCollections}
                        gseaParams={gseaParams}
                        onUpdateGseaParams={setGseaParams}
                        onRemoveGeneSetCollection={handleRemoveGeneSetCollection}
                    />
                )}
            </div>
            <DatabasePopup setIsVisible={setIsVisible} isVisible={isVisible} onDatasetSelect={handleDatasetSelect} />
            <div id="analysis_visualization_section">
                <div id="analysis_tab_nav">
                    <TabButton
                        label="Data Exploration"
                        onClick={() => handleStageChange('exploration')}
                        isActive={currentStage === 'exploration'}
                        isLocked={lockedStates.exploration}
                    />
                    <TabButton
                        label="Differential Expression Analysis"
                        onClick={() => handleStageChange('differential')}
                        isActive={currentStage === 'differential'}
                        isLocked={lockedStates.differential}
                    />
                    <TabButton
                        label="Gene Set Enrichment Analysis"
                        onClick={() => handleStageChange('enrichment')}
                        isActive={currentStage === 'enrichment'}
                        isLocked={lockedStates.enrichment}
                    />
                </div>
                <div id="analysis_content">
                    <div id="view_content">
                        {currentStage === 'exploration' && (
                            <ExplorationContent
                                data={edaData}
                                activeTab={activeTab}
                                setActiveTab={setActiveTab}
                                isLoading={isLoading}
                                progress={progress}
                                renderTable={renderTable}
                            />
                        )}
                        {currentStage === 'differential' && (
                            <DifferentialExpressionContent
                                data={deData}
                                activeTab={activeTab}
                                setActiveTab={setActiveTab}
                                isLoading={isLoading}
                                progress={progress}
                                renderTable={renderTable}
                                currentTable={currentTable}
                                setCurrentTable={setCurrentTable}
                                currentPlot={currentPlot}
                                setCurrentPlot={setCurrentPlot}
                            />
                        )}
                        {currentStage === 'enrichment' && (
                            <GSEAContent
                                data={gseaData}
                                activeTab={activeTab}
                                setActiveTab={setActiveTab}
                                onAddSamplesToGroup={handleAddSamplesToGroup}
                                onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                                contrastGroup={contrastGroup}
                                referenceGroup={referenceGroup}
                                isLoading={isLoading}
                                progress={progress}
                            />
                        )}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default Analysis;
