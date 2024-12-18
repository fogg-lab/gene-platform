import React, { useState, useEffect, useCallback } from 'react';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
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
import { formatNumber } from '../lib/utils';

const Analysis = () => {
    const [activeTab, setActiveTab] = useState('table');
    const [error, setError] = useState(null);
    const [isVisible, setIsVisible] = useState(false);
    const [dataset, setDataset] = useState(null);
    const [edaData, setEdaData] = useState(null);
    const [deData, setDeData] = useState(null);
    const [enrichData, setEnrichData] = useState(null);
    const [currentTable, setCurrentTable] = useState(null);
    const [currentPlot, setCurrentPlot] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [tableColumns, setTableColumns] = useState([]);

    const [contrastGroup, setContrastGroup] = useState({ samples: [] });
    const [referenceGroup, setReferenceGroup] = useState({ samples: [] });

    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);

    const [geneSetCollections, setGeneSetCollections] = useState([]);
    const [enrichParams, setEnrichParams] = useState({
        weight: 1,
        minSize: 15,
        maxSize: 500,
        nperm: 1000
    });

    const [currentStage, setCurrentStage] = useState('exploration');
    const [lockedStates, setLockedStates] = useState({
        exploration: false,
        differential: true,
        enrichment: true
    });
    const [transformMethod, setTransformMethod] = useState('log2');

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
        if (type === 'external') {
            setDataset(data);
            dset = data;
        } else if (type === 'upload') {
            setDataset(prev => ({
                ...prev,
                ...data
            }));
            if (data.tableData && data.tableColumns) {
                setTableData(data.tableData);
                setTableColumns(data.tableColumns);
            }
            return;
        } else if (type === 'example') {
            setIsLoading(true);
            try {
                const exampleData = await getExternalDataset('GDC', 'TCGA-CESC');
                dset = exampleData;
                setDataset(exampleData);

                handleClearGroup(true);
                handleClearGroup(false);

                // Get all samples and split them by FIGO stage
                const allSamples = exampleData.coldataTable.data.map(row => ({
                    id: row[exampleData.coldataTable.cols.indexOf('sample_id')],
                    [exampleData.coldataTable.cols[0]]: row[0],
                    figoStage: row[exampleData.coldataTable.cols.indexOf('figo_stage')]
                }));

                // Reference group: Stage I samples
                setReferenceGroup({
                    samples: allSamples.filter(sample => {
                        const stage = sample.figoStage;
                        return stage && 
                               stage.startsWith('Stage I') && 
                               !stage.startsWith('Stage II') &&
                               !stage.startsWith('Stage IV');
                    })
                });

                // Contrast group: Stage III samples
                setContrastGroup({
                    samples: allSamples.filter(sample => {
                        const stage = sample.figoStage;
                        return stage && stage.startsWith('Stage III');
                    })
                });

                const collectionSpecies = "human";
                const collectionId = "C2:CP:REACTOME";
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
                if (!currentPlot) setCurrentPlot('pca');
                break;
            case 'differential':
                setActiveTab('table');
                if (deData) {
                    setCurrentTable('de_results');
                    if (!currentPlot) setCurrentPlot('volcano_plot');
                } else {
                    setCurrentTable('coldata');
                }
                break;
            case 'enrichment':
                if (enrichData) {
                    setCurrentTable('gsea_results');
                    //if (!currentPlot) setCurrentPlot('gene_concept_network');
                }
                break;
        }
    };

    const runAnalysis = async (processedData) => {
        setIsLoading(true);
        setProgress(0);

        if (currentStage === 'enrichment') {
            if (!dataset || !dataset.deResults) {
                setError("Run differential expression analysis first");
                setIsLoading(false);
                return;
            }

            try {
                // Filter gene sets by minimum size
                const filteredGeneSets = {
                    sets: [],
                    symbols: []
                };

                geneSetCollections.forEach(collection => {
                    collection.data.geneSets.forEach((set, idx) => {
                        if (set.length >= (enrichParams.minSize || 15)) {
                            filteredGeneSets.sets.push(set);
                            filteredGeneSets.symbols.push(collection.data.geneSetSymbols[idx]);
                        }
                    });
                });

                // Create group vector
                const group = dataset.includedColdata.data.map(row => {
                    const sampleId = row[dataset.includedColdata.cols.indexOf('sample_id')];
                    return dataset.contrastSampleIds.has(sampleId) ? "contrastGrp" : "referenceGrp";
                });

                // Create design matrix
                const design = dataset.includedColdata.data.map((_, i) => [
                    group[i] === "referenceGrp" ? 1 : 0,
                    group[i] === "contrastGrp" ? 1 : 0
                ]);

                const gseaResults = await WorkerManager.runTask('r', 'run_camera', {
                    counts: dataset.filteredCounts,
                    numSamples: dataset.includedSampleIds.length,
                    numGenes: dataset.geneSymbols.length,
                    geneSymbols: dataset.geneSymbols,
                    geneSetNames: filteredGeneSets.sets,
                    geneSetSymbols: filteredGeneSets.symbols,
                    design: design,
                    contrast: [-1, 1],
                    minSetSize: enrichParams.minSize || 15
                });

                setEnrichData({
                    tables: {
                        results: {
                            cols: ['Name', 'NGenes', 'Direction', 'PValue', 'FDR'],
                            data: gseaResults.Name.map((name, i) => [
                                name,
                                gseaResults.NGenes[i],
                                gseaResults.Direction[i],
                                gseaResults.PValue[i],
                                gseaResults.FDR[i]
                            ])
                        }
                    }
                });

            } catch (error) {
                console.error("Error in GSEA:", error);
                setError("An error occurred during the gene set enrichment analysis.");
            }
            setIsLoading(false);
            return;
        }

        const analysisData = processedData || dataset;
        const numInputGenes = analysisData.countsTable.rows.length;
        const numSamples = analysisData.coldataTable.rows.length;
        const expressionData = analysisData.expression;

        const geneInfo = await WorkerManager.runTask('py', 'check_genes', {
            genesQuery: analysisData.countsTable.rows,
            humanGeneReference: humanGenes.genes,
            mouseGeneReference: mouseGenes.genes
        });

        const geneReference = geneInfo.detectedSpecies === "mouse" ? mouseGenes.genes : humanGenes.genes;

        // Filter and map genes
        const geneIndex = Array.from({ length: numInputGenes }, (_, i) => i);
        const geneMask = geneIndex.map((i) => geneInfo.queryResultCode[i] == 0);
        const mappedGeneIdx = geneIndex.filter((i) => geneMask[i]);
        const mappedGeneSymbols = mappedGeneIdx.map((i) => geneReference[geneInfo.queryResultGeneIndices[i]][2]);
        const numGenes = mappedGeneSymbols.length;
        const mappedGenesExpression = new (expressionData || dataset.expression).constructor(numGenes * numSamples);
        let mappedGenesExpressionIndex = 0;
        let colIndex = 0;
        for (let i = 0; i < (expressionData || dataset.expression).length; i++) {
            if (geneMask[colIndex++]) {
                mappedGenesExpression[mappedGenesExpressionIndex++] = (expressionData || dataset.expression)[i];
            }
            if (colIndex === numInputGenes) colIndex = 0;
        }
        const mappedGenesCountsTable = {
            cols: ["Symbol", ...(analysisData.countsTable || dataset.countsTable).cols],
            rows: mappedGeneSymbols,
            data: mappedGeneIdx.map((i) => [mappedGeneSymbols[i], ...(analysisData.countsTable || dataset.countsTable).data[i]])
        }
        setProgress(20);

        // EDA
        if (currentStage === 'exploration') {
            try {
                const transformedCounts = await WorkerManager.runTask('py', `transform_${transformMethod}`, {
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
                        coldata: (analysisData.coldataTable || dataset.coldataTable),
                        counts: mappedGenesCountsTable
                    },
                    plots: { pca: pcaPlot, tsne: tsnePlot, heatmap: heatmap }
                });
                setCurrentTable('coldata');
                updateLockedState('differential', false);
            } catch (error) {
                console.error("Error in exploration analysis:", error);
                setError("An error occurred during the exploration analysis.");
            }
        }

        // DE Analysis
        if (currentStage === 'differential') {
            try {
                const effectiveLibSizes = Array.from(await WorkerManager.runTask('py', 'compute_tmm_effective_library_sizes', {
                    expression: mappedGenesExpression,
                    numSamples: numSamples,
                    numGenes: numGenes
                }));

                const contrastSampleIds = new Set(contrastGroup.samples.map(s => s.id));
                const referenceSampleIds = new Set(referenceGroup.samples.map(s => s.id));
                const includedSampleIds = new Set([...contrastSampleIds, ...referenceSampleIds]);

                const sampleIds = mappedGenesCountsTable.cols.slice(1).filter(col => includedSampleIds.has(col));
                const countsColMask = new Array(mappedGenesCountsTable.cols.length).fill(false);
                countsColMask[0] = true;
                sampleIds.forEach(sampleId => {
                    const index = mappedGenesCountsTable.cols.indexOf(sampleId);
                    countsColMask[index] = true;
                });

                const filteredCountsTable = {
                    cols: ["Symbol", ...sampleIds],
                    data: mappedGenesCountsTable.data.map(row => 
                        row.filter((_, index) => countsColMask[index])
                    )
                };

                const filteredCounts = Int32Array.from(
                    filteredCountsTable.data.flatMap(row => row.slice(1))
                );
                const filteredEffectiveLibSizes = effectiveLibSizes.filter((_, index) => countsColMask[index + 1]);
                const filteredNumSamples = filteredEffectiveLibSizes.length;

                const includedColdata = {
                    cols: dataset.coldataTable.cols,
                    data: dataset.coldataTable.data.filter(
                        row => includedSampleIds.has(row[dataset.coldataTable.cols.indexOf('sample_id')])
                    )
                };

                setProgress(40);
                const deResults = await WorkerManager.runTask('r', 'run_de_analysis', {
                    counts: filteredCounts,
                    geneSymbols: mappedGeneSymbols,
                    coldata: includedColdata,
                    contrastGroup: contrastGroup.samples.map(sample => sample.id),
                    referenceGroup: referenceGroup.samples.map(sample => sample.id),
                    libSizes: filteredEffectiveLibSizes,
                    numSamples: filteredNumSamples,
                    numGenes: numGenes
                });
                setDataset(prevDataset => {
                    return {
                        ...prevDataset,
                        includedColdata: includedColdata,
                        deResults: deResults,
                        geneSymbols: mappedGeneSymbols,
                        filteredCounts: filteredCounts,
                        filteredCountsTable: filteredCountsTable,
                        contrastSampleIds: contrastSampleIds,
                        referenceSampleIds: referenceSampleIds,
                        includedSampleIds: includedSampleIds
                    };
                });
                const deTable = {
                    cols: ['symbol', 'logFC', 't', 'p_value', 'p_value_adj'],
                    data: [...Array(deResults.row_names.length)].map((_, i) => [deResults.row_names[i], deResults.logFC[i], deResults.t[i], deResults.p_value[i], deResults.p_value_adj[i]])
                };

                setProgress(60);
                const volcanoPlot = await WorkerManager.runTask('py', 'create_volcano_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    pval_thresh: 0.05,
                    lfc_thresh: 1.5,
                    cohort_name: 'cohort name goes here 🚧'
                });

                setProgress(80);
                const meanDifferencePlot = await WorkerManager.runTask('py', 'create_mean_difference_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    fdr: 0.05,
                    cohort_name: 'cohort name goes here 🚧'
                });

                setDeData({ table: deTable, plots: { meanDifference: meanDifferencePlot, volcano: volcanoPlot } });
                updateLockedState('enrichment', false);
                setCurrentTable('de_results');
            } catch (error) {
                console.error("Error in differential expression analysis:", error);
                setError("An error occurred during the differential expression analysis.");
            }
        }

        setProgress(100);
        setIsLoading(false);
    };

    const renderTable = (filterProps) => {
        if ((!tableData || !tableColumns.length) && currentStage !== 'differential' && currentStage !== 'enrichment') {
            return (
                <div className='analysisContentGuide'>
                    <h1>Run analysis to view results</h1>
                </div>
            );
        }

        const enableSampleSelection = currentStage === 'differential' && currentTable === 'coldata';
        let currentTableData = tableData;
        let currentTableColumns = tableColumns;

        if (currentStage === 'exploration') {
            switch (currentTable) {
                case 'coldata':
                    currentTableData = dataset.coldataTable.data.map((row, index) => {
                        const obj = {};
                        dataset.coldataTable.cols.forEach((col, colIndex) => {
                            obj[col] = row[colIndex];
                        });
                        obj.id = index;
                        return obj;
                    });
                    currentTableColumns = dataset.coldataTable.cols.map(col => ({ key: col, name: col }));
                    break;
                case 'counts':
                    currentTableData = dataset.countsTable.data.map((row, index) => {
                        const obj = {};
                        dataset.countsTable.cols.forEach((col, colIndex) => {
                            obj[col] = row[colIndex];
                        });
                        obj.id = index;
                        return obj;
                    });
                    currentTableColumns = dataset.countsTable.cols.map(col => ({ key: col, name: col }));
                    break;
                case 'transformed_counts':
                    if (edaData && edaData.tables && edaData.tables.transformed_counts) {
                        currentTableData = edaData.tables.transformed_counts.data.map((row, index) => {
                            const obj = {};
                            edaData.tables.transformed_counts.cols.forEach((col, colIndex) => {
                                obj[col] = row[colIndex];
                            });
                            obj.id = index;
                            return obj;
                        });
                        currentTableColumns = edaData.tables.transformed_counts.cols.map(col => ({ key: col, name: col }));
                    }
                    break;
            }
        } else if (currentStage === 'differential') {
            switch (currentTable) {
                case 'coldata':
                    currentTableData = dataset.coldataTable.data.map((row, index) => {
                        const obj = {};
                        dataset.coldataTable.cols.forEach((col, colIndex) => {
                            obj[col] = row[colIndex];
                        });
                        obj.id = index;
                        return obj;
                    });
                    currentTableColumns = dataset.coldataTable.cols.map(col => ({ key: col, name: col }));
                    break;
                case 'counts':
                    currentTableData = dataset.countsTable.data.map((row, index) => {
                        const obj = {};
                        dataset.countsTable.cols.forEach((col, colIndex) => {
                            obj[col] = row[colIndex];
                        });
                        obj.id = index;
                        return obj;
                    });
                    currentTableColumns = dataset.countsTable.cols.map(col => ({ key: col, name: col }));
                    break;
                case 'de_results':
                    if (deData && deData.table) {
                        currentTableData = deData.table.data.map((row, index) => {
                            const obj = {};
                            deData.table.cols.forEach((col, colIndex) => {
                                obj[col] = row[colIndex];
                            });
                            obj.id = index;
                            return obj;
                        });
                        currentTableColumns = deData.table.cols.map(col => ({ key: col, name: col }));
                    }
                    break;
            }
        }

        return (
            <DataTable
                data={currentTableData}
                columns={currentTableColumns}
                contrastGroup={contrastGroup}
                referenceGroup={referenceGroup}
                onAddSamplesToGroup={handleAddSamplesToGroup}
                onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                requiresToolbar={enableSampleSelection}
                activeFilter={filterProps?.activeFilter}
            />
        );
    };

    const handleRemoveGeneSetCollection = (index) => {
        setGeneSetCollections(prev => prev.filter((_, i) => i !== index));
    };

    const handleTransformMethodChange = (method) => {
        setTransformMethod(method);
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
                        onTransformMethodChange={handleTransformMethodChange}
                        dataset={dataset}
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
                        onAddGeneSetCollection={(species, collectionId) =>
                            getExternalGeneSetCollection(species, collectionId)
                                .then(newCollection =>
                                    setGeneSetCollections(prev => [...prev, { name: collectionId, data: newCollection }])
                                )
                        }
                        geneSetCollections={geneSetCollections}
                        enrichParams={enrichParams}
                        onUpdateEnrichParams={setEnrichParams}
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
                                currentTable={currentTable}
                                setCurrentTable={setCurrentTable}
                                dataset={dataset}
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
                                dataset={dataset}
                            />
                        )}
                        {currentStage === 'enrichment' && (
                            <GSEAContent
                                data={enrichData}
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
