import React, { useState, useEffect, useCallback, useRef } from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
import WorkerManager from '../utils/Workers';
import PlotArea from '../components/ui/PlotArea';
import ProgressBar from '../components/ui/ProgressBar';
import { getExternalDataset } from '../services/api';
import DifferentialExpressionContent from '../components/ui/DifferentialExpressionContent';
import ExplorationContent from '../components/ui/ExplorationContent';
import GSEAContent from '../components/ui/GSEAContent';
import GSEAInputForm from '../components/form/GSEAInputForm';
import DEAInputForm from '../components/form/DEAInputForm';
import EDAInputForm from '../components/form/EDAInputForm';


const Analysis = () => {
    const [activeTab, setActiveTab] = useState('table');
    const [tableScrollPosition, setTableScrollPosition] = useState(0);
    const [error, setError] = useState(null);
    const [isVisible, setIsVisible] = useState(false);
    const tableContainerRef = useRef(null);
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

    const [currentStage, setCurrentStage] = useState('exploration');

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

    useEffect(() => {
        if (activeTab === 'table' && tableContainerRef.current) {
            tableContainerRef.current.scrollTop = tableScrollPosition;
        }
    }, [activeTab, tableScrollPosition]);

    const handleTableScroll = (event) => {
        setTableScrollPosition(event.target.scrollTop);
    };

    const handleDatasetSelect = async (type, data) => {
        if (type === 'external' || type === 'upload') {
            setDataset(data);
        } else if (type === 'example') {
            setIsLoading(true);
            try {
                const exampleData = await getExternalDataset('GDC', 'CDDP_EAGLE-1');
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
            } catch (error) {
                console.error('Error loading example dataset:', error);
                setError('Failed to load example dataset');
            } finally {
                setIsLoading(false);
            }
        }

        if (data && data.coldataTable) {
            const tableData = data.coldataTable.data.map((row, index) => {
                const obj = {};
                data.coldataTable.cols.forEach((col, colIndex) => {
                    obj[col] = row[colIndex];
                });
                obj.id = index;
                return obj;
            });

            const columns = data.coldataTable.cols.map(col => ({ key: col, name: col }));

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

    const renderTableButtons = () => {
        let tables = [];
        if (edaData && edaData.tables) {
            tables = Object.keys(edaData.tables);
        }
        if (deData && deData.table) {
            tables.push('de_results');
        }
        if (gseaData && gseaData.table) {
            tables.push('gsea_results');
        }
        return tables.map(table => (
            <button
                key={table}
                onClick={() => setCurrentTable(table)}
                className={`view-toggle-btn ${currentTable === table ? 'active' : ''}`}
            >
                {table}
            </button>
        ));
    };

    const renderPlotButtons = () => {
        let plots = [];
        if (edaData && edaData.plots) {
            plots = Object.keys(edaData.plots);
        }
        if (deData && deData.plots) {
            plots = [...plots, ...Object.keys(deData.plots)];
        }
        if (gseaData && gseaData.plots) {
            plots = [...plots, ...Object.keys(gseaData.plots)];
        }
        return plots.map(plot => (
            <button
                key={plot}
                onClick={() => setCurrentPlot(plot)}
                className={`view-toggle-btn ${currentPlot === plot ? 'active' : ''}`}
            >
                {plot}
            </button>
        ));
    };

    // const renderTable = () => {
    //     let tableData, tableColumns;
    //     if (currentTable === 'coldata' || currentTable === 'counts') {
    //         if (edaData && edaData.tables && edaData.tables[currentTable]) {
    //             const rawData = edaData.tables[currentTable].data;
    //             const cols = edaData.tables[currentTable].cols;
    //             tableData = rawData.map(row => {
    //                 const obj = {};
    //                 cols.forEach((col, index) => {
    //                     obj[col] = row[index];
    //                 });
    //                 return obj;
    //             });
    //             tableColumns = cols.map(col => ({ key: col, name: col }));
    //         }
    //     } else if (currentTable === 'de_results') {
    //         if (deData && deData.table) {
    //             const rawData = deData.table.data;
    //             const cols = deData.table.cols;
    //             tableData = rawData.map(row => {
    //                 const obj = {};
    //                 cols.forEach((col, index) => {
    //                     obj[col] = row[index];
    //                 });
    //                 return obj;
    //             });
    //             tableColumns = cols.map(col => ({ key: col, name: col }));
    //         }
    //     } else if (currentTable === 'gsea_results') {
    //         if (gseaData && gseaData.table) {
    //             const rawData = gseaData.table.data;
    //             const cols = gseaData.table.cols;
    //             tableData = rawData.map(row => {
    //                 const obj = {};
    //                 cols.forEach((col, index) => {
    //                     obj[col] = row[index];
    //                 });
    //                 return obj;
    //             });
    //             tableColumns = cols.map(col => ({ key: col, name: col }));
    //         }
    //     }

    //     if (tableData && tableColumns) {
    //         return (
    //             <DataTable
    //                 data={tableData}
    //                 columns={tableColumns}
    //                 contrastGroup={contrastGroup}
    //                 referenceGroup={referenceGroup}
    //                 onAddSamplesToGroup={handleAddSamplesToGroup}
    //                 onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
    //             />
    //         );
    //     } else {
    //         return <p>No data available</p>;
    //     }
    // };

    const renderPlot = () => {
        let plotHtml;
        if (edaData && edaData.plots && edaData.plots[currentPlot]) {
            plotHtml = edaData.plots[currentPlot];
        } else if (deData && deData.plots && deData.plots[currentPlot]) {
            plotHtml = deData.plots[currentPlot];
        } else if (gseaData && gseaData.plots && gseaData.plots[currentPlot]) {
            plotHtml = gseaData.plots[currentPlot];
        }
        console.log('Plot HTML:', plotHtml);
        return plotHtml ? (
            <PlotArea htmlContent={plotHtml} />
        ) : (
            <p>No plot available</p>
        );
    };

    const runAnalysis = async () => {
        setIsLoading(true);
        setProgress(0);

        const numGenes = dataset.countsTable.rows.length;
        const numSamples = dataset.coldataTable.rows.length;

        // EDA
        if (currentStage === 'exploration') {
            try {
                setProgress(10);
                const transformedCounts = await WorkerManager.runTask('py', 'transform_vst', {
                    expression: dataset.expression,
                    numSamples: numSamples,
                    numGenes: numGenes
                });

                setProgress(20);
                const pcaPlot = await WorkerManager.runTask('py', 'create_pca', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: dataset.countsTable.cols.slice(2)
                });

                setProgress(27);
                const tsnePlot = await WorkerManager.runTask('py', 'create_tsne', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: dataset.countsTable.cols.slice(2)
                });

                setProgress(34);
                const heatmap = await WorkerManager.runTask('py', 'create_heatmap', {
                    counts: transformedCounts,
                    numSamples: numSamples,
                    numGenes: numGenes,
                    sample_ids: dataset.countsTable.cols.slice(2)
                });

                setEdaData({
                    tables: {
                        coldata: dataset.coldataTable,
                        counts: dataset.countsTable
                    },
                    plots: { pca: pcaPlot, tsne: tsnePlot, heatmap: heatmap }
                });
            } catch (error) {
                console.error("Error in exploration analysis:", error);
                setError("An error occurred during the exploration analysis. Please try again.");
            }
        }
        // DE Analysis
        if (currentStage === 'differential') {
            try {
                setProgress(10);
                // Calculate effective library sizes
                const effectiveLibSizes = Array.from(await WorkerManager.runTask('py', 'compute_tmm_effective_library_sizes', {
                    expression: dataset.expression,
                    numSamples: numSamples,
                    numGenes: numGenes
                }));

                // Get the sample IDs from contrast and reference groups
                const selectedSampleIds = new Set([
                    ...contrastGroup.samples.map(s => s.id),
                    ...referenceGroup.samples.map(s => s.id)
                ]);

                const sampleIds = dataset.countsTable.cols.slice(2).filter(col => selectedSampleIds.has(col));
                const filteredColdata = {
                    cols: dataset.coldataTable.cols,
                    data: sampleIds.map(sampleId => {
                        const rowIndex = dataset.coldataTable.data.findIndex(row => row[0] === sampleId);
                        return dataset.coldataTable.data[rowIndex];
                    })
                }

                const countsColMask = new Array(dataset.countsTable.cols.length).fill(false);
                countsColMask[0] = true;
                countsColMask[1] = true;
                sampleIds.forEach(sampleId => {
                    const index = dataset.countsTable.cols.indexOf(sampleId);
                    countsColMask[index] = true;
                });
                const filteredCountsTable = {
                    cols: [dataset.countsTable.cols[0], dataset.countsTable.cols[1], ...sampleIds],
                    data: dataset.countsTable.data.map(row => row.filter((_, index) => countsColMask[index]))
                }
                const filteredCounts = Int32Array.from(filteredCountsTable.data.flatMap(row => row.slice(2)));
                const filteredEffectiveLibSizes = effectiveLibSizes.filter((_, index) => countsColMask[index + 2]);
                const filteredNumSamples = filteredColdata.data.length;

                setProgress(50);
                const deResults = await WorkerManager.runTask('r', 'run_de_analysis', {
                    counts: filteredCounts,
                    ensemblIds: dataset.countsTable.rows,
                    coldata: filteredColdata,
                    contrastGroup: contrastGroup.samples.map(sample => sample.id),
                    referenceGroup: referenceGroup.samples.map(sample => sample.id),
                    libSizes: filteredEffectiveLibSizes,
                    numSamples: filteredNumSamples,
                    numGenes: numGenes
                });

                const deTable = {
                    cols: ['ensembl_id', 'logFC', 't', 'p_value', 'p_value_adj'],
                    data: [...Array(deResults.row_names.length)].map((_, i) => [deResults.row_names[i], deResults.logFC[i], deResults.t[i], deResults.p_value[i], deResults.p_value_adj[i]])
                }

                setProgress(70);
                const volcanoPlot = await WorkerManager.runTask('py', 'create_volcano_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    pval_thresh: 0.05,
                    lfc_thresh: 1.5,
                    cohort_name: 'cohort name goes here 🚧'
                });

                setProgress(90);
                const meanDifferencePlot = await WorkerManager.runTask('py', 'create_mean_difference_plot', {
                    data: deResults.data,
                    row_names: deResults.row_names,
                    column_names: deResults.column_names,
                    fdr: 0.05,
                    cohort_name: 'cohort name goes here 🚧'
                });

                setDeData({ table: deTable, plots: { meanDifference: meanDifferencePlot, volcano: volcanoPlot } });
            } catch (error) {
                console.error("Error in differential expression analysis:", error);
                setError("An error occurred during the differential expression analysis. Please try again.");
            }
        }
        // GSEA
        // setProgress(70);
        // const gseaResults = await WorkerManager.runTask('rust', 'run_gsea', {
        //     genes: deResults.row_names,
        //     metric: deResults.t,
        //     geneSets: geneSetCollections,
        //     weight: 1,
        //     minSize: 15,
        //     maxSize: 500,
        //     nperm: 1000,
        //     seed: Date.now()
        // });
        // setProgress(80);
        // const geneConceptNetwork = await WorkerManager.runTask('py', 'create_gene_concept_network', {
        //     gsea_res: gseaResults,
        //     de_res: deResults,
        //     ensembl_to_symbol: ensemblToSymbol,
        //     color_metric: 'P.Value',
        //     pvalue_threshold: 0.05,
        //     layout_seed: Date.now(),
        //     color_seed: Date.now()
        // });
        // setProgress(90);

        // Update state with results
        // setGseaData({ table: gseaResults, plots: { geneConceptNetwork: geneConceptNetwork } });
        setIsLoading(false);
    };

    const renderTable = () => {
        if ((!tableData || !tableColumns.length) && currentStage !== 'differential') {
            return (
                <div className='analysisContentGuide'>
                    <h1>Select dataset then run analysis to see results</h1>
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
                    />
                )}
                {currentStage === 'differential' && (
                    <DEAInputForm
                        contrastGroup={contrastGroup}
                        referenceGroup={referenceGroup}
                        onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                        runAnalysis={runAnalysis}
                    />
                )}
                {currentStage === 'enrichment' && (
                    <GSEAInputForm
                        setIsVisible={setIsVisible}
                        onDatasetSelect={handleDatasetSelect}
                        runAnalysis={runAnalysis}
                        isLoading={isLoading}
                    />
                )}
            </div>
            <DatabasePopup setIsVisible={setIsVisible} isVisible={isVisible} onDatasetSelect={handleDatasetSelect} />

            <div id="analysis_visualization_section">
                <div id="analysis_tab_nav">
                    <TabButton label="1. Data Exploration" onClick={() => handleStageChange('exploration')} />
                    <TabButton label="2. Differential Expression Analysis" onClick={() => handleStageChange('differential')} />
                    <TabButton label="3. Gene Set Enrichment Analysis" onClick={() => handleStageChange('enrichment')} />
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
