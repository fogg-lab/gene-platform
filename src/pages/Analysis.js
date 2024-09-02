import React, { useState, useEffect, useCallback, useRef } from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
import WorkerManager from '../utils/Workers';
import PlotArea from '../components/ui/PlotArea';
import ProgressBar from '../components/ui/ProgressBar';
import { getExternalDataset } from '../services/api';

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

    const [contrastGroup, setContrastGroup] = useState({ samples: [] });
    const [referenceGroup, setReferenceGroup] = useState({ samples: [] });

    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);

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
    };

    const handleStageChange = (stage) => {
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

    const renderTable = () => {
        let tableData, tableColumns;
        if (currentTable === 'coldata' || currentTable === 'counts') {
            if (edaData && edaData.tables && edaData.tables[currentTable]) {
                const rawData = edaData.tables[currentTable].data;
                const cols = edaData.tables[currentTable].cols;
                tableData = rawData.map(row => {
                    const obj = {};
                    cols.forEach((col, index) => {
                        obj[col] = row[index];
                    });
                    return obj;
                });
                tableColumns = cols.map(col => ({ key: col, name: col }));
            }
        } else if (currentTable === 'de_results') {
            if (deData && deData.table) {
                const rawData = deData.table.data;
                const cols = deData.table.cols;
                tableData = rawData.map(row => {
                    const obj = {};
                    cols.forEach((col, index) => {
                        obj[col] = row[index];
                    });
                    return obj;
                });
                tableColumns = cols.map(col => ({ key: col, name: col }));
            }
        } else if (currentTable === 'gsea_results') {
            if (gseaData && gseaData.table) {
                const rawData = gseaData.table.data;
                const cols = gseaData.table.cols;
                tableData = rawData.map(row => {
                    const obj = {};
                    cols.forEach((col, index) => {
                        obj[col] = row[index];
                    });
                    return obj;
                });
                tableColumns = cols.map(col => ({ key: col, name: col }));
            }
        }

        if (tableData && tableColumns) {
            return <DataTable
                data={tableData}
                columns={tableColumns}
                contrastGroup={contrastGroup}
                referenceGroup={referenceGroup}
                onAddSamplesToGroup={handleAddSamplesToGroup}
            />;
        } else {
            return <p>No data available</p>;
        }
    };

    const renderPlot = () => {
        let plotHtml;
        if (edaData && edaData.plots && edaData.plots[currentPlot]) {
            plotHtml = edaData.plots[currentPlot];
        } else if (deData && deData.plots && deData.plots[currentPlot]) {
            plotHtml = deData.plots[currentPlot];
        } else if (gseaData && gseaData.plots && gseaData.plots[currentPlot]) {
            plotHtml = gseaData.plots[currentPlot];
        }
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
        setProgress(10);
        const transformedCounts = await WorkerManager.runTask('py', 'transform_log2', {
            counts: dataset.counts,
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

        // DE Analysis
        setProgress(40);
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

        // enforced order of samples
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

        setProgress(50);
        const volcanoPlot = await WorkerManager.runTask('py', 'create_volcano_plot', {
            data: deResults.data,
            row_names: deResults.row_names,
            column_names: deResults.column_names,
            fdr: 0.05,
            cohort_name: 'DE Analysis'
        });

        setProgress(60);
        const meanDifferencePlot = await WorkerManager.runTask('py', 'create_mean_difference_plot', {
            data: deResults.data,
            row_names: deResults.row_names,
            column_names: deResults.column_names,
            fdr: 0.05,
            cohort_name: 'DE Analysis'
        });

        // GSEA
        setProgress(70);
        const gseaResults = await WorkerManager.runTask('rust', 'run_gsea', {
            genes: deResults.genes,
            metric: deResults.logFC,
            geneSets: geneSetCollections,
            weight: 1,
            minSize: 15,
            maxSize: 500,
            nperm: 1000,
            seed: Date.now()
        });
        setProgress(80);
        const geneConceptNetwork = await WorkerManager.runTask('py', 'create_gene_concept_network', {
            gsea_res: gseaResults,
            de_res: deResults,
            ensembl_to_symbol: ensemblToSymbol,
            color_metric: 'P.Value',
            pvalue_threshold: 0.05,
            layout_seed: Date.now(),
            color_seed: Date.now()
        });
        setProgress(90);

        // Update state with results
        setEdaData({ plots: { pca: pcaPlot, tsne: tsnePlot, heatmap: heatmap } });
        setDeData({ table: deResults, plots: { meanDifference: meanDifferencePlot, volcano: volcanoPlot } });
        setGseaData({ table: gseaResults, plots: { geneConceptNetwork: geneConceptNetwork } });
        setProgress(100);

        setIsLoading(false);
    };

    if (error) {
        return <div>Error: {error}</div>;
    }

    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                <AnalysisInputForm
                    setIsVisible={setIsVisible}
                    onDatasetSelect={handleDatasetSelect}
                    contrastGroup={contrastGroup}
                    referenceGroup={referenceGroup}
                    onRemoveSamplesFromGroup={handleRemoveSamplesFromGroup}
                    runAnalysis={runAnalysis}
                    isLoading={isLoading}
                />
            </div>
            <DatabasePopup setIsVisible={setIsVisible} isVisible={isVisible} onDatasetSelect={handleDatasetSelect} />
            <div id="analysis_visualization_section">
                <div id="analysis_tab_nav">
                    <TabButton label="Data Exploration" onClick={() => handleStageChange('exploration')} />
                    <TabButton label="Differential Expression Analysis" onClick={() => handleStageChange('differential')} />
                    <TabButton label="Gene Set Enrichment Analysis" onClick={() => handleStageChange('enrichment')} />
                </div>
                <div id="analysis_content">
                    <div id="view_toggle">
                        <button
                            className={`view-toggle-btn ${activeTab === 'table' ? 'active' : ''}`}
                            onClick={() => setActiveTab('table')}
                        >
                            Table View
                        </button>
                        <button
                            className={`view-toggle-btn ${activeTab === 'plot' ? 'active' : ''}`}
                            onClick={() => setActiveTab('plot')}
                        >
                            Plot View
                        </button>
                    </div>
                    <div id="view_content">
                        <div
                            style={{ display: activeTab === 'table' ? 'block' : 'none', height: '100%', overflow: 'auto' }}
                            ref={tableContainerRef}
                            onScroll={handleTableScroll}
                        >
                            <div id="table_toggle">
                                {renderTableButtons()}
                            </div>
                            {renderTable()}
                        </div>
                        <div
                            style={{
                                display: activeTab === 'plot' ? 'block' : 'none',
                                height: 'calc(100vh)',
                                width: '100%',
                                overflow: 'hidden'
                            }}
                        >
                            <div id="plot_toggle">
                                {renderPlotButtons()}
                            </div>
                            {renderPlot()}
                        </div>
                    </div>
                </div>
            </div>
            <button onClick={runAnalysis} disabled={isLoading}>
                Run Analysis
            </button>
            {isLoading && <ProgressBar progress={progress} />}
        </div>
    );
};

export default Analysis;
