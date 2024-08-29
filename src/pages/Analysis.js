import React, { useState, useEffect, useCallback, useRef } from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
import WorkerManager from '../utils/Workers';
import PlotArea from '../components/ui/PlotArea';
import ProgressBar from '../components/ui/ProgressBar';

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

    const [selectedSamples, setSelectedSamples] = useState([]);
    const [contrastGroup, setContrastGroup] = useState({ samples: [] });
    const [referenceGroup, setReferenceGroup] = useState({ samples: [] });
    const [groupCounter, setGroupCounter] = useState(1);

    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);

    const handleSelectionChange = useCallback((newSelectedRows) => {
        console.log('Analysis - New selected rows:', newSelectedRows);
        setSelectedSamples(newSelectedRows);
    }, []);

    const handleAddGroup = useCallback((isContrast) => {
        const newGroup = {
            id: groupCounter,
            name: isContrast ? "Contrast Group" : "Reference Group",
            samples: []
        };
        if (isContrast) {
            setContrastGroup(newGroup);
        } else {
            setReferenceGroup(newGroup);
        }
        setGroupCounter(prevCounter => prevCounter + 1);
    }, [groupCounter]);

    const handleUpdateGroup = useCallback((groupId, updates, isContrast) => {
        const updateGroup = (group) => ({ ...group, ...updates });
        if (isContrast) {
            setContrastGroup(updateGroup);
        } else {
            setReferenceGroup(updateGroup);
        }
    }, []);

    const handleAddSamplesToGroup = useCallback((groupId, isContrast, samplesToAdd) => {
        console.log("Adding samples to group:", { groupId, isContrast, samplesToAdd });

        if (isContrast) {
            setContrastGroup(prevGroup => ({
                ...prevGroup,
                samples: [...prevGroup.samples, ...samplesToAdd]
            }));
        } else {
            setReferenceGroup(prevGroup => ({
                ...prevGroup,
                samples: [...prevGroup.samples, ...samplesToAdd]
            }));
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

    useEffect(() => {
        console.log('Analysis - Updated selectedSamples:', selectedSamples);
    }, [selectedSamples]);

    const handleTableScroll = (event) => {
        setTableScrollPosition(event.target.scrollTop);
    };

    const handleDatasetSelect = (type, data) => {
        if (type === 'external' || type === 'upload') {
            setDataset(data);
        } else if (type === 'example') {
            // Load example dataset (you may need to implement this)
            loadExampleDataset();
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
                onSelectionChange={handleSelectionChange}
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

        try {
            // EDA
            setProgress(10);
            const transformedCounts = await WorkerManager.runTask('py', 'transform_vst', {
                counts: dataset.counts
            });
            setProgress(20);
            const pcaPlot = await WorkerManager.runTask('py', 'create_pca', {
                counts: transformedCounts,
                sample_ids: dataset.coldataTable.rows
            });
            setProgress(27);
            const tsnePlot = await WorkerManager.runTask('py', 'create_tsne', {
                counts: transformedCounts,
                sample_ids: dataset.coldataTable.rows
            });
            setProgress(34);
            const heatmap = await WorkerManager.runTask('py', 'create_heatmap', {
                counts: transformedCounts,
                sample_ids: dataset.coldataTable.rows
            });

            // DE Analysis
            setProgress(40);
            const deResults = await WorkerManager.runTask('r', 'run_de_analysis', {
                counts: transformedCounts,
                coldata: dataset.coldataTable,
                contrastGroup,
                referenceGroup
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
        } catch (error) {
            console.error('Analysis error:', error);
            setError(`An error occurred during analysis: ${error.message}`);
        } finally {
            setIsLoading(false);
        }
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
                    onAddGroup={handleAddGroup}
                    onUpdateGroup={handleUpdateGroup}
                    selectedSamples={selectedSamples}
                    onAddSamplesToGroup={handleAddSamplesToGroup}
                    runAnalysis={runAnalysis}
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
