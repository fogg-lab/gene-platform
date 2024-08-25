import React, { useState, useEffect, useRef } from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import DataTable from '../components/ui/DataTable';
import DatabasePopup from '../components/ui/DatabasePopup';
import Papa from 'papaparse';
import { getPublicUrl } from '../utils/environment';

const Analysis = () => {
    const [activeTab, setActiveTab] = useState('table');
    const [tableScrollPosition, setTableScrollPosition] = useState(0);
    const [currentStage, setCurrentStage] = useState('exploration');
    const [error, setError] = useState(null);
    const [isVisible, setIsVisible] = useState(false);
    const tableContainerRef = useRef(null);
    const [dataset, setDataset] = useState(null);
    const [edaData, setEdaData] = useState(null);
    const [deData, setDeData] = useState(null);
    const [gseaData, setGseaData] = useState(null);
    const [currentTable, setCurrentTable] = useState(null);
    const [currentPlot, setCurrentPlot] = useState(null);
    const [shouldDisplayPlot, setShouldDisplayPlot] = useState(false);

    const file_to_display_name = {
        'coldata.csv': 'Sample Metadata',
        'DE_results.csv': 'Differential Expression',
        'GSEA_results.csv': 'Gene Set Enrichment',
        'pca_3d.html': 'PCA 3D Embedding',
        'umap_3d.html': 'UMAP 3D Embedding',
        'sample_correlation_heatmap.html': 'Sample Correlation Heatmap',
        'mean_difference.html': 'Mean Difference Plot',
        'volcano_plot.html': 'Volcano Plot',
        'gene_concept_network.html': 'Gene Concept Network'
    };

    const stages = {
        exploration: {
            tables: ['coldata.csv'],
            plots: ['pca_3d.html', 'umap_3d.html', 'sample_correlation_heatmap.html']
        },
        differential: {
            tables: ['DE_results.csv'],
            plots: ['mean_difference.html', 'volcano_plot.html']
        },
        enrichment: {
            tables: ['GSEA_results.csv'],
            plots: ['gene_concept_network.html']
        }
    };

    useEffect(() => {
        loadTableData(currentTable);
    }, [currentTable]);

    useEffect(() => {
        if (activeTab === 'plot' && currentPlot) {
            setShouldDisplayPlot(true);
        } else {
            setShouldDisplayPlot(false);
        }
    }, [activeTab, currentPlot]);

    useEffect(() => {
        if (activeTab === 'table' && tableContainerRef.current) {
            tableContainerRef.current.scrollTop = tableScrollPosition;
        }
    }, [activeTab, tableScrollPosition]);
    const handleTableScroll = (event) => {
        setTableScrollPosition(event.target.scrollTop);
    };

    const loadTableData = (filename) => {
        if (dataset && (filename === 'coldata' || filename === 'counts')) {
            setEdaData(prevData => ({
                ...prevData,
                tables: {
                    ...(prevData?.tables || {}),
                    [filename]: dataset[`${filename}Table`]
                }
            }));
            setCurrentTable(filename);
        } else {
            // Existing logic for loading example data
            fetch(`${getPublicUrl()}/data/${filename}`)
                .then(response => {
                    if (!response.ok) {
                        throw new Error(`HTTP error! status: ${response.status}`);
                    }
                    return response.text();
                })
                .then(csvString => {
                    const result = Papa.parse(csvString, { header: true });
                    if (result.data.length > 0) {
                        setEdaData(prevData => ({
                            ...prevData,
                            tables: {
                                ...(prevData?.tables || {}),
                                [filename]: {
                                    data: result.data,
                                    cols: Object.keys(result.data[0])
                                }
                            }
                        }));
                        setCurrentTable(filename);
                    } else {
                        setError('CSV file is empty or could not be parsed correctly.');
                    }
                })
                .catch(error => {
                    console.error('Error loading CSV:', error);
                    setError(`Failed to load CSV: ${error.message}`);
                });
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

    const runDifferentialExpression = async () => {
        // Call API to run DE analysis
        // Update deData state with results
        // setDeData(results);
    };

    const runGSEA = async () => {
        // Call API to run GSEA
        // Update gseaData state with results
        // setGseaData(results);
    };

    const handleTabChange = (tab) => {
        setActiveTab(tab);
    };

    const renderTableButtons = () => {
        let tables = [];
        switch (currentStage) {
            case 'exploration':
                tables = edaData && edaData.tables ? Object.keys(edaData.tables) : [];
                break;
            case 'differential':
                tables = deData && deData.table ? ['de_results'] : [];
                break;
            case 'enrichment':
                tables = gseaData && gseaData.table ? ['gsea_results'] : [];
                break;
        }
        return tables.map(table => (
            <button
                key={table}
                onClick={() => setCurrentTable(table)}
                className={`view-toggle-btn ${currentTable === table ? 'active' : ''}`}
            >
                {file_to_display_name[table] || table}
            </button>
        ));
    };

    const renderPlotButtons = () => {
        let plots = [];
        switch (currentStage) {
            case 'exploration':
                plots = edaData && edaData.plots ? Object.keys(edaData.plots) : [];
                break;
            case 'differential':
                plots = deData && deData.plots ? Object.keys(deData.plots) : [];
                break;
            case 'enrichment':
                plots = gseaData && gseaData.plots ? Object.keys(gseaData.plots) : [];
                break;
        }
        return plots.map(plot => (
            <button
                key={plot}
                onClick={() => setCurrentPlot(plot)}
                className={`view-toggle-btn ${currentPlot === plot ? 'active' : ''}`}
            >
                {file_to_display_name[plot] || plot}
            </button>
        ));
    };

    const handleDatasetSelect = (type, data) => {
        if (type === 'external' && data) {
            setDataset(data);
            setEdaData({
                tables: {
                    coldata: {
                        data: data.coldataTable.data,
                        cols: data.coldataTable.cols
                    },
                    counts: {
                        data: data.countsTable.data,
                        cols: data.countsTable.cols
                    }
                },
                plots: data.plots || {}
            });
            setCurrentTable('coldata');
            setCurrentPlot(null);
        } else if (type === 'example') {
            // Load example dataset
            loadTableData('coldata.csv');
        }
    };

    const renderTable = () => {
        let tableData, tableColumns;
        switch (currentStage) {
            case 'exploration':
                if (edaData && edaData.tables && currentTable && edaData.tables[currentTable]) {
                    tableData = edaData.tables[currentTable].data;
                    tableColumns = edaData.tables[currentTable].cols.map(col => ({ key: col, name: col }));
                }
                break;
            case 'differential':
                if (deData && deData.table) {
                    tableData = deData.table.data;
                    tableColumns = deData.table.cols.map(col => ({ key: col, name: col }));
                }
                break;
            case 'enrichment':
                if (gseaData && gseaData.table) {
                    tableData = gseaData.table.data;
                    tableColumns = gseaData.table.cols.map(col => ({ key: col, name: col }));
                }
                break;
        }
        
        if (tableData && tableColumns) {
            // Ensure tableData is an array of objects
            if (Array.isArray(tableData[0])) {
                tableData = tableData.map(row => {
                    let obj = {};
                    tableColumns.forEach((col, index) => {
                        obj[col.key] = row[index];
                    });
                    return obj;
                });
            }
            return <DataTable data={tableData} columns={tableColumns} />;
        } else {
            return <p>No data available</p>;
        }
    };

    const renderPlot = () => {
        let plotData;
        switch (currentStage) {
            case 'exploration':
                plotData = edaData && currentPlot ? edaData.plots[currentPlot] : null;
                break;
            case 'differential':
                plotData = deData && currentPlot ? deData.plots[currentPlot] : null;
                break;
            case 'enrichment':
                plotData = gseaData && currentPlot ? gseaData.plots[currentPlot] : null;
                break;
        }
        return plotData ? (
            <PlotArea data={plotData.data} layout={plotData.layout} config={plotData.config} />
        ) : (
            <p>No plot available</p>
        );
    };

    if (error) {
        return <div>Error: {error}</div>;
    }

    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                <AnalysisInputForm setIsVisible={setIsVisible} onDatasetSelect={handleDatasetSelect} />
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
        </div >
    );
};

export default Analysis;
