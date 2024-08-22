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
    const [tableData, setTableData] = useState([]);
    const [tableColumns, setTableColumns] = useState([]);
    const [error, setError] = useState(null);
    const [currentTable, setCurrentTable] = useState('coldata.csv');
    const [currentPlot, setCurrentPlot] = useState('pca_3d.html');
    const [shouldDisplayPlot, setShouldDisplayPlot] = useState(false);
    const [isVisible, setIsVisible] = useState(false);
    const [setSelectedRadio] = useState(null);
    const tableContainerRef = useRef(null);

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
        // Load plot when switching to 'plot' tab or when changing stage while in 'plot' tab
        if (activeTab === 'plot') {
            setShouldDisplayPlot(true);
        } else {
            setShouldDisplayPlot(false);
        }
    }, [activeTab, currentStage]);

    useEffect(() => {
        if (activeTab === 'table' && tableContainerRef.current) {
            tableContainerRef.current.scrollTop = tableScrollPosition;
        }
    }, [activeTab, tableScrollPosition]);
    const handleTableScroll = (event) => {
        setTableScrollPosition(event.target.scrollTop);
    };

    const loadTableData = (filename) => {
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
                    setTableData(result.data);
                    setTableColumns(Object.keys(result.data[0]).map(field => ({ key: field, name: field })));
                } else {
                    setError('CSV file is empty or could not be parsed correctly.');
                }
            })
            .catch(error => {
                console.error('Error loading CSV:', error);
                setError(`Failed to load CSV: ${error.message}`);
            });
    };

    const handleStageChange = (stage) => {
        setCurrentStage(stage);
        setCurrentTable(stages[stage].tables[0]);
        setCurrentPlot(stages[stage].plots[0]);
    };

    const handleTabChange = (tab) => {
        setActiveTab(tab);
    };

    const renderTableButtons = () => {
        return stages[currentStage].tables.map(table => (
            <button
                key={table}
                onClick={() => setCurrentTable(table)}
                className={`view-toggle-btn ${currentTable === table ? 'active' : ''}`}
            >
                {file_to_display_name[table]}
            </button>
        ));
    };

    const renderPlotButtons = () => {
        return stages[currentStage].plots.map(plot => (
            <button
                key={plot}
                onClick={() => setCurrentPlot(plot)}
                className={`view-toggle-btn ${currentPlot === plot ? 'active' : ''}`}
            >
                {file_to_display_name[plot]}
            </button>
        ));
    };

    if (error) {
        return <div>Error: {error}</div>;
    }

    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                <AnalysisInputForm setIsVisible={setIsVisible} />
            </div>
            <DatabasePopup setIsVisible={setIsVisible} isVisible={isVisible} setSelectedRadio={setSelectedRadio} />
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
                            onClick={() => handleTabChange('table')}
                        >
                            Table View
                        </button>
                        <button
                            className={`view-toggle-btn ${activeTab === 'plot' ? 'active' : ''}`}
                            onClick={() => handleTabChange('plot')}
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
                            <div
                                style={{
                                    position: 'relative',
                                    paddingBottom: '90%',
                                    height: 0,
                                    overflow: 'hidden'
                                }}
                            >
                                <div
                                    style={{
                                        position: 'absolute',
                                        top: 0,
                                        left: 0,
                                        width: '100%',
                                        height: '100%',
                                        border: 'none',
                                        padding: '10px'
                                    }}
                                >
                                    {tableData.length > 0 && tableColumns.length > 0 ? (
                                        <DataTable data={tableData} columns={tableColumns} />
                                    ) : (
                                        <p>Loading table data...</p>
                                    )}
                                </div>
                            </div>
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
                            <div style={{
                                position: 'relative',
                                paddingBottom: '56.25%', // 16:9 aspect ratio
                                height: 0,
                                overflow: 'hidden'
                            }}>
                                {shouldDisplayPlot && (
                                    <iframe
                                        src={`${getPublicUrl()}/plots/${currentPlot}`}
                                        style={{
                                            position: 'absolute',
                                            top: 0,
                                            left: 0,
                                            width: '100%',
                                            height: '100%',
                                            border: 'none'
                                        }}
                                        title="Analysis Plot"
                                    />
                                )}
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div >
    );
};

export default Analysis;
