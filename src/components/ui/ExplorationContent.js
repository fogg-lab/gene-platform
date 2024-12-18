import React, { useState } from 'react';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';

const ExplorationContent = ({
    data,
    activeTab,
    setActiveTab,
    isLoading,
    progress,
    renderTable,
    currentTable,
    setCurrentTable,
    dataset
}) => {
    const [currentPlot, setCurrentPlot] = useState('pca');

    const renderInitialGuide = () => (
        <div className='analysisContentGuide'>
            <h2>First, load a dataset via 1 of 3 options:</h2>
            <div className="guide-options">
                <div className="guide-option">
                    <h3>1. Example dataset</h3>
                    <p>TCGA-CESC Stage I vs Stage III cervical cancer</p>
                </div>
                <div className="guide-option">
                    <h3>2. GEO/GDC</h3>
                    <p>Load bulk RNA-Seq data from GEO or GDC</p>
                </div>
                <div className="guide-option">
                    <h3>3. Uploaded CSV files</h3>
                    <div className="tip">
                        <span className="tip-label">Tip:</span>
                        <p>
                            If you're unsure how the tables should be structured, 
                            use the example dataset and then use the export button on the 
                            top toolbar above the displayed coldata and counts tables
                        </p>
                    </div>
                </div>
            </div>
        </div>
    );

    const renderPlotTabs = () => {
        if (!data || !data.plots) {
            return <p>Run analysis to see plots</p>;
        }

        const availablePlots = Object.keys(data.plots);

        return (
            <div className="plot-container">
                <div className="plot-tabs">
                    {availablePlots.map(plot => (
                        <button
                            key={plot}
                            className={`plot-tab ${currentPlot === plot ? 'active' : ''}`}
                            onClick={() => setCurrentPlot(plot)}
                        >
                            {plot.toUpperCase()}
                        </button>
                    ))}
                </div>
                <div className="plot-content">
                    {data.plots[currentPlot] ? (
                        <PlotArea htmlContent={data.plots[currentPlot]} />
                    ) : (
                        <p>Run analysis to see plots</p>
                    )}
                </div>
            </div>
        );
    };

    const renderTableTabs = () => {
        const availableTables = ['coldata', 'counts'];
        if (data && data.tables && data.tables.transformed_counts) {
            availableTables.push('transformed_counts');
        }
        return (
            <div id="table_subtabs">
                {availableTables.map(tbl => (
                    <button
                        key={tbl}
                        className={`view-toggle-btn ${currentTable === tbl ? 'active' : ''}`}
                        onClick={() => setCurrentTable(tbl)}
                    >
                        {tbl === 'transformed_counts' ? 'TRANSFORMED COUNTS' : tbl.toUpperCase()}
                    </button>
                ))}
            </div>
        );
    };

    return (
        <div className="exploration-content">
            <div>
                {!dataset ? (
                    renderInitialGuide()
                ) : (
                    <>
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
                        {activeTab === 'table' && renderTableTabs()}
                        {isLoading && <ProgressBar progress={progress} />}
                        <div className={`table-view ${activeTab === 'table' ? 'active' : ''}`}>
                            {renderTable()}
                        </div>
                        <div className={`plot-view ${activeTab === 'plot' ? 'active' : ''}`}>
                            {renderPlotTabs()}
                        </div>
                    </>
                )}
            </div>
        </div>
    );
};

export default ExplorationContent;
