import React, { useEffect, useState } from 'react';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';
import GeneSetFilterPopup from './GeneSetFilterPopup';

const DifferentialExpressionContent = ({
    data,
    activeTab,
    setActiveTab,
    isLoading,
    progress,
    renderTable,
    currentTable,
    setCurrentTable,
    currentPlot,
    setCurrentPlot,
    dataset,
}) => {
    useEffect(() => {
        if (data && !activeTab) {
            setActiveTab('table');
        }
    }, [data, activeTab, setActiveTab]);

    useEffect(() => {
        if (data?.plots && Object.keys(data.plots).length > 0 && !currentPlot) {
            setCurrentPlot(Object.keys(data.plots)[0]);
        }
    }, [data, currentPlot, setCurrentPlot]);

    const [showGeneSetFilter, setShowGeneSetFilter] = useState(false);
    const [activeFilter, setActiveFilter] = useState(null);

    const handleGeneSetFilter = (geneSetId, genes) => {
        setActiveFilter({ id: geneSetId, genes });
    };

    const renderTableTabs = () => {
        const availableTables = ['coldata', 'counts'];
        if (data && data.table) {
            availableTables.push('de_results');
        }
        return (
            <div id="table_subtabs">
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                    <div>
                        {availableTables.map(tbl => (
                            <button
                                key={tbl}
                                className={`view-toggle-btn ${currentTable === tbl ? 'active' : ''}`}
                                onClick={() => setCurrentTable(tbl)}
                            >
                                {tbl.toUpperCase()}
                            </button>
                        ))}
                    </div>
                    {currentTable === 'de_results' && (
                        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                            {activeFilter && (
                                <span style={{ fontSize: '0.9em', color: '#666' }}>
                                    Filtered by: {activeFilter.id}
                                    <button
                                        onClick={() => setActiveFilter(null)}
                                        style={{
                                            marginLeft: '8px',
                                            border: 'none',
                                            background: 'none',
                                            cursor: 'pointer',
                                            color: '#999'
                                        }}
                                    >
                                        ×
                                    </button>
                                </span>
                            )}
                            <button
                                className="filter-btn"
                                onClick={() => setShowGeneSetFilter(true)}
                                title="Filter table rows using gene sets from MSigDB"
                            >
                                Filter by Gene Set
                            </button>
                        </div>
                    )}
                </div>
            </div>
        );
    };

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
                        <p>Select plot. If no plot is available, run analysis.</p>
                    )}
                </div>
            </div>
        );
    };

    return (
        <div className="exploration-content">
            <div>
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
                    {renderTable({ activeFilter })}
                </div>
                <div className={`plot-view ${activeTab === 'plot' ? 'active' : ''}`}>
                    {renderPlotTabs()}
                </div>
            </div>
            <GeneSetFilterPopup
                isVisible={showGeneSetFilter}
                setIsVisible={setShowGeneSetFilter}
                onFilterApply={handleGeneSetFilter}
                species={dataset?.species || 'human'}
            />
        </div>
    );
}

export default DifferentialExpressionContent;