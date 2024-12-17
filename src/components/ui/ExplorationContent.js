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
}) => {
    const [currentPlot, setCurrentPlot] = useState('pca');

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

    return (
        <div className="exploration-content">
            <div id="view_toggle">
                <button
                    className={`view-toggle-btn ${activeTab === 'table' ? 'active' : ''}`}
                    onClick={() => setActiveTab('table')}
                >
                    Table View
                </button>
                <button
                    className={`view-toggle-btn ${activeTab === 'plot' ? 'active' : ''}`}
                    onClick={() => {
                        setActiveTab('plot');
                        // showError("Huge error lay ahead");
                    }}
                >
                    Plot View
                </button>
            </div>
            {isLoading && <ProgressBar progress={progress} />}
            <div className={`table-view ${activeTab === 'table' ? 'active' : ''}`}>
                {renderTable()}
            </div>
            <div className={`plot-view ${activeTab === 'plot' ? 'active' : ''}`}>
                {renderPlotTabs()}
            </div>
        </div>
    );
}

export default ExplorationContent;
