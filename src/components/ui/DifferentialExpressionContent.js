import React, { useEffect } from 'react';
import DataTable from './DataTable';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';

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
}) => {
    const renderPlotTabs = () => {
        if (!data || !data.plots) {
            return <p>No plots available</p>;
        }

        const availablePlots = Object.keys(data.plots);

        useEffect(() => {
            if (data && data.plots && Object.keys(data.plots).length > 0) {
                const firstPlot = Object.keys(data.plots)[0];
                if (currentPlot !== firstPlot) {
                    setCurrentPlot(firstPlot);
                }
            }
        }, [data, currentPlot, setCurrentPlot]);

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
                {isLoading && <ProgressBar progress={progress} />}
                <div className={`table-view ${activeTab === 'table' ? 'active' : ''}`}>
                    {renderTable()}
                </div>
                <div className={`plot-view ${activeTab === 'plot' ? 'active' : ''}`}>
                    {renderPlotTabs()}
                </div>
            </div>
        </div>
    );
}

export default DifferentialExpressionContent;