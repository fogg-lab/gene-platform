import React, { useState } from 'react';
import DataTable from './DataTable';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';
import { formatNumber } from '../../lib/utils';

const GSEAContent = ({
    data,
    activeTab,
    setActiveTab,
    onAddSamplesToGroup,
    onRemoveSamplesFromGroup,
    contrastGroup,
    referenceGroup,
    isLoading,
    progress
}) => {
    const [currentPlot, setCurrentPlot] = useState('pca');

    const renderTable = () => {
        if (data && data.tables && data.tables.results) {
            const rawData = data.tables.results.data;
            const cols = data.tables.results.cols;
            const numericColumns = ['PValue', 'FDR'];
            const tableData = rawData.map((row, index) => {
                const obj = {};
                cols.forEach((col, colIndex) => {
                    obj[col] = numericColumns.includes(col) ? formatNumber(row[colIndex]) : row[colIndex];
                });
                obj.id = index;
                return obj;
            });

            const columns = data.tables.results.cols.map(col => ({ key: col, name: col }));

            return (
                <DataTable
                    data={tableData}
                    columns={columns}
                    contrastGroup={contrastGroup}
                    referenceGroup={referenceGroup}
                    onAddSamplesToGroup={onAddSamplesToGroup}
                    onRemoveSamplesFromGroup={onRemoveSamplesFromGroup}
                />
            );
        }

        return (
            <div className='analysisContentGuide'>
                <h1>To run Gene Set Enrichment Analysis:</h1>
                <p>Step 1. Add a gene set</p>
                <p>(Optional) Change standard configuration.</p>
                <p>Step 2. Run</p>
            </div>
        );
    };

    const renderPlotTabs = () => {
        if (!data || !data.plots) {
            return <p>No plots available 🚧</p>;
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
                        <p>No {currentPlot} plot available</p>
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

export default GSEAContent;