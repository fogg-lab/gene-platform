import React, { useState } from 'react';
import DataTable from './DataTable';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';

const ExplorationContent = ({
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
        if (!data || !data.tables || !data.tables.coldata) {
            return (
                <div className='analysisContentGuide'>
                    <h1>To run Exploratory Data Analysis:</h1>
                    <p>Step 1. Choose an example dataset, or select one from GDC / GEO, or input a custom dataset.</p>
                    <p>Step 2. Select desired data transform.</p>
                    <p>Step 3. Run</p>
                </div>);
        }

        const tableData = data.tables.coldata.data.map((row, index) => {
            const obj = {};
            data.tables.coldata.cols.forEach((col, colIndex) => {
                obj[col] = row[colIndex];
            });
            obj.id = index;
            return obj;
        });

        const columns = data.tables.coldata.cols.map(col => ({ key: col, name: col }));

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
    };

    const renderPlotTabs = () => {
        if (!data || !data.plots) {
            return <p>No plots available</p>;
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
    );
}

export default ExplorationContent;