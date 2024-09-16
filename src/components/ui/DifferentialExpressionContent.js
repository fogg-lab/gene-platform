import React, { useRef } from 'react';
import DataTable from './DataTable';
import ProgressBar from './ProgressBar';
import PlotArea from './PlotArea';

const DifferentialExpressionContent = ({
    data,
    activeTab,
    onAddSamplesToGroup,
    onRemoveSamplesFromGroup,
    contrastGroup,
    referenceGroup,
    isLoading,
    progress
}) => {
    const tableContainerRef = useRef(null);

    const renderTable = () => {
        if (!data || !data.tables || !data.tables.coldata) {
            return <p>No data available</p>;
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

    const renderPlot = () => {
        if (!data || !data.plots || !data.plots.pca) {
            return <p>No plot available</p>;
        }
        return <PlotArea htmlContent={data.plots.pca} />;
    };

    return (
        <div id="view_content" style={{ height: '100%', width: '100%', overflow: 'hidden' }}>
            {isLoading && <ProgressBar progress={progress} />}
            <div
                style={{
                    display: activeTab === 'table' ? 'block' : 'none',
                    height: '100%',
                    overflow: 'auto'
                }}
                ref={tableContainerRef}
            >
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
                {renderPlot()}
            </div>
        </div>
    );
}

export default DifferentialExpressionContent;