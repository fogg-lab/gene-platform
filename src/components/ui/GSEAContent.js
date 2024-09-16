import React from 'react';
import DataTable from './DataTable';
import ProgressBar from './ProgressBar';

const GSEAContent = ({
    data,
    activeTab,
    onAddSamplesToGroup,
    onRemoveSamplesFromGroup,
    contrastGroup,
    referenceGroup,
    isLoading,
    progress
}) => {
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
        return <div dangerouslySetInnerHTML={{ __html: data.plots.pca }} />;
    };

    return (
        <div>
            <p>Gene Set Enrichment Analysis</p>
            {isLoading && <ProgressBar progress={progress} />}
            <div style={{ display: activeTab === 'table' ? 'block' : 'none' }}>
                {renderTable()}
            </div>
            <div style={{ display: activeTab === 'plot' ? 'block' : 'none' }}>
                {renderPlot()}
            </div>
        </div>
    );
}

export default GSEAContent;