import React, { useState, useMemo } from 'react';
import { DataGrid, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector } from '@mui/x-data-grid';
import PropTypes from 'prop-types';

const CustomToolbar = ({
  selectedSamples,
  contrastGroups,
  referenceGroups,
  onAddSamplesToGroup
}) => {
  const [selectedGroupType, setSelectedGroupType] = useState('contrast');
  const [selectedGroupId, setSelectedGroupId] = useState('');

  const handleAddSamples = () => {
    if (selectedGroupId) {
      onAddSamplesToGroup(parseInt(selectedGroupId), selectedGroupType === 'contrast');
    }
  };

  return (
    <GridToolbarContainer>
      <GridToolbarColumnsButton />
      <GridToolbarFilterButton />
      <GridToolbarDensitySelector />
      <GridToolbarExport />
      <div style={{ marginLeft: 'auto', display: 'flex', alignItems: 'center', gap: '10px' }}>
        <select
          value={selectedGroupType}
          onChange={(e) => {
            setSelectedGroupType(e.target.value);
            setSelectedGroupId('');
          }}
          style={{ padding: '5px' }}
        >
          <option value="contrast">Contrast Groups</option>
          <option value="reference">Reference Groups</option>
        </select>
        <select
          value={selectedGroupId}
          onChange={(e) => setSelectedGroupId(e.target.value)}
          disabled={!contrastGroups.length && !referenceGroups.length}
          style={{ padding: '5px' }}
        >
          <option value="">Select a group</option>
          {(selectedGroupType === 'contrast' ? contrastGroups : referenceGroups).map(group => (
            <option key={group.id} value={group.id}>{group.name}</option>
          ))}
        </select>
        <button
          onClick={handleAddSamples}
          disabled={!selectedGroupId || !selectedSamples.length}
          style={{ padding: '5px 10px', cursor: 'pointer' }}
        >
          Add Selected ({selectedSamples.length})
        </button>
      </div>
    </GridToolbarContainer>
  );
};

CustomToolbar.propTypes = {
  selectedSamples: PropTypes.array.isRequired,
  contrastGroups: PropTypes.array.isRequired,
  referenceGroups: PropTypes.array.isRequired,
  onAddSamplesToGroup: PropTypes.func.isRequired,
};

const DataTable = ({
  data,
  columns,
  onSelectionChange,
  contrastGroups,
  referenceGroups,
  onAddSamplesToGroup
}) => {
  const [sortModel, setSortModel] = useState([]);
  const [selectionModel, setSelectionModel] = useState([]);
  const [filterModel, setFilterModel] = useState({
    items: [],
  });

  const gridColumns = useMemo(() => {
    return columns.map(col => ({
      field: col.key,
      headerName: col.name,
      width: 150,
      sortable: true,
      resizable: true,
      filterable: true,
    }));
  }, [columns]);

  const filteredRows = useMemo(() => {
    return data.filter(row => {
      return Object.values(row).some(value =>
        value !== null && value !== undefined && value !== ''
      );
    }).map((row, index) => ({
      ...row,
      id: row.id || index, // Ensure each row has a unique id
    }));
  }, [data]);

  const handleSelectionModelChange = (newSelectionModel) => {
    setSelectionModel(newSelectionModel);
    const selectedRows = filteredRows.filter(row => newSelectionModel.includes(row.id));
    onSelectionChange(selectedRows);
  };

  return (
    <div style={{
      height: '80vh',
      width: '100%',
      maxHeight: 'calc(100vh - 20px)',
      overflow: 'hidden'
    }}>
      <DataGrid
        rows={filteredRows}
        columns={gridColumns}
        sortModel={sortModel}
        onSortModelChange={(newSortModel) => setSortModel(newSortModel)}
        filterModel={filterModel}
        onFilterModelChange={(newFilterModel) => setFilterModel(newFilterModel)}
        checkboxSelection
        disableColumnMenu={false}
        disableSelectionOnClick
        density="compact"
        getRowId={(row) => row.id}
        selectionModel={selectionModel}
        onSelectionModelChange={handleSelectionModelChange}
        pagination
        paginationMode="client"
        rowsPerPageOptions={[25, 50, 100]}
        initialState={{
          pagination: {
            pageSize: 100,
          },
        }}
        components={{
          Toolbar: CustomToolbar,
        }}
        componentsProps={{
          toolbar: {
            selectedSamples: selectionModel,
            contrastGroups,
            referenceGroups,
            onAddSamplesToGroup,
          },
        }}
        sx={{
          height: '100%',
          width: '100%',
          '& .MuiDataGrid-main': { overflow: 'hidden' },
          '& .MuiDataGrid-virtualScroller': { overflow: 'auto' },
          '& .MuiDataGrid-footerContainer': { overflow: 'hidden' },
          '& .MuiDataGrid-toolbarContainer': {
            padding: '8px',
            justifyContent: 'space-between',
          },
          '& .MuiDataGrid-toolbarContainer .MuiButton-root': {
            color: '#D73F09',
          },
          '& .MuiDataGrid-toolbarContainer .MuiButton-root .MuiSvgIcon-root': {
            color: '#D73F09',
          },
          '& .MuiDataGrid-toolbarContainer .MuiButton-root .MuiButton-startIcon + .MuiButton-text': {
            color: '#D73F09',
          },
          '& .MuiDataGrid-columnHeader .MuiDataGrid-iconButtonContainer .MuiIconButton-root': {
            color: '#D73F09',
          },
          '& .MuiDataGrid-columnHeader .MuiDataGrid-filterIcon': {
            color: '#D73F09',
          },
        }}
      />
    </div>
  );
};

DataTable.propTypes = {
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  columns: PropTypes.arrayOf(
    PropTypes.shape({
      key: PropTypes.string.isRequired,
      name: PropTypes.string.isRequired,
    })
  ).isRequired,
  onSelectionChange: PropTypes.func.isRequired,
  contrastGroups: PropTypes.array.isRequired,
  referenceGroups: PropTypes.array.isRequired,
  onAddSamplesToGroup: PropTypes.func.isRequired,
};

export default DataTable;