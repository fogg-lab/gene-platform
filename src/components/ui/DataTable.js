import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { DataGrid, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector } from '@mui/x-data-grid';
import PropTypes from 'prop-types';

const CustomToolbar = ({ selectedSamples, contrastGroups, referenceGroups, onAddSamplesToGroup }) => {
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
      <select
        value={selectedGroupType}
        onChange={(e) => {
          setSelectedGroupType(e.target.value);
          setSelectedGroupId('');
        }}
      >
        <option value="contrast">Contrast Groups</option>
        <option value="reference">Reference Groups</option>
      </select>
      <select
        value={selectedGroupId}
        onChange={(e) => setSelectedGroupId(e.target.value)}
        disabled={!contrastGroups.length && !referenceGroups.length}
      >
        <option value="">Select a group</option>
        {(selectedGroupType === 'contrast' ? contrastGroups : referenceGroups).map(group => (
          <option key={group.id} value={group.id}>{group.name}</option>
        ))}
      </select>
      <button
        onClick={handleAddSamples}
        disabled={!selectedGroupId || selectedSamples.length === 0}
      >
        Add Selected ({selectedSamples.length})
      </button>
    </GridToolbarContainer>
  );
};

const DataTable = ({ data, columns, onSelectionChange, contrastGroups, referenceGroups, onAddSamplesToGroup }) => {
  const [sortModel, setSortModel] = useState([]);
  const [selectionModel, setSelectionModel] = useState([]);
  const [filterModel, setFilterModel] = useState({ items: [] });

  const gridColumns = useMemo(() => {
    return columns.map(col => ({
      field: col.key,
      headerName: col.name,
      width: 150,
      sortable: true,
      filterable: true,
      renderCell: (params) => {
        if (typeof params.value === 'object' && params.value !== null) {
          return JSON.stringify(params.value);
        }
        return params.value;
      }
    }));
  }, [columns]);

  const filteredRows = useMemo(() => {
    return data.filter(row => {
      return Object.values(row).some(value =>
        value !== null && value !== undefined && value !== ''
      );
    }).map((row, index) => ({
      ...row,
      id: row.id || index,
    }));
  }, [data]);

  const updateSelectedRows = useCallback((newSelectionModel) => {
    const newSelectedRows = filteredRows.filter(row => newSelectionModel.includes(row.id));
    onSelectionChange(newSelectedRows);
  }, [filteredRows, onSelectionChange]);

  const handleHeaderClick = useCallback((params) => {
    if (params.field === '__check__') {
      console.log("Select-all header clicked");
      const newSelectionModel = selectionModel.length === filteredRows.length ? [] : filteredRows.map(row => row.id);
      setSelectionModel(newSelectionModel);
      updateSelectedRows(newSelectionModel);
    }
  }, [selectionModel, filteredRows, updateSelectedRows]);

  const handleSelectionModelChange = useCallback((newSelectionModel) => {
    setSelectionModel(newSelectionModel);
    updateSelectedRows(newSelectionModel);
    console.log('Selection changed. New selection:', newSelectionModel);
    console.log('Number of selected rows:', newSelectionModel.length);
  }, [updateSelectedRows]);

  const handleCellClick = useCallback((params) => {
    const { id } = params;
    const newSelectionModel = selectionModel.includes(id)
      ? selectionModel.filter(rowId => rowId !== id)
      : [...selectionModel, id];

    setSelectionModel(newSelectionModel);
    updateSelectedRows(newSelectionModel);

    console.log('Clicked row data:', params.row);
    console.log('Clicked cell value:', params.value);
    console.log('Clicked column field:', params.field);
  }, [selectionModel, updateSelectedRows]);

  return (
    <div style={{ height: '80vh', width: '100%', maxHeight: 'calc(100vh - 20px)', overflow: 'hidden' }}>
      <DataGrid
        rows={filteredRows}
        columns={gridColumns}
        sortModel={sortModel}
        onSortModelChange={setSortModel}
        filterModel={filterModel}
        onFilterModelChange={setFilterModel}
        checkboxSelection
        disableColumnMenu={false}
        disableSelectionOnClick={true}
        density="compact"
        getRowId={(row) => row.id}
        onColumnHeaderClick={handleHeaderClick}
        selectionModel={selectionModel}
        onSelectionModelChange={handleSelectionModelChange}
        onCellClick={handleCellClick}
        pagination
        paginationMode="client"
        rowsPerPageOptions={[25, 50, 100]}
        initialState={{ pagination: { pageSize: 100 } }}
        components={{ Toolbar: CustomToolbar }}
        componentsProps={{
          toolbar: {
            selectedSamples: filteredRows.filter(row => selectionModel.includes(row.id)),
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