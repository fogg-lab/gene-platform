import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { DataGrid, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector } from '@mui/x-data-grid';
import PropTypes from 'prop-types';

const CustomToolbar = ({ selectedSamples, contrastGroup, referenceGroup, onAddSamplesToGroup }) => {
  const [selectedGroupType, setSelectedGroupType] = useState('reference');

  const handleAddSamples = () => {
    const groupId = selectedGroupType === 'contrast' ? 1 : 2; // Assuming 1 for contrast, 2 for reference
    onAddSamplesToGroup(groupId, selectedGroupType === 'contrast', selectedSamples);
  };

  return (
    <GridToolbarContainer>
      <GridToolbarColumnsButton />
      <GridToolbarFilterButton />
      <GridToolbarDensitySelector />
      <GridToolbarExport />
      <select
        value={selectedGroupType}
        onChange={(e) => setSelectedGroupType(e.target.value)}
      >
        <option value="contrast">Contrast Group</option>
        <option value="reference">Reference Group</option>
      </select>
      <button
        onClick={handleAddSamples}
        disabled={selectedSamples.length === 0}
      >
        Add Selected ({selectedSamples.length})
      </button>
    </GridToolbarContainer>
  );
};

const DataTable = ({ data, columns, onSelectionChange, contrastGroup, referenceGroup, onAddSamplesToGroup }) => {
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
      const newSelectionModel = selectionModel.length > 0 ? [] : filteredRows.map(row => row.id);
      setSelectionModel(newSelectionModel);
      updateSelectedRows(newSelectionModel);
    }
  }, [selectionModel, filteredRows, updateSelectedRows]);

  const handleSelectionModelChange = useCallback((newSelectionModel) => {
    setSelectionModel(newSelectionModel);
    const selectedRows = data.filter(row => newSelectionModel.includes(row.id));
    onSelectionChange(selectedRows);
  }, [data, onSelectionChange]);


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
            contrastGroup,
            referenceGroup,
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
          // New styles for checkboxes
          '& .MuiCheckbox-root': {
            color: '#D73F09',
          },
          '& .MuiCheckbox-root.Mui-checked': {
            color: '#D73F09',
          },
          '& .MuiCheckbox-root:hover': {
            backgroundColor: 'rgba(215, 63, 9, 0.04)',
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
  contrastGroup: PropTypes.object.isRequired,
  referenceGroup: PropTypes.object.isRequired,
  onAddSamplesToGroup: PropTypes.func.isRequired,
};

export default DataTable;
