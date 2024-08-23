import React, { useState, useMemo } from 'react';
import { DataGridPro, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector } from '@mui/x-data-grid-pro';
import PropTypes from 'prop-types';

const CustomToolbar = () => {
  return (
    <GridToolbarContainer>
      <GridToolbarColumnsButton />
      <GridToolbarFilterButton />
      <GridToolbarDensitySelector />
      <GridToolbarExport />
    </GridToolbarContainer>
  );
};

const DataTable = ({ data, columns }) => {
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
    // You can add additional logic here, e.g., callback to parent component
  };

  return (
    <div style={{
      height: '80vh',
      width: '100%',
      maxHeight: 'calc(100vh - 20px)',
      overflow: 'hidden'
    }}>
      <DataGridPro
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
};

export default DataTable;