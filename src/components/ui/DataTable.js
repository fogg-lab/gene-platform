import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { DataGridPro, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector, gridExpandedSortedRowIdsSelector, useGridApiContext, } from '@mui/x-data-grid-pro';
import PropTypes from 'prop-types';

const getFilteredRows = ({ apiRef }) => gridExpandedSortedRowIdsSelector(apiRef);
const getIntersectingRows = (filteredRows, selectionModel) => {
  return selectionModel.filter(id => filteredRows.includes(id));
};

const CustomToolbar = ({ onAddSamplesToGroup, selectionModel, rows, clearSelection, requiresToolbar, contrastGroup, referenceGroup }) => {
  const apiRef = useGridApiContext();

  const basicControls = (
    <div style={{ display: 'flex', gap: '8px' }}>
      <GridToolbarColumnsButton />
      <GridToolbarFilterButton />
      <GridToolbarDensitySelector />
      <GridToolbarExport
        printOptions={{ disableToolbarButton: false }}
        csvOptions={{ disableToolbarButton: false }}
      />
    </div>
  );

  const handleAddSamples = (isContrast) => {
    const filteredSelectionModel = getIntersectingRows(getFilteredRows({ apiRef }), selectionModel);
    const selectedSamples = filteredSelectionModel.map(id => {
      const row = rows.find(r => r.id === id);
      return {
        id: row.sample_id || row.id,
        [Object.keys(row)[0]]: row[Object.keys(row)[0]]
      };
    });
    onAddSamplesToGroup(isContrast, selectedSamples);
    clearSelection();
  };

  const numSelectedFiltered = getIntersectingRows(getFilteredRows({ apiRef }), selectionModel).length;
  const referenceCount = referenceGroup.samples.length;
  const contrastCount = contrastGroup.samples.length;

  return (
    <GridToolbarContainer sx={{ 
      display: 'flex', 
      justifyContent: 'space-between', 
      padding: '8px 16px',
      borderBottom: '1px solid rgba(224, 224, 224, 1)'
    }}>
      {basicControls}
      {requiresToolbar && (
        <div style={{ 
          display: 'flex', 
          gap: '8px', 
          alignItems: 'center'
        }}>
          <span style={{ 
            fontFamily: 'KievitOffc',
            marginRight: '8px'
          }}>
            {numSelectedFiltered} samples selected
          </span>
          <button
            onClick={() => handleAddSamples(false)}
            disabled={numSelectedFiltered === 0}
            style={{
              backgroundColor: '#D73F09',
              color: 'white',
              border: 'none',
              padding: '8px 16px',
              borderRadius: '4px',
              cursor: numSelectedFiltered === 0 ? 'not-allowed' : 'pointer',
              opacity: numSelectedFiltered === 0 ? 0.6 : 1,
              fontFamily: 'KievitOffc',
              display: 'flex',
              flexDirection: 'column',
              alignItems: 'center',
              gap: '2px'
            }}
          >
            <span>Assign to Reference Group</span>
            <span style={{ fontSize: '0.8em' }}>({referenceCount} samples)</span>
          </button>
          <button
            onClick={() => handleAddSamples(true)}
            disabled={numSelectedFiltered === 0}
            style={{
              backgroundColor: '#D73F09',
              color: 'white',
              border: 'none',
              padding: '8px 16px',
              borderRadius: '4px',
              cursor: numSelectedFiltered === 0 ? 'not-allowed' : 'pointer',
              opacity: numSelectedFiltered === 0 ? 0.6 : 1,
              fontFamily: 'KievitOffc',
              display: 'flex',
              flexDirection: 'column',
              alignItems: 'center',
              gap: '2px'
            }}
          >
            <span>Assign to Contrast Group</span>
            <span style={{ fontSize: '0.8em' }}>({contrastCount} samples)</span>
          </button>
        </div>
      )}
    </GridToolbarContainer>
  );
};

CustomToolbar.propTypes = {
  onAddSamplesToGroup: PropTypes.func.isRequired,
  selectionModel: PropTypes.array.isRequired,
  rows: PropTypes.array.isRequired,
  clearSelection: PropTypes.func.isRequired,
  requiresToolbar: PropTypes.bool.isRequired,
  contrastGroup: PropTypes.object.isRequired,
  referenceGroup: PropTypes.object.isRequired,
};

const DataTable = ({
  data,
  columns,
  contrastGroup,
  referenceGroup,
  onAddSamplesToGroup,
  onRemoveSamplesFromGroup,
  requiresToolbar,
  activeFilter,
}) => {
  const [sortModel, setSortModel] = useState([]);
  const [selectionModel, setSelectionModel] = useState([]);
  const [filterModel, setFilterModel] = useState({ items: [] });
  const [key, setKey] = useState(0);
  const [rowGroups, setRowGroups] = useState({});

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
      id: row.id || row.sample_id || index,
      groupType: rowGroups[row.id || row.sample_id || index] || ''
    }));
  }, [data, rowGroups]);

  useEffect(() => {
    if (activeFilter) {
      const upperGenes = activeFilter.genes.map(gene => gene.toUpperCase());

      setFilterModel({
        items: [{
          field: 'symbol',
          operator: 'isAnyOf',
          value: upperGenes
        }]
      });

      // Debug check if genes exist in the data
      const foundGenes = data.filter(row => 
        upperGenes.includes(row.symbol?.toUpperCase())
      );
      console.log('Found matching genes:', foundGenes.length);
      console.log('Sample matching genes:', foundGenes.slice(0, 5));
    } else {
      setFilterModel({ items: [] });
    }
  }, [activeFilter, data]);

  useEffect(() => {
    if (activeFilter) {
      console.log('Active filter genes:', activeFilter.genes);
      console.log('Available rows:', filteredRows);
      console.log('Current filter model:', filterModel);
    }
  }, [activeFilter, filteredRows, filterModel]);

  const handleHeaderClick = useCallback((params) => {
    if (params.field === '__check__') {
      const newSelectionModel = selectionModel.length > 0 ? [] : filteredRows.map(row => row.id);
      setSelectionModel(newSelectionModel);
    }
  }, [selectionModel, filteredRows]);

  const handleCellClick = useCallback((params) => {
    const { id } = params;
    const newSelectionModel = selectionModel.includes(id)
      ? selectionModel.filter(rowId => rowId !== id)
      : [...selectionModel, id];

    setSelectionModel(newSelectionModel);
  }, [selectionModel]);

  const getRowClassName = useCallback((params) => {
    if (!requiresToolbar) return '';
    
    const sampleId = params.row.sample_id || params.id;
    if (contrastGroup.samples.some(sample => sample.id === sampleId)) {
      return 'contrast-row';
    } else if (referenceGroup.samples.some(sample => sample.id === sampleId)) {
      return 'reference-row';
    }
    return '';
  }, [contrastGroup, referenceGroup, requiresToolbar]);

  const clearSelection = useCallback(() => {
    console.log("clearSelection function called");
    setSelectionModel([]);
    setKey(prevKey => prevKey + 1);
    console.log("Selection cleared and DataGrid re-render triggered");
  }, []);

  return (
    <div style={{ height: '80vh', width: '100%', maxHeight: 'calc(100vh - 20px)', overflow: 'hidden' }}>
      <DataGridPro
        key={key}
        getRowClassName={getRowClassName}
        rows={filteredRows}
        columns={gridColumns}
        sortModel={sortModel}
        onSortModelChange={setSortModel}
        filterModel={filterModel}
        onFilterModelChange={(model) => {
          setFilterModel(model);
        }}
        checkboxSelection={requiresToolbar}
        disableColumnMenu={false}
        disableSelectionOnClick={false}
        density="compact"
        getRowId={(row) => row.id}
        onColumnHeaderClick={handleHeaderClick}
        selectionModel={selectionModel}
        onSelectionModelChange={(newSelectionModel) => {
          setSelectionModel(newSelectionModel);
        }}
        onCellClick={handleCellClick}
        pagination
        paginationMode="client"
        rowsPerPageOptions={[25, 50, 100]}
        initialState={{ pagination: { pageSize: 100 } }}
        slots={{
          toolbar: CustomToolbar,
        }}
        slotProps={{
          toolbar: {
            onAddSamplesToGroup,
            selectionModel,
            rows: filteredRows,
            clearSelection,
            requiresToolbar,
            contrastGroup,
            referenceGroup,
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
          '& .MuiCheckbox-root': {
            color: '#D73F09',
          },
          '& .MuiCheckbox-root.Mui-checked': {
            color: '#D73F09',
          },
          '& .MuiCheckbox-root:hover': {
            backgroundColor: 'rgba(215, 63, 9, 0.04)',
          },
          '& .reference-row': {
            backgroundColor: 'rgb(235, 253, 255) !important',
            '&:hover': {
              backgroundColor: 'rgb(235, 253, 255) !important',
            },
          },
          '& .contrast-row': {
            backgroundColor: 'rgb(237, 255, 235) !important',
            '&:hover': {
              backgroundColor: 'rgb(237, 255, 235) !important',
            },
          },
          '& .MuiDataGrid-cell': {
            userSelect: requiresToolbar ? 'none' : 'default',
            cursor: requiresToolbar ? 'pointer' : 'default',
          },
        }}
        isRowSelectable={(params) => requiresToolbar}
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
  contrastGroup: PropTypes.object.isRequired,
  referenceGroup: PropTypes.object.isRequired,
  onAddSamplesToGroup: PropTypes.func.isRequired,
};

export default DataTable;
