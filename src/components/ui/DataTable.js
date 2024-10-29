import React, { useState, useMemo, useCallback } from 'react';
import { DataGridPro, GridToolbarContainer, GridToolbarColumnsButton, GridToolbarFilterButton, GridToolbarExport, GridToolbarDensitySelector, gridExpandedSortedRowIdsSelector, useGridApiContext } from '@mui/x-data-grid-pro';
import PropTypes from 'prop-types';

const getFilteredRows = ({ apiRef }) => gridExpandedSortedRowIdsSelector(apiRef);
const getIntersectingRows = (filteredRows, selectionModel) => {
  return selectionModel.filter(id => filteredRows.includes(id));
};

const CustomToolbar = ({ onAddSamplesToGroup, selectionModel, rows, clearSelection }) => {
  const [selectedGroupType, setSelectedGroupType] = useState('reference');
  const apiRef = useGridApiContext();

  const handleAddSamples = () => {
    console.log("handleAddSamples called");
    console.log("Current selectionModel:", selectionModel);
    const filteredSelectionModel = getIntersectingRows(getFilteredRows({ apiRef }), selectionModel);
    console.log("Filtered selectionModel:", filteredSelectionModel);
    const selectedSamples = filteredSelectionModel.map(id => {
      const row = rows.find(r => r.id === id);
      return {
        id: row.sample_id || row.id,
        [Object.keys(row)[0]]: row[Object.keys(row)[0]]
      };
    });
    console.log("Selected samples:", selectedSamples);
    onAddSamplesToGroup(selectedGroupType === 'contrast', selectedSamples);
    console.log("Calling clearSelection");
    clearSelection();
    console.log("clearSelection called");
  };

  const numSelectedFiltered = getIntersectingRows(getFilteredRows({ apiRef }), selectionModel).length;

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
        <option value="reference">Reference Group</option>
        <option value="contrast">Contrast Group</option>
      </select>
      <button
        onClick={handleAddSamples}
        disabled={numSelectedFiltered === 0}
      >
        Add Selected ({numSelectedFiltered})
      </button>
    </GridToolbarContainer>
  );
};
const DataTable = ({
  data,
  columns,
  contrastGroup,
  referenceGroup,
  onAddSamplesToGroup,
}) => {
  const [sortModel, setSortModel] = useState([]);
  const [selectionModel, setSelectionModel] = useState([]);
  const [filterModel, setFilterModel] = useState({ items: [] });
  const [key, setKey] = useState(0);
  const [rowGroups, setRowGroups] = useState({});

  console.log("DataTable rendering, rowGroups:", rowGroups);

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


  const handleHeaderClick = useCallback((params) => {
    if (params.field === '__check__') {
      console.log("Select-all header clicked");
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

  const handleAddSamplesToGroup = useCallback((isContrast, samples) => {
    console.log("handleAddSamplesToGroup called", { isContrast, samples });
    const groupType = isContrast ? 'contrast' : 'reference';
    setRowGroups(prevGroups => {
      const newGroups = { ...prevGroups };
      samples.forEach(sample => {
        newGroups[sample.id || sample.sample_id] = groupType;
      });
      console.log("New rowGroups:", newGroups);
      return newGroups;
    });
    onAddSamplesToGroup(isContrast, samples);
    setKey(prevKey => {
      console.log("Incrementing key to force re-render");
      return prevKey + 1;
    });
  }, [onAddSamplesToGroup]);

  const getRowClassName = useCallback((params) => {
    const sampleId = params.row.sample_id || params.id;
    if (contrastGroup.samples.some(sample => sample.id === sampleId)) {
      return 'contrast-row';
    } else if (referenceGroup.samples.some(sample => sample.id === sampleId)) {
      return 'reference-row';
    }
    return '';
  }, [contrastGroup, referenceGroup]);

  const clearSelection = useCallback(() => {
    console.log("clearSelection function called");
    setSelectionModel([]);
    setKey(prevKey => prevKey + 1);
    console.log("Selection cleared and DataGrid re-render triggered");
  }, []);

  console.log("Rendering DataTable, current selectionModel:", selectionModel);

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
        onFilterModelChange={setFilterModel}
        checkboxSelection
        disableColumnMenu={false}
        disableSelectionOnClick={true}
        density="compact"
        getRowId={(row) => row.id}
        onColumnHeaderClick={handleHeaderClick}
        selectionModel={selectionModel}
        onSelectionModelChange={(newSelectionModel) => {
          console.log("onSelectionModelChange called", newSelectionModel);
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
            onAddSamplesToGroup: handleAddSamplesToGroup,
            selectionModel: selectionModel,
            rows: filteredRows,
            clearSelection: clearSelection,
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
  contrastGroup: PropTypes.object.isRequired,
  referenceGroup: PropTypes.object.isRequired,
  onAddSamplesToGroup: PropTypes.func.isRequired,
};

export default DataTable;
