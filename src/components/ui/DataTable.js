import React, { useState, useMemo } from 'react';
import DataGrid, { SelectColumn, textEditor } from 'react-data-grid';
import 'react-data-grid/lib/styles.css';
import PropTypes from 'prop-types';

const DataTable = ({ data, columns }) => {
  const [rows, setRows] = useState(data);
  const [selectedRows, setSelectedRows] = useState(new Set());
  const [sortColumns, setSortColumns] = useState([]);

  const gridColumns = useMemo(() => {
    return [
      SelectColumn,
      ...columns.map(col => ({
        ...col,
        sortable: true,
        resizable: true,
        editor: textEditor
      }))
    ];
  }, [columns]);

  const sortedRows = useMemo(() => {
    if (sortColumns.length === 0) return rows;

    return [...rows].sort((a, b) => {
      for (const sort of sortColumns) {
        const comparator = (a, b) => {
          if (a[sort.columnKey] === b[sort.columnKey]) return 0;
          return a[sort.columnKey] > b[sort.columnKey] ? 1 : -1;
        };
        const compResult = comparator(a, b);
        if (compResult !== 0) {
          return sort.direction === 'ASC' ? compResult : -compResult;
        }
      }
      return 0;
    });
  }, [rows, sortColumns]);

  function rowKeyGetter(row) {
    return row.id;
  }

  return (
    <DataGrid
      columns={gridColumns}
      rows={sortedRows}
      rowKeyGetter={rowKeyGetter}
      onRowsChange={setRows}
      selectedRows={selectedRows}
      onSelectedRowsChange={setSelectedRows}
      sortColumns={sortColumns}
      onSortColumnsChange={setSortColumns}
      className="fill-grid"
    />
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
