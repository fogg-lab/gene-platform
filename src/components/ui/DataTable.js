import React, { useState, useMemo, useEffect } from 'react';
import DataGrid from 'react-data-grid';
import 'react-data-grid/lib/styles.css';
import PropTypes from 'prop-types';

const DataTable = ({ data, columns }) => {
  const [sortColumns, setSortColumns] = useState([]);

  const gridColumns = useMemo(() => {
    return columns.map(col => ({
      ...col,
      sortable: true,
      resizable: true,
      width: 150,
      minwidth: 50,
    }));
  }, [columns]);

  const filteredData = useMemo(() => {
    return data.filter(row => {
      // Check if the row has at least one non-empty value
      return Object.values(row).some(value => 
        value !== null && value !== undefined && value !== ''
      );
    });
  }, [data]);

  const sortedRows = useMemo(() => {
    if (sortColumns.length === 0) return filteredData;

    return [...filteredData].sort((a, b) => {
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
  }, [filteredData, sortColumns]);

  function rowKeyGetter(row) {
    return row.id;
  }

  // ResizeObserver error workaround
  useEffect(() => {
    const resizeObserverError = error => {
      if (error.message.includes('ResizeObserver loop')) {
        const resizeObserverError = new Error('Resize observer error');
        console.log('ResizeObserver error caught and suppressed:', resizeObserverError);
      }
    };
    window.addEventListener('error', resizeObserverError);
    return () => window.removeEventListener('error', resizeObserverError);
  }, []);

  return (
    <DataGrid
      columns={gridColumns}
      rows={sortedRows}
      rowKeyGetter={rowKeyGetter}
      sortColumns={sortColumns}
      onSortColumnsChange={setSortColumns}
      className="fill-grid"
      defaultColumnOptions={{
        resizable: true,
        sortable: true,
      }}
      style={{
        height: 'calc(100% - 40px)',
      }}
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
