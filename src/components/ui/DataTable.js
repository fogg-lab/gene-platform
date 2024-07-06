import React, { useState, useEffect } from 'react';
import { useVirtualizer } from '@tanstack/react-virtual';
import { Card, CardContent } from '@/components/ui/card';
import PropTypes from 'prop-types';

const DataTable = ({ data, columns }) => {
  const [tableData, setTableData] = useState([]);
  const [tableColumns, setTableColumns] = useState([]);

  useEffect(() => {
    setTableData(data);
    setTableColumns(columns);
  }, [data, columns]);

  const parentRef = React.useRef();

  const rowVirtualizer = useVirtualizer({
    count: tableData.length + 1, // +1 for header row
    getScrollElement: () => parentRef.current,
    estimateSize: () => 35,
    overscan: 5,
  });

  const columnVirtualizer = useVirtualizer({
    horizontal: true,
    count: tableColumns.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 150,
    overscan: 5,
  });

  return (
    <Card className="w-full">
      <CardContent>
        <div
          ref={parentRef}
          className="overflow-auto"
          style={{
            height: '400px',
            width: '100%',
          }}
        >
          <div
            style={{
              height: `${rowVirtualizer.getTotalSize()}px`,
              width: `${columnVirtualizer.getTotalSize()}px`,
              position: 'relative',
            }}
          >
            {rowVirtualizer.getVirtualItems().map((virtualRow) => (
              <React.Fragment key={virtualRow.index}>
                {columnVirtualizer.getVirtualItems().map((virtualColumn) => (
                  <div
                    key={`${virtualRow.index}-${virtualColumn.index}`}
                    style={{
                      position: 'absolute',
                      top: 0,
                      left: 0,
                      width: `${virtualColumn.size}px`,
                      height: `${virtualRow.size}px`,
                      transform: `translateX(${virtualColumn.start}px) translateY(${virtualRow.start}px)`,
                    }}
                  >
                    <div
                      style={{
                        padding: '4px',
                        borderBottom: '1px solid #eee',
                        borderRight: '1px solid #eee',
                        background: virtualRow.index === 0 ? '#f0f0f0' : 'white',
                      }}
                    >
                      {virtualRow.index === 0
                        ? tableColumns[virtualColumn.index]
                        : tableData[virtualRow.index - 1][virtualColumn.index]}
                    </div>
                  </div>
                ))}
              </React.Fragment>
            ))}
          </div>
        </div>
      </CardContent>
    </Card>
  );
};

DataTable.propTypes = {
  data: PropTypes.arrayOf(PropTypes.arrayOf(PropTypes.any)).isRequired,
  columns: PropTypes.arrayOf(PropTypes.string).isRequired,
};

export default DataTable;

