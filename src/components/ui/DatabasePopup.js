import React, { useState, useEffect, useCallback, useRef } from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';
import { getExternalDataset } from '../../services/api';

const DatabasePopup = ({ setIsVisible, isVisible, onDatasetSelect }) => {
  const [dataSrc, setDataSrc] = useState('GDC');
  const [datasets, setDatasets] = useState([]);
  const [searchTerm, setSearchTerm] = useState('');
  const [suggestions, setSuggestions] = useState([]);
  const [selectedDataset, setSelectedDataset] = useState(null);
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 50;
  const wrapperRef = useRef(null);

  useEffect(() => {
    import(`../../assets/external_data_index/${dataSrc}.json`)
      .then(data => setDatasets(Object.entries(data.datasets)))
      .catch(error => console.error('Error loading datasets:', error));
  }, [dataSrc]);

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (suggestions.length > 0 && !event.target.closest('.inputWrapper')) {
        setSuggestions([]);
      }
    };

    document.addEventListener('mousedown', handleClickOutside);

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [suggestions]);

  const handleSearch = useCallback((value) => {
    setSearchTerm(value);
    setCurrentPage(1);
    if (value.length > 0) {
      const filtered = datasets.filter(([id, info]) => 
        id.toLowerCase().includes(value.toLowerCase()) || 
        info.title.toLowerCase().includes(value.toLowerCase())
      );
      setSuggestions(filtered);
    } else {
      setSuggestions(datasets);
    }
  }, [datasets]);

  const paginatedSuggestions = suggestions.slice(
    (currentPage - 1) * itemsPerPage,
    currentPage * itemsPerPage
  );

  const totalPages = Math.ceil(suggestions.length / itemsPerPage);

  const handlePageChange = (newPage) => {
    setCurrentPage(newPage);
  };

  const handleSelectDataset = (dataset) => {
    setSelectedDataset(dataset);
    setSearchTerm(dataset[1].title);
    setSuggestions([]);
  };

  const handleLoadDataset = async () => {
    if (selectedDataset) {
      try {
        const data = await getExternalDataset(dataSrc, selectedDataset[0]);
        onDatasetSelect('external', data);
        setIsVisible(false);
      } catch (error) {
        console.error('Error fetching external dataset:', error);
      }
    }
  };

  const handleClose = () => {
    setIsVisible(false);
    setSearchTerm('');
    setSelectedDataset(null);
    setSuggestions([]);
  };

  const truncateTitle = (title) => {
    return title.length > 49 ? title.substring(0, 46) + '...' : title;
  };

  return (
    <div>
      {isVisible && (
        <div className="databasePopup">
          <div
            className="closeButton"
            onClick={handleClose}
            role="button"
            tabIndex="0"
          >
            <img src={closeIcon} alt="Close" />
          </div>
          <div>
            <h1>Search GEO and GDC Datasets</h1>
            <div className='inputWrapper' ref={wrapperRef}>
              <select
                value={dataSrc}
                onChange={(e) => {
                  setDataSrc(e.target.value);
                  setSuggestions([]);
                }}
                style={{ marginRight: '10px' }}
              >
                <option value="GDC">GDC</option>
                <option value="GEO">GEO</option>
              </select>
              <input
                type="text"
                value={searchTerm}
                onChange={(e) => handleSearch(e.target.value)}
                onFocus={(e) => handleSearch(e.target.value)}
                placeholder="Search datasets..."
              />
              {suggestions.length > 0 && (
                <>
                  <ul className="suggestions">
                    {paginatedSuggestions.map(([id, info]) => (
                      <li key={id} onClick={() => handleSelectDataset([id, info])}>
                        {id}: {truncateTitle(info.title)}
                      </li>
                    ))}
                  </ul>
                  {totalPages > 1 && (
                    <div className="pagination">
                      <button 
                        onClick={() => handlePageChange(currentPage - 1)}
                        disabled={currentPage === 1}
                      >
                        Previous
                      </button>
                      <span>{currentPage} / {totalPages}</span>
                      <button 
                        onClick={() => handlePageChange(currentPage + 1)}
                        disabled={currentPage === totalPages}
                      >
                        Next
                      </button>
                    </div>
                  )}
                </>
              )}
            </div>
          </div>
          <div className="selectedDataset">
            {selectedDataset && (
              <p>{selectedDataset[0]}: {truncateTitle(selectedDataset[1].title)}</p>
            )}
          </div>
          <div className='databasePopupButton'>
            <IconButton 
              icon={terminal}
              label="Select dataset" 
              onClick={handleLoadDataset}
              disabled={!selectedDataset}
            />
          </div>
        </div>
      )}
    </div>
  );
};

DatabasePopup.propTypes = {
  setIsVisible: PropTypes.func.isRequired,
  isVisible: PropTypes.bool.isRequired,
  onDatasetSelect: PropTypes.func.isRequired,
};

export default DatabasePopup;
