import React, { useState, useEffect, useCallback, useRef } from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';
import IconButton from './IconButton';
import terminal from '../../assets/icons/terminal.png';
import { getExternalGeneSet } from '../../services/api';

const GeneSetFilterPopup = ({ isVisible, setIsVisible, onFilterApply, species }) => {
  const [searchTerm, setSearchTerm] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [suggestions, setSuggestions] = useState([]);
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 50;
  const wrapperRef = useRef(null);

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
    if (value.length > 2) {
      try {
        import(`../../assets/external_data_index/MSigDBGeneSets.json`).then(response => {
          const geneSets = response[`${species}Sets`];
          const filtered = Object.keys(geneSets)
            .filter(name => name.toLowerCase().includes(value.toLowerCase()));
          setSuggestions(filtered);
        });
      } catch (error) {
        console.error('Error loading gene sets:', error);
        setSuggestions([]);
      }
    } else {
      setSuggestions([]);
    }
  }, [species]);

  const paginatedSuggestions = suggestions.slice(
    (currentPage - 1) * itemsPerPage,
    currentPage * itemsPerPage
  );

  const totalPages = Math.ceil(suggestions.length / itemsPerPage);

  const handlePageChange = (newPage) => {
    setCurrentPage(newPage);
  };

  const handleSelectGeneSet = (geneSet) => {
    setSearchTerm(geneSet);
    setSuggestions([]);
  };

  const handleApplyFilter = async () => {
    if (!searchTerm) return;
    
    setIsLoading(true);
    setError(null);
    try {
      const data = await getExternalGeneSet(species, searchTerm);
      onFilterApply(searchTerm, data.geneSymbols);
      setIsVisible(false);
    } catch (error) {
      console.error('Error applying gene set filter:', error);
      setError('Gene set not found. Please check the name and try again.');
    } finally {
      setIsLoading(false);
    }
  };

  const handleClose = () => {
    setIsVisible(false);
    setSearchTerm('');
    setError(null);
    setSuggestions([]);
  };

  return (
    <div>
      {isVisible && (
        <div className="databasePopup">
          <div className="closeButton" onClick={handleClose}>
            <img src={closeIcon} alt="Close" />
          </div>
          <div>
            <h1>Filter by Gene Set</h1>
            <div className="inputWrapper" ref={wrapperRef}>
              <input
                type="text"
                value={searchTerm}
                onChange={(e) => handleSearch(e.target.value)}
                onFocus={(e) => handleSearch(e.target.value)}
                placeholder="Search gene sets (e.g., HALLMARK_HYPOXIA)"
              />
              {suggestions.length > 0 && (
                <>
                  <ul className="suggestions">
                    {paginatedSuggestions.map(name => (
                      <li key={name} onClick={() => handleSelectGeneSet(name)}>
                        {name}
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
              {error && (
                <div className="error-message" style={{ color: 'red', marginTop: '8px', fontSize: '14px' }}>
                  {error}
                </div>
              )}
            </div>
          </div>
          <div className="databasePopupButton">
            {isLoading ? (
              <div className="loader"></div>
            ) : (
              <IconButton
                icon={terminal}
                label="Apply Filter"
                onClick={handleApplyFilter}
                disabled={!searchTerm}
              />
            )}
          </div>
        </div>
      )}
    </div>
  );
};

GeneSetFilterPopup.propTypes = {
  isVisible: PropTypes.bool.isRequired,
  setIsVisible: PropTypes.func.isRequired,
  onFilterApply: PropTypes.func.isRequired,
  species: PropTypes.string.isRequired,
};

export default GeneSetFilterPopup;
