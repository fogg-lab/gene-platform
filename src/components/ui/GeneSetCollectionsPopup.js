import React, { useState, useEffect, useCallback, useRef } from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';

const GeneSetCollectionsPopup = ({ setIsVisible, isVisible, onCollectionSelect }) => {
  const [species, setSpecies] = useState('human');
  const [collections, setCollections] = useState([]);
  const [searchTerm, setSearchTerm] = useState('');
  const [suggestions, setSuggestions] = useState([]);
  const [selectedCollection, setSelectedCollection] = useState(null);
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 50;
  const wrapperRef = useRef(null);

  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    import(`../../assets/external_data_index/MSigDBCollections.json`)
      .then(data => setCollections(Object.entries(data[`${species}Collections`])))
      .catch(error => console.error('Error loading collections:', error));
  }, [species]);

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
      const filtered = collections.filter(([id, info]) =>
        id.toLowerCase().includes(value.toLowerCase()) ||
        info.title.toLowerCase().includes(value.toLowerCase()) ||
        info.description.toLowerCase().includes(value.toLowerCase())
      );
      setSuggestions(filtered);
    } else {
      setSuggestions(collections);
    }
  }, [collections]);

  const paginatedSuggestions = suggestions.slice(
    (currentPage - 1) * itemsPerPage,
    currentPage * itemsPerPage
  );

  const totalPages = Math.ceil(suggestions.length / itemsPerPage);

  const handlePageChange = (newPage) => {
    setCurrentPage(newPage);
  };

  const handleSelectCollection = (collection) => {
    setSelectedCollection(collection);
    setSearchTerm(collection[1].title);
    setSuggestions([]);
  };

  const handleLoadCollection = async () => {
    if (selectedCollection) {
      setIsLoading(true);
      console.log("Attempting to search collection")
      try {
        onCollectionSelect(species, selectedCollection[0]);
        console.log("Collection received.")
        setIsVisible(false);
      } catch (error) {
        console.error('Error fetching external collection:', error);
      } finally {
        setIsLoading(false);
      }
    }
  };

  const handleClose = () => {
    setIsVisible(false);
    setSearchTerm('');
    setSelectedCollection(null);
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
            <h1>Search human and mouse gene sets</h1>
            <div className='inputWrapper' ref={wrapperRef}>
              <select
                value={species}
                onChange={(e) => {
                  setSpecies(e.target.value);
                  setSuggestions([]);
                }}
                style={{ marginRight: '10px' }}
              >
                <option value="human">Human</option>
                <option value="mouse">Mouse</option>
              </select>
              <input
                type="text"
                value={searchTerm}
                onChange={(e) => handleSearch(e.target.value)}
                onFocus={(e) => handleSearch(e.target.value)}
                placeholder="Search collections..."
              />
              {suggestions.length > 0 && (
                <>
                  <ul className="suggestions">
                    {paginatedSuggestions.map(([id, info]) => (
                      <li key={id} onClick={() => handleSelectCollection([id, info])}>
                        {id}: {truncateTitle(`${info.title} | ${info.description}`)}
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
          <div className="selectedCollection">
            {selectedCollection && (
              <p>{selectedCollection[0]}: {truncateTitle(selectedCollection[1].title)}</p>
            )}
          </div>
          <div className='databasePopupButton'>
            {isLoading ? (
              <div className="loader"></div>
            ) : (
              <IconButton
                icon={terminal}
                label="Select collection"
                onClick={handleLoadCollection}
                disabled={!selectedCollection}
              />
            )}
          </div>
        </div>
      )}
    </div>
  );
};

GeneSetCollectionsPopup.propTypes = {
  setIsVisible: PropTypes.func.isRequired,
  isVisible: PropTypes.bool.isRequired,
  onCollectionSelect: PropTypes.func.isRequired,
};

export default GeneSetCollectionsPopup;
