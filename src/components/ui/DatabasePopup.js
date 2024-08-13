import React from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';
import IconButton from '../ui/IconButton';

const DatabasePopup = ({ setIsVisible, isVisible, setSelectedRadio }) => {

  const handleClose = () => {
    setIsVisible(false);
    setSelectedRadio(null); // Uncheck the radio button
  };

  return (
    <div>
      {isVisible ? (
        <div className="databasePopup">
          <div className="closeButton" onClick={handleClose} role="button" tabIndex={0} onKeyPress={handleClose}>
            <img src={closeIcon} alt="Close" />
          </div>
          <div>
            <h1>Search GEO and GDC Datasets</h1>
            <form className='inputWrapper'>
              <input type="text" id="geoInput" name="geoInput" placeholder="GSE12345..."></input>
            </form>
          </div>
          <div className="selectedDataset">
            <p>Selected datasets will go here</p>
          </div>
          <div className='databasePopupButton'>
            <IconButton iconFilename="terminal.png" label="Select dataset" />
          </div>
        </div>
      ) : null}
    </div>
  );
};

DatabasePopup.propTypes = {
  setIsVisible: PropTypes.func.isRequired,
  isVisible: PropTypes.bool.isRequired,
  setSelectedRadio: PropTypes.func.isRequired,
  selectedRadio: PropTypes.string,
};

export default DatabasePopup;
