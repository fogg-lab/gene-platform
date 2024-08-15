import React, { useState, useEffect, useCallback, useRef } from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';

const DatabasePopup = ({ setIsVisible, isVisible, isCheckedRadioButton, setSelectedRadio, selectedRadio }) => {
  let iconSrc;

  try {
    iconSrc = require(`../../assets/icons/close.svg`).default;
  } catch (error) {
    console.error(`Icon not found in assets/icons`);
    return null;
  }

  const handleClose = () => {
    setIsVisible(false);
    setSelectedRadio(null); // Uncheck the radio button
  };

  return (
    <div>
      {isVisible ? (
        <div className="databasePopup">
          <div className="closeButton" onClick={handleClose}>
            <img src={iconSrc} alt="Close" />
          </div>
          <div>
            <h1>Search GEO and GDC Datasets</h1>
            <form className='inputWrapper'>
              <input type="text" id="geoInput" name="geoInput" placeholder="GSE12345..."></input>
            </form>
          </div>
          <div className="selectedDataset">
            <p>selected datasets will go here</p>
          </div>
          <div className='databasePopupButton'>
            <IconButton icon={terminal} label="Select dataset" />
          </div>
        </div>
      ) : (
        null
      )}
    </div>
  );
};

DatabasePopup.propTypes = {
  setIsVisible: PropTypes.func.isRequired,
  isVisible: PropTypes.bool.isRequired,
  isCheckedRadioButton: PropTypes.bool.isRequired,
  setSelectedRadio: PropTypes.func.isRequired,
  selectedRadio: PropTypes.string,
};

export default DatabasePopup;
