import React from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';

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
    //setIsCheckedRadiobutton(false); // Uncheck the radio button
  };

  return (
    <div>
      {isVisible ? (
        <div className="databasePopup">
          <div class="closeButton" onClick={handleClose}>
            <img src={iconSrc} alt="" />
          </div>
          <div>
            <h1>Search GEO datasets</h1>
            <form className='inputWrapper'>
              <input type="text" id="textInput" name="textInput" placeholder="GSE12345..."></input>
            </form>
          </div>
          <div>
            <h1>Search GDC datasets</h1>
            <form className='inputWrapper'>
              <input type="text" id="textInput" name="textInput" placeholder="TCGA-LUAD..."></input>
            </form>
          </div>
          <div className="selectedDataset">
            <p>selected datasets will go here</p>
          </div>
        </div>
      ) : (
        null // Or you can display an alternative message or nothing at all
      )}
    </div>
  );
};

export default DatabasePopup;
