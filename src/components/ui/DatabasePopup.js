import React from 'react';
import ExternalDatasetSelector from './ExternalDatasetSelector';

const DatabasePopup = ({ isVisible, setIsVisible, onDatasetSelect }) => {
    if (!isVisible) return null;

    return (
        <div className="popup">
            <div className="popup-content">
                <h2>Select External Dataset</h2>
                <ExternalDatasetSelector onDatasetSelect={onDatasetSelect} />
                <button onClick={() => setIsVisible(false)}>Close</button>
            </div>
        </div>
    );
};

export default DatabasePopup;
