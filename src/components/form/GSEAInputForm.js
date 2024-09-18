import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';
import pako from 'pako';
import SampleField from '../ui/SampleField';

const GSEAInputForm = ({
    setIsVisible,
    onDatasetSelect,
    contrastGroup,
    referenceGroup,
    onRemoveSamplesFromGroup,
    runAnalysis,
    isLoading
}) => {
    const [libraries, setLibraries] = useState([]);
    const [minSize, setMinSize] = useState(15);
    const [maxSize, setMaxSize] = useState(500);
    const [permutations, setPermutations] = useState(1000);
    const [seed, setSeed] = useState(123);
    const [showLibrariesPopup, setShowLibrariesPopup] = useState(false);

    const handleButtonClick = (datasetType) => {
        if (datasetType === 'external') {
            setIsVisible(true); // Show plot area when 'Use External Dataset' is selected
        } else {
            setIsVisible(false); // Hide plot area when 'Use Example Dataset' is selected
            // Load example dataset
            onDatasetSelect('example', null);
        }
    };

    const handleRunAnalysis = () => {
        // Pass all GSEA options to the runAnalysis function
        runAnalysis({ libraries, minSize, maxSize, permutations, seed });
    };

    const handleAddLibrary = (library) => {
        setLibraries([...libraries, library]);
        setShowLibrariesPopup(false);
    };

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data - GSEAInputForm</h3>
            <div>
                {isLoading ? (
                    <div className="loader"></div>
                ) : (
                    <>
                        <button
                            className="analysisInputButton"
                            onClick={() => handleButtonClick('external')}
                        >
                            Add gene set
                        </button>
                    </>
                )}
            </div>
            <h3>Configuration</h3>
            <div className="form-container">
                <form id="gseaOptionsForm">
                    <div className="form-group">
                        <label className="form-label" htmlFor="minSize">Minimum Size:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="minSize"
                            name="minSize"
                            value={minSize}
                            onChange={(e) => setMinSize(Number(e.target.value))}
                            required
                        />
                    </div>
                    <div className="form-group">
                        <label className="form-label" htmlFor="maxSize">Maximum Size:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="maxSize"
                            name="maxSize"
                            value={maxSize}
                            onChange={(e) => setMaxSize(Number(e.target.value))}
                            required
                        />
                    </div>
                    <div className="form-group">
                        <label className="form-label" htmlFor="permutations">Permutations:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="permutations"
                            name="permutations"
                            value={permutations}
                            onChange={(e) => setPermutations(Number(e.target.value))}
                            required
                        />
                    </div>
                    <div className="form-group">
                        <label className="form-label" htmlFor="seed">Seed:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="seed"
                            name="seed"
                            value={seed}
                            onChange={(e) => setSeed(Number(e.target.value))}
                            required
                        />
                    </div>
                </form>
            </div>
            <div id="runAnalysisContainer">
                <IconButton icon={terminal} label="Run Analysis" onClick={handleRunAnalysis} />
            </div>
        </div>
    );
};

GSEAInputForm.propTypes = {
    setIsVisible: PropTypes.func.isRequired,
    onDatasetSelect: PropTypes.func.isRequired,
    contrastGroup: PropTypes.object.isRequired,
    referenceGroup: PropTypes.object.isRequired,
    onRemoveSamplesFromGroup: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    isLoading: PropTypes.bool.isRequired,
};

export default GSEAInputForm;