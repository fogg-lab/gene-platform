import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import IconButton from '../ui/IconButton';
import ToolTip from '../ui/ToolTip';
import terminal from '../../assets/icons/terminal.png';

const GSEAInputForm = ({
    setIsVisible,
    onDatasetSelect,
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
            setIsVisible(true);
        } else {
            setIsVisible(false);
            onDatasetSelect('example', null);
        }
    };

    const handleRunAnalysis = () => {
        runAnalysis({ libraries, minSize, maxSize, permutations, seed });
    };

    const handleAddLibrary = (library) => {
        setLibraries([...libraries, library]);
        setShowLibrariesPopup(false);
    };

    return (
        <div id="analysisInputContainer_comp">
            <div className="form-with-tooltips">
                <div className="form-content">
                    <div className="form-fields">
                        <h3>Data</h3>
                        <div className="form-field-row">
                            <div>
                                {isLoading ? (
                                    <div className="loader"></div>
                                ) : (
                                    <button
                                        className="analysisInputButton"
                                        onClick={() => handleButtonClick('external')}
                                    >
                                        Add gene set
                                    </button>
                                )}
                            </div>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Add a gene set library for enrichment analysis" />
                            </div>
                        </div>
                        <h3>Configuration</h3>
                        <div className="form-container">
                            <form id="gseaOptionsForm">
                                <div className="form-field-row">
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
                                    <div className="tooltip-wrapper">
                                        <ToolTip content="Minimum number of genes required in a gene set for it to be considered in the analysis" />
                                    </div>
                                </div>

                                <div className="form-field-row">
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
                                    <div className="tooltip-wrapper">
                                        <ToolTip content="Maximum number of genes allowed in a gene set for it to be included in the analysis" />
                                    </div>
                                </div>

                                <div className="form-field-row">
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
                                    <div className="tooltip-wrapper">
                                        <ToolTip content="Number of permutations to perform for statistical significance testing" />
                                    </div>
                                </div>

                                <div className="form-field-row">
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
                                    <div className="tooltip-wrapper">
                                        <ToolTip content="Random seed for reproducible results" />
                                    </div>
                                </div>
                            </form>
                        </div>

                        <div id="runAnalysisContainer">
                            <IconButton icon={terminal} label="Run Analysis" onClick={handleRunAnalysis} />
                        </div>
                    </div>
                </div>
            </div>

            <style jsx>{`
                .form-with-tooltips {
                    display: flex;
                    width: 100%;
                    gap: 1rem;
                }

                .form-content {
                    flex: 1;
                }

                .form-fields {
                    display: flex;
                    flex-direction: column;
                    gap: 1rem;
                }

                .form-field-row {
                    display: flex;
                    align-items: flex-start;
                    gap: 1rem;
                    min-height: 2rem;
                }

                .form-field-row > *:first-child {
                    flex: 1;
                }

                .tooltip-wrapper {
                    width: 24px;
                    display: flex;
                    align-items: flex-start;
                    padding-top: 0.25rem;
                }

                .form-group {
                    display: flex;
                    flex-direction: column;
                    gap: 0.5rem;
                    width: 100%;
                }

                .form-label {
                    font-weight: 500;
                }

                .form-input {
                    padding: 0.5rem;
                    border: 1px solid #ccc;
                    border-radius: 4px;
                    width: 100%;
                }

                .analysisInputButton {
                    padding: 0.5rem 1rem;
                    border: 1px solid #ccc;
                    border-radius: 4px;
                }

                #runAnalysisContainer {
                    margin-top: 2rem;
                }

                .loader {
                    width: 24px;
                    height: 24px;
                    border: 2px solid #f3f3f3;
                    border-top: 2px solid #3498db;
                    border-radius: 50%;
                    animation: spin 1s linear infinite;
                }

                @keyframes spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
            `}</style>
        </div>
    );
};

GSEAInputForm.propTypes = {
    setIsVisible: PropTypes.func.isRequired,
    onDatasetSelect: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    isLoading: PropTypes.bool.isRequired,
};

export default GSEAInputForm;