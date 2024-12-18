import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import IconButton from '../ui/IconButton';
import ToolTip from '../ui/ToolTip';
import terminal from '../../assets/icons/terminal.png';
import GeneSetCollectionsPopup from '../ui/GeneSetCollectionsPopup';

const GSEAInputForm = ({
    setIsVisible,
    onDatasetSelect,
    runAnalysis,
    isLoading,
    onAddGeneSetCollection,
    geneSetCollections,
    enrichParams,
    onUpdateEnrichParams,
    onRemoveGeneSetCollection
}) => {
    const [showGeneSetPopup, setShowGeneSetPopup] = useState(false);

    const handleRunAnalysis = () => {
        runAnalysis();
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
                                        onClick={() => setShowGeneSetPopup(true)}
                                    >
                                        Add gene set
                                    </button>
                                )}
                            </div>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Add a gene set library for enrichment analysis" />
                            </div>
                        </div>

                        <div className="form-field-row">
                            <div className="form-group">
                                <label className="form-label">Minimum Gene Set Size:</label>
                                <input
                                    type="number"
                                    className="form-input"
                                    value={enrichParams.minSize || 15}
                                    onChange={(e) => onUpdateEnrichParams({
                                        ...enrichParams,
                                        minSize: parseInt(e.target.value) || 15
                                    })}
                                    min="1"
                                />
                            </div>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Minimum number of genes required in a gene set" />
                            </div>
                        </div>

                        {geneSetCollections.length > 0 && (
                            <div className="selected-gene-sets">
                                <h4>Selected Gene Sets:</h4>
                                {geneSetCollections.map((collection, index) => (
                                    <div key={index} className="gene-set-item">
                                        <span>{collection.name}</span>
                                        <button
                                            className="remove-button"
                                            onClick={() => onRemoveGeneSetCollection(index)}
                                            aria-label="Remove gene set"
                                        >
                                            ×
                                        </button>
                                    </div>
                                ))}
                            </div>
                        )}

                        <div id="runAnalysisContainer">
                            <IconButton
                                icon={terminal}
                                label="Run Analysis"
                                onClick={handleRunAnalysis}
                                disabled={geneSetCollections.length === 0}
                            />
                        </div>
                    </div>
                </div>
            </div>

            <GeneSetCollectionsPopup
                isVisible={showGeneSetPopup}
                setIsVisible={setShowGeneSetPopup}
                onCollectionSelect={onAddGeneSetCollection}
            />

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

                .selected-gene-sets {
                    margin-top: 1rem;
                    padding: 1rem;
                    border: 1px solid #eee;
                    border-radius: 4px;
                }

                .selected-gene-sets h4 {
                    margin: 0 0 0.5rem 0;
                    font-size: 0.9rem;
                    color: #666;
                }

                .gene-set-item {
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    padding: 0.5rem;
                    margin-bottom: 0.5rem;
                    background-color: #f5f5f5;
                    border-radius: 4px;
                }

                .gene-set-item:last-child {
                    margin-bottom: 0;
                }

                .remove-button {
                    background: none;
                    border: none;
                    color: #666;
                    font-size: 1.2rem;
                    cursor: pointer;
                    padding: 0 0.5rem;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                }

                .remove-button:hover {
                    color: #ff4444;
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
    onAddGeneSetCollection: PropTypes.func.isRequired,
    geneSetCollections: PropTypes.array.isRequired,
    enrichParams: PropTypes.object.isRequired,
    onUpdateEnrichParams: PropTypes.func.isRequired,
    onRemoveGeneSetCollection: PropTypes.func.isRequired,
};

export default GSEAInputForm;
