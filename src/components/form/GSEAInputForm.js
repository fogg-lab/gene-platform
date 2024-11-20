import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';
import GeneSetCollectionsPopup from '../ui/GeneSetCollectionsPopup';

const GSEAInputForm = ({
    setIsVisible,
    onDatasetSelect,
    runAnalysis,
    isLoading,
    onAddGeneSetCollection,
    geneSetCollections,
    gseaParams,
    onUpdateGseaParams,
    onRemoveGeneSetCollection
}) => {
    const [isGeneSetPopupVisible, setIsGeneSetPopupVisible] = useState(false);

    const handleRunAnalysis = () => {
        runAnalysis();
    };

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data</h3>
            <div>
                {isLoading ? (
                    <div className="loader"></div>
                ) : (
                    <>
                        <button
                            className="analysisInputButton"
                            onClick={() => setIsGeneSetPopupVisible(true)}
                        >
                            Add Gene Set Collection
                        </button>
                        <div className="geneSetCollections">
                            {geneSetCollections.map((collection, index) => (
                                <div key={index} className="geneSetCollection">
                                    <span>{collection.name}</span>
                                    <button 
                                        className="removeButton"
                                        onClick={() => onRemoveGeneSetCollection(index)}
                                        aria-label={`Remove ${collection.name}`}
                                    >
                                        ×
                                    </button>
                                </div>
                            ))}
                        </div>
                    </>
                )}
            </div>
            <h3>Configuration</h3>
            <div className="form-container">
                <form id="gseaOptionsForm">
                    <div className="form-group">
                        <label className="form-label" htmlFor="weight">Weight:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="weight"
                            name="weight"
                            value={gseaParams.weight}
                            onChange={(e) => onUpdateGseaParams('weight', parseFloat(e.target.value))}
                            required
                        />
                    </div>
                    <div className="form-group">
                        <label className="form-label" htmlFor="minSize">Minimum Size:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="minSize"
                            name="minSize"
                            value={gseaParams.minSize}
                            onChange={(e) => onUpdateGseaParams('minSize', parseInt(e.target.value))}
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
                            value={gseaParams.maxSize}
                            onChange={(e) => onUpdateGseaParams('maxSize', parseInt(e.target.value))}
                            required
                        />
                    </div>
                    <div className="form-group">
                        <label className="form-label" htmlFor="nperm">Permutations:</label>
                        <input
                            className="form-input"
                            type="number"
                            id="nperm"
                            name="nperm"
                            value={gseaParams.nperm}
                            onChange={(e) => onUpdateGseaParams('nperm', parseInt(e.target.value))}
                            required
                        />
                    </div>
                </form>
            </div>
            <GeneSetCollectionsPopup
                isVisible={isGeneSetPopupVisible}
                setIsVisible={setIsGeneSetPopupVisible}
                onCollectionSelect={onAddGeneSetCollection}
            />
            <div id="runAnalysisContainer">
                <IconButton icon={terminal} label="Run Analysis" onClick={handleRunAnalysis} />
            </div>
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
    gseaParams: PropTypes.object.isRequired,
    onUpdateGseaParams: PropTypes.func.isRequired,
    onRemoveGeneSetCollection: PropTypes.func.isRequired,
};

export default GSEAInputForm;