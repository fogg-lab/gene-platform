import React from 'react';
import PropTypes from 'prop-types';
import { useErrorPopup } from '../ui/ErrorPopup';
import IconButton from '../ui/IconButton';
import ToolTip from '../ui/ToolTip';
import SampleField from '../ui/SampleField';
import terminal from '../../assets/icons/terminal.png';
import next from '../../assets/icons/next.svg';

const DEAInputForm = ({
    contrastGroup,
    referenceGroup,
    onRemoveSamplesFromGroup,
    runAnalysis,
    handleStageChange,
    currentStage,
    deData,
}) => {
    const { showError } = useErrorPopup();

    return (
        <div id="analysisInputContainer_comp">
            <div className="form-with-tooltips">
                <div className="form-content">
                    <div className="form-fields">
                        <div className="form-field-row">
                            <label className="radioLabel opacity-50 pointer-events-none">
                                <input className="radioInput" type="radio" name="covariates" />
                                <span>Add covariates (🚧)</span>
                            </label>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Include additional variables that might affect gene expression" />
                            </div>
                        </div>

                        <div className="form-field-row">
                            <label className="radioLabel opacity-50 pointer-events-none">
                                <span id="adjustmentSubfield">Adjustment method (🚧):</span>
                                <select id="adjustmentMethod" name="adjustmentMethod">
                                    <option value="bh">Benjamini and Hochberg</option>
                                    <option value="bonferroni">Bonferroni</option>
                                </select>
                            </label>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Choose the method to control for multiple testing" />
                            </div>
                        </div>

                        <div className="form-field-row">
                            <div className="dataSubfieldSampleField">
                                <SampleField
                                    headerName="Reference Group"
                                    samples={referenceGroup.samples}
                                    onRemoveSample={(sampleId) => onRemoveSamplesFromGroup(false, [sampleId])}
                                    isContrast={false}
                                />
                                <SampleField
                                    headerName="Contrast Group"
                                    samples={contrastGroup.samples}
                                    onRemoveSample={(sampleId) => onRemoveSamplesFromGroup(true, [sampleId])}
                                    isContrast={true}
                                />
                            </div>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Select samples for your reference (control) and contrast (treatment) groups" />
                            </div>
                        </div>

                        <div className="form-field-row">
                            <label className="radioLabel opacity-50 pointer-events-none">
                                <input className="radioInput" type="radio" name="name" />
                                <span>Batch correction (🚧)</span>
                            </label>
                            <div className="tooltip-wrapper">
                                <ToolTip content="Correct for systematic variations between sample batches" />
                            </div>
                        </div>
                    </div>

                    <div id="runAnalysisContainer">
                        <IconButton
                            icon={deData?.plots?.volcano ? next : terminal}
                            label={deData?.plots?.volcano ? "Next Stage" : "Run Analysis"}
                            onClick={() => {
                                if (deData?.plots?.volcano) {
                                    handleStageChange('enrichment');
                                } else {
                                    if (referenceGroup.samples.length < 1 || contrastGroup.samples.length < 1) {
                                        showError("Each group must have at least one sample to run the analysis.");
                                        return;
                                    }
                                    runAnalysis();
                                }
                            }}
                        />
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

                .radioLabel {
                    display: flex;
                    align-items: center;
                    gap: 0.5rem;
                }

                #runAnalysisContainer {
                    margin-top: 1rem;
                }
            `}</style>
        </div>
    );
};

DEAInputForm.propTypes = {
    contrastGroup: PropTypes.object.isRequired,
    referenceGroup: PropTypes.object.isRequired,
    onRemoveSamplesFromGroup: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    handleStageChange: PropTypes.func.isRequired,
    currentStage: PropTypes.string.isRequired,
    deData: PropTypes.object,
};

export default DEAInputForm;