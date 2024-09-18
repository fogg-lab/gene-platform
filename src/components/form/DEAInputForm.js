import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';
import pako from 'pako';
import SampleField from '../ui/SampleField';

const DEAInputForm = ({
    setIsVisible,
    onDatasetSelect,
    contrastGroup,
    referenceGroup,
    onRemoveSamplesFromGroup,
    runAnalysis,
    isLoading
}) => {
    const [countsFileName, setCountsFileName] = useState('');
    const [coldataFileName, setColdataFileName] = useState('');

    const handleRemoveContrastSample = (sampleId) => {
        console.log("Removing contrast sample:", sampleId);
        onRemoveSamplesFromGroup(true, [sampleId]);
    };

    const handleRemoveReferenceSample = (sampleId) => {
        console.log("Removing reference sample:", sampleId);
        onRemoveSamplesFromGroup(false, [sampleId]);
    };

    const onDropCounts = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setCountsFileName(acceptedFiles[0].name);
        }
        console.log("AnalysisInputForm.js:onDropCounts requires further implementation");
    }, []);

    const onDropColdata = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setColdataFileName(acceptedFiles[0].name);
        }
        console.log("AnalysisInputForm.js:onDropColdata requires further implementation");
    }, []);

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data - DEAInputForm</h3>
            <h3>Configuration</h3>
            <div>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="covariates" />
                    <span>Add covariates (🚧)</span>
                </label>
                <label className="radioLabel">
                    <span id="adjustmentSubfield">Adjustment method (🚧):</span>
                    <select id="adjustmentMethod" name="adjustmentMethod">
                        <option value="option1">Bonferroni</option>
                        <option value="option2">Benjamini and Hochberg</option>
                    </select>
                </label>
                <div className='dataSubfieldSampleField'>
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
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="name" />
                    <span>Batch correction (🚧)</span>
                </label>
            </div>
            <div id="runAnalysisContainer">
                <IconButton icon={terminal} label="Run Analysis" onClick={runAnalysis} />
            </div>
        </div>
    );
};

DEAInputForm.propTypes = {
    setIsVisible: PropTypes.func.isRequired,
    onDatasetSelect: PropTypes.func.isRequired,
    contrastGroup: PropTypes.object.isRequired,
    referenceGroup: PropTypes.object.isRequired,
    onRemoveSamplesFromGroup: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    isLoading: PropTypes.bool.isRequired,
};

export default DEAInputForm;
