import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
// import Papa from 'papaparse';
import terminal from '../../assets/icons/terminal.png';
import SampleField from '../ui/SampleField';

const FileDropArea = ({ title, onDrop, fileName }) => {
    const [isFileTypeValid, setIsFileTypeValid] = useState(true);

    const onDragEnter = useCallback((event) => {
        const fileType = event.dataTransfer.items[0].type;
        if (fileType !== 'text/tab-separated-values' && fileType !== 'text/csv') {
            setIsFileTypeValid(false);
        } else {
            setIsFileTypeValid(true);
        }
    }, []);

    const onDragLeave = useCallback(() => {
        setIsFileTypeValid(true);
    }, []);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: (acceptedFiles) => {
            const validFiles = acceptedFiles.filter(file => file.type === 'text/tab-separated-values' || file.type === 'text/csv');
            onDrop(validFiles);
        },
        onDragEnter,
        onDragLeave,
        accept: '.tsv, .csv'
    });

    return (
        <div {...getRootProps()} className="filedropArea">
            <input {...getInputProps()} className="fileDrop" />
            <h4>{title}</h4>
            <span>Drop .tsv/.csv file here or</span>
            <button className="openFilesystemButton">
                <span>Browse</span>
            </button>
            {isDragActive ? (
                <p>{isFileTypeValid ? '' : <span style={{ color: 'red' }}>Invalid file type.</span>}</p>
            ) : (
                <p></p>
            )}
            {fileName && <p className="fileName">{fileName}</p>}
        </div>
    );
};

const AnalysisInputForm = ({
    setIsVisible,
    contrastGroups,
    referenceGroups,
    onAddGroup,
    onUpdateGroup,
    selectedSamples,
    onAddSamplesToGroup
}) => {
    // const [countsFileName, setCountsFileName] = useState('');
    // const [coldataFileName, setColdataFileName] = useState('');

    // const cleanData = (data) => {
    //     return data
    //         .map(item => item && item.trim()) // Trim whitespace
    //         .filter(item => item); // Filter out empty strings or null values
    // };

    // const onDropCounts = useCallback((acceptedFiles) => {
    //     if (acceptedFiles.length > 0) {
    //         setCountsFileName(acceptedFiles[0].name);
    //     }
    // }, []);

    const handleButtonClick = (datasetType) => {
        if (datasetType === 'external') {
            setIsVisible(true); // Show plot area when 'Use External Dataset' is selected
        } else {
            setIsVisible(false); // Hide plot area when 'Use Example Dataset' is selected
        }
    };

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data</h3>
            <div className='dataSubfield'>
                {/* Load example data / External dataset */}
                <button
                    className="analysisInputButton"
                    onClick={() => handleButtonClick('example')}
                >
                    Use Example Dataset
                </button>
                <button
                    className="analysisInputButton"
                    onClick={() => handleButtonClick('external')}
                >
                    Use External Dataset
                </button>
            </div>
            <div id="filedropContainer">
                <FileDropArea title="Counts" />
                <FileDropArea title="Coldata" />
            </div>
            <h3>Configuration</h3>
            <div>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="covariates" />
                    <span>Add covariates</span>
                </label>
                <label className="radioLabel">
                    <span>Data type:</span>
                    <select id="exampleDropdown" name="dataType">
                        <option value="option1">Microarray</option>
                        <option value="option2">RNA-Seq</option>
                    </select>
                </label>
                <label className="radioLabel">
                    <span id="adjustmentSubfield">Adjustment method:</span>
                    <select id="exampleDropdown" name="exampleDropdown">
                        <option value="option1">Bonferroni</option>
                        <option value="option2">Benjamini and Hochberg</option>
                    </select>
                </label>
                <div className='dataSubfield'>
                    <div className="contrastReferenceFields">
                        <SampleField
                            headerName='Contrast Groups'
                            groups={contrastGroups}
                            onAddGroup={() => onAddGroup(true)}
                            onUpdateGroup={(id, updates) => onUpdateGroup(id, updates, true)}
                        />
                        <button
                            onClick={() => onAddSamplesToGroup(contrastGroups[0]?.id, true)}
                            disabled={!contrastGroups.length || !selectedSamples.length}
                        >
                            Add to Contrast Group
                        </button>
                    </div>
                    <div className="contrastReferenceFields">
                        <SampleField
                            headerName='Reference Groups'
                            groups={referenceGroups}
                            onAddGroup={() => onAddGroup(false)}
                            onUpdateGroup={(id, updates) => onUpdateGroup(id, updates, false)}
                        />
                        <button
                            onClick={() => onAddSamplesToGroup(referenceGroups[0]?.id, false)}
                            disabled={!referenceGroups.length || !selectedSamples.length}
                        >
                            Add to Reference Group
                        </button>
                    </div>
                </div>
                <label className="radioLabel">
                    <span>Data transformation:</span>
                    <select id="exampleDropdown" name="exampleDropdown">
                        <option value="option1">None</option>
                        <option value="option1">VST</option>
                        <option value="option2">rlog</option>
                    </select>
                </label>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="name" />
                    <span>Batch correction</span>
                </label>
            </div>
            <div id="runAnalysisContainer">
                <IconButton icon={terminal} label="Run Analysis" />
            </div>
        </div>
    );
};

FileDropArea.propTypes = {
    title: PropTypes.string.isRequired,
    onDrop: PropTypes.func.isRequired,
    fileName: PropTypes.string.isRequired,
};

AnalysisInputForm.propTypes = {
    setIsVisible: PropTypes.func.isRequired,
    contrastGroups: PropTypes.array.isRequired,
    referenceGroups: PropTypes.array.isRequired,
    onAddGroup: PropTypes.func.isRequired,
    onUpdateGroup: PropTypes.func.isRequired,
    selectedSamples: PropTypes.array.isRequired,
    onAddSamplesToGroup: PropTypes.func.isRequired,
};


export default AnalysisInputForm;
