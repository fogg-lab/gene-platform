import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
// import Papa from 'papaparse';
import terminal from '../../assets/icons/terminal.png';
import pako from 'pako'; // Import pako for gzip decompression
import SampleField from '../ui/SampleField';

function validFileType(filetype) {
    return filetype.startsWith("text/") || filetype == "application/gzip" || filetype == "application/x-gzip";
}

const FileDropArea = ({ title, onDrop, fileName }) => {
    const [isFileTypeValid, setIsFileTypeValid] = useState(true);

    const onDragEnter = useCallback((event) => {
        const fileType = event.dataTransfer.items[0].type;
        setIsFileTypeValid(validFileType(fileType));
    }, []);

    const onDragLeave = useCallback(() => {
        setIsFileTypeValid(true);
    }, []);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: (acceptedFiles) => {
            const validFiles = acceptedFiles.filter(file => validFileType(file.type));
            onDrop(validFiles);
        },
        onDragEnter,
        onDragLeave,
        accept: {
            "text/*": [".csv", ".tsv", ".txt"],
            "application/gzip": [".csv.gz", ".tsv.gz", ".txt.gz"],
            "application/x-gzip": [".csv.gz", ".tsv.gz", ".txt.gz"]
        }
    });

    return (
        <div {...getRootProps()} className="filedropArea">
            <input {...getInputProps()} className="fileDrop" />
            <h4>{title}</h4>
            <span>Drop file here or</span>
            <button className="openFilesystemButton">
                <span>Browse</span>
            </button>
            {isDragActive ? (
                <p>{isFileTypeValid ? <span style={{ color: 'green' }}>Drop here</span> : <span style={{ color: 'red' }}>Invalid file type</span>}</p>
            ) : (
                <p></p>
            )}
            {fileName && <p className="fileName">{fileName}</p>}
        </div>
    );
};

const AnalysisInputForm = ({
    setIsVisible,
    onDatasetSelect,
    contrastGroup,
    referenceGroup,
    onAddGroup,
    onUpdateGroup,
    selectedSamples,
    onAddSamplesToGroup,
    runAnalysis,
    isLoading
}) => {
    const [countsFileName, setCountsFileName] = useState('');
    const [coldataFileName, setColdataFileName] = useState('');

    // const cleanData = (data) => {
    //     return data
    //         .map(item => item && item.trim()) // Trim whitespace
    //         .filter(item => item); // Filter out empty strings or null values
    // };

    const handleRemoveContrastSample = (sampleId) => {
        const updatedContrastGroup = {
            ...contrastGroup,
            samples: contrastGroup.samples.filter((sample) => sample.id !== sampleId),
        };
        setContrastGroup(updatedContrastGroup);
    };

    const handleRemoveReferenceSample = (sampleId) => {
        const updatedReferenceGroup = {
            ...referenceGroup,
            samples: referenceGroup.samples.filter((sample) => sample.id !== sampleId),
        };
        setReferenceGroup(updatedReferenceGroup);
    };

    const decompressAndParseFile = useCallback((file, onParsed) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            let content = event.target.result;
            if (file.name.endsWith('.gz')) {
                // Decompress gzip content
                const compressed = new Uint8Array(content);
                content = pako.inflate(compressed, { to: 'string' });
            }
            const data = Papa.parse(content, { header: true, delimiter: '\t' }).data;
            onParsed(data);
        };
        reader.readAsArrayBuffer(file);
    }, []);

    const onDropCounts = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setCountsFileName(acceptedFiles[0].name);
        }
    }, []);

    const onDropColdata = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setColdataFileName(acceptedFiles[0].name);
            const file = acceptedFiles[0];

            decompressAndParseFile(file, (data) => {
                const conditions = cleanData([...new Set(data.map(item => item.condition))]);
                const phases = cleanData([...new Set(data.map(item => item.phase))]);
                const contrasts = [...conditions, ...phases];

                setReferenceLevels({ conditions, phases });
                setContrastLevels(cleanData(contrasts));
            });
        }
    }, [decompressAndParseFile]);

    const handleButtonClick = (datasetType) => {
        if (datasetType === 'external') {
            setIsVisible(true); // Show plot area when 'Use External Dataset' is selected
        } else {
            setIsVisible(false); // Hide plot area when 'Use Example Dataset' is selected
            // Load example dataset
            onDatasetSelect('example', null);
        }
    };

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data</h3>
            <div className='dataSubfield'>
                {isLoading ? (
                    <div className="loader"></div>
                ) : (
                    <>
                        <button
                            className="analysisInputButton"
                            onClick={() => onDatasetSelect('example')}
                        >
                            Use Example Dataset
                        </button>
                        <button
                            className="analysisInputButton"
                            onClick={() => handleButtonClick('external')}
                        >
                            Use External Dataset
                        </button>
                    </>
                )}
            </div>
            <div id="filedropContainer">
                <FileDropArea
                    title="Counts"
                    onDrop={onDropCounts}
                    fileName={countsFileName}
                />
                <FileDropArea
                    title="Coldata"
                    onDrop={onDropColdata}
                    fileName={coldataFileName}
                />
            </div>
            <h3>Configuration</h3>
            <div>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="covariates" />
                    <span>Add covariates</span>
                </label>
                <label className="radioLabel">
                    <span id="adjustmentSubfield">Adjustment method:</span>
                    <select id="adjustmentMethod" name="adjustmentMethod">
                        <option value="option1">Bonferroni</option>
                        <option value="option2">Benjamini and Hochberg</option>
                    </select>
                </label>
                <div className='dataSubfieldSampleField'>
                    <SampleField
                        headerName="Reference Group"
                        samples={referenceGroup.samples}
                        onRemoveSample={handleRemoveReferenceSample}
                    />
                    <SampleField
                        headerName="Contrast Group"
                        samples={contrastGroup.samples}
                        onRemoveSample={handleRemoveContrastSample}
                    />
                </div>
                <label className="radioLabel">
                    <span>Data Exploration Transform:</span>
                    <select id="transformationMethod" name="transformationMethod">
                        <option value="option1">VST</option>
                        <option value="option2">log2(counts + 1)</option>
                        <option value="option3">ln(counts + 1)</option>
                        <option value="option4">log10(counts + 1)</option>
                    </select>
                </label>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="name" />
                    <span>Batch correction</span>
                </label>
            </div>
            <div id="runAnalysisContainer">
                <IconButton icon={terminal} label="Run Analysis" onClick={runAnalysis} />
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
    onDatasetSelect: PropTypes.func.isRequired,
    contrastGroup: PropTypes.object.isRequired,
    referenceGroup: PropTypes.object.isRequired,
    onAddGroup: PropTypes.func.isRequired,
    onUpdateGroup: PropTypes.func.isRequired,
    selectedSamples: PropTypes.array.isRequired,
    onAddSamplesToGroup: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    isLoading: PropTypes.bool.isRequired,
};

export default AnalysisInputForm;
