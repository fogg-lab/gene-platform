import React from 'react'
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
import { useCallback, useState } from 'react';
import Papa from 'papaparse';



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
        onDrop: (acceptedFiles, rejectedFiles) => {
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


const AnalysisInput = ({/* Put parameters for different inputs here later (ie which tab and change accordingly) */ }) => {

    const [countsFileName, setCountsFileName] = useState('');
    const [coldataFileName, setColdataFileName] = useState('');
    const [referenceLevels, setReferenceLevels] = useState({
        conditions: [],
        phases: []
    });
    const [contrastLevels, setContrastLevels] = useState([]);

    const cleanData = (data) => {
        return data
            .map(item => item && item.trim()) // Trim whitespace
            .filter(item => item); // Filter out empty strings or null values
    };

    const onDropCounts = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setCountsFileName(acceptedFiles[0].name);
        }
    }, []);

    const onDropColdata = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setColdataFileName(acceptedFiles[0].name);
            const file = acceptedFiles[0];
            const reader = new FileReader();
            reader.onload = (event) => {
                const text = event.target.result;
                const data = Papa.parse(text, { header: true, delimiter: '\t' }).data;

                const conditions = cleanData([...new Set(data.map(item => item.condition))]);
                const phases = cleanData([...new Set(data.map(item => item.phase))]);
                const contrasts = [...conditions, ...phases];

                setReferenceLevels({ conditions, phases });
                setContrastLevels(cleanData(contrasts));
            };
            reader.readAsText(file);
        }
    }, []);

    return (
        <div id="analysisInputContainer_comp">
            <h3>Data</h3>
            <div className='dataSubfield'>
                {/* Load example data / External dataset */}
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="exampleDataset" />
                    <span>Use Example Dataset</span>
                </label>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="exampleDataset" />
                    <span>Use External Dataset</span>
                </label>
            </div>
            <div id="filedropContainer">
                <FileDropArea title="Counts" onDrop={onDropCounts} fileName={countsFileName} />
                <FileDropArea title="Coldata" onDrop={onDropColdata} fileName={coldataFileName} />
                {/* <div className="filedropArea">
                    <h4>Counts</h4>
                    <span>Drop .tsv/.csv file here or</span>
                    <input type="file" title=" " name="myFile" />
                    <button className="openFilesystemButton">
                        <span>Browse</span>
                    </button>
                </div>
                <div className="filedropArea">
                    <h4>Coldata</h4>
                    <span>Drop .tsv/.csv file here or</span>
                    <input className="fileDrop" type="file" title=" " name="myFile" />
                    <button className="openFilesystemButton">
                        <span>Browse</span>
                    </button>
                </div> */}
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
                <div className='dataSubfield'>
                    <label className="radioLabel">
                        <span id="adjustmentSubfield">Adjustment method:</span>
                        <select id="exampleDropdown" name="exampleDropdown">
                            <option value="option1">Bonferroni</option>
                            <option value="option2">Benjamini and Hochberg</option>
                        </select>
                    </label>
                    <label className="radioLabel">
                        <span id="adjustmentSubfield">Contrast level:</span>
                        <select id="exampleDropdown" name="contrastLevel">
                            {contrastLevels.map(level => (
                                <option key={level} value={level}>{level}</option>
                            ))}
                        </select>
                    </label>
                    {/* <label className="radioLabel">
                        <span id="adjustmentSubfield">Reference level:</span>
                        <select id="exampleDropdown" name="exampleDropdown">
                            <option value="option1">Example</option>
                            <option value="option2">Example</option>
                        </select>
                    </label> */}
                    <label className="radioLabel">
                        <span id="adjustmentSubfield">Reference level:</span>
                        <select id="exampleDropdown" name="exampleDropdown">
                            {referenceLevels.conditions.map(condition => (
                                <option key={condition} value={condition}>{condition}</option>
                            ))}
                            {referenceLevels.phases.map(phase => (
                                <option key={phase} value={phase}>{phase}</option>
                            ))}
                        </select>
                    </label>
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
                <IconButton iconFilename="terminal.png" label="Run Analysis" />
            </div>
        </div >
    );
};



export default AnalysisInput;