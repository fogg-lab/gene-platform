import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
import ToolTip from '../ui/ToolTip';
import terminal from '../../assets/icons/terminal.png';
import next from '../../assets/icons/next.svg';
import Papa from 'papaparse';

function validFileType(filetype) {
    return filetype.startsWith("text/") || filetype == "application/gzip" || filetype == "application/x-gzip";
}

const readFileAsText = (file) => {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => resolve(event.target.result);
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
};


async function processUploadedFiles(countsFile, coldataFile) {
    // Parse coldata file
    const coldataText = await readFileAsText(coldataFile);
    let coldataData = Papa.parse(coldataText, { header: true }).data;
    if (Object.keys(coldataData.at(-1)).length === 1) {
        coldataData = coldataData.slice(0, -1);
    }
    const coldataTable = structureColdataTable(coldataData);

    // Parse counts file
    const countsText = await readFileAsText(countsFile);
    let countsData = Papa.parse(countsText, { header: true }).data;
    if (Object.keys(countsData.at(-1)).length === 1) {
        countsData = countsData.slice(0, -1);
    }
    const { countsTable, expression, counts } = structureCountsData(countsData);

    return {
        expression,
        counts,
        countsTable,
        coldataTable,
    };
}

function structureColdataTable(coldataData) {
    const cols = Object.keys(coldataData[0]);
    return {
        cols,
        rows: coldataData.map(row => row.sample_id),
        data: coldataData.map(row => cols.map(col => row[col]))
    };
}

function structureCountsData(countsData) {
    const sampleIds = Object.keys(countsData[0]).slice(1);
    const expressionData = countsData.map(row => sampleIds.map(id => parseInt(row[id])));

    const expression = new Int32Array(expressionData.flat());
    const counts = new Int32Array(expressionData.map((row, i) => row.map((val, j) => expressionData[j][i])).flat());
    const geneIdType = Object.keys(countsData[0])[0];
    const countsTable = {
        cols: [...sampleIds],
        rows: countsData.map(row => row[geneIdType]),
        data: countsData.map(row => sampleIds.map(id => parseInt(row[id])))
    };

    return { countsTable, expression, counts };
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

const EDAInputForm = ({
    setIsVisible,
    onDatasetSelect,
    onRemoveSamplesFromGroup,
    runAnalysis,
    isLoading,
    handleStageChange,
    currentStage,
    edaData
}) => {
    const [countsFile, setCountsFile] = useState(null);
    const [coldataFile, setColdataFile] = useState(null);
    const [countsFileName, setCountsFileName] = useState('');
    const [coldataFileName, setColdataFileName] = useState('');

    const onDropCounts = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setCountsFile(acceptedFiles[0]);
            setCountsFileName(acceptedFiles[0].name);
        }
    }, []);

    const onDropColdata = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setColdataFile(acceptedFiles[0]);
            setColdataFileName(acceptedFiles[0].name);
        }
    }, []);

    const handleButtonClick = (datasetType) => {
        if (datasetType === 'external') {
            setIsVisible(true);
        } else {
            setIsVisible(false);
            onDatasetSelect('example', null);
        }
    };

    const handleRunAnalysis = async () => {
        if (edaData?.plots?.pca) {
            handleStageChange('differential');
        }
        else {
            if (countsFile && coldataFile) {
                try {
                    const processedData = await processUploadedFiles(countsFile, coldataFile);
                    runAnalysis(processedData);
                } catch (error) {
                    console.error("Error processing uploaded files:", error);
                }
            } else {
                runAnalysis();
            }
        }
    };

    return (
        <div id="analysisInputContainer_comp">
            <div className="form-with-tooltips">
                <div className="form-content">
                    {/* Left side: All the form fields */}
                    <h3>Use GDC/GEO data or run an example dataset</h3>
                    <div className='dataSubfield flex flex-col gap-2 border border-black rounded p-4'>
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
                                <div className="border-b border-black w-full"></div>
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
                    <div className="opacity-50 pointer-events-none">
                        <label className="radioLabel">
                            <span>Data Exploration Transform (🚧):</span>
                            <select id="transformationMethod" name="transformationMethod">
                                <option value="option1">VST</option>
                                <option value="option2">log2(counts + 1)</option>
                                <option value="option3">ln(counts + 1)</option>
                                <option value="option4">log10(counts + 1)</option>
                            </select>
                        </label>
                    </div>
                    <div id="runAnalysisContainer" style={{ marginTop: '50px' }}>
                        <IconButton
                            icon={edaData?.plots?.pca ? next : terminal}
                            label={edaData?.plots?.pca ? "Next Stage" : "Run Analysis"}
                            onClick={handleRunAnalysis}
                        />
                    </div>
                </div>
                <div className="tooltips-column">
                    <div className="tooltip-row" style={{ marginTop: '40px' }}>
                        <ToolTip content="Choose between example and external datasets" />
                    </div>
                    <div className="tooltip-row" style={{ marginTop: '80px' }}>
                        <ToolTip content="Upload your gene expression matrix and sample metadata files" />
                    </div>
                    <div className="tooltip-row" style={{ marginTop: '150px' }}>
                        <ToolTip content="Select the transformation method to normalize your count data" />
                    </div>
                </div>
            </div>
        </div>
    );
};

FileDropArea.propTypes = {
    title: PropTypes.string.isRequired,
    onDrop: PropTypes.func.isRequired,
    fileName: PropTypes.string.isRequired,
};

EDAInputForm.propTypes = {
    setIsVisible: PropTypes.func.isRequired,
    onDatasetSelect: PropTypes.func.isRequired,
    onRemoveSamplesFromGroup: PropTypes.func.isRequired,
    runAnalysis: PropTypes.func.isRequired,
    isLoading: PropTypes.bool.isRequired,
};

export default EDAInputForm;
