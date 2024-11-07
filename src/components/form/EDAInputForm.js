import React, { useCallback, useState } from 'react';
import PropTypes from 'prop-types';
import { useDropzone } from 'react-dropzone';
import IconButton from '../ui/IconButton';
import terminal from '../../assets/icons/terminal.png';
import next from '../../assets/icons/next.svg';
import pako from 'pako';
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
    const coldataData = Papa.parse(coldataText, { header: true }).data;
    const coldataTable = structureColdataTable(coldataData);

    // Parse counts file
    const countsText = await readFileAsText(countsFile);
    const countsData = Papa.parse(countsText, { header: true }).data;
    const { countsTable, genesTable, expression, counts } = structureCountsData(countsData);

    return {
        expression,
        counts,
        countsTable,
        coldataTable,
        genesTable,
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
    const sampleIds = Object.keys(countsData[0]).slice(2);
    const genesData = countsData.map(row => [row['Ensembl gene'], row['Symbol']]);
    const expressionData = countsData.map(row => sampleIds.map(id => parseInt(row[id])));

    const expression = new Int32Array(expressionData.flat());
    const counts = new Int32Array(expressionData.map((row, i) => row.map((val, j) => expressionData[j][i])).flat());

    const countsTable = {
        cols: ['Ensembl gene', 'Symbol', ...sampleIds],
        rows: genesData.map(gene => gene[0]),
        data: countsData.map(row => [row['Ensembl gene'], row['Symbol'], ...sampleIds.map(id => parseInt(row[id]))])
    };

    const genesTable = {
        cols: ['ensembl_gene', 'symbol'],
        rows: genesData.map(gene => gene[0]),
        data: genesData
    };

    return { countsTable, genesTable, expression, counts };
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
        <div {...getRootProps()} className="filedropArea filedropArea-disabled">
            <input {...getInputProps()} className="fileDrop" disabled />
            <h4>{title}</h4>
            <span>Drop file here or</span>
            <button className="openFilesystemButton" disabled>
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
            // If analysis is complete, move to next stage
            handleStageChange('differential');
        }
        else {
            if (countsFile && coldataFile) {
                try {
                    const processedData = await processUploadedFiles(countsFile, coldataFile);
                    runAnalysis(processedData);
                } catch (error) {
                    console.error("Error processing uploaded files:", error);
                    // Handle error (e.g., show error message to user)
                }
            } else {
                runAnalysis();
            }
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
                    <span>Data Exploration Transform (🚧):</span>
                    <select id="transformationMethod" name="transformationMethod">
                        <option value="option1">VST</option>
                        <option value="option2">log2(counts + 1)</option>
                        <option value="option3">ln(counts + 1)</option>
                        <option value="option4">log10(counts + 1)</option>
                    </select>
                </label>
            </div>
            <div id="runAnalysisContainer">
                <IconButton
                    icon={edaData?.plots?.pca ? next : terminal}
                    label={edaData?.plots?.pca ? "Next Stage" : "Run Analysis"}
                    onClick={handleRunAnalysis} />
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
