import React from 'react'
import IconButton from '../ui/IconButton';

const AnalysisInput = ({/* Put parameters for different inputs here later(ie, which tab and change accordingly) */ }) => {

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
                <div className="filedropArea">
                    <h4>Counts</h4>
                    <span>Drop .tsv file here or</span>
                    <input type="file" title=" " name="myFile" />
                    <button className="openFilesystemButton">
                        <span>Browse</span>
                    </button>
                </div>
                <div className="filedropArea">
                    <h4>Coldata</h4>
                    <span>Drop .tsv file here or</span>
                    <input className="fileDrop" type="file" title=" " name="myFile" />
                    <button className="openFilesystemButton">
                        <span>Browse</span>
                    </button>
                </div>
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
                        {/* &nbsp; to add whitespace is not professional. I will figure out a better way to do this
                        vertical alignment later. -Lucas */}
                        <span id="adjustmentSubfield">Contrast level:&nbsp;&nbsp;</span>
                        <input type="text" id="textInput" name="textInput" />
                    </label>
                    <label className="radioLabel">
                        <span id="adjustmentSubfield">Reference level:</span>
                        <input type="text" id="textInput" name="textInput" />
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