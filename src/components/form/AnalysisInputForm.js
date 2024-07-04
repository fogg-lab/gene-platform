import React from 'react'
const AnalysisInput = ({/* Put parameters for different inputs here later(ie, which tab and change accordingly) */ }) => {

    return (
        <div>
            <h3>Data</h3>
            <div className='dataSubfield'>
                {/* Load example data / External dataset */}
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="name" />
                    <span>Use Example Dataset</span>
                </label>
                <label className="radioLabel">
                    <input className="radioInput" type="radio" name="name" />
                    <span>Use External Dataset</span>
                </label>
            </div>
            <div id="filedropContainer">
                <div className="filedropArea">
                    <h2>Counts</h2>
                    <span>Drop .tsv file here or</span>
                    <input type="file" title=" " name="myFile" />
                    <button className="openFilesystemButton">
                        <span>Browse</span>
                    </button>
                </div>
                <div className="filedropArea">
                    <h2>Coldata</h2>
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
                    <span>Add covariates</span>
                    <input className="radioInput" type="radio" name="name" />
                </label>
                <label className="radioLabel">
                    <span>Select datatype:</span>
                    <select id="exampleDropdown" name="exampleDropdown">
                        <option value="option1">Option 1</option>
                        <option value="option2">Option 2</option>
                        <option value="option3">Option 3</option>
                    </select>
                </label>
                <div>
                    {/* Add adjustment method, etc. */}
                </div>
            </div>
        </div>
    );
};

export default AnalysisInput;