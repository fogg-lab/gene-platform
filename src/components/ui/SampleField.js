import React, { useState } from 'react';
import PropTypes from 'prop-types';

function SampleField({ headerName, samples, onRemoveSample }) {
    const [expandedSamples, setExpandedSamples] = useState({});

    const handleRemoveSample = (sampleId) => {
        onRemoveSample(sampleId);
    };

    const toggleSampleDetails = (sampleId) => {
        setExpandedSamples(prev => ({
            ...prev,
            [sampleId]: !prev[sampleId]
        }));
    };


    const renderSampleDetails = (sample) => {
        const excludeKeys = ['id', 'name', 'sample'];
        const isExpanded = expandedSamples[sample.id];

        return (
            <div className="sampleDetails">
                <button onClick={() => toggleSampleDetails(sample.id)}>
                    {isExpanded ? 'Collapse' : 'Expand'} Details
                </button>
                {isExpanded && (
                    <div className="">
                        {Object.entries(sample).map(([key, value]) => {
                            if (!excludeKeys.includes(key)) {
                                return (
                                    <div key={key}>
                                        <strong>{key}:</strong> {
                                            typeof value === 'object' ? JSON.stringify(value) : value
                                        }
                                    </div>
                                );
                            }
                            return null;
                        })}
                    </div>
                )}
            </div>
        );
    };

    return (
        <div className='sampleFieldContainer'>
            <div className='sampleGroupHeader'>
                <h3>{headerName}</h3>
            </div>
            <div className='groupsList'>
                {samples && (
                    <ul>

                        {samples.map((sample) => {
                            const firstEntry = Object.entries(sample)[0];
                            const [firstKey, firstValue] = firstEntry || [];

                            return (
                                <li key={sample.id} className='sampleListElement'>
                                    <span className="sampleText">
                                        <strong>{firstValue}</strong>
                                    </span>
                                    <button className="removeButton" onClick={() => handleRemoveSample(sample.id)}>
                                        X
                                    </button>
                                    {/* {renderSampleDetails(sample)} */}
                                </li>
                            );
                        })}
                    </ul>
                )}
            </div>
        </div>
    );
}

SampleField.propTypes = {
    headerName: PropTypes.string.isRequired,
    samples: PropTypes.arrayOf(PropTypes.object).isRequired,
    onRemoveSample: PropTypes.func.isRequired,
};

export default SampleField;