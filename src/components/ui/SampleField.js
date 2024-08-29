import React, { useState, useRef, useEffect } from 'react';
import PropTypes from 'prop-types';

function SampleField({ headerName, groups, onAddGroup, onUpdateGroup }) {
    const [editingId, setEditingId] = useState(null);
    const [expandedSamples, setExpandedSamples] = useState({});
    const editInputRef = useRef(null);

    useEffect(() => {
        if (editingId !== null && editInputRef.current) {
            editInputRef.current.focus();
        }
    }, [editingId]);

    useEffect(() => {
        // Log the full structure of groups
        console.log('Groups data:', groups);

        // Log each group's samples for more detail
        groups.forEach(group => {
            console.log(`Group ${group.id} - ${group.name}:`, group.samples);
        });
    }, [groups]);

    const startEditing = (id) => {
        setEditingId(id);
    };

    const handleNameChange = (id, newName) => {
        onUpdateGroup(id, { name: newName });
    };

    const stopEditing = () => {
        setEditingId(null);
    };

    const handleRemoveSample = (groupId, sampleId) => {
        const group = groups.find(g => g.id === groupId);
        if (group) {
            const updatedSamples = group.samples.filter(sample => sample.id !== sampleId);
            onUpdateGroup(groupId, { samples: updatedSamples });
        }
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
                <button className='squareButton' onClick={onAddGroup}>
                    +
                </button>
            </div>
            <div className='groupsList'>
                {groups.map(group => (
                    <div key={group.id} className='groupItem'>
                        {editingId === group.id ? (
                            <input
                                ref={editInputRef}
                                type="text"
                                value={group.name}
                                onChange={(e) => handleNameChange(group.id, e.target.value)}
                                onBlur={stopEditing}
                                onKeyPress={(e) => {
                                    if (e.key === 'Enter') {
                                        stopEditing();
                                    }
                                }}
                            />
                        ) : (
                            <button
                                className="groupNameButton"
                                onClick={() => startEditing(group.id)}
                            >
                                {group.name}
                            </button>
                        )}
                        <ul>
                            {group.samples.map((sample) => {
                                return (
                                    <li key={sample.id} className='sampleListElement'>
                                        <span className="sampleText">
                                            <strong>{sample.name}</strong>
                                        </span>
                                        <button className="removeButton" onClick={() => handleRemoveSample(group.id, sample.id)}>
                                            X
                                        </button>
                                        {renderSampleDetails(sample)}
                                    </li>
                                );
                            })}
                        </ul>
                    </div>
                ))}
            </div>
        </div >
    );
}

SampleField.propTypes = {
    headerName: PropTypes.string.isRequired,
    groups: PropTypes.arrayOf(PropTypes.shape({
        id: PropTypes.number.isRequired,
        name: PropTypes.string.isRequired,
        samples: PropTypes.arrayOf(PropTypes.object).isRequired,
    })).isRequired,
    onAddGroup: PropTypes.func.isRequired,
    onUpdateGroup: PropTypes.func.isRequired,
};

export default SampleField;
