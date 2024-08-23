import React, { useState, useRef, useEffect } from 'react';
import PropTypes from 'prop-types';

function SampleField({ headerName, groups, onAddGroup, onUpdateGroup }) {
    const [editingId, setEditingId] = useState(null);
    const editInputRef = useRef(null);

    useEffect(() => {
        if (editingId !== null && editInputRef.current) {
            editInputRef.current.focus();
        }
    }, [editingId]);

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
                            {group.samples.map((sample) => (
                                <li key={sample.id}>
                                    {sample.name}
                                    <button onClick={() => handleRemoveSample(group.id, sample.id)}>
                                        X
                                    </button>
                                </li>
                            ))}
                        </ul>
                    </div>
                ))}
            </div>
        </div>
    );
}

SampleField.propTypes = {
    headerName: PropTypes.string.isRequired,
    groups: PropTypes.arrayOf(PropTypes.shape({
        id: PropTypes.number.isRequired,
        name: PropTypes.string.isRequired,
        samples: PropTypes.arrayOf(PropTypes.shape({
            id: PropTypes.oneOfType([PropTypes.string, PropTypes.number]).isRequired,
            name: PropTypes.string.isRequired,
        })).isRequired,
    })).isRequired,
    onAddGroup: PropTypes.func.isRequired,
    onUpdateGroup: PropTypes.func.isRequired,
};

export default SampleField;