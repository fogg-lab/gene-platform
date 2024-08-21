import React, { useState, useRef, useEffect } from 'react';

function SampleField() {
    const [groups, setGroups] = useState([]);
    const [groupCounter, setGroupCounter] = useState(1);
    const [editingId, setEditingId] = useState(null);
    const editInputRef = useRef(null);

    useEffect(() => {
        if (editingId !== null && editInputRef.current) {
            editInputRef.current.focus();
        }
    }, [editingId]);

    const handleAddGroup = () => {
        const newGroup = {
            id: groupCounter,
            name: `Group ${groupCounter}`,
            samples: []
        };
        setGroups([...groups, newGroup]);
        setGroupCounter(groupCounter + 1);
    };

    const startEditing = (id) => {
        setEditingId(id);
    };

    const handleNameChange = (id, newName) => {
        setGroups(groups.map(group =>
            group.id === id ? { ...group, name: newName } : group
        ));
    };

    const stopEditing = () => {
        setEditingId(null);
    };

    return (
        <div className='sampleFieldContainer'>
            <div className='sampleGroupHeader'>
                <h3>Groups</h3>
                <button className='squareButton' onClick={handleAddGroup}>
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
                            {group.samples.map((sample, index) => (
                                <li key={index}>{sample}</li>
                            ))}
                        </ul>
                    </div>
                ))}
            </div>
        </div>
    );
}

export default SampleField;