import React from 'react';
import PropTypes from 'prop-types';

const TabButton = ({ label, onClick }) => {
    const buttonStyle = {
        backgroundColor: '#D73F09',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '5px 20px 5px 20px',
        margin: '3px',
        border: 'none',
        color: 'white',
        cursor: 'pointer',
        borderRadius: '20px',
        fontSize: '1em'
    };

    return (
        <button className="tab-button" style={buttonStyle} onClick={onClick}>
            <span>{label}</span>
        </button>
    );
};

TabButton.propTypes = {
    label: PropTypes.string.isRequired,
    onClick: PropTypes.func.isRequired,
};


export default TabButton;
