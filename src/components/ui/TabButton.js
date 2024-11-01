import React from 'react';
import PropTypes from 'prop-types';

const TabButton = ({ label, onClick, isActive }) => {
    const buttonStyle = {
        backgroundColor: isActive ? '#D73F09' : '#A62E06',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '5px 20px 5px 20px',
        margin: '10px',
        border: 'none',
        color: 'white',
        cursor: 'pointer',
        borderRadius: '5px',
        fontSize: '1em',
        transition: 'all 0.3s ease',
        position: 'relative',
        zIndex: isActive ? 1 : 0,
        transform: isActive ? 'translateY(-6px)' : 'translateY(0)',
    };

    const shadowStyle = {
        content: '""',
        position: 'absolute',
        top: '4px',
        left: '2px',
        right: '-5px',
        bottom: '-5px',
        background: 'rgba(215, 63, 9, 0.5)', // Lighter version of #A62E06
        borderRadius: '5px',
        zIndex: -1,
        opacity: isActive ? 1 : 0,
        transition: 'opacity 0.5s ease',
    };

    return (
        <button
            className={`tab-button ${isActive ? 'active' : ''}`}
            style={buttonStyle}
            onClick={onClick}
        >
            <span>{label}</span>
            <div style={shadowStyle}></div>
        </button>
    );
};

TabButton.propTypes = {
    label: PropTypes.string.isRequired,
    onClick: PropTypes.func.isRequired,
    isActive: PropTypes.bool.isRequired,
};

export default TabButton;