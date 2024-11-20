import React from 'react';
import PropTypes from 'prop-types';
import lockIcon from '../../assets/icons/lock.svg';

const TabButton = ({ label, onClick, isActive, isLocked }) => {

    const iconStyle = {
        margin: '0 8px 0 0',
        padding: '0px',
        width: '20px',
        height: '20px',
    };

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
        background: 'rgba(215, 63, 9, 0.5)',
        borderRadius: '5px',
        zIndex: -1,
        opacity: isActive ? 1 : 0,
        transition: 'opacity 0.5s ease',
    };

    const handleClick = (e) => {
        if (!isLocked) {
            onClick(e);
        }
    };

    return (
        <button
            className={`tab-button ${isActive ? 'active' : ''} ${isLocked ? 'locked' : ''}`}
            style={buttonStyle}
            onClick={handleClick}
            disabled={isLocked}
        >
            {isLocked && <img src={lockIcon} style={iconStyle} alt="lock" />}
            <span>{label}</span>
            <div style={shadowStyle}></div>
        </button>
    );
};

TabButton.propTypes = {
    label: PropTypes.string.isRequired,
    onClick: PropTypes.func.isRequired,
    isActive: PropTypes.bool.isRequired,
    isLocked: PropTypes.bool,
};

TabButton.defaultProps = {
    isLocked: false,
};

export default TabButton;