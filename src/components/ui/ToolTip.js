import React from 'react';

const ToolTip = ({ }) => {
    let iconSrc = require(`../../assets/icons/help.png`).default;

    const iconStyle = {
        margin: '0px',
        padding: '0px',
        width: '20px',
        height: '20px',
    };

    const tooltipStyle = {
        backgroundColor: 'grey',
        padding: '0px',
        margin: '0px',
        width: '24px',
        height: '24px',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        border: 'none',
        color: 'grey',
        cursor: 'pointer',
        borderRadius: '200px',
    };

    return (
        <button className="tooltip" style={tooltipStyle}>
            <img src={iconSrc} style={iconStyle} />
        </button>
    );
};

export default ToolTip;