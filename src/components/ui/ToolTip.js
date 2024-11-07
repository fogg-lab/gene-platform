import React from 'react';

const ToolTip = ({ content }) => {
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
        position: 'relative',
    };

    const tooltipTextStyle = {
        visibility: 'hidden',
        width: '200px',
        backgroundColor: '#555',
        color: 'white',
        textAlign: 'center',
        padding: '1px 1px',
        borderRadius: '6px',
        position: 'fixed',
        zIndex: 9999,
        opacity: 0,
        transition: 'opacity 0.3s',
    };

    const handleMouseMove = (e) => {
        const tooltip = e.currentTarget.querySelector('.tooltip-text');
        if (tooltip) {
            tooltip.style.left = `${e.clientX + 20}px`;
            tooltip.style.top = `${e.clientY}px`;
        }
    };

    return (
        <div className="tooltip-wrapper" style={{ position: 'relative' }}>
            <button
                className="tooltip"
                style={tooltipStyle}
                onMouseMove={handleMouseMove}
            >
                <img src={iconSrc} style={iconStyle} alt="help" />
                <span className="tooltip-text" style={tooltipTextStyle}>{content}</span>
            </button>
            <style>
                {`
                    .tooltip:hover .tooltip-text {
                        visibility: visible !important;
                        opacity: 1 !important;
                    }
                `}
            </style>
        </div>
    );
};

export default ToolTip;