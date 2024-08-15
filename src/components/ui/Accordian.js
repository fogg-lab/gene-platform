import React, { useState } from 'react';
import PropTypes from 'prop-types';

const Accordion = ({ data }) => {
    const [activeIndex, setActiveIndex] = useState(null);

    const handleAccordionClick = (index) => {
        setActiveIndex(activeIndex === index ? null : index);
    };

    const accordionStyle = {
        width: '80%',
        margin: '0 auto',
        border: '1px solid #ccc',
        borderRadius: '5px',
    };

    const itemStyle = {
        borderBottom: '1px solid #ccc',
    };

    const headerStyle = {
        cursor: 'pointer',
    };

    const questionBarStyle = {
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        backgroundColor: '#D73F09',
        padding: '15px',
        borderRadius: '5px 5px 0 0',
        transition: 'background-color 0.3s ease',
    };

    const questionBarHoverStyle = {
        ...questionBarStyle,
        backgroundColor: '#D73F09',
    };

    const questionStyle = {
        color: 'white',
        fontSize: '1.2em',
    };

    const iconStyle = {
        height: '20px',
        width: '20px',
        transition: 'transform 0.3s ease',
    };

    const contentStyle = {
        padding: '15px',
        backgroundColor: '#fff',
    };

    const pStyle = {
        margin: '0',
    };

    const iconSrc = require(`../../assets/icons/add_circle.png`).default;

    return (
        <div style={accordionStyle}>
            {data.map((item, index) => (
                <div key={index} style={itemStyle}>
                    <div
                        style={headerStyle}
                        onClick={() => handleAccordionClick(index)}
                        onKeyDown={(e) => {
                            if (e.key === 'Enter' || e.key === ' ') {
                                handleAccordionClick(index);
                            }
                        }}
                        role="button"
                        tabIndex={0}
                    >
                        <div
                            style={activeIndex !== index ? questionBarStyle : questionBarHoverStyle}
                        >
                            {activeIndex !== index && <span style={questionStyle}>{item.question}</span>}
                            {activeIndex !== index && (
                                <img
                                    src={iconSrc}
                                    alt="Add Circle"
                                    style={iconStyle}
                                />
                            )}
                        </div>
                    </div>
                    {activeIndex === index && (
                        <div style={contentStyle}>
                            <p style={pStyle}>{item.answer}</p>
                        </div>
                    )}
                </div>
            ))}
        </div>
    );
};

Accordion.propTypes = {
    data: PropTypes.arrayOf(
        PropTypes.shape({
            question: PropTypes.string.isRequired,
            answer: PropTypes.string.isRequired,
        })
    ).isRequired,
};

export default Accordion;
