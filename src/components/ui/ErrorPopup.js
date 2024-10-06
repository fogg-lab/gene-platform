import React, { createContext, useState, useContext, useCallback } from 'react';
import PropTypes from 'prop-types';
import closeIcon from '../../assets/icons/close.svg';

const ErrorPopupContext = createContext();

export const useErrorPopup = () => {
    const context = useContext(ErrorPopupContext);
    if (!context) {
        throw new Error('useErrorPopup must be used within an ErrorPopupProvider');
    }
    return context;
};

export const ErrorPopupProvider = ({ children }) => {
    const [isVisible, setIsVisible] = useState(false);
    const [errorMessage, setErrorMessage] = useState('');

    const showError = useCallback((message) => {
        setErrorMessage(message);
        setIsVisible(true);
    }, []);

    const hideError = useCallback(() => {
        setIsVisible(false);
    }, []);

    return (
        <ErrorPopupContext.Provider value={{ showError, hideError }}>
            {children}
            <ErrorPopup isVisible={isVisible} errorMessage={errorMessage} onClose={hideError} />
        </ErrorPopupContext.Provider>
    );
};

ErrorPopupProvider.propTypes = {
    children: PropTypes.node.isRequired,
};

const ErrorPopup = ({ isVisible, errorMessage, onClose }) => {
    let iconSrc = require(`../../assets/icons/warning.svg`).default;

    const iconStyle = {
        margin: '0px',
        padding: '0px',
        width: '40%',
        height: '40%',
    };

    if (!isVisible) return null;

    return (
        <div className="errorPopup">
            <div
                className="closeButton"
                onClick={onClose}
                role="button"
                tabIndex="0"
            >
                <img src={closeIcon} alt="Close" />
            </div>
            <div className="errorPopupContent">
                <img src={iconSrc} style={iconStyle} alt="Warning" />
                <p>{errorMessage}</p>
            </div>
        </div>
    );
};

ErrorPopup.propTypes = {
    isVisible: PropTypes.bool.isRequired,
    errorMessage: PropTypes.string.isRequired,
    onClose: PropTypes.func.isRequired,
};

export default ErrorPopup;