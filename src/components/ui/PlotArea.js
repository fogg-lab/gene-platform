import React, { useEffect, useRef } from 'react';
import PropTypes from 'prop-types';

const PlotArea = ({ htmlContent }) => {
    const iframeRef = useRef(null);

    useEffect(() => {
        if (iframeRef.current && htmlContent) {
            const iframe = iframeRef.current;
            const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
            iframeDoc.open();
            iframeDoc.write(htmlContent);
            iframeDoc.close();
        }
    }, [htmlContent]);

    return (
        <iframe
            ref={iframeRef}
            style={{ width: '100%', height: '100%', border: 'none' }}
            title="Plot"
        />
    );
};

PlotArea.propTypes = {
    htmlContent: PropTypes.string.isRequired,
};

export default PlotArea;
