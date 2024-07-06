import React, { useEffect, useRef, useCallback } from 'react';
import PropTypes from 'prop-types';

const PlotArea = ({ data = [], layout, config }) => {
    const plotRef = useRef(null);

    const initializePlot = useCallback(() => {
        if (plotRef.current) {
            import('plotly.js').then(Plotly => {
                Plotly.newPlot(plotRef.current, data, layout, config);
            });
        }
    }, [data, layout, config]);

    useEffect(() => {
        // Delay initialization to ensure DOM is ready
        const timer = setTimeout(() => {
            initializePlot();
        }, 0);

        return () => clearTimeout(timer);
    }, [initializePlot]);

    return <div ref={plotRef} style={{ width: '100%', height: '400px' }} />;
};

PlotArea.propTypes = {
    data: PropTypes.array,
    layout: PropTypes.object,
    config: PropTypes.object
};

PlotArea.defaultProps = {
    layout: {},
    config: {}
};

export default PlotArea;
