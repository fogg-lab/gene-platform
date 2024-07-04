import React, { useEffect, useRef } from 'react';
import Plotly from 'plotly.js';

const PlotArea = ({ data, layout, config }) => {
    const plotRef = useRef(null);

    useEffect(() => {
        Plotly.newPlot(plotRef.current, data, layout, config);

        return () => {
            Plotly.purge(plotRef.current);
        };
    }, [data, layout, config]);

    return <div ref={plotRef} style={{ width: '500px', height: '400' }} />;
};

export default PlotArea;
