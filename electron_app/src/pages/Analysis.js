import React from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';

const Analysis = () => {
    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                <AnalysisInputForm />
            </div>
            <div id="analysis_visualization_section">
            </div>

        </div>
    );
};

export default Analysis;