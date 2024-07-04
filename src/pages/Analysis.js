import React from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import IconButton from '../components/ui/IconButton';
import ToolTip from '../components/ui/ToolTip';
import PlotArea from '../components/ui/PlotArea';

const Analysis = () => {
    return (
        <div id="analysis_container">
            <div id="analysis_user_input">
                <AnalysisInputForm />
            </div>
            <div id="analysis_visualization_section">
                <div id="analysis_tab_nav">
                    <TabButton label="Exploratory" />
                    <TabButton label="Differential Expression Analysis" />
                    <TabButton label="Differential Expression Results" />
                    <TabButton label="Gene Set Enrichment Analysis" />
                    <TabButton label="Analysis Log" />
                    {/* <ToolTip /> */}
                </div>
                <div id="analysis_run_section">
                    <PlotArea />
                </div>
            </div>
        </div>
    );
};

export default Analysis;