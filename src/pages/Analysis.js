import React, { useState, useMemo, useRef, useEffect } from 'react';
import AnalysisInputForm from '../components/form/AnalysisInputForm';
import TabButton from '../components/ui/TabButton';
import IconButton from '../components/ui/IconButton';
import ToolTip from '../components/ui/ToolTip';
import PlotArea from '../components/ui/PlotArea';
import DataTable from '../components/ui/DataTable';

const Analysis = () => {
    const [activeTab, setActiveTab] = useState('table');
    const [tableScrollPosition, setTableScrollPosition] = useState(0);
    const tableContainerRef = useRef(null);

    const data = useMemo(() => {
        const rows = [];
        for (let i = 0; i < 1000; i++) {
            rows.push({
                id: i + 1,
                gene: `ENSG${String(i + 1).padStart(11, '0')}`,
                symbol: `GENE${i + 1}`,
                baseMean: Math.floor(Math.random() * 10000),
                log2FoldChange: (Math.random() * 4 - 2).toFixed(6)
            });
        }
        return rows;
    }, []);

    const columns = [
        { key: 'id', name: 'ID' },
        { key: 'gene', name: 'Gene' },
        { key: 'symbol', name: 'Symbol' },
        { key: 'baseMean', name: 'Base Mean' },
        { key: 'log2FoldChange', name: 'Log2 Fold Change' },
    ];

    useEffect(() => {
        if (activeTab === 'table' && tableContainerRef.current) {
            tableContainerRef.current.scrollTop = tableScrollPosition;
        }
    }, [activeTab, tableScrollPosition]);

    const handleTableScroll = (event) => {
        setTableScrollPosition(event.target.scrollTop);
    };

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
                </div>
                <div id="analysis_content">
                    <div id="view_toggle">
                        <button 
                            className={`view-toggle-btn ${activeTab === 'table' ? 'active' : ''}`}
                            onClick={() => setActiveTab('table')}
                        >
                            Table View
                        </button>
                        <button 
                            className={`view-toggle-btn ${activeTab === 'plot' ? 'active' : ''}`}
                            onClick={() => setActiveTab('plot')}
                        >
                            Plot View
                        </button>
                    </div>
                    <div id="view_content">
                        <div 
                            style={{ display: activeTab === 'table' ? 'block' : 'none', height: '600px', overflow: 'auto' }}
                            ref={tableContainerRef}
                            onScroll={handleTableScroll}
                        >
                            <DataTable data={data} columns={columns} />
                        </div>
                        <div style={{ display: activeTab === 'plot' ? 'block' : 'none' }}>
                            <PlotArea />
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default Analysis;
