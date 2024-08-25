import React, { useState, useEffect } from 'react';
import { getExternalDataset } from '../../services/api';

const ExternalDatasetSelector = ({ onDatasetSelect }) => {
    const [dataSrc, setDataSrc] = useState('GDC');
    const [datasets, setDatasets] = useState([]);
    const [selectedDataset, setSelectedDataset] = useState('');

    useEffect(() => {
        // Load datasets from GDC.json or GEO.json
        import(`../../assets/external_data_index/${dataSrc}.json`)
            .then(data => setDatasets(Object.entries(data.datasets)))
            .catch(error => console.error('Error loading datasets:', error));
    }, [dataSrc]);

    const handleDatasetSelect = async () => {
        if (selectedDataset) {
            try {
                const data = await getExternalDataset(dataSrc, selectedDataset);
                onDatasetSelect('external', data);
            } catch (error) {
                console.error('Error fetching external dataset:', error);
            }
        }
    };

    return (
        <div>
            <select value={dataSrc} onChange={(e) => setDataSrc(e.target.value)}>
                <option value="GDC">GDC</option>
                <option value="GEO">GEO</option>
            </select>
            <select value={selectedDataset} onChange={(e) => setSelectedDataset(e.target.value)}>
                <option value="">Select a dataset</option>
                {datasets.map(([id, info]) => (
                    <option key={id} value={id}>{info.title}</option>
                ))}
            </select>
            <button onClick={handleDatasetSelect}>Load Dataset</button>
        </div>
    );
};

export default ExternalDatasetSelector;
