import React from 'react';
import { useNavigate } from 'react-router-dom';
import IconButton from '../components/ui/IconButton';

const Home = () => {
	const navigate = useNavigate();

	const handleAnalysisClick = () => {
		navigate('/Analysis');
	};

	const handleGuideClick = () => {
		navigate('/Guide');
	};

	return (
		<div id="homepage_main_div">
			<h1>GENE Platform</h1>
			<h2>Gene Expression Explorer</h2>
			<div id="homepage_buttons">
				<IconButton iconFilename="biotech.png" label="Start Analysis" onClick={handleAnalysisClick} />
				<IconButton iconFilename="quick_reference.png" label="Guide / FAQ" onClick={handleGuideClick} />
			</div>
		</div>
	);
};

export default Home;
