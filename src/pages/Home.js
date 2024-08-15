import React from 'react';
import { useNavigate } from 'react-router-dom';
import IconButton from '../components/ui/IconButton';
import logo from '../assets/fogg_logo.svg';
import biotech from '../assets/icons/biotech.png';
import quick_reference from '../assets/icons/quick_reference.png';

const Home = () => {
	const navigate = useNavigate();

	const handleAnalysisClick = () => {
		console.log("HERE. handleAnalysisClick");
		navigate('/Analysis');
	};

	const handleGuideClick = () => {
		navigate('/Guide');
	};

	return (
		<div id="homepage_main_div">
			<img src={logo} id="fogg_logo_home" alt="GENE Logo" />
			<h2>Gene Expression Explorer</h2>
			<div id="homepage_buttons">
				<IconButton icon={biotech} label="Start Analysis" onClick={handleAnalysisClick} />
				<IconButton icon={quick_reference} label="Guide / FAQ" onClick={handleGuideClick} />
			</div>
		</div>
	);
};

export default Home;
