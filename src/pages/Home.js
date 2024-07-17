import React from 'react';
import { useNavigate } from 'react-router-dom';
import IconButton from '../components/ui/IconButton';
import logo from '../assets/fogg_logo.svg';

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
			{/* <h1>GENE Platform</h1> */}
			<img src={logo} id="fogg_logo_home" alt="Fogg Labs Logo" />
			<h2>Gene Expression Explorer</h2>
			<div id="homepage_buttons">
				<IconButton iconFilename="biotech.png" label="Start Analysis" onClick={handleAnalysisClick} />
				<IconButton iconFilename="quick_reference.png" label="Guide / FAQ" onClick={handleGuideClick} />
			</div>
		</div>
	);
};

export default Home;
