import React from 'react';
import { useRouter } from 'next/router';
import IconButton from '../components/ui/IconButton';
import logo from '../assets/fogg_logo.svg';

const Home = () => {
	const router = useRouter();

	const handleAnalysisClick = () => {
		console.log("HERE. handleAnalysisClick");
		router.push('/Analysis'); // Navigate using Next.js router
	};

	const handleGuideClick = () => {
		router.push('/Guide'); // Navigate using Next.js router
	};

	return (
		<div id="homepage_main_div">
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
