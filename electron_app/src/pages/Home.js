import React from 'react';
import IconButton from '../components/ui/IconButton';

const Home = () => {
	return (
		<div id="homepage_main_div">
			<h1>GENE Platform</h1>
			<h2>Gene Expression Explorer</h2>
			<div id="homepage_buttons">
				<IconButton iconFilename="biotech.png" label="Start Analysis" />
				<IconButton iconFilename="quick_reference.png" label="Guide / FAQ" />
			</div>
		</div>
	);
};

export default Home;
