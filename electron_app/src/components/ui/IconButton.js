import React from 'react';
import PropTypes from 'prop-types';

const IconButton = ({ iconFilename, label }) => {
	let iconSrc;

	try {
		iconSrc = require(`../../assets/icons/${iconFilename}`).default;
	} catch (error) {
		console.error(`Icon ${iconFilename} not found in assets/icons`);
		return null;
	}

	console.log(`Icon source for ${iconFilename}:`, iconSrc);

	const iconStyle = {
		width: '24px',
		height: '24px',
		marginRight: '8px',
	};

	const buttonStyle = {
		backgroundColor: '#D73F09',
		display: 'flex', // Use flexbox for alignment
		alignItems: 'center', // Center vertically
		justifyContent: 'center', // Center horizontally if needed
		padding: '5px',
		border: 'none',
		color: 'white',
		cursor: 'pointer',
		borderRadius: '3px',
		fontSize: '16px'
	};

	return (
		<button className="icon-button" style={buttonStyle}>
			<img src={iconSrc} alt={label} style={iconStyle} />
			<span>{label}</span>
		</button>
	);
};

IconButton.propTypes = {
	iconFilename: PropTypes.string.isRequired,
	label: PropTypes.string.isRequired,
};

export default IconButton;
