import React from 'react';
import PropTypes from 'prop-types';

const IconButton = ({ iconFilename, label, onClick }) => {
	let iconSrc;

	try {
		iconSrc = require(`../../assets/icons/${iconFilename}`).default;
	} catch (error) {
		console.error(`Icon ${iconFilename} not found in assets/icons`);
		return null;
	}

	const iconStyle = {
		width: '24px',
		height: '24px',
		marginRight: '8px',
	};

	const buttonStyle = {
		backgroundColor: '#D73F09',
		display: 'flex',
		alignItems: 'center',
		justifyContent: 'center',
		padding: '5px',
		border: 'none',
		color: 'white',
		cursor: 'pointer',
		borderRadius: '3px',
		fontSize: '16px'
	};

	return (
		<button className="icon-button" style={buttonStyle} onClick={onClick}>
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
