import React from 'react';
import PropTypes from 'prop-types';

const IconButton = ({ icon, label, onClick }) => {
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
			<img src={icon} alt={label} style={iconStyle} />
			<span>{label}</span>
		</button>
	);
};

IconButton.propTypes = {
	iconFilename: PropTypes.string.isRequired,
	label: PropTypes.string.isRequired,
};

export default IconButton;
