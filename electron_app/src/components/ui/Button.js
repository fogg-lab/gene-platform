// src/components/Button.js
import React, { useState } from 'react';

function Button()
{
	const [count, setCount] = useState(0);

	const handleClick = () =>{
		setCount(count + 1);
	};

	return (
		<div>
			<button onClick={handleClick}>
				Clicked {count} times
			</button>
		</div>
	)
}

export default Button;
