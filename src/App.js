import React from 'react';
import Link from 'next/link';
import Image from 'next/image';

const App = () => {
	return (
		<div>
			<nav>
				<ul>
					<li>
						<Link href="/">
							<Image src="/fogg_logo.png" alt="Fogg Labs Logo" id="fogg_logo" width={100} height={50} />
						</Link>
					</li>
				</ul>
				<ul id="primary_navbar">
					<li>
						<Link href="/analysis">Analysis</Link>
					</li>
					<li>
						<Link href="/guide">Guide</Link>
					</li>
					<li>
						<Link href="/view-runs">View Runs</Link>
					</li>
				</ul>
			</nav>
		</div>
	);
};

export default App;
