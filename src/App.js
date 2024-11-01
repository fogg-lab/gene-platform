import React from "react";
import { HashRouter as Router, Route, Routes, Link } from 'react-router-dom';
import Home from './pages/Home'
import Analysis from './pages/Analysis'
import Guide from './pages/Guide'
import ViewRuns from './pages/ViewRuns'
import './assets/styles.css';
import logo from './assets/fogg_logo.png';
import { ErrorPopupProvider } from "./components/ui/ErrorPopup";

const App = () => {
	return (
		<ErrorPopupProvider>
			<Router>
				<nav>
					<ul>
						<li>
							<Link to="/">
								<img src={logo} id="fogg_logo" alt="Fogg Labs Logo" />
							</Link>
						</li>
					</ul>
					<ul id="primary_navbar">
						<li>
							<Link to="/Analysis">Analysis</Link>
						</li>
						<li>
							<Link to="/Guide">Guide</Link>
						</li>
						<li>
							<Link to="/ViewRuns">View Runs</Link>
						</li>
					</ul>
				</nav>
				<Routes>
					<Route path="/" element={<Home />} />
					<Route path="/Analysis" element={<Analysis />} />
					<Route path="/Guide" element={<Guide />} />
					<Route path="/ViewRuns" element={<ViewRuns />} />
				</Routes>
			</Router>
		</ErrorPopupProvider>
	);
};

export default App;
