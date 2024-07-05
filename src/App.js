import React from 'react';
import { BrowserRouter as Router, Route, Routes, Link } from 'react-router-dom';
import Home from './pages/Home'
import Analysis from './pages/Analysis'
import Guide from './pages/Guide'
import ViewRuns from './pages/ViewRuns'
import './assets/styles.css';
import logo from './assets/fogg_logo.png';

class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false };
  }

  static getDerivedStateFromError(error) {
    return { hasError: true };
  }

  componentDidCatch(error, errorInfo) {
    console.log('Error caught by ErrorBoundary:', error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      return <h1>Something went wrong.</h1>;
    }

    return this.props.children;
  }
}

const App = () => {
  console.log('App component rendering');
  return (
    <ErrorBoundary>
      <Router>
        <div>
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
        </div>
      </Router>
    </ErrorBoundary>
  );
};

export default App;