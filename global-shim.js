// Set up global object if it doesn't exist
if (typeof global === 'undefined') {
    window.global = window;
}

// Load environment variables from electron when available
if (window.electron?.getEnvVar) {
    window.electron.getEnvVar('REACT_APP__MUI_KEY').then(key => {
        if (key) {
            process.env.REACT_APP__MUI_KEY = key;
        }
    });
}
