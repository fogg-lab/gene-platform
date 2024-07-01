Note: This is a rough readme and stand-in for the future readme.md

Prerequisites: 
1. NPM installed
2. Node.js installed

Install and run steps:
1. install NPM dependencies: npm install
2. start webpack dev server: npm start
3. run electron application: npm run electron

This should launch the electron/react window.

Project Structure

```
electron_app/
│
├── src/
│   ├── index.js // Main entry for React / where root gets rendered
│   ├── App.js //routing / state management
│   └── styles.css (not implemented)
│
├── public/
│   └── index.html //root React element
│
├── dist/
│   └── (output files after build)
│
├── main.js //react + electron bundler
├── package.json
├── webpack.config.js
├── .babelrc
└── README.md
```
